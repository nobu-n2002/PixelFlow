! ibm_3d_omp_cpu.f90
! link: lib/global.f90 lib/output.f90 lib/utils.f90

program main
  !$ use omp_lib
  use global_3d
  use valiables
  use output_3d
  use grid_3d
  use utils
  implicit none
  integer:: istep
  real,dimension(0:md,0:nd,0:ld):: u, v, w, p, u_old, v_old, w_old
  real,dimension(0:md,0:nd,0:ld):: porosity
  real:: dx, dy, dz, dt
  real,dimension(0:md):: xp
  real,dimension(0:nd):: yp
  real,dimension(0:ld):: zp
  integer:: m, n, l
  integer:: i, j, k
  !--- namelist valiables Declaration in valiables module
  ! real:: xnue, xlambda, density, width, height, depth, time
  ! real:: inlet_velocity, outlet_pressure, AoA
  ! integer:: istep_max, istep_out
  ! real:: thickness, threshold, radius, center_x, center_y, center_z
  ! logical:: nonslip
  ! character(len=50) :: output_folder
  ! character(len=50) :: csv_file
  ! integer:: iter_max
  ! real:: relux_factor
  ! ----------------

  call get_now_time()
  
  ! ----------------
  ! calculation start
  ! ----------------
  ! set up condistions
  m=0         ! setup switch for grid conditions
  density=0.  ! setup switch for physical conditions

  call read_settings(&
  xnue, xlambda, density, width, height, depth, time,&
  inlet_velocity, outlet_pressure, AoA,&
  istep_max, istep_out,&
  thickness, threshold, radius, center_x, center_y, center_z,&
  nonslip,&
  output_folder,csv_file,&
  iter_max, relux_factor)

  call system('mkdir -p '//trim(output_folder))
  call system('mkdir -p etc')

  call grid_conditions_yz_periodic (&
  xp, yp, zp, dx, dy, dz, dt, xnue, xlambda, density, width, height, depth, &
  thickness, threshold, radius,&
  center_x, center_y, center_z, time,&
  inlet_velocity, AoA, porosity, m, n, l, istep_max, &
  csv_file)

  call output_grid_3d (xp, yp, zp, m, n, l)  
  istep = 0
  time = istep * dt
  ! ----------------
  
  write(*,*) '# istep_max= ', istep_max,'   istep_out= ', istep_out
  
  call initial_conditions (p, u, v, w, xp, yp, zp, width, height, depth &
                        , inlet_velocity, outlet_pressure, AoA, m, n, l)
  call boundary (p, u, v, w, xp, yp, zp, width, height, depth &
                         , inlet_velocity, outlet_pressure, AoA, porosity, m, n, l)
  
  ! print initial conditions
  ! call  output_solution (p, u, v, w, m, n, l)
  
  ! ----------------
  ! MAC algorithm start
  call get_now_time()
  write(*,*) '# --- MAC algorithm start'

  do istep = 1, istep_max
    
    time = istep * dt
    
    !$omp parallel private(i, j, k) &
    !$omp & shared(m, n, l) &
    !$omp & shared(u_old, v_old, w_old, u, v, w) &
    !$omp & default(none)
    !$omp do
    do k = 0, l+1
      do j = 0, n+1
        do i = 0, m+1
          u_old(i,j,k) = u(i,j,k)
          v_old(i,j,k) = v(i,j,k)
          w_old(i,j,k) = w(i,j,k)
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
  
    write(*,*)'--- time_steps= ',istep, ' --  time = ',time
  
    call solve_p (p, u, v, w, u_old, v_old, w_old, porosity, &
    xnue, xlambda, density, height, thickness, &
    yp, dx, dy, dz, dt, m, n, l, &
    nonslip, iter_max, relux_factor)
  
    !-- solve u, v, w (fix u, v, w)
    !$omp parallel private(i, j, k) &
    !$omp & shared(m, n, l, dt, dx, dy, dz, density) &
    !$omp & shared(p, u, v, w) &
    !$omp & default(none)
    !$omp do
    do k = 1, l
      do j = 1, n
        do i = 1, m
          u(i,j,k) = u(i,j,k) - dt/density*(p(i+1,j,k)-p(i-1,j,k))/dx*0.5
          v(i,j,k) = v(i,j,k) - dt/density*(p(i,j+1,k)-p(i,j-1,k))/dy*0.5
          w(i,j,k) = w(i,j,k) - dt/density*(p(i,j,k+1)-p(i,j,k-1))/dz*0.5
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel
  
    call boundary(p, u, v, w, xp, yp, zp, width, height, depth  &
                      , inlet_velocity, outlet_pressure, AoA, porosity, m, n, l)
    
    if(mod(istep,istep_out)==0) call  output_paraview_temp_3d (p, u, v, w, porosity, xp, yp, zp, m, n, l, istep)
    
  end do
  call get_now_time()
  ! MAC algorithm end
  ! ----------------
  
  ! print solutions 
  call  output_solution_post_3d (p, u, v, w, xp, yp, zp, porosity, m, n, l)
  call  output_divergent_3d (p, u, v, w, porosity, dx, dy, dz, m, n, l)
  call  output_paraview_3d (p, u, v, w, porosity, xp, yp, zp, m, n, l)
  
  write(*,*) 'program finished'
  call get_now_time()
  
end program main
!******************
  
!  solve variables  

!******************
subroutine  solve_p (p, u, v, w, u_old, v_old, w_old, porosity, &
                     xnue, xlambda, density, height, thickness, &
                     yp, dx, dy, dz, dt, m, n, l, &
                     nonslip, iter_max, relux_factor)
  use global_3d
  implicit none
  real,intent(in):: dx, dy, dz, dt
  real,intent(in):: xnue, xlambda, density, height, thickness
  real,intent(inout),dimension(0:md,0:nd,0:ld):: u, v, w, p, u_old, v_old, w_old 
  real,intent(in),dimension(0:md,0:nd,0:ld):: porosity
  real,intent(in),dimension(0:nd) :: yp
  integer,intent(in):: m, n, l
  logical,intent(in):: nonslip
  integer,intent(in):: iter_max
  real,intent(in):: relux_factor
 
  !-----------------
  ! local variables 
  real, parameter :: small = 1.e-6
  real, parameter :: alpha = 32.0
 
  real,dimension(0:md,0:nd,0:ld) :: ap, ae, aw, an, as, at, ab, bb, div
  integer :: i, j, k
  
  !$omp parallel private(i, j, k) &
  !$omp & shared(dx, dy, dz, dt, xnue, density, height, thickness, m, n, l) &
  !$omp & shared(p, u, v, w, u_old, v_old, w_old, porosity, yp) &
  !$omp & shared(xlambda, nonslip) &
  !$omp & shared(ap, ae, aw, an, as, at, ab, bb, div) &
  !$omp & default(none)

  !-----------------
  !  divergence term  div(u)
  !-----------------
  !$omp do
  do k = 1, l
    do j = 1, n
      do i = 1, m
        div(i,j,k)= (u_old(i+1,j,k)-u_old(i-1,j,k))/dx*0.5 &
                  + (v_old(i,j+1,k)-v_old(i,j-1,k))/dy*0.5 &
                  + (w_old(i,j,k+1)-w_old(i,j,k-1))/dz*0.5 
      end do
    end do
  end do
  !$omp end do
  
  !$omp do
  do k = 1, l
    do j = 1, n
      div(0,j,k)  = 0.  ! inlet yz
      div(m+1,j,k)= 0.  ! outlet yz
    end do
  end do
  !$omp end do
 
  !$omp do 
  do k = 1, l
    do i = 1, m
      div(i,0,k)  = div(i,n,k)  ! periodic condition xz
      div(i,n+1,k)= div(i,1,k)
    end do
  end do
  !$omp end do
  
  !$omp do
  do j = 1, n
    do i = 1, m
      div(i,j,0)  = div(i,j,l)  ! periodic condition xy
      div(i,j,l+1)= div(i,j,1)
    end do
  end do
  !$omp end do
  

  ! ----------------
  !   velocity u
  ! ----------------
  !$omp do
  do k = 1, l
    do j = 1, n
      do i = 1, m
        !-- convection_x  (2nd central scheme)
        u(i,j,k) = u_old(i,j,k) &
                   - dt*u_old(i,j,k)*(u_old(i+1,j,k)-u_old(i-1,j,k))/dx*0.5
        
        !-- convection_y  (2nd central scheme)
        u(i,j,k) = u(i,j,k) &
                   - dt*v_old(i,j,k)*(u_old(i,j+1,k)-u_old(i,j-1,k))/dy*0.5
    
        !-- convection_w  (2nd central scheme)
        u(i,j,k) = u(i,j,k) &
                   - dt*w_old(i,j,k)*(u_old(i,j,k+1)-u_old(i,j,k-1))/dz*0.5
        
        !-- diffusion_x
        u(i,j,k)=u(i,j,k) + dt*xnue*(u_old(i+1,j,k)-2.*u_old(i,j,k)+u_old(i-1,j,k))/dx/dx 
                            
        !-- diffusion_y     
        u(i,j,k)=u(i,j,k) + dt*xnue*(u_old(i,j+1,k)-2.*u_old(i,j,k)+u_old(i,j-1,k))/dy/dy
                            
        !-- diffusion_z     
        u(i,j,k)=u(i,j,k) + dt*xnue*(u_old(i,j,k+1)-2.*u_old(i,j,k)+u_old(i,j,k-1))/dz/dz
                            
        !-- divergence term 
        u(i,j,k)=u(i,j,k) + dt*(xnue + xlambda)*(div(i+1,j,k)-div(i-1,j,k))/dx*0.5
    
        !-- additional terms by porosity profile
        u(i,j,k) = u(i,j,k) &
                   + dt*( &
                          ((u_old(i+1,j,k)-u_old(i-1,j,k))/dx*0.5+(u_old(i+1,j,k)-u_old(i-1,j,k))/dx*0.5)  &
                            *xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*0.5                             &
                          +((u_old(i,j+1,k)-u_old(i,j-1,k))/dy*0.5+(v_old(i+1,j,k)-v_old(i-1,j,k))/dx*0.5) &
                            *xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*0.5                             &
                          +((u_old(i,j,k+1)-u_old(i,j,k-1))/dz*0.5+(w_old(i+1,j,k)-w_old(i-1,j,k))/dx*0.5) &
                            *xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*0.5                             &
                          + div(i,j,k)*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*0.5*xlambda                 &
                     )/porosity(i,j,k)
  
        !-- force on wall
        if(nonslip) then
          u(i,j,k)=u(i,j,k)-dt*xnue*u_old(i,j,k)/(thickness*dx)**2*alpha*porosity(i,j,k)*(1.-porosity(i,j,k))*(1.-porosity(i,j,k))
        end if

      end do
    end do
  end do
  !$omp end do

  ! ----------------
  !   velocity v
  ! ----------------
  !$omp do
  do k = 1, l
    do j = 1, n
      do i = 1, m
        !-- convection_x  (2nd central scheme)
        v(i,j,k) = v_old(i,j,k) &
                   - dt*u_old(i,j,k)*(v_old(i+1,j,k)-v_old(i-1,j,k))/dx*0.5
    
        !-- convection_y (2nd central scheme)
        v(i,j,k) = v(i,j,k) &
                   - dt*v_old(i,j,k)*(v_old(i,j+1,k)-v_old(i,j-1,k))/dy*0.5
    
        !-- convection_z (2nd central scheme)
        v(i,j,k) = v(i,j,k) &
                   - dt*w_old(i,j,k)*(v_old(i,j,k+1)-v_old(i,j,k-1))/dz*0.5
        
        !-- diffusion_x
        v(i,j,k)=v(i,j,k) + dt*xnue*(v_old(i+1,j,k)-2.*v_old(i,j,k)+v_old(i-1,j,k))/dx/dx
                            
        !-- diffusion_y     
        v(i,j,k)=v(i,j,k) + dt*xnue*(v_old(i,j+1,k)-2.*v_old(i,j,k)+v_old(i,j-1,k))/dy/dy
                            
        !-- diffusion_z     
        v(i,j,k)=v(i,j,k) + dt*xnue*(v_old(i,j,k+1)-2.*v_old(i,j,k)+v_old(i,j,k-1))/dz/dz

        !-- divergence term
        v(i,j,k)=v(i,j,k) + dt*(xnue + xlambda)*(div(i,j+1,k)-div(i,j-1,k))/dy*0.5

        !-- additional terms by porosity profile
        v(i,j,k) = v(i,j,k) &
                   + dt*( &
                         ((v_old(i+1,j,k)-v_old(i-1,j,k))/dx*0.5+(u_old(i,j+1,k)-u_old(i,j-1,k))/dy*0.5) &
                           *xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*0.5                            &
                         +((v_old(i,j+1,k)-v_old(i,j-1,k))/dy*.5+(v_old(i,j+1,k)-v_old(i,j-1,k))/dy*0.5) &
                           *xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*0.5                            &
                         +((v_old(i,j,k+1)-v_old(i,j,k-1))/dz*.5+(w_old(i,j+1,k)-w_old(i,j-1,k))/dy*0.5) &
                           *xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*0.5                            &
                         + div(i,j,k)*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*0.5*xlambda                &
                     )/porosity(i,j,k)
        !-- force on wall
        if(nonslip) then
          v(i,j,k)=v(i,j,k)-dt*xnue*v_old(i,j,k)/(thickness*dy)**2*alpha*porosity(i,j,k)*(1.-porosity(i,j,k))*(1.-porosity(i,j,k))
        end if

      end do
    end do
  end do
  !$omp end do
  

  ! ----------------
  !   velocity w
  ! ----------------
  !$omp do
  do k = 1, l
    do j = 1, n
      do i = 1, m
        !-- convection_x  (2nd central scheme)
        w(i,j,k) = w_old(i,j,k) &
                   - dt*u_old(i,j,k)*(w_old(i+1,j,k)-w_old(i-1,j,k))/dx*0.5
        
        !-- convection_y  (2nd central scheme)
        w(i,j,k) = w(i,j,k) &
                   - dt*v_old(i,j,k)*(w_old(i,j+1,k)-w_old(i,j-1,k))/dy*0.5
    
        !-- convection_z  (2nd central scheme)
        w(i,j,k) = w(i,j,k) &
                   - dt*w_old(i,j,k)*(w_old(i,j,k+1)-w_old(i,j,k-1))/dz*0.5
        
        !-- diffusion_x
        w(i,j,k)=w(i,j,k) +dt*xnue*(w_old(i+1,j,k)-2.*w_old(i,j,k)+w_old(i-1,j,k))/dx/dx

        !-- diffusion_y
        w(i,j,k)=w(i,j,k) +dt*xnue*(w_old(i,j+1,k)-2.*w_old(i,j,k)+w_old(i,j-1,k))/dy/dy

        !-- diffusion_z
        w(i,j,k)=w(i,j,k) +dt*xnue*(w_old(i,j,k+1)-2.*w_old(i,j,k)+w_old(i,j,k-1))/dz/dz

        !-- divergence term
        w(i,j,k)=w(i,j,k) +dt*(xnue + xlambda)*(div(i,j,k+1)-div(i,j,k-1))/dz*0.5

        !-- additional terms by porosity profile
        w(i,j,k) = w(i,j,k) &
                   + dt*( &
                          ((w_old(i+1,j,k)-w_old(i-1,j,k))/dx*0.5+(u_old(i,j,k+1)-u_old(i,j,k-1))/dz*0.5)  &
                            *xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*0.5                             &
                          +((w_old(i,j+1,k)-w_old(i,j-1,k))/dy*0.5+(v_old(i,j,k+1)-v_old(i,j,k-1))/dz*0.5) &
                            *xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*0.5                             &
                          +((w_old(i,j,k+1)-w_old(i,j,k-1))/dz*0.5+(w_old(i,j,k+1)-w_old(i,j,k-1))/dz*0.5) &
                            *xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*0.5                             &
                          + div(i,j,k)*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*0.5*xlambda                 &
                     )/porosity(i,j,k)
        ! force on wall
        if (nonslip) then
          w(i,j,k)=w(i,j,k)- dt*xnue*w_old(i,j,k)/(thickness*dz)**2*alpha*porosity(i,j,k)*(1.-porosity(i,j,k))*(1.-porosity(i,j,k))
        end if
    
      end do
    end do
  end do
  !$omp end do
  

  ! ----------------
  ! matrix solution  !  formulation of porous media
  !$omp do
  do k = 1, l
    do j = 1, n
      do i = 1, m
        ae(i,j,k) = dt*max(small,(porosity(i+1,j,k)+porosity(i,j,k))*0.5)/dx/dx
    
        aw(i,j,k) = dt*max(small,(porosity(i,j,k)+porosity(i-1,j,k))*0.5)/dx/dx
    
        an(i,j,k) = dt*max(small,(porosity(i,j+1,k)+porosity(i,j,k))*0.5)/dy/dy
    
        as(i,j,k) = dt*max(small,(porosity(i,j,k)+porosity(i,j-1,k))*0.5)/dy/dy
    
        at(i,j,k) = dt*max(small,(porosity(i,j,k+1)+porosity(i,j,k))*0.5)/dz/dz
    
        ab(i,j,k) = dt*max(small,(porosity(i,j,k)+porosity(i,j,k-1))*0.5)/dz/dz
    
        ap(i,j,k) = -ae(i,j,k)-aw(i,j,k)-an(i,j,k)-as(i,j,k)-at(i,j,k)-ab(i,j,k)
    
        bb(i,j,k) = ((porosity(i+1,j,k)*u(i,j,k)+porosity(i,j,k)*u(i+1,j,k))*0.5              &
                     -(porosity(i-1,j,k)*u(i,j,k)+porosity(i,j,k)*u(i-1,j,k))*0.5)*density/dx &
                    +((porosity(i,j+1,k)*v(i,j,k)+porosity(i,j,k)*v(i,j+1,k))*0.5             &
                     -(porosity(i,j-1,k)*v(i,j,k)+porosity(i,j,k)*v(i,j-1,k))*0.5)*density/dy &
                    +((porosity(i,j,k+1)*w(i,j,k)+porosity(i,j,k)*w(i,j,k+1))*0.5             &
                     -(porosity(i,j,k-1)*w(i,j,k)+porosity(i,j,k)*w(i,j,k-1))*0.5)*density/dz 
    
      end do
    end do
  end do
  !$omp end do
  
  !$omp end parallel

  call boundrary_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, l, height, yp)
  
  ! call solve_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, l)
  call solve_matrix_vec_omp (p, ap, ae, aw, an, as, at, ab, bb, &
  m, n, l, iter_max, relux_factor)
  ! call solve_matrix_vec_oacc (p, ap, ae, aw, an, as, at, ab, bb, &
  ! m, n, l, iter_max, relux_factor)
  
end subroutine solve_p
!******************

!******************
! OpenMP Parallelized
! Written only for CPU machine
! No efficiency ensured on GPU machine 
subroutine  solve_matrix_vec_omp (p, ap, ae, aw, an, as, at, ab, bb, &
  m, n, l, iter_max, relux_factor)
  use global_3d
  implicit none
  real,intent(inout),dimension(0:md,0:nd,0:ld):: p
  real,intent(in),dimension(0:md,0:nd,0:ld):: ap, ae, aw, an, as, at, ab, bb
  integer,intent(in):: m, n, l
  integer,intent(in):: iter_max
  real,intent(in):: relux_factor

  ! local variables
  real:: error
  real,dimension(0:md,0:nd,0:ld):: p_old
  integer::i, j, k, iter, ii

  !$omp parallel private(iter, i, j, k, ii) &
  !$omp & shared(iter_max, relux_factor, m, n, l) &
  !$omp & shared(error, p_old, p, ap, ae, aw, an, as, at, ab, bb) &
  !$omp & default(none)

  ! ----------------
  !   SOR algorithm
  ! ----------------
  !$omp single
  error = 0.0
  !$omp end single

  do iter = 1, iter_max

  ! default periodic condition in xz-direction
  !$omp do
  do i = 1, m
    do k = 1, l
      p(i,0,k)   = p(i,n,k)
      p(i,n+1,k) = p(i,1,k)
    end do
  end do
  !$omp end do

  ! default periodic condition in xy-direction
  !$omp do
  do i = 1, m
    do j = 1, n
      p(i,j,0)   = p(i,j,l)
      p(i,j,l+1) = p(i,j,1)
    end do
  end do
  !$omp end do

  !$omp do
  do i = 0, m+1
    do j = 0, n+1
      do k = 0, l+1
        p_old(i,j,k) = p(i,j,k)
      end do
    end do
  end do
  !$omp end do

  !$omp do
  do ii = 2, m*n*l, 2 ! evenspace
    k = (ii - 1) / (m * n)  + 1
    j = ((ii - 1) / m + 1) - (k - 1) * n
    i = (ii - (j - 1) * m) - (k - 1) * m * n

    !--- IF (m, n) is (ODD, EVEN) (Based on Column-Major Order; FORTRAN)
    if (mod(m,2) /= 0 .and. mod(n,2)==0 .and. mod(k,2)==0) then
      if (mod(i,2) /= 0) then 
          j = j - 1
      else 
          j = j + 1
      end if
    ! --- IF (m, n) is (EVEN, EVEN) or (EVEN, ODD) (Based on Column-Major Order; FORTRAN)
    else if (mod(m,2) == 0 .and. (mod(j,2) + mod(k,2) == 1)) then
        i = i - 1
    end if

    p(i,j,k) = (bb(i,j,k) &
              - ae(i,j,k)*p_old(i+1,j,k)-aw(i,j,k)*p_old(i-1,j,k)  &
              - an(i,j,k)*p_old(i,j+1,k)-as(i,j,k)*p_old(i,j-1,k)  &
              - at(i,j,k)*p_old(i,j,k+1)-ab(i,j,k)*p_old(i,j,k-1)) &
              / ap(i,j,k)*relux_factor &
              + p_old(i,j,k)*(1.-relux_factor)
  end do
  !$omp end do

  ! default periodic condition in xz-direction
  !$omp do
  do i = 1, m
    do k = 1, l
      p(i,0,k)   = p(i,n,k)
      p(i,n+1,k) = p(i,1,k)
    end do
  end do
  !$omp end do

  ! default periodic condition in xy-direction
  !$omp do
  do i = 1, m
    do j = 1, n
      p(i,j,0)   = p(i,j,l)
      p(i,j,l+1) = p(i,j,1)
    end do
  end do
  !$omp end do

  !$omp do
  do i = 0, m+1
    do j = 0, n+1
      do k = 0, l+1
        p_old(i,j,k) = p(i,j,k)
      end do
    end do
  end do
  !$omp end do

  !$omp do
  do ii = 1, m*n*l, 2 ! odd space
    k = (ii - 1) / (m * n)  + 1
    j = ((ii - 1) / m + 1) - (k - 1) * n
    i = (ii - (j - 1) * m) - (k - 1) * m * n

    !--- IF (m, n) is (ODD, EVEN) (Based on Column-Major Order; FORTRAN)
    if (mod(m,2) /= 0 .and. mod(n,2)==0 .and. mod(k,2)==0) then
      if (mod(i,2) /= 0) then 
          j = j + 1
      else 
          j = j - 1
      end if
    ! --- IF (m, n) is (EVEN, EVEN) or (EVEN, ODD) (Based on Column-Major Order; FORTRAN)
    else if (mod(m,2) == 0 .and. (mod(j,2) + mod(k,2) == 1)) then
        i = i + 1
    end if
    p(i,j,k) = (bb(i,j,k) &
              - ae(i,j,k)*p_old(i+1,j,k)-aw(i,j,k)*p_old(i-1,j,k)  &
              - an(i,j,k)*p_old(i,j+1,k)-as(i,j,k)*p_old(i,j-1,k)  &
              - at(i,j,k)*p_old(i,j,k+1)-ab(i,j,k)*p_old(i,j,k-1)) &
              / ap(i,j,k)*relux_factor &
              + p_old(i,j,k)*(1.-relux_factor)
  end do
  !$omp end do

  !$omp do reduction(max:error)
  do i = 1, m
    do j = 1, n
      do k = 1, l
        error = max(error, abs(p(i,j,k)-p_old(i,j,k)))
      end do
    end do
  end do
  !$omp end do

  end do

  ! default periodic condition in xz-direction
  !$omp do
  do i = 1, m
    do k = 1, l
      p(i,0,k)  =p(i,n,k)
      p(i,n+1,k)=p(i,1,k)
    end do
  end do
  !$omp end do

  ! default periodic condition in xy-direction
  !$omp do
  do i = 1, m
    do j = 1, n
      p(i,j,0)  =p(i,j,l)
      p(i,j,l+1)=p(i,j,1)
    end do
  end do
  !$omp end do

  !$omp master
  write(*,*) 'SOR iteration no.', iter_max, '-- p error:', error
  !$omp end master

  !$omp end parallel

  return
end subroutine solve_matrix_vec_omp
!******************

!******************
subroutine  boundrary_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, l, height, yp)

  use global_3d
  implicit none
  real,intent(in)::height
  real,intent(in),dimension(0:md,0:nd,0:ld)::p
  real,intent(inout),dimension(0:md,0:nd,0:ld)::ap, ae, aw, an, as, at, ab, bb
  real,intent(in),dimension(0:nd)::yp
  integer,intent(in)::m, n, l

  ! local variables
  integer::i, j, k
  
  !$omp parallel private(i, j, k) &
  !$omp & shared(m, n, l, height) &
  !$omp & shared(p, ap, ae, aw, an, as, at, ab, bb, yp) &
  !$omp & default(none)

  ! inlet (dp/x=0 at i=1)
  !$omp do
  do k = 1, l
    do j = 1, n
      ae(1,j,k) =ae(1,j,k)+aw(1,j,k)
      aw(1,j,k) =0.
    end do
  end do
  !$omp end do
  
  ! outlet (p=outlet_pressure at i=m)
  !$omp do
  do k = 1, l
    do j = 1, n
      bb(m,j,k) = bb(m,j,k)+ae(m,j,k)*p(m+1,j,k)
      ae(m,j,k) = 0.
      aw(m,j,k) = 0.
      an(m,j,k) = 0.
      as(m,j,k) = 0.
      at(m,j,k) = 0.
      ab(m,j,k) = 0.
    end do
  end do
  !$omp end do

  !$omp end parallel
  
end subroutine  boundrary_matrix
!******************

!  conditions  

!******************
subroutine  boundary(p, u, v, w, xp, yp, zp, width, height, depth    &
                     , inlet_velocity, outlet_pressure, AoA, porosity, m, n, l)
  use global_3d
  implicit none
  real,intent(in)::width, height, depth, inlet_velocity, outlet_pressure, AoA
  real,intent(inout),dimension(0:md,0:nd,0:ld)::u, v, w, p
  real,intent(in),dimension(0:md,0:nd,0:ld)::porosity
  real,intent(in),dimension(0:md)::xp
  real,intent(in),dimension(0:nd)::yp
  real,intent(in),dimension(0:ld)::zp
  integer,intent(in)::m, n, l

  ! local variables
  real, parameter :: pi=atan(1.)*4.
  integer::i, j, k
  
  !$omp parallel private(i, j, k) &
  !$omp & shared(p, u, v, w, xp, yp, zp) &
  !$omp & shared(width, height, depth, inlet_velocity, outlet_pressure) &
  !$omp & shared(AoA, porosity, m, n, l) &
  !$omp & default(none)

  ! inlet (u=inlet_velocity, v=0., dp/dx=0 at i=1)
  !$omp do
  do k = 1, l
    do j = 1, n
      u(1,j,k) = inlet_velocity*cos(AoA/1300.*pi)
      v(1,j,k) = inlet_velocity*sin(AoA/1300.*pi)
      w(1,j,k) = 0.
      u(0,j,k) = u(1,j,k)    ! dummy
      v(0,j,k) = v(1,j,k)    ! dummy
      w(0,j,k) = w(1,j,k)    ! dummy
      p(0,j,k) = p(2,j,k)
    end do
  end do
  !$omp end do
  
  ! outlet (du/dx=0., dv/dx=0., p=outlet_pressure at i=m)
  !$omp do
  do k = 1, l
    do j = 1, n
      u(m+1,j,k) = u(m-1,j,k)
      v(m+1,j,k) = v(m-1,j,k)
      w(m+1,j,k) = w(m-1,j,k)
      p(m+1,j,k) = outlet_pressure   ! dummy
    end do
  end do
  !$omp end do

  ! default: periodic condition (xz-direction at j=1 & n)
  !$omp do
  do k = 0, l+1
    do i = 0, m+1
      u(i,0,k)   = u(i,n,k)
      v(i,0,k)   = v(i,n,k)
      w(i,0,k)   = w(i,n,k)
      p(i,0,k)   = p(i,n,k)
      u(i,n+1,k) = u(i,1,k)
      v(i,n+1,k) = v(i,1,k)
      w(i,n+1,k) = w(i,1,k)
      p(i,n+1,k) = p(i,1,k)
    end do
  end do
  !$omp end do

  ! default: periodic condition (xy-direction at k=1 & l)
  !$omp do
  do j = 0, n+1
    do i = 0, m+1
      u(i,j,0)   = u(i,j,l)
      v(i,j,0)   = v(i,j,l)
      w(i,j,0)   = w(i,j,l)
      p(i,j,0)   = p(i,j,l)
      u(i,j,l+1) = u(i,j,1)
      v(i,j,l+1) = v(i,j,1)
      w(i,j,l+1) = w(i,j,1)
      p(i,j,l+1) = p(i,j,1)
    end do
  end do
  !$omp end do

  !$omp end parallel
  
end subroutine boundary
!*****************************

!******************
subroutine  initial_conditions (p, u, v, w, xp, yp, zp, width, height, depth  &
                               , inlet_velocity, outlet_pressure, AoA, m, n, l)

  use global_3d
  implicit none
  real,intent(in)::width, height, depth, inlet_velocity, outlet_pressure, AoA
  real,intent(out),dimension(0:md,0:nd,0:ld)::u, v, w, p 
  real,intent(in),dimension(0:md)::xp
  real,intent(in),dimension(0:nd)::yp
  real,intent(in),dimension(0:ld)::zp
  integer,intent(in)::m, n, l

  ! local variables
  integer::i, j, k
  real,parameter:: pi=atan(1.)*4.   
  
  !$omp parallel private(i, j, k) &
  !$omp & shared(p, u, v, w, xp, yp, zp) &
  !$omp & shared(inlet_velocity, outlet_pressure, AoA, m, n, l) &
  !$omp & shared(width, height, depth) &
  !$omp & default(none)

  !$omp do
  do k = 1, l
    do j = 1, n
      do i = 1, m
        u(i,j,k) = inlet_velocity * cos(AoA/360*pi)
        v(i,j,k) = inlet_velocity * sin(AoA/360*pi)
        w(i,j,k) = 0.
        p(i,j,k) = outlet_pressure
      end do
    end do
  end do
  !$omp end do

  !$omp end parallel

end subroutine initial_conditions
!******************