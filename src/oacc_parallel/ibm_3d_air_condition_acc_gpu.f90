! ibm_3d_air_condition_acc_gpu.f90
! link: lib/global.f90 lib/output.f90 lib/utils.f90

module wall_conditions
  implicit none
  ! wall condition:
  ! 0: wall 
  ! 1: inlet (@poro>=0.9 then inlet, @poro<0.9 then wall)
  ! 2: outlet (@poro>=0.9 then p = outlet_pressure, @poro<0.9 then wall)
  integer, parameter:: top_wall     = 1
  integer, parameter:: bottom_wall  = 0
  integer, parameter:: east_wall    = 0
  integer, parameter:: west_wall    = 0
  integer, parameter:: south_wall   = 2
  integer, parameter:: north_wall   = 0
end module wall_conditions

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

  call grid_conditions_wall (&
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
  call  output_paraview_temp_3d (p, u, v, w, porosity, xp, yp, zp, m, n, l, 0)  
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
  ! ----------------
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

  ! --- wall condition xy
  !$omp do
  do i = 1, m
    do j = 1, n
      div(i,j,0)  = 0.
      div(i,j,l+1)= 0.
    end do
  end do
  !$omp end do
  ! ---
  ! --- wall condition yz
  !$omp do
  do j = 1, n
    do k = 1, l
      div(0,j,k)  = 0.
      div(m+1,j,k)= 0.
    end do
  end do
  !$omp end do
  ! --- wall condition zx
  !$omp do
  do i = 1, m
    do k = 1, l
      div(i,0,k)  = 0.
      div(i,n+1,k)= 0.
    end do
  end do
  !$omp end do
  ! ---

  ! --- periodic condition xy
  ! !$omp do
  ! do i = 1, m
  !   do j = 1, n
  !     div(i,j,0)  = div(i,j,l)
  !     div(i,j,l+1)= div(i,j,1)
  !   end do
  ! end do
  ! !$omp end do
  ! ---
  ! --- periodic condition yz
  ! !$omp do
  ! do j = 1, n
  !   do k = 1, l
  !     div(0,j,k)  = div(m,j,k)
  !     div(m+1,j,k)= div(1,j,k)
  !   end do
  ! end do
  ! !$omp end do
  ! ---  
  ! --- periodic condition zx
  ! !$omp do
  ! do i = 1, m
  !   do k = 1, l
  !     div(i,0,k)  = div(i,n,k)
  !     div(i,n+1,k)= div(i,1,k)
  !   end do
  ! end do
  ! !$omp end do
  ! ---


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
        u(i,j,k)=u(i,j,k) + dt*(xnue + xlambda/density)*(div(i+1,j,k)-div(i-1,j,k))/dx*0.5
    
        !-- additional terms by porosity profile
        u(i,j,k) = u(i,j,k) &
                   + dt*( &
                          ((u_old(i+1,j,k)-u_old(i-1,j,k))/dx*0.5+(u_old(i+1,j,k)-u_old(i-1,j,k))/dx*0.5)  &
                            *xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*0.5                             &
                          +((u_old(i,j+1,k)-u_old(i,j-1,k))/dy*0.5+(v_old(i+1,j,k)-v_old(i-1,j,k))/dx*0.5) &
                            *xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*0.5                             &
                          +((u_old(i,j,k+1)-u_old(i,j,k-1))/dz*0.5+(w_old(i+1,j,k)-w_old(i-1,j,k))/dx*0.5) &
                            *xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*0.5                             &
                          + div(i,j,k)*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*0.5*xlambda/density                 &
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
        v(i,j,k)=v(i,j,k) + dt*(xnue + xlambda/density)*(div(i,j+1,k)-div(i,j-1,k))/dy*0.5

        !-- additional terms by porosity profile
        v(i,j,k) = v(i,j,k) &
                   + dt*( &
                         ((v_old(i+1,j,k)-v_old(i-1,j,k))/dx*0.5+(u_old(i,j+1,k)-u_old(i,j-1,k))/dy*0.5) &
                           *xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*0.5                            &
                         +((v_old(i,j+1,k)-v_old(i,j-1,k))/dy*.5+(v_old(i,j+1,k)-v_old(i,j-1,k))/dy*0.5) &
                           *xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*0.5                            &
                         +((v_old(i,j,k+1)-v_old(i,j,k-1))/dz*.5+(w_old(i,j+1,k)-w_old(i,j-1,k))/dy*0.5) &
                           *xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*0.5                            &
                         + div(i,j,k)*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*0.5*xlambda/density                &
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
        w(i,j,k)=w(i,j,k) +dt*(xnue + xlambda/density)*(div(i,j,k+1)-div(i,j,k-1))/dz*0.5

        !-- additional terms by porosity profile
        w(i,j,k) = w(i,j,k) &
                   + dt*( &
                          ((w_old(i+1,j,k)-w_old(i-1,j,k))/dx*0.5+(u_old(i,j,k+1)-u_old(i,j,k-1))/dz*0.5)  &
                            *xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*0.5                             &
                          +((w_old(i,j+1,k)-w_old(i,j-1,k))/dy*0.5+(v_old(i,j,k+1)-v_old(i,j,k-1))/dz*0.5) &
                            *xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*0.5                             &
                          +((w_old(i,j,k+1)-w_old(i,j,k-1))/dz*0.5+(w_old(i,j,k+1)-w_old(i,j,k-1))/dz*0.5) &
                            *xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*0.5                             &
                          + div(i,j,k)*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*0.5*xlambda/density                 &
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

  call boundary_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, l, height, yp, porosity)

  ! call solve_matrix_vec_omp (p, ap, ae, aw, an, as, at, ab, bb, &
  ! m, n, l, iter_max, relux_factor)
  call solve_matrix_vec_oacc (p, ap, ae, aw, an, as, at, ab, bb, &
  m, n, l, iter_max, relux_factor)
  ! ----------------
  ! ----------------
  return
end subroutine solve_p
!******************


!******************
! OpenACC Parallelized
! Written only for GPU machine
! No efficiency ensured on CPU machine 
subroutine  solve_matrix_vec_oacc (p, ap, ae, aw, an, as, at, ab, bb, &
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

  ! ----------------
  !   SOR algorithm
  ! ----------------

  !$acc data copy(p) &
  !$acc & copyin(ap, ae, aw, an, as, at, ab, bb) &
  !$acc & create(p_old)

  error = 0.0

  do iter = 1, iter_max
  ! write(*,*)'CHECK iteration no.'


  ! ! default periodic condition in yz-direction
  ! !$acc kernels
  ! !$acc loop independent
  ! do i = 1, m
  ! !$acc loop independent
  ! do k = 1, l
  !   p(i,0,k)  =p(i,n,k)
  !   p(i,n+1,k)=p(i,1,k)
  ! end do
  ! end do
  ! !$acc end kernels

  ! ! default periodic condition in xy-direction
  ! !$acc kernels
  ! !$acc loop independent
  ! do i = 1, m
  ! !$acc loop independent
  ! do j = 1, n
  !   p(i,j,0)  =p(i,j,l)
  !   p(i,j,l+1)=p(i,j,1)
  ! end do
  ! end do
  ! !$acc end kernels

  !$acc kernels
  !$acc loop independent
  do i = 0, m+1
  !$acc loop independent
  do j = 0, n+1
  !$acc loop independent
  do k = 0, l+1
    p_old(i,j,k) = p(i,j,k)
  end do
  end do
  end do
  !$acc end kernels

  !-- EVEN SPACE process
  !$acc kernels
  !$acc loop independent
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
  !$acc end kernels

  ! ! default periodic condition in yz-direction
  ! !$acc kernels
  ! !$acc loop independent
  ! do i = 1, m
  ! !$acc loop independent
  ! do k = 1, l
  !   p(i,0,k)  =p(i,n,k)
  !   p(i,n+1,k)=p(i,1,k)
  ! end do
  ! end do
  ! !$acc end kernels

  ! ! default periodic condition in xy-direction
  ! !$acc kernels
  ! !$acc loop independent
  ! do i = 1, m
  ! !$acc loop independent
  ! do j = 1, n
  !   p(i,j,0)  =p(i,j,l)
  !   p(i,j,l+1)=p(i,j,1)
  ! end do
  ! end do
  ! !$acc end kernels

  !$acc kernels
  !$acc loop independent
  do i = 0, m+1
  !$acc loop independent
  do j = 0, n+1
  !$acc loop independent
    do k = 0, l+1
      p_old(i,j,k) = p(i,j,k)
    end do
  end do
  end do
  !$acc end kernels


  !-- ODD SPACE process
  !$acc kernels
  !$acc loop independent
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
  !$acc end kernels

  !$acc kernels
  !$acc loop independent reduction(max:error)
  do k = 1, l
  !$acc loop independent reduction(max:error)
    do j = 1, n
  !$acc loop independent reduction(max:error)
      do i = 1, m
        error = max(error, abs(p_old(i,j,k)-p(i,j,k)))
      end do 
    end do
  end do
  !$acc end kernels

  end do

  ! ! default periodic condition in yz-direction
  ! !$acc kernels
  ! !$acc loop independent
  ! do i = 1, m
  ! !$acc loop independent
  !   do k = 1, l
  !     p(i,0,k)  =p(i,n,k)
  !     p(i,n+1,k)=p(i,1,k)
  !   end do
  ! end do
  ! !$acc end kernels

  ! ! default periodic condition in xy-direction
  ! !$acc kernels
  ! !$acc loop independent
  ! do i = 1, m
  ! !$acc loop independent
  !   do j = 1, n
  !     p(i,j,0)  =p(i,j,l)
  !     p(i,j,l+1)=p(i,j,1)
  !   end do
  ! end do
  ! !$acc end kernels

  write(*,*)'SOR iteration no.', iter-1,'  -- error=', error
  
  !$acc end data 

  return
end subroutine solve_matrix_vec_oacc
!******************


!******************
subroutine  boundary_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, l, height, yp, porosity)
  use global_3d
  use wall_conditions
  implicit none
  real,intent(in)::height
  real,intent(in),dimension(0:md,0:nd,0:ld)::p
  real,intent(inout),dimension(0:md,0:nd,0:ld)::ap, ae, aw, an, as, at, ab, bb
  real,intent(in),dimension(0:nd)::yp
  real,intent(in),dimension(0:md,0:nd,0:ld)::porosity
  integer,intent(in)::m, n, l
  
  ! local variables
  integer::i, j, k
  
  !$omp parallel private(i, j, k) &
  !$omp & shared(m, n, l, height) &
  !$omp & shared(p, ap, ae, aw, an, as, at, ab, bb, yp, porosity) &
  ! !$omp & shared(top_wall, bottom_wall, east_wall, west_wall, north_wall, south_wall) &
  !$omp & default(none)

  !--- top
  !$omp do
  do i = 0, m+1
    do j = 0, n+1
      if (top_wall == 0 .or. top_wall == 1) then
        ! --- wall (dp/dz=0. at k=l)
        ab(i,j,l) =ab(i,j,l)+at(i,j,l)
        at(i,j,l) =0.
        ! ---
      else if (top_wall == 2) then
        if (porosity(i,j,l) < 0.9) then
          ! --- wall (dp/dz=0. at k=1)
          ab(i,j,l) =ab(i,j,l)+at(i,j,l)
          at(i,j,l) =0.
          ! ---          
        else
          ! --- outlet (p=outlet_pressure at i=m)
          bb(i,j,l) = bb(i,j,1) + at(i,j,l) * p(i,j,l+1)
          ae(i,j,l) = 0.
          aw(i,j,l) = 0.
          an(i,j,l) = 0.
          as(i,j,l) = 0.
          at(i,j,l) = 0.
          ab(i,j,l) = 0.
          ! ---
        end if
      end if
    end do
  end do
  !$omp end do

  ! bottom
  !$omp do
  do i = 0, m+1
    do j = 0, n+1
      if (bottom_wall == 0 .or. bottom_wall == 1) then
        ! --- wall (dp/dz=0. at k=1)
        at(i,j,1) = at(i,j,1)+ab(i,j,1)
        ab(i,j,1) = 0.
        ! ---
      else if (bottom_wall == 2) then
        if (porosity(i,j,1) < 0.9) then
          ! --- wall (dp/dz=0. at k=1)
          at(i,j,1) = at(i,j,1)+ab(i,j,1)
          ab(i,j,1) = 0.
          ! ---          
        else
          ! --- outlet (p=outlet_pressure at i=m)
          bb(i,j,1) = bb(i,j,1) + ab(i,j,1) * p(i,j,0)
          ae(i,j,1) = 0.
          aw(i,j,1) = 0.
          an(i,j,1) = 0.
          as(i,j,1) = 0.
          at(i,j,1) = 0.
          ab(i,j,1) = 0.
          ! ---
        end if
      end if
    end do
  end do
  !$omp end do

  ! east
  !$omp do
  do j = 0, n+1
    do k = 0, l+1
      if (east_wall == 0 .or. east_wall == 1) then
      ! --- wall (dp/dx=0. at i=m)
      aw(m,j,k) = aw(m,j,k)+ae(m,j,k)
      ae(m,j,k) = 0.
      ! ---
      else if (east_wall == 2) then
        if (porosity(m,j,k) < 0.9) then
          ! --- wall (dp/dx=0. at i=m)
          aw(m,j,k) = aw(m,j,k)+ae(m,j,k)
          ae(m,j,k) = 0.
          ! ---          
        else
          ! --- outlet (p=outlet_pressure at i=m)
          bb(m,j,k) = bb(m,j,k) + ae(m,j,k) * p(m+1,j,k)
          ae(m,j,k) = 0.
          aw(m,j,k) = 0.
          an(m,j,k) = 0.
          as(m,j,k) = 0.
          at(m,j,k) = 0.
          ab(m,j,k) = 0.
          ! ---
        end if
      end if
    end do
  end do
  !$omp end do

  ! west 
  !$omp do
  do j = 0, n+1
    do k = 0, l+1
      if (west_wall == 0 .or. west_wall == 1) then
        ! --- wall (dp/dx=0. at i=1)
        ae(1,j,k) = ae(1,j,k)+aw(1,j,k)
        aw(1,j,k) = 0.
        ! ---
      else if (west_wall == 2) then
        if (porosity(1,j,k) < 0.9) then
          ! --- wall (dp/dx=0. at i=1)
          ae(1,j,k) = ae(1,j,k)+aw(1,j,k)
          aw(1,j,k) = 0.
          ! ---
        else
          ! --- outlet (p=outlet_pressure at i=1)
          bb(1,j,k) = bb(1,j,k) + aw(1,j,k) * p(0,j,k)
          ae(1,j,k) = 0.
          aw(1,j,k) = 0.
          an(1,j,k) = 0.
          as(1,j,k) = 0.
          at(1,j,k) = 0.
          ab(1,j,k) = 0.
          ! ---
        end if
      end if
    end do
  end do
  !$omp end do

  ! north 
  !$omp do
  do i = 0, m+1
    do k = 0, l+1
      if (north_wall == 0 .or. north_wall == 1) then
        ! --- wall (dp/dy=0. at j=n)
        as(i,n,k) = as(i,n,k)+an(i,n,k)
        an(i,n,k) = 0.
        ! ---
      else if (north_wall == 2) then
        if (porosity(i,n,k) < 0.9) then
          ! --- wall (dp/dy=0. at j=n)
          as(i,n,k) = as(i,n,k)+an(i,n,k)
          an(i,n,k) = 0.
          ! ---
        else
          ! --- outlet (p=outlet_pressure at j=n)
          bb(i,n,k) = bb(i,n,k) + an(i,n,k) * p(i,n+1,k)
          ae(i,n,k) = 0.
          aw(i,n,k) = 0.
          an(i,n,k) = 0.
          as(i,n,k) = 0.
          at(i,n,k) = 0.
          ab(i,n,k) = 0.
        end if
      end if
    end do
  end do
  !$omp end do

  ! south
  !$omp do
  do i = 0, m+1
    do k = 0, l+1
      if (south_wall == 0 .or. south_wall == 1) then
        an(i,1,k) = an(i,1,k)+as(i,1,k)
        as(i,1,k) = 0.        
      else if (south_wall == 2) then
        ! --- wall (dp/dy=0. at j=1)
        if (porosity(i,1,k) < 0.9) then
          an(i,1,k) = an(i,1,k)+as(i,1,k)
          as(i,1,k) = 0.
        else
        ! --- outlet (p=outlet_pressure at j=1)          
          bb(i,1,k) = bb(i,1,k) + as(i,1,k) * p(i,0,k)
          ae(i,1,k) = 0.
          aw(i,1,k) = 0.
          an(i,1,k) = 0.
          as(i,1,k) = 0.
          at(i,1,k) = 0.
          ab(i,1,k) = 0.
        end if
      end if
    end do
  end do
  !$omp end do
  !$omp end parallel
  return
end subroutine  boundary_matrix
!******************

!  conditions  

!******************
subroutine  boundary(p, u, v, w, xp, yp, zp, width, height, depth    &
                      , inlet_velocity, outlet_pressure, AoA, porosity, m, n, l)
  use global_3d
  use wall_conditions
  implicit none
  real,intent(in)::width, height, depth, inlet_velocity, outlet_pressure, AoA
  real,intent(inout),dimension(0:md,0:nd,0:ld)::u, v, w, p
  real,intent(in),dimension(0:md,0:nd,0:ld)::porosity
  real,intent(in),dimension(0:md)::xp
  real,intent(in),dimension(0:nd)::yp
  real,intent(in),dimension(0:ld)::zp
  integer,intent(in)::m, n, l
  
  ! local variables
  real, parameter::small=1.e-6, big=1.e6, zero=0., pai=atan(1.)*4.
  integer::i, j, k

  ! top
  do i = 0, m+1
    do j = 0, n+1
      if (top_wall == 0) then
        ! --- wall (u=0., v=0., w=0., dp/dz=0. at k=l, @poro<0.9 (solid))
        u(i,j,l) = 0.
        v(i,j,l) = 0.
        w(i,j,l) = 0.
        w(i,j,l+1) =-w(i,j,l-1) ! dummy
        p(i,j,l+1) = p(i,j,l-1) ! dummy
      else if (top_wall == 1) then
        if (porosity(i,j,l) >= 0.9) then
        ! --- inlet (w=inlet_velocity, u,v=0., dp/dz=0 at i=l, @poro>=0.9 (fluid))
          u(i,j,l) =0.
          v(i,j,l) =0.
          w(i,j,l) =-inlet_velocity
          u(i,j,l+1) =u(i,j,l)  ! dummy
          v(i,j,l+1) =v(i,j,l)  ! dummy
          w(i,j,l+1) =w(i,j,l)  ! dummy
          p(i,j,l+1) =p(i,j,l-1)  ! dummy
        else
        ! --- wall (u=0., v=0., w=0., dp/dz=0. at k=l, @poro<0.9 (solid))
          u(i,j,l) = 0.
          v(i,j,l) = 0.
          w(i,j,l) = 0.
          w(i,j,l+1) =-w(i,j,l-1) ! dummy
          p(i,j,l+1) = p(i,j,l-1) ! dummy
        end if
      else if (top_wall == 2) then
        if (porosity(i,j,l) >= 0.9) then
        ! --- outlet (du/dx=0., dv/dx=0., dw/dz=0., p=outlet_pressure at k=l)
          u(i,j,l+1) =u(i,j,l-1)  ! dummy
          v(i,j,l+1) =v(i,j,l-1)  ! dummy
          w(i,j,l+1) =w(i,j,l-1)  ! dummy
          p(i,j,l+1) =outlet_pressure ! dummy
        else
        ! --- wall (u=0., v=0., w=0., dp/dz=0. at k=l, @poro<0.9 (solid))
          u(i,j,l) = 0.
          v(i,j,l) = 0.
          w(i,j,l) = 0.
          w(i,j,l+1) =-w(i,j,l-1) ! dummy
          p(i,j,l+1) = p(i,j,l-1) ! dummy
        end if
      end if
    end do
  end do

  ! bottom 
  do i = 0, m+1
    do j = 0, n+1
      if (bottom_wall == 0) then
        ! --- wall (u=0., v=0., w=0., dp/dz=0. at k=1, @poro<0.9 (solid))
        u(i,j,1) = 0.
        v(i,j,1) = 0.
        w(i,j,1) = 0.
        w(i,j,0) =-w(i,j,2) ! dummy
        p(i,j,0) = p(i,j,2) ! dummy
      else if (bottom_wall == 1) then
        if (porosity(i,j,l) >= 0.9) then
        ! --- inlet (w=inlet_velocity, u,v=0., dp/dz=0 at i=1, @poro>=0.9 (fluid))
          u(i,j,1) =0.
          v(i,j,1) =0.
          w(i,j,1) =inlet_velocity
          u(i,j,0) =u(i,j,1)  ! dummy
          v(i,j,0) =v(i,j,1)  ! dummy
          w(i,j,0) =w(i,j,1)  ! dummy
          p(i,j,0) =p(i,j,2)  ! dummy
        else
        ! --- wall (u=0., v=0., w=0., dp/dz=0. at k=1, @poro<0.9 (solid))
          u(i,j,1) = 0.
          v(i,j,1) = 0.
          w(i,j,1) = 0.
          w(i,j,0) =-w(i,j,2) ! dummy
          p(i,j,0) = p(i,j,2) ! dummy
        end if
      else if (bottom_wall == 2) then
        if (porosity(i,j,1) >= 0.9) then
        ! --- outlet (du/dx=0., dv/dx=0., dw/dz=0., p=outlet_pressure at k=1)
          u(i,j,0) =u(i,j,1)  ! dummy
          v(i,j,0) =v(i,j,1)  ! dummy
          w(i,j,0) =w(i,j,1)  ! dummy
          p(i,j,0) =outlet_pressure ! dummy
        else
        ! --- wall (u=0., v=0., w=0., dp/dz=0. at k=1, @poro<0.9 (solid))
          u(i,j,1) = 0.
          v(i,j,1) = 0.
          w(i,j,1) = 0.
          w(i,j,0) =-w(i,j,2) ! dummy
          p(i,j,0) = p(i,j,2) ! dummy 
        end if
      end if
    end do
  end do

  ! west 
  do j = 0, n+1
    do k = 0, l+1
      ! --- wall (u=0., v=0., w=0., dp/dx=0. at i=1, @poro<0.9 (solid))
      if (west_wall == 0) then
        u(1,j,k) = 0.
        v(1,j,k) = 0.
        w(1,j,k) = 0.
        u(0,j,k) =-u(2,j,k) ! dummy
        p(0,j,k) = p(2,j,k) ! dummy
      else if (west_wall == 1) then
        ! --- inlet (u=inlet_velocity, v,w=0., dp/dz=0 at i=1, @poro>=0.9 (fluid))
        if (porosity(1,j,k) >= 0.9) then
          u(1,j,k) =inlet_velocity
          v(1,j,k) =0.
          w(1,j,k) =0.
          u(0,j,k) =u(1,j,k)  ! dummy
          v(0,j,k) =v(1,j,k)  ! dummy
          w(0,j,k) =w(1,j,k)  ! dummy
          p(0,j,k) =p(2,j,k)  ! dummy
        else
        ! --- wall (u=0., v=0., w=0., dp/dx=0. at i=1, @poro<0.9 (solid))
          u(1,j,k) = 0.
          v(1,j,k) = 0.
          w(1,j,k) = 0.
          u(0,j,k) =-u(2,j,k) ! dummy
          p(0,j,k) = p(2,j,k) ! dummy          
        end if
      else if (west_wall == 2) then
        if (porosity(1,j,k) >= 0.9) then
        ! --- outlet (du/dx=0., dv/dx=0., dw/dz=0., p=outlet_pressure at i=1)
          u(0,j,k) =u(1,j,k)  ! dummy
          v(0,j,k) =v(1,j,k)  ! dummy
          w(0,j,k) =w(1,j,k)  ! dummy
          p(0,j,k) =outlet_pressure ! dummy
        else
        ! --- wall (u=0., v=0., w=0., dp/dx=0. at i=1, @poro<0.9 (solid))
          u(1,j,k) = 0.
          v(1,j,k) = 0.
          w(1,j,k) = 0.
          u(0,j,k) =-u(2,j,k) ! dummy
          p(0,j,k) = p(2,j,k) ! dummy            
        end if
      end if
    end do
  end do

  ! east
  do j = 0, n+1
    do k = 0, l+1
      if (east_wall == 0) then
        ! --- wall (u=0., v=0., w=0., dp/dx=0. at i=m, @poro<0.9 (solid))
          u(m,j,k) = 0.
          v(m,j,k) = 0.
          w(m,j,k) = 0.
          u(m+1,j,k) =-u(m-1,j,k) ! dummy
          p(m+1,j,k) = p(m-1,j,k) ! dummy
        ! ---
      else if (east_wall == 1) then
        if (porosity(m,j,k) >= 0.9) then
        ! --- inlet (u=inlet_velocity, v,w=0., dp/dz=0 at i=m, @poro>=0.9 (fluid))
          u(m,j,k) =-inlet_velocity
          v(m,j,k) =0.
          w(m,j,k) =0.
          u(m+1,j,k) =u(m,j,k)  ! dummy
          v(m+1,j,k) =v(m,j,k)  ! dummy
          w(m+1,j,k) =w(m,j,k)  ! dummy
          p(m+1,j,k) =p(m-1,j,k)  ! dummy
        else
        ! --- wall (u=0., v=0., w=0., dp/dx=0. at i=m, @poro<0.9 (solid))
          u(m,j,k) = 0.
          v(m,j,k) = 0.
          w(m,j,k) = 0.
          u(m+1,j,k) =-u(m-1,j,k) ! dummy
          p(m+1,j,k) = p(m-1,j,k) ! dummy
        end if
      ! ---
      else if (east_wall == 2) then
        if (porosity(m,j,k) >= 0.9) then
        ! --- outlet (du/dx=0., dv/dx=0., dw/dz=0., p=outlet_pressure at i=1)
          u(m+1,j,k) =u(m,j,k)  ! dummy
          v(m+1,j,k) =v(m,j,k)  ! dummy
          w(m+1,j,k) =w(m,j,k)  ! dummy
          p(m+1,j,k) =outlet_pressure ! dummy
        else
        ! --- wall (u=0., v=0., w=0., dp/dx=0. at i=m, @poro<0.9 (solid))
          u(m,j,k) = 0.
          v(m,j,k) = 0.
          w(m,j,k) = 0.
          u(m+1,j,k) =-u(m-1,j,k) ! dummy
          p(m+1,j,k) = p(m-1,j,k) ! dummy
        end if
      end if
    end do
  end do

  ! north
  do i = 0, m+1
    do k = 0, l+1
      if (north_wall == 0) then
        u(i,n,k) = 0.
        v(i,n,k) = 0.
        w(i,n,k) = 0.
        u(i,n+1,k) =-u(i,n-1,k) ! dummy
        p(i,n+1,k) = p(i,n-1,k) ! dummy
      else if (north_wall == 1) then
        if (porosity(i,n,k) >= 0.9) then
          u(i,n,k) =-inlet_velocity
          v(i,n,k) =0.
          w(i,n,k) =0.
          u(i,n+1,k) =u(i,n,k)  ! dummy
          v(i,n+1,k) =v(i,n,k)  ! dummy
          w(i,n+1,k) =w(i,n,k)  ! dummy
          p(i,n+1,k) =p(i,n-1,k)  ! dummy          
        else
          u(i,n,k) = 0.
          v(i,n,k) = 0.
          w(i,n,k) = 0.
          u(i,n+1,k) =-u(i,n-1,k) ! dummy
          p(i,n+1,k) = p(i,n-1,k) ! dummy
        end if        
      else if (north_wall == 2) then
        if (porosity(i,n,k) >= 0.9) then
          u(i,n+1,k) =u(i,n,k)  ! dummy
          v(i,n+1,k) =v(i,n,k)  ! dummy
          w(i,n+1,k) =w(i,n,k)  ! dummy
          p(i,n+1,k) =outlet_pressure ! dummy         
        else
          u(i,n,k) = 0.
          v(i,n,k) = 0.
          w(i,n,k) = 0.
          u(i,n+1,k) =-u(i,n-1,k) ! dummy
          p(i,n+1,k) = p(i,n-1,k) ! dummy
        end if                
      end if
    end do
  end do

  ! south 
  do i = 0, m+1
    do k = 0, l+1
      if (south_wall == 0) then
      ! --- wall (u=0., v=0., w=0., dp/dy=0. at j=1, @poro<0.9 (solid))
        u(i,1,k) = 0.
        v(i,1,k) = 0.
        w(i,1,k) = 0.
        v(i,0,k) =-v(i,2,k) ! dummy
        p(i,0,k) = p(i,2,k) ! dummy
      else if (south_wall == 1) then
        if (porosity(i,1,k) >= 0.9) then
        ! --- inlet (u=inlet_velocity, v,w=0., dp/dz=0 at j=1, @poro>=0.9 (fluid))
          u(i,1,k) =inlet_velocity
          v(i,1,k) =0.
          w(i,1,k) =0.
          u(i,0,k) =u(i,1,k)  ! dummy
          v(i,0,k) =v(i,1,k)  ! dummy
          w(i,0,k) =w(i,1,k)  ! dummy
          p(i,0,k) =p(i,2,k)  ! dummy          
        else
        ! --- wall (u=0., v=0., w=0., dp/dy=0. at j=1, @poro<0.9 (solid))
          u(i,1,k) = 0.
          v(i,1,k) = 0.
          w(i,1,k) = 0.
          v(i,0,k) =-v(i,2,k) ! dummy
          p(i,0,k) = p(i,2,k) ! dummy
        end if
      else if (south_wall == 2) then
        ! --- south outlet (du/dx=0., dv/dx=0., dw/dz=0., p=outlet_pressure at j=1)
          if (porosity(i,1,k) >= 0.9) then
            u(i,0,k) =u(i,2,k)  ! dummy
            v(i,0,k) =v(i,2,k)  ! dummy
            w(i,0,k) =w(i,2,k)  ! dummy
            p(i,0,k) =outlet_pressure   ! dummy
          else
          ! --- wall (u=0., v=0., w=0., dp/dy=0. at j=1, @poro<0.9 (solid))
            u(i,1,k) = 0.
            v(i,1,k) = 0.
            w(i,1,k) = 0.
            v(i,0,k) =-v(i,2,k) ! dummy
            p(i,0,k) = p(i,2,k) ! dummy
          end if
      end if
    end do
  end do
  
  return
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
  u(i,j,k) = 0.
  v(i,j,k) = 0.
  w(i,j,k) = 0.
  p(i,j,k) = outlet_pressure
  end do
  end do
  end do
  !$omp end do

  !$omp end parallel
  return
end subroutine initial_conditions
!******************