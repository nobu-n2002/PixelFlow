! ibm_3d_acc_gpu.f90
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
  real::t_0,t_1,t_2,t_3
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
  call cpu_time(t_0)
  
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

  !$acc enter data copyin (p,u,v,w,xp,yp,zp,porosity) create(u_old,v_old,w_old)
  
  call boundary (p, u, v, w, xp, yp, zp, width, height, depth &
                         , inlet_velocity, outlet_pressure, AoA, porosity, m, n, l)
  
  call  output_paraview_temp_3d (p, u, v, w, porosity, xp, yp, zp, m, n, l, 0)  
  ! print initial conditions
  ! call  output_solution_3d (p, u, v, w, m, n, l)
  
  ! ----------------
  ! MAC algorithm start
  call get_now_time()
  call cpu_time(t_1)
  write(*,*) '# --- MAC algorithm start'

  do istep = 1, istep_max
    
    time = istep * dt

    !$acc parallel loop present(u,v,w,u_old,v_old,w_old)
    do k = 0, l+1
      do j = 0, n+1
        do i = 0, m+1
          u_old(i,j,k) = u(i,j,k)
          v_old(i,j,k) = v(i,j,k)
          w_old(i,j,k) = w(i,j,k)
        end do
      end do
    end do
    !$acc end parallel
  
    write(*,*)'--- time_steps= ',istep, ' --  time = ',time
  
    call solve_p (p, u, v, w, u_old, v_old, w_old, porosity, &
    xnue, xlambda, density, height, thickness, &
    yp, dx, dy, dz, dt, m, n, l, &
    nonslip, iter_max, relux_factor)
  
    !-- solve u, v, w (fix u, v, w)
    !$acc parallel loop present(p,u,v,w)
    do k = 1, l
      do j = 1, n
        do i = 1, m
          u(i,j,k) = u(i,j,k) - dt/density*(p(i+1,j,k)-p(i-1,j,k))/dx*0.5
          v(i,j,k) = v(i,j,k) - dt/density*(p(i,j+1,k)-p(i,j-1,k))/dy*0.5
          w(i,j,k) = w(i,j,k) - dt/density*(p(i,j,k+1)-p(i,j,k-1))/dz*0.5
        end do
      end do
    end do
    !$acc end parallel
  
    call boundary(p, u, v, w, xp, yp, zp, width, height, depth  &
                      , inlet_velocity, outlet_pressure, AoA, porosity, m, n, l)
    
    if(mod(istep,istep_out)==0) then
      !$acc update host(p,u,v,w,porosity)
      call  output_paraview_temp_3d (p, u, v, w, porosity, xp, yp, zp, m, n, l, istep)
    end if
    
  end do
  call get_now_time()
  ! MAC algorithm end
  ! ----------------
  
  ! print solutions 
  !$acc exit data copyout(p,u,v,w,porosity)
  ! call  output_solution_post_3d (p, u, v, w, xp, yp, zp, porosity, m, n, l)
  ! call  output_divergent_3d (p, u, v, w, porosity, dx, dy, dz, m, n, l)
  call  output_paraview_3d (p, u, v, w, porosity, xp, yp, zp, m, n, l)
  
  call get_now_time()
  call cpu_time(t_2)

  write(*,*) "################################"
  write(*,*) "Initialization time:", t_1-t_0, "s"
  write(*,*) "Total Elapsed time:", t_2-t_0, "s"
  write(*,*) "################################"
  write(*,*) 'program finished'
  
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
  real:: alpha
 
  real,dimension(0:md,0:nd,0:ld) :: ap, ae, aw, an, as, at, ab, bb, div
  integer :: i, j, k

  if (nonslip) then
    alpha = 32.0
  else
    alpha = 0.0
  end if
  
  !$acc enter data create(ap,ae,aw,an,as,at,ab,bb,div)
  !-----------------
  !  divergence term  div(u)
  !-----------------
  !$acc kernels present(u_old,v_old,w_old,div)
  !$acc loop
  do k = 1, l
    do j = 1, n
      do i = 1, m
        div(i,j,k)= (u_old(i+1,j,k)-u_old(i-1,j,k))/dx*0.5 &
                  + (v_old(i,j+1,k)-v_old(i,j-1,k))/dy*0.5 &
                  + (w_old(i,j,k+1)-w_old(i,j,k-1))/dz*0.5 
      end do
    end do
  end do
  !$acc end kernels
  
  !$acc kernels present(div)
  !$acc loop
  do k = 1, l
    do j = 1, n
      div(0,j,k)  = 0.  ! inlet yz
      div(m+1,j,k)= 0.  ! outlet yz
    end do
  end do
  !$acc end kernels
 
  !$acc kernels present(div)
  !$acc loop
  do k = 1, l
    do i = 1, m
      div(i,0,k)  = div(i,n,k)  ! periodic condition xz
      div(i,n+1,k)= div(i,1,k)
    end do
  end do
  !$acc end kernels
  
  !$acc kernels present(div)
  !$acc loop
  do j = 1, n
    do i = 1, m
      div(i,j,0)  = div(i,j,l)  ! periodic condition xy
      div(i,j,l+1)= div(i,j,1)
    end do
  end do
  !$acc end kernels

  ! ----------------
  !   velocity u
  ! ----------------
  !$acc kernels present(u,v,w,u_old,v_old,w_old,div,porosity)
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
        u(i,j,k)=u(i,j,k)-dt*xnue*u_old(i,j,k)/(thickness*dx)**2*alpha*porosity(i,j,k)*(1.-porosity(i,j,k))*(1.-porosity(i,j,k))

      end do
    end do
  end do
  !$acc end kernels

  ! ----------------
  !   velocity v
  ! ----------------
  !$acc kernels present(u,v,w,u_old,v_old,w_old,div,porosity)
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
        v(i,j,k)=v(i,j,k)-dt*xnue*v_old(i,j,k)/(thickness*dy)**2*alpha*porosity(i,j,k)*(1.-porosity(i,j,k))*(1.-porosity(i,j,k))

      end do
    end do
  end do
  !$acc end kernels
  

  ! ----------------
  !   velocity w
  ! ----------------
  !$acc kernels present(u,v,w,u_old,v_old,w_old,div,porosity)
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
        w(i,j,k)=w(i,j,k)- dt*xnue*w_old(i,j,k)/(thickness*dz)**2*alpha*porosity(i,j,k)*(1.-porosity(i,j,k))*(1.-porosity(i,j,k))
    
      end do
    end do
  end do
  !$acc end kernels
  

  ! ----------------
  ! matrix solution  !  formulation of porous media
  !$acc kernels present(u,v,w,ap,ae,aw,an,as,at,ab,bb,porosity)
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
  !$acc end kernels

  call boundrary_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, l, height, yp)
  
  call solve_matrix_3x_async_1x8 (p, ap, ae, aw, an, as, at, ab, bb, &
  m, n, l, iter_max, relux_factor)

  ! call solve_matrix_RB (p, ap, ae, aw, an, as, at, ab, bb, &
  ! m, n, l, iter_max, relux_factor)

  !$acc exit data delete (ap,ae,aw,an,as,at,ab,bb,div)
  
  return
end subroutine solve_p
!******************

!******************
! OpenACC Parallelized
! Written only for GPU machine
! No efficiency ensured on CPU machine 

subroutine  solve_matrix_3x_async_1x8 (p, ap, ae, aw, an, as, at, ab, bb, &
  m, n, l, iter_max, relux_factor)
  use global_3d
  implicit none
  real,intent(inout),dimension(0:md,0:nd,0:ld):: p
  real,intent(in),dimension(0:md,0:nd,0:ld):: ap, ae, aw, an, as, at, ab, bb
  integer,intent(in):: m, n, l
  integer,intent(in):: iter_max
  real,intent(in):: relux_factor

  ! local variables
  real:: error, correct
  real,dimension(0:(l+3)/2,0:(n+3)/2,0:(m+3)/2):: &
                f000,f001,f010,f011,f100,f101,f110,f111
  real,dimension(0:(l+3)/2,0:(n+3)/2,0:(m+3)/2):: &
                ap000,ap001,ap010,ap011,ap100,ap101,ap110,ap111, & 
                ae000,ae001,ae010,ae011,ae100,ae101,ae110,ae111, & 
                aw000,aw001,aw010,aw011,aw100,aw101,aw110,aw111, & 
                an000,an001,an010,an011,an100,an101,an110,an111, & 
                as000,as001,as010,as011,as100,as101,as110,as111, & 
                at000,at001,at010,at011,at100,at101,at110,at111, & 
                ab000,ab001,ab010,ab011,ab100,ab101,ab110,ab111, & 
                bb000,bb001,bb010,bb011,bb100,bb101,bb110,bb111
  integer::i, j, k, iter
  integer::i_o, j_o, k_o, i_e, j_e, k_e, ii, jj, kk, mm, nn, ll
  integer::i0,i1,j0,j1,k0,k1
  !check counter

  ! ----------------
  !   SOR algorithm
  ! ----------------

  !Nakamichi 01Sep2025
  !$acc data create(f000,f001,f010,f011,f100,f101,f110,f111,&
  !$acc &        ap000,ap001,ap010,ap011,ap100,ap101,ap110,ap111,&
  !$acc &        ae000,ae001,ae010,ae011,ae100,ae101,ae110,ae111,&
  !$acc &        aw000,aw001,aw010,aw011,aw100,aw101,aw110,aw111,&
  !$acc &        an000,an001,an010,an011,an100,an101,an110,an111,&
  !$acc &        as000,as001,as010,as011,as100,as101,as110,as111,&
  !$acc &        at000,at001,at010,at011,at100,at101,at110,at111,&
  !$acc &        ab000,ab001,ab010,ab011,ab100,ab101,ab110,ab111,&
  !$acc &        bb000,bb001,bb010,bb011,bb100,bb101,bb110,bb111)

  !convert local matrix  p to (f000,...,f111)
  i_e= mod(m,2)  ! 0 for m=even / 1 for m=odd
  j_e= mod(n,2)  ! 0 for n=even / 1 for n=odd
  k_e= mod(l,2)  ! 0 for l=even / 1 for l=odd
  i_o= 1-i_e     ! 0 for m=odd / 1 for m=even  
  j_o= 1-j_e     ! 0 for n=odd / 1 for n=even
  k_o= 1-k_e     ! 0 for l=odd / 1 for l=even
  mm =(m+i_e)/2 ! p(1,*,*)=f**1(*,*,1) for m=odd / p(1,*,*)=f**0(*,*,1) for m=even
  nn =(n+j_e)/2 ! p(*,1,*)=f*1*(*,1,*) for n=odd / p(*,1,*)=f*0*(*,1,*) for n=even    
  ll =(l+k_e)/2 ! p(*,*,1)=f1**(1,*,*) for l=odd / p(*,*,1)=f0**(1,*,*) for l=even    
                    ! (0,0,0) <= f***(kk,jj,ii) <= (ll+1,nn+1,mm+1)


  ! convert data  p => (f000,...,f111)
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111,p)
  !$acc loop gang, vector(128) collapse(3) !Nakamichi 01Sep2025
  do ii = 1, mm
  do jj = 1, nn
  do kk = 1, ll
  i0= (ii-1)*2 + i_o;  i1 = i0 + 1
  j0= (jj-1)*2 + j_o;  j1 = j0 + 1
  k0= (kk-1)*2 + k_o;  k1 = k0 + 1

  f000(kk,jj,ii)=p(i0,j0,k0)
  f001(kk,jj,ii)=p(i1,j0,k0)
  f010(kk,jj,ii)=p(i0,j1,k0)
  f011(kk,jj,ii)=p(i1,j1,k0)
  
  f100(kk,jj,ii)=p(i0,j0,k1)
  f101(kk,jj,ii)=p(i1,j0,k1)
  f110(kk,jj,ii)=p(i0,j1,k1)
  f111(kk,jj,ii)=p(i1,j1,k1)
  end do
  end do
  end do
  !$acc end kernels
  
  ! convert matrix components  ap,ae,aw,an,an,at,ab,bb
  !$acc kernels present(ap000,ap001,ap010,ap011,ap100,ap101,ap110,ap111,&
  !$acc &               ae000,ae001,ae010,ae011,ae100,ae101,ae110,ae111,&
  !$acc &               aw000,aw001,aw010,aw011,aw100,aw101,aw110,aw111,&
  !$acc &               an000,an001,an010,an011,an100,an101,an110,an111,&
  !$acc &               as000,as001,as010,as011,as100,as101,as110,as111,&
  !$acc &               at000,at001,at010,at011,at100,at101,at110,at111,&
  !$acc &               ab000,ab001,ab010,ab011,ab100,ab101,ab110,ab111,&
  !$acc &               bb000,bb001,bb010,bb011,bb100,bb101,bb110,bb111)
  !$acc loop gang, vector(128) collapse(3) !Nakamichi 01Sep2025
  do ii = 1, mm
  do jj = 1, nn
  do kk = 1, ll
  i0= (ii-1)*2 + i_o;  i1 = i0 + 1
  j0= (jj-1)*2 + j_o;  j1 = j0 + 1
  k0= (kk-1)*2 + k_o;  k1 = k0 + 1

  ap000(kk,jj,ii)=ap(i0,j0,k0)
  ap001(kk,jj,ii)=ap(i1,j0,k0)
  ap010(kk,jj,ii)=ap(i0,j1,k0)
  ap011(kk,jj,ii)=ap(i1,j1,k0)
  ap100(kk,jj,ii)=ap(i0,j0,k1)
  ap101(kk,jj,ii)=ap(i1,j0,k1)
  ap110(kk,jj,ii)=ap(i0,j1,k1)
  ap111(kk,jj,ii)=ap(i1,j1,k1)

  ae000(kk,jj,ii)=ae(i0,j0,k0)
  ae001(kk,jj,ii)=ae(i1,j0,k0)
  ae010(kk,jj,ii)=ae(i0,j1,k0)
  ae011(kk,jj,ii)=ae(i1,j1,k0)
  ae100(kk,jj,ii)=ae(i0,j0,k1)
  ae101(kk,jj,ii)=ae(i1,j0,k1)
  ae110(kk,jj,ii)=ae(i0,j1,k1)
  ae111(kk,jj,ii)=ae(i1,j1,k1)

  aw000(kk,jj,ii)=aw(i0,j0,k0)
  aw001(kk,jj,ii)=aw(i1,j0,k0)
  aw010(kk,jj,ii)=aw(i0,j1,k0)
  aw011(kk,jj,ii)=aw(i1,j1,k0)
  aw100(kk,jj,ii)=aw(i0,j0,k1)
  aw101(kk,jj,ii)=aw(i1,j0,k1)
  aw110(kk,jj,ii)=aw(i0,j1,k1)
  aw111(kk,jj,ii)=aw(i1,j1,k1)

  an000(kk,jj,ii)=an(i0,j0,k0)
  an001(kk,jj,ii)=an(i1,j0,k0)
  an010(kk,jj,ii)=an(i0,j1,k0)
  an011(kk,jj,ii)=an(i1,j1,k0)
  an100(kk,jj,ii)=an(i0,j0,k1)
  an101(kk,jj,ii)=an(i1,j0,k1)
  an110(kk,jj,ii)=an(i0,j1,k1)
  an111(kk,jj,ii)=an(i1,j1,k1)

  as000(kk,jj,ii)=as(i0,j0,k0)
  as001(kk,jj,ii)=as(i1,j0,k0)
  as010(kk,jj,ii)=as(i0,j1,k0)
  as011(kk,jj,ii)=as(i1,j1,k0)
  as100(kk,jj,ii)=as(i0,j0,k1)
  as101(kk,jj,ii)=as(i1,j0,k1)
  as110(kk,jj,ii)=as(i0,j1,k1)
  as111(kk,jj,ii)=as(i1,j1,k1)

  at000(kk,jj,ii)=at(i0,j0,k0)
  at001(kk,jj,ii)=at(i1,j0,k0)
  at010(kk,jj,ii)=at(i0,j1,k0)
  at011(kk,jj,ii)=at(i1,j1,k0)
  at100(kk,jj,ii)=at(i0,j0,k1)
  at101(kk,jj,ii)=at(i1,j0,k1)
  at110(kk,jj,ii)=at(i0,j1,k1)
  at111(kk,jj,ii)=at(i1,j1,k1)

  ab000(kk,jj,ii)=ab(i0,j0,k0)
  ab001(kk,jj,ii)=ab(i1,j0,k0)
  ab010(kk,jj,ii)=ab(i0,j1,k0)
  ab011(kk,jj,ii)=ab(i1,j1,k0)
  ab100(kk,jj,ii)=ab(i0,j0,k1)
  ab101(kk,jj,ii)=ab(i1,j0,k1)
  ab110(kk,jj,ii)=ab(i0,j1,k1)
  ab111(kk,jj,ii)=ab(i1,j1,k1)

  bb000(kk,jj,ii)=bb(i0,j0,k0)
  bb001(kk,jj,ii)=bb(i1,j0,k0)
  bb010(kk,jj,ii)=bb(i0,j1,k0)
  bb011(kk,jj,ii)=bb(i1,j1,k0)
  bb100(kk,jj,ii)=bb(i0,j0,k1)
  bb101(kk,jj,ii)=bb(i1,j0,k1)
  bb110(kk,jj,ii)=bb(i0,j1,k1)
  bb111(kk,jj,ii)=bb(i1,j1,k1)
  end do
  end do
  end do
  !$acc end kernels

  do iter = 1, iter_max   ! -- SOR iteration start
 
  ! default periodic condition in y-direction
  if(j_e.eq.0)then  ! n=even j_e=0 & j_o=1
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111)
  !$acc loop gang, vector(128) collapse(2)
  do ii = 1, mm
  do kk = 1, ll
  f010(kk,   0,ii)=f010(kk,nn,ii)
  f011(kk,   0,ii)=f011(kk,nn,ii)
  f000(kk,nn+1,ii)=f000(kk,1 ,ii)
  f001(kk,nn+1,ii)=f001(kk,1 ,ii)
  f110(kk,   0,ii)=f110(kk,nn,ii)
  f111(kk,   0,ii)=f111(kk,nn,ii)
  f100(kk,nn+1,ii)=f100(kk,1 ,ii)
  f101(kk,nn+1,ii)=f101(kk,1 ,ii)
  end do
  end do
  !$acc end kernels
  else  ! n=odd j_e=1 & j_o=0  
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111)
  !$acc loop gang, vector(128) collapse(2)
  do ii = 1, mm
  do kk = 1, ll
  f000(kk,   1,ii)=f010(kk,nn,ii)
  f001(kk,   1,ii)=f011(kk,nn,ii)
  f000(kk,nn+1,ii)=f010(kk,1 ,ii)
  f001(kk,nn+1,ii)=f011(kk,1 ,ii)
  f100(kk,   1,ii)=f110(kk,nn,ii)
  f101(kk,   1,ii)=f111(kk,nn,ii)
  f100(kk,nn+1,ii)=f110(kk,1 ,ii)
  f101(kk,nn+1,ii)=f111(kk,1 ,ii)
  end do
  end do
  !$acc end kernels
  end if

  ! default periodic condition in z-direction
  if(k_e.eq.0)then  ! l=even k_e=0 & k_o=1
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111)
  !$acc loop gang, vector(128) collapse(2)
  do ii = 1, mm
  do jj = 1, nn
  f010(ll+1,jj,ii)=f010(1,jj,ii)
  f011(ll+1,jj,ii)=f011(1,jj,ii)
  f000(ll+1,jj,ii)=f000(1,jj,ii)
  f001(ll+1,jj,ii)=f001(1,jj,ii)
  f110(0,jj,ii)=f110(ll,jj,ii)
  f111(0,jj,ii)=f111(ll,jj,ii)
  f100(0,jj,ii)=f100(ll,jj,ii)
  f101(0,jj,ii)=f101(ll,jj,ii)
  end do
  end do
  !$acc end kernels
  else  ! l=odd k_e=1 & k_o=0  
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111)
  !$acc loop gang, vector(128) collapse(2)
  do ii = 1, mm
  do jj = 1, nn
  f010(ll+1,jj,ii)=f110(1,jj,ii)
  f011(ll+1,jj,ii)=f111(1,jj,ii)
  f000(ll+1,jj,ii)=f100(1,jj,ii)
  f001(ll+1,jj,ii)=f101(1,jj,ii)
  f010(1,jj,ii)=f110(ll,jj,ii)
  f011(1,jj,ii)=f111(ll,jj,ii)
  f000(1,jj,ii)=f100(ll,jj,ii)
  f001(1,jj,ii)=f101(ll,jj,ii)
  end do
  end do
  !$acc end kernels
  end if

  ! default Dirichlet(i=m+1) / Neumann(i=0) condition in x-direction
  if(i_e.eq.0)then  ! m=even i_e=0 & i_o=1
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111)
  !$acc loop gang, vector(128) collapse(2)
  do jj = 1, nn
  do kk = 1, ll
  f001(kk,jj,0)=f000(kk,jj,1)  ! Neumann
  f011(kk,jj,0)=f010(kk,jj,1)  ! Neumann
  f101(kk,jj,0)=f100(kk,jj,1)  ! Neumann
  f111(kk,jj,0)=f110(kk,jj,1)  ! Neumann
  f000(kk,jj,mm+1)=0.          ! Dirichlet
  f010(kk,jj,mm+1)=0.          ! Dirichlet
  f100(kk,jj,mm+1)=0.          ! Dirichlet
  f110(kk,jj,mm+1)=0.          ! Dirichlet
  end do
  end do
  !$acc end kernels
  else  ! m=odd i_e=1 & i_o=0  
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111)
  !$acc loop gang, vector(128) collapse(2)
  do jj = 1, nn
  do kk = 1, ll
  f000(kk,jj,1)=f001(kk,jj,1)  ! Neumann
  f010(kk,jj,1)=f011(kk,jj,1)  ! Neumann
  f100(kk,jj,1)=f101(kk,jj,1)  ! Neumann
  f110(kk,jj,1)=f111(kk,jj,1)  ! Neumann
  f000(kk,jj,mm+1)=0.          ! Dirichlet
  f010(kk,jj,mm+1)=0.          ! Dirichlet
  f100(kk,jj,mm+1)=0.          ! Dirichlet
  f110(kk,jj,mm+1)=0.          ! Dirichlet
  end do
  end do
  !$acc end kernels
  end if

  !
  ! NO need old data replace
  !

  ! f000 corretion
  !$acc kernels async(1) present(f000,f001,f010,f100,bb000,ae000,aw000,an000,as000,at000,ab000,ap000)
  !$acc loop gang, collapse(3)
  do ii = 1+i_e, mm
    do jj = 1+j_e, nn
      do kk = 1+k_e, ll
        correct = ( bb000(kk,jj,ii)                                 &
                  - ae000(kk,jj,ii)*f001(kk,  jj,  ii  )            &
                  - aw000(kk,jj,ii)*f001(kk,  jj,  ii-1)            &
                  - an000(kk,jj,ii)*f010(kk,  jj,  ii  )            &
                  - as000(kk,jj,ii)*f010(kk,  jj-1,ii  )            &
                  - at000(kk,jj,ii)*f100(kk,  jj,  ii  )            &
                  - ab000(kk,jj,ii)*f100(kk-1,jj,  ii  ) )          &
                  / ap000(kk,jj,ii)

        f000(kk,jj,ii) = f000(kk,jj,ii) * (1.0 - relux_factor) &
                         + correct * relux_factor
      end do
    end do
  end do
  !$acc end kernels

  ! f011 corretion
  !$acc kernels async(1) present(f011,f010,f001,f111,bb011,ae011,aw011,an011,as011,at011,ab011,ap011)
  !$acc loop gang, collapse(3)
  do ii = 1, mm
    do jj = 1, nn
      do kk = 1+k_e, ll
        correct = ( bb011(kk,jj,ii)                                  &
                  - ae011(kk,jj,ii)*f010(kk,  jj,  ii+1)             &
                  - aw011(kk,jj,ii)*f010(kk,  jj,  ii  )             &
                  - an011(kk,jj,ii)*f001(kk,  jj+1,ii  )             &
                  - as011(kk,jj,ii)*f001(kk,  jj,  ii  )             &
                  - at011(kk,jj,ii)*f111(kk,  jj,  ii  )             &
                  - ab011(kk,jj,ii)*f111(kk-1,jj,  ii  ) )           &
                  / ap011(kk,jj,ii)

        f011(kk,jj,ii) = f011(kk,jj,ii) * (1.0 - relux_factor) &
                          + correct * relux_factor
      end do
    end do
  end do
  !$acc end kernels

  ! f101 corretion
  !$acc kernels async(1) present(f101,f100,f111,f001,bb101,ae101,aw101,an101,as101,at101,ab101,ap101)
  !$acc loop gang, collapse(3)
  do ii = 1, mm
    do jj = 1+j_e, nn
      do kk = 1, ll
        correct = ( bb101(kk,jj,ii)                                  &
                  - ae101(kk,jj,ii)*f100(kk,  jj,  ii+1)             &
                  - aw101(kk,jj,ii)*f100(kk,  jj,  ii  )             &
                  - an101(kk,jj,ii)*f111(kk,  jj,  ii  )             &
                  - as101(kk,jj,ii)*f111(kk,  jj-1,ii  )             &
                  - at101(kk,jj,ii)*f001(kk+1,jj,  ii  )             &
                  - ab101(kk,jj,ii)*f001(kk,  jj,  ii  ) )           &
                  / ap101(kk,jj,ii)

        f101(kk,jj,ii) = f101(kk,jj,ii) * (1.0 - relux_factor) &
                          + correct * relux_factor
      end do
    end do
  end do
  !$acc end kernels

  ! f110 corretion
  !$acc kernels async(1) present(f110,f111,f100,f010,bb110,ae110,aw110,an110,as110,at110,ab110,ap110)
  !$acc loop gang, collapse(3)
  do ii = 1+i_e, mm
    do jj = 1, nn
      do kk = 1, ll
        correct = ( bb110(kk,jj,ii)                                  &
                  - ae110(kk,jj,ii)*f111(kk,  jj,  ii  )             &
                  - aw110(kk,jj,ii)*f111(kk,  jj,  ii-1)             &
                  - an110(kk,jj,ii)*f100(kk,  jj+1,ii  )             &
                  - as110(kk,jj,ii)*f100(kk,  jj,  ii  )             &
                  - at110(kk,jj,ii)*f010(kk+1,jj,  ii  )             &
                  - ab110(kk,jj,ii)*f010(kk,  jj,  ii  ) )           &
                  / ap110(kk,jj,ii)

        f110(kk,jj,ii) = f110(kk,jj,ii) * (1.0 - relux_factor) &
                          + correct * relux_factor
      end do
    end do
  end do
  !$acc end kernels

  !$acc wait

  ! default periodic condition in y-direction
  if(j_e.eq.0)then  ! n=even j_e=0 & j_o=1
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111)
  !$acc loop gang, vector(128) collapse(2)
  do ii = 1, mm
  do kk = 1, ll
  f010(kk,   0,ii)=f010(kk,nn,ii)
  f011(kk,   0,ii)=f011(kk,nn,ii)
  f000(kk,nn+1,ii)=f000(kk,1 ,ii)
  f001(kk,nn+1,ii)=f001(kk,1 ,ii)
  f110(kk,   0,ii)=f110(kk,nn,ii)
  f111(kk,   0,ii)=f111(kk,nn,ii)
  f100(kk,nn+1,ii)=f100(kk,1 ,ii)
  f101(kk,nn+1,ii)=f101(kk,1 ,ii)
  end do
  end do
  !$acc end kernels
  else  ! n=odd j_e=1 & j_o=0  
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111)
  !$acc loop gang, vector(128) collapse(2)
  do ii = 1, mm
  do kk = 1, ll
  f000(kk,   1,ii)=f010(kk,nn,ii)
  f001(kk,   1,ii)=f011(kk,nn,ii)
  f000(kk,nn+1,ii)=f010(kk,1 ,ii)
  f001(kk,nn+1,ii)=f011(kk,1 ,ii)
  f100(kk,   1,ii)=f110(kk,nn,ii)
  f101(kk,   1,ii)=f111(kk,nn,ii)
  f100(kk,nn+1,ii)=f110(kk,1 ,ii)
  f101(kk,nn+1,ii)=f111(kk,1 ,ii)
  end do
  end do
  !$acc end kernels
  end if

  ! default periodic condition in z-direction
  if(k_e.eq.0)then  ! l=even k_e=0 & k_o=1
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111)
  !$acc loop gang, vector(128) collapse(2)
  do ii = 1, mm
  do jj = 1, nn
  f010(ll+1,jj,ii)=f010(1,jj,ii)
  f011(ll+1,jj,ii)=f011(1,jj,ii)
  f000(ll+1,jj,ii)=f000(1,jj,ii)
  f001(ll+1,jj,ii)=f001(1,jj,ii)
  f110(0,jj,ii)=f110(ll,jj,ii)
  f111(0,jj,ii)=f111(ll,jj,ii)
  f100(0,jj,ii)=f100(ll,jj,ii)
  f101(0,jj,ii)=f101(ll,jj,ii)
  end do
  end do
  !$acc end kernels
  else  ! l=odd k_e=1 & k_o=0  
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111)
  !$acc loop gang, vector(128) collapse(2)
  do ii = 1, mm
  do jj = 1, nn
  f010(ll+1,jj,ii)=f110(1,jj,ii)
  f011(ll+1,jj,ii)=f111(1,jj,ii)
  f000(ll+1,jj,ii)=f100(1,jj,ii)
  f001(ll+1,jj,ii)=f101(1,jj,ii)
  f010(1,jj,ii)=f110(ll,jj,ii)
  f011(1,jj,ii)=f111(ll,jj,ii)
  f000(1,jj,ii)=f100(ll,jj,ii)
  f001(1,jj,ii)=f101(ll,jj,ii)
  end do
  end do
  !$acc end kernels
  end if

  ! default Dirichlet(i=m+1) / Neumann(i=0) condition in x-direction
  if(i_e.eq.0)then  ! m=even i_e=0 & i_o=1
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111)
  !$acc loop gang, vector(128) collapse(2)
  do jj = 1, nn
  do kk = 1, ll
  f001(kk,jj,0   )=f000(kk,jj,1 )  ! Neumann
  f011(kk,jj,0   )=f010(kk,jj,1 )  ! Neumann
  f101(kk,jj,0   )=f100(kk,jj,1 )  ! Neumann
  f111(kk,jj,0   )=f110(kk,jj,1 )  ! Neumann
  f000(kk,jj,mm+1)=0.          ! Dirichlet
  f010(kk,jj,mm+1)=0.          ! Dirichlet
  f100(kk,jj,mm+1)=0.          ! Dirichlet
  f110(kk,jj,mm+1)=0.          ! Dirichlet
  end do
  end do
  !$acc end kernels
  else  ! m=odd i_e=1 & i_o=0  
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111)
  !$acc loop gang, vector(128) collapse(2)
  do jj = 1, nn
  do kk = 1, ll
  f000(kk,jj,1   )=f001(kk,jj,1 )  ! Neumann
  f010(kk,jj,1   )=f011(kk,jj,1 )  ! Neumann
  f100(kk,jj,1   )=f101(kk,jj,1 )  ! Neumann
  f110(kk,jj,1   )=f111(kk,jj,1 )  ! Neumann
  f000(kk,jj,mm+1)=0.          ! Dirichlet
  f010(kk,jj,mm+1)=0.          ! Dirichlet
  f100(kk,jj,mm+1)=0.          ! Dirichlet
  f110(kk,jj,mm+1)=0.          ! Dirichlet
  end do
  end do
  !$acc end kernels
  end if

  ! f100 corretion
  !$acc kernels async(1) present(f100,f101,f110,f000,bb100,ae100,aw100,an100,as100,at100,ab100,ap100)
  !$acc loop gang, collapse(3)
  do ii = 1+i_e, mm
    do jj = 1+j_e, nn
      do kk = 1, ll
        correct = ( bb100(kk,jj,ii)                                  &
                  - ae100(kk,jj,ii)*f101(kk,  jj,  ii  )             &
                  - aw100(kk,jj,ii)*f101(kk,  jj,  ii-1)             &
                  - an100(kk,jj,ii)*f110(kk,  jj,  ii  )             &
                  - as100(kk,jj,ii)*f110(kk,  jj-1,ii  )             &
                  - at100(kk,jj,ii)*f000(kk+1,jj,  ii  )             &
                  - ab100(kk,jj,ii)*f000(kk,  jj,  ii  ) )           &
                  / ap100(kk,jj,ii)

        f100(kk,jj,ii) = f100(kk,jj,ii) * (1.0 - relux_factor) &
                          + correct * relux_factor
      end do
    end do
  end do
  !$acc end kernels

  ! f111 corretion
  !$acc kernels async(1) present(f111,f010,f001,f011, bb111,ae111,aw111,an111,as111,at111,ab111,ap111)
  !$acc loop gang, collapse(3)
  do ii = 1, mm
    do jj = 1, nn
      do kk = 1, ll
        correct = ( bb111(kk,jj,ii)                                  &
                  - ae111(kk,jj,ii)*f110(kk,  jj,  ii+1)             &
                  - aw111(kk,jj,ii)*f110(kk,  jj,  ii  )             &
                  - an111(kk,jj,ii)*f101(kk,  jj+1,ii  )             &
                  - as111(kk,jj,ii)*f101(kk,  jj,  ii  )             &
                  - at111(kk,jj,ii)*f011(kk+1,jj,  ii  )             &
                  - ab111(kk,jj,ii)*f011(kk,  jj,  ii  ) )           &
                  / ap111(kk,jj,ii)

        f111(kk,jj,ii) = f111(kk,jj,ii) * (1.0 - relux_factor) &
                          + correct * relux_factor
      end do
    end do
  end do
  !$acc end kernels

  ! f001 corretion
  !$acc kernels async(1) present(f001,f000,f011,f101,bb001,ae001,aw001,an001,as001,at001,ab001,ap001)
  !$acc loop gang, collapse(3)
  do ii = 1, mm
    do jj = 1+j_e, nn
      do kk = 1+k_e, ll
        correct = ( bb001(kk,jj,ii)                                  &
                  - ae001(kk,jj,ii)*f000(kk,  jj,  ii+1)             &
                  - aw001(kk,jj,ii)*f000(kk,  jj,  ii  )             &
                  - an001(kk,jj,ii)*f011(kk,  jj,  ii  )             &
                  - as001(kk,jj,ii)*f011(kk,  jj-1,ii  )             &
                  - at001(kk,jj,ii)*f101(kk,  jj,  ii  )             &
                  - ab001(kk,jj,ii)*f101(kk-1,jj,  ii  ) )           &
                  / ap001(kk,jj,ii)

        f001(kk,jj,ii) = f001(kk,jj,ii) * (1.0 - relux_factor) &
                          + correct * relux_factor
      end do
    end do
  end do
  !$acc end kernels

  ! f010 corretion
  !$acc kernels async(1) present(f010,f011,f000,f110,bb010,ae010,aw010,an010,as010,at010,ab010,ap010)
  !$acc loop gang, collapse(3)
  do ii = 1+i_e, mm
    do jj = 1, nn
      do kk = 1+k_e, ll
        correct = ( bb010(kk,jj,ii)                                  &
                  - ae010(kk,jj,ii)*f011(kk,  jj,  ii  )             &
                  - aw010(kk,jj,ii)*f011(kk,  jj,  ii-1)             &
                  - an010(kk,jj,ii)*f000(kk,  jj+1,ii  )             &
                  - as010(kk,jj,ii)*f000(kk,  jj,  ii  )             &
                  - at010(kk,jj,ii)*f110(kk,  jj,  ii  )             &
                  - ab010(kk,jj,ii)*f110(kk-1,jj,  ii  ) )           &
                  / ap010(kk,jj,ii)

        f010(kk,jj,ii) = f010(kk,jj,ii) * (1.0 - relux_factor) &
                          + correct * relux_factor
      end do
    end do
  end do
  !$acc end kernels

  !$acc wait

  end do   ! SOR intereation finish  

  ! reconvert data  p <= f00,f01,f10,f11
  !$acc kernels present(f000,f001,f010,f011,f100,f101,f110,f111,p)
  !$acc loop gang, vector(128) collapse(3)
  do ii = 1, mm
  do jj = 1, nn
  do kk = 1, ll
  i0= (ii-1)*2 + i_o;  i1 = i0 + 1
  j0= (jj-1)*2 + j_o;  j1 = j0 + 1
  k0= (kk-1)*2 + k_o;  k1 = k0 + 1
  p(i0,j0,k0)=f000(kk,jj,ii)
  p(i1,j0,k0)=f001(kk,jj,ii)
  p(i0,j1,k0)=f010(kk,jj,ii)
  p(i1,j1,k0)=f011(kk,jj,ii)
  p(i0,j0,k1)=f100(kk,jj,ii)
  p(i1,j0,k1)=f101(kk,jj,ii)
  p(i0,j1,k1)=f110(kk,jj,ii)
  p(i1,j1,k1)=f111(kk,jj,ii)
  end do
  end do
  end do
  !$acc end kernels

  ! default Dirichlet(i=m+1) / Neumann(i=0) condition in x-direction
  !$acc kernels present(p)
  !$acc loop gang, vector(128) collapse(2)
  do j=1,n
  do k=1,l
  p(0,j,k)  = p(1,j,k)
  p(m+1,j,k)= 0.
  end do
  end do
  !$acc end kernels

  ! default periodic condition in y-direction
  !$acc kernels present(p)
  !$acc loop gang, vector(128) collapse(2)
  do i=1,m
  do k=1,l
  p(i,  0,k)= p(i,n,k)
  p(i,n+1,k)= p(i,1,k)
  end do
  end do
  !$acc end kernels

  ! default periodic condition in z-direction
  !$acc kernels present(p)
  !$acc loop gang, vector(128) collapse(2)
  do i=0,m
  do j=0,n
  p(i,j,  0)= p(i,j,l)
  p(i,j,l+1)= p(i,j,1)
  end do
  end do
  !$acc end kernels

  !$acc end data !Nakamichi 01Sep2025

  ! ----------------

  return
end subroutine solve_matrix_3x_async_1x8
!******************

subroutine  solve_matrix_RB (p, ap, ae, aw, an, as, at, ab, bb, &
  m, n, l, iter_max, relux_factor)
  use global_3d
  implicit none
  real,intent(inout),dimension(0:md,0:nd,0:ld):: p
  real,intent(in),dimension(0:md,0:nd,0:ld):: ap, ae, aw, an, as, at, ab, bb
  integer,intent(in):: m, n, l
  integer,intent(in):: iter_max
  real,intent(in):: relux_factor

  ! local variables
  real:: error, correct
  real,dimension(0:md,0:nd,0:ld):: p_old
  integer::i, j, k, iter, ii
  integer:: phase, istart

  ! ----------------
  !   SOR algorithm
  ! ----------------

  !$acc data create(p_old)

  !$acc kernels present(p,p_old)
  !$acc loop collapse(3)
  do i = 0, m+1
  do j = 0, n+1
  do k = 0, l+1
  p_old(i,j,k) = p(i,j,k)
  end do
  end do
  end do
  !$acc end kernels

  do iter = 1, iter_max
  error = 0.0

  ! default periodic condition in y-direction
  !$acc kernels present(p)
  !$acc loop collapse(2)
  do i = 1, m
  do k = 1, l
    p(i,0,k)  =p(i,n,k)
    p(i,n+1,k)=p(i,1,k)
  end do
  end do
  !$acc end kernels

  ! default periodic condition in z-direction
  !$acc kernels present(p)
  !$acc loop collapse(2)
  do i = 1, m
  do j = 1, n
    p(i,j,0)  =p(i,j,l)
    p(i,j,l+1)=p(i,j,1)
  end do
  end do
  !$acc end kernels

  ! default Dirichlet(i=m+1) / Neumann(i=0) condition in x-direction
  !$acc kernels present(p)
  !$acc loop collapse(2)
  do j=1,n
  do k=1,l
  p(0,j,k)  = p(1,j,k)
  p(m+1,j,k)= 0. !dummy
  end do
  end do
  !$acc end kernels

  !$acc kernels present(p,p_old)
  !$acc loop collapse(3)
  do i = 0, m+1
  do j = 0, n+1
  do k = 0, l+1
  p_old(i,j,k) = p(i,j,k)
  end do
  end do
  end do
  !$acc end kernels

  !-- EVEN SPACE process
  phase = 0
  !$acc kernels present(ae,aw,an,as,at,ab,ap,bb,p,p_old)
    !$acc loop gang collapse(2)
    do k = 1, l
      do j = 1, n
        istart = 1 + iand(j + k + phase, 1)   ! =1 or 2
        !$acc loop, vector(128)
        do i = istart, m, 2
          correct = ( bb(i,j,k)                                         &
                  - ae(i,j,k)*p_old(i+1,j,k) - aw(i,j,k)*p_old(i-1,j,k)&
                  - an(i,j,k)*p_old(i,j+1,k) - as(i,j,k)*p_old(i,j-1,k)&
                  - at(i,j,k)*p_old(i,j,k+1) - ab(i,j,k)*p_old(i,j,k-1)) &
                  / ap(i,j,k)

          p(i,j,k) = p_old(i,j,k)*(1.0 - relux_factor) + relux_factor*correct
        end do
      end do
    end do
  !$acc end kernels

  ! default periodic condition in y-direction
  !$acc kernels present(p)
  !$acc loop collapse(2)
  do i = 1, m
  do k = 1, l
    p(i,0,k)  =p(i,n,k)
    p(i,n+1,k)=p(i,1,k)
  end do
  end do
  !$acc end kernels

  ! default periodic condition in z-direction
  !$acc kernels present(p)
  !$acc loop collapse(2)
  do i = 1, m
  do j = 1, n
    p(i,j,0)  =p(i,j,l)
    p(i,j,l+1)=p(i,j,1)
  end do
  end do
  !$acc end kernels

  ! default Dirichlet(i=m+1) / Neumann(i=0) condition in x-direction
  !$acc kernels present(p)
  !$acc loop collapse(2)
  do j=1,n
  do k=1,l
  p(0,j,k)  = p(1,j,k)
  p(m+1,j,k)= 0. !dummy
  end do
  end do
  !$acc end kernels

  !$acc kernels present(p,p_old)
  !$acc loop collapse(3)
  do i = 0, m+1
  do j = 0, n+1
  do k = 0, l+1
  p_old(i,j,k) = p(i,j,k)
  end do
  end do
  end do
  !$acc end kernels

  !-- ODD SPACE process
  phase = 1
  !$acc kernels present(ae,aw,an,as,at,ab,ap,bb,p,p_old)
    !$acc loop gang collapse(2)
    do k = 1, l
      do j = 1, n
        istart = 1 + iand(j + k + phase, 1)   ! =1 or 2
        !$acc loop, vector(128)
        do i = istart, m, 2
          correct = ( bb(i,j,k)                                         &
                  - ae(i,j,k)*p_old(i+1,j,k) - aw(i,j,k)*p_old(i-1,j,k)&
                  - an(i,j,k)*p_old(i,j+1,k) - as(i,j,k)*p_old(i,j-1,k)&
                  - at(i,j,k)*p_old(i,j,k+1) - ab(i,j,k)*p_old(i,j,k-1)) &
                  / ap(i,j,k)

          p(i,j,k) = p_old(i,j,k)*(1.0 - relux_factor) + relux_factor*correct
        end do
      end do
    end do
  !$acc end kernels

  end do

  ! default periodic condition in y-direction
  !$acc kernels present(p)
  !$acc loop collapse(2)
  do i = 1, m
  do k = 1, l
    p(i,0,k)  =p(i,n,k)
    p(i,n+1,k)=p(i,1,k)
  end do
  end do
  !$acc end kernels

  ! default periodic condition in z-direction
  !$acc kernels present(p)
  !$acc loop collapse(2)
  do i = 1, m
  do j = 1, n
    p(i,j,0)  =p(i,j,l)
    p(i,j,l+1)=p(i,j,1)
  end do
  end do
  !$acc end kernels

  ! default Dirichlet(i=m+1) / Neumann(i=0) condition in x-direction
  !$acc kernels present(p)
  !$acc loop collapse(2)
  do j=1,n
  do k=1,l
  p(0,j,k)  = p(1,j,k)
  p(m+1,j,k)= 0. !dummy
  end do
  end do
  !$acc end kernels

  !$acc end data

  ! write(*,*)'SOR iteration no.', iter-1,'  -- error=', error
  

  return
end subroutine solve_matrix_RB


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

  ! inlet (dp/x=0 at i=1)
  !$acc kernels present(ae,aw)
  do k = 1, l
    do j = 1, n
      ae(1,j,k) =ae(1,j,k)+aw(1,j,k)
      aw(1,j,k) =0.
    end do
  end do
  !$acc end kernels
  
  ! outlet (p=outlet_pressure at i=m)
  !$acc kernels present(ae,aw,an,as,at,ab,bb,p)
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
  !$acc end kernels

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
  
  ! inlet (u=inlet_velocity, v=0., dp/dx=0 at i=1)
  !$acc kernels present(u,v,w,p)
  !$acc loop
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
  !$acc end kernels
  
  ! outlet (du/dx=0., dv/dx=0., p=outlet_pressure at i=m)
  !$acc kernels present(u,v,w,p)
  !$acc loop
  do k = 1, l
    do j = 1, n
      u(m+1,j,k) = u(m-1,j,k)
      v(m+1,j,k) = v(m-1,j,k)
      w(m+1,j,k) = w(m-1,j,k)
      p(m+1,j,k) = outlet_pressure   ! dummy
    end do
  end do
  !$acc end kernels

  ! default: periodic condition (xz-direction at j=1 & n)
  !$acc kernels present(u,v,w,p)
  !$acc loop
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
  !$acc end kernels

  ! default: periodic condition (xy-direction at k=1 & l)
  !$acc kernels present(u,v,w,p)
  !$acc loop
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
  !$acc end kernels

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

end subroutine initial_conditions
!******************
