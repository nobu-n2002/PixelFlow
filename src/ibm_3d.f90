program main
  !$ use omp_lib
   implicit none
   integer,parameter:: md = 300, nd = 300, ld = 300   ! md, nd, ld > grid size (m,n,l)
   real:: dx, dy, dz, dt
   real:: xnue, xlamda, density, width, height, depth, time
   real:: inlet_velocity, outlet_pressure, AoA, thickness
   real,dimension(0:md,0:nd,0:ld):: u, v, w, p, u_old, v_old, w_old
   real,dimension(0:md,0:nd,0:ld):: porosity
   real,dimension(0:md):: xp
   real,dimension(0:nd):: yp
   real,dimension(0:ld):: zp
   integer:: m, n, l, istep, istep_max, iset, istep_out
   integer:: i, j, k
   character(len=50) :: output_folder
   character(len=50) :: csv_file
  ! ----------------
  ! read input data by using namelist 
  namelist /file_control/istep_out
  namelist /grid_control/istep_max
  namelist /directory_control/output_folder, csv_file
  open(11,file="config/controlDict.txt",status="old",action="read")
  read(11,nml=file_control)
  read(11,nml=grid_control)
  read(11,nml=directory_control)
  close(11)
  ! ----------------
  ! write(*,*)'porosity setting:0 or calculation start:1 ?'
  ! read(*,*) iset
  ! make output directory
  
  call system('mkdir -p '//trim(output_folder))
  call system('mkdir -p etc')
  ! -----------------
  
  iset=1
  
  !-----------------
  ! porosity setting
  
  if (iset==0)then
    m=0         ! setup switch for grid conditions
    density=0.  ! setup switch for physical conditions
  
    call  physical_conditions (xnue, xlamda, density, width, height, depth, time &
                          , inlet_velocity, outlet_pressure, AoA, m, n, l)
    call  grid_conditions (xp, yp, zp, dx, dy, dz, dt, xnue, density, width, height, depth, thickness, time &
                          , inlet_velocity, AoA, porosity, m, n, l, istep_max, iset)
    ! call  output_grid_list (xp, yp, zp, m, n, l, angle_of_attack)
    stop
  end if
  
  ! ----------------
  ! calculation start  (if iest=!0)
  ! ----------------
  ! set up condistions
  m=0         ! setup switch for grid conditions
  density=0.  ! setup switch for physical conditions
  
  call  physical_conditions (xnue, xlamda, density, width, height, depth, time &
                         , inlet_velocity, outlet_pressure, AoA, m, n, l)
  call  grid_conditions (xp, yp, zp, dx, dy, dz, dt, xnue, density, width, height, depth, thickness, time &
                        , inlet_velocity, AoA, porosity, m, n, l, istep_max, iset)
  call  output_grid (xp, yp, zp, m, n, l)  
  istep = 0
  time = istep * dt
  ! ----------------
  
  write(*,*) 'istep_max= ', istep_max,'   istep_out= ', istep_out
  
  call  initial_conditions (p, u, v, w, xp, yp, zp, width, height, depth &
                         , inlet_velocity, outlet_pressure, AoA, m, n, l)
  call  boundary (p, u, v, w, xp, yp, zp, width, height, depth &
                         , inlet_velocity, outlet_pressure, AoA, porosity, m, n, l)
  
  ! print initial conditions
  ! call  output_solution (p, u, v, w, m, n, l)
  
  ! ----------------
  ! MAC algorithm start
  do istep = 1, istep_max
  
    time=istep* dt
  
    do i = 0, m+1
      do j = 0, n+1
        do k = 0, l+1
          u_old(i,j,k) = u(i,j,k)
          v_old(i,j,k) = v(i,j,k)
          w_old(i,j,k) = w(i,j,k)
        end do
      end do
    end do

    write(*,*)'--- time_steps= ',istep, ' --  time = ',time

    call solve_p (p, u, v, w, u_old, v_old, w_old, porosity, xnue, xlamda, &
                        density, height, thickness, yp, dx, dy, dz, dt, m, n, l)
    call solve_u (p, u, v, w, u_old, v_old, w_old, porosity, xnue, xlamda, density, dx, dy, dz, dt, m, n, l)
    call solve_v (p, u, v, w, u_old, v_old, w_old, porosity, xnue, xlamda, density, dx, dy, dz, dt, m, n, l)
    call solve_w (p, u, v, w, u_old, v_old, w_old, porosity, xnue, xlamda, density, dx, dy, dz, dt, m, n, l)
    call boundary(p, u, v, w, xp, yp, zp, width, height, depth  &
                      , inlet_velocity, outlet_pressure, AoA, porosity, m, n, l)
  
    if(mod(istep,istep_out)==0) call  output_paraview_temp (p, u, v, w, porosity, xp, yp, zp, m, n, l, istep)
    
    end do
    ! MAC algorithm end
    ! ----------------
    
    ! print conditions (recall)
    ! call  physical_conditions (xnue, density, width, height, depth, time, inlet_velocity, outlet_pressure, AoA, m, n, l)
    ! call  grid_conditions (xp, yp, zp, dx, dy, dz, dt, xnue, density, width, height, depth, thickness, time, inlet_velocity, AoA, porosity, m, n, l, istep_max, iset)
    
    ! print solutions 
    call  output_solution_post (p, u, v, w, xp, yp, zp, porosity, m, n, l)
    call  output_divergent (p, u, v, w, porosity, dx, dy, dz, m, n, l)
    call  output_paraview (p, u, v, w, porosity, xp, yp, zp, m, n, l)
    
    write(*,*) 'program finished'
  
  end program main
  !******************
  
  !  solve variables  
  
  !******************
  subroutine  solve_p (p, u, v, w, u_old, v_old, w_old, porosity, &
    xnue, xlamda, density, height, thickness, yp, dx, dy, dz, dt, m, n, l)
   implicit none
   integer,parameter:: md = 300, nd = 300, ld = 300
   real,intent(in):: dx, dy, dz, dt
   real,intent(in):: xnue, density, height, thickness
   real,intent(inout),dimension(0:md,0:nd,0:ld):: u, v, w, p, u_old, v_old, w_old 
   real,intent(in),dimension(0:md,0:nd,0:ld):: porosity
   real,intent(in),dimension(0:nd):: yp
   integer,intent(in):: m, n, l
  
  !-----------------
  ! local variables 
   real, parameter:: small = 1.e-6, big = 1.e6, zero = 0.0
   real, parameter::alpha = 32.0

   real:: u_stg, v_stg
   real,dimension(0:md,0:nd,0:ld):: ap, ae, aw, an, as, at, ab, bb, div
   integer:: i, j, k
   real:: fc
   real:: xlamda
  
  !-----------------
  !  divergence term  div(u)
  !-----------------
  ! ----------------
  ! read input data by using namelist 
  ! by Nobuto Nakamichi 27/7/2023
  logical::nonslip
  namelist /calculation_method/nonslip
  open(11,file="config/controlDict.txt",status="old",action="read")
  read(11,nml=calculation_method)
  close(11)
  ! ----------------
  do i = 1, m
    do j = 1, n
      do k = 1, l
        div(i,j,k)= (u_old(i+1,j,k)-u_old(i-1,j,k))/dx*0.5 &
                  + (v_old(i,j+1,k)-v_old(i,j-1,k))/dy*0.5 &
                  + (w_old(i,j,k+1)-w_old(i,j,k-1))/dz*0.5 
      end do
    end do
  end do
  
  do j = 1, n
    do k = 1, l
      div(0,j,k)  = 0.  ! inlet yz
      div(m+1,j,k)= 0.  ! outlet yz
    end do
  end do
  
  do i = 1, m
    do k = 1, l
      div(i,0,k)  = div(i,n,k)  ! periodic condition xz
      div(i,n+1,k)= div(i,1,k)
    end do
  end do

  do i = 1, m
    do j = 1, n
      div(i,j,0)  = div(i,j,l)  ! periodic condition xy
      div(i,j,l+1)= div(i,j,1)
    end do
  end do
  
  ! ----------------
  fc=0.
  
  do i = 1, m
  do j = 1, n
  do k = 1, l
  ! ----------------
  !   velocity u
  ! ----------------
  ! convection_x  (1st upwind scheme)
  u(i,j,k)=u_old(i,j,k)-dt*( &
          fc*( max(u_old(i,j,k),0.)*(u_old(i,j,k)-u_old(i-1,j,k))/dx   &      ! u>0 1st upwind scheme 
              +min(u_old(i,j,k),0.)*(u_old(i+1,j,k)-u_old(i,j,k))/dx ) &      ! u<0 1st upwind scheme
  +(1.-fc)* u_old(i,j,k)*(u_old(i+1,j,k)-u_old(i-1,j,k))/dx*0.5)    ! 2nd central scheme
  
  ! convection_y
  u(i,j,k)=u(i,j,k)-dt*( &
         fc*( max(v_old(i,j,k),0.)*(u_old(i,j,k)-u_old(i,j-1,k))/dy   &   ! v>0 1st upwind scheme 
              +min(v_old(i,j,k),0.)*(u_old(i,j+1,k)-u_old(i,j,k))/dy) &   ! v<0 1st upwind scheme
  +(1.-fc)* v_old(i,j,k)*(u_old(i,j+1,k)-u_old(i,j-1,k))/dy*0.5) ! 2nd central scheme

  ! convection_w
  u(i,j,k)=u(i,j,k)-dt*( &
          fc*( max(w_old(i,j,k),0.)*(u_old(i,j,k)-u_old(i,j,k-1))/dz   &      ! w>0 1st upwind scheme 
              +min(w_old(i,j,k),0.)*(u_old(i,j,k-1)-u_old(i,j,k))/dz ) &      ! w<0 1st upwind scheme
  +(1.-fc)* w_old(i,j,k)*(u_old(i,j,k+1)-u_old(i,j,k-1))/dz*0.5)    ! 2nd central scheme
  
  ! diffusion_x
  u(i,j,k)=u(i,j,k) +dt*xnue*(u_old(i+1,j,k)-2.*u_old(i,j,k)+u_old(i-1,j,k))/dx/dx 
  !      +dt*xnue/(small+porosity(i,j,k))*(u_old(i+1,j,k)-u_old(i-1,j,k))*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx/dx*0.25 ! non-conseved term
  ! diffusion_y
  u(i,j,k)=u(i,j,k) +dt*xnue*(u_old(i,j+1,k)-2.*u_old(i,j,k)+u_old(i,j-1,k))/dy/dy
  !      +dt*xnue/(small+porosity(i,j,k))*(u_old(i,j+1,k)-u_old(i,j-1,k))*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy/dy*0.25 ! non-conseved term
  ! diffusion_z
  u(i,j,k)=u(i,j,k) +dt*xnue*(u_old(i,j,k+1)-2.*u_old(i,j,k)+u_old(i,j,k-1))/dz/dz
  !      +dt*xnue/(small+porosity(i,j,k))*(u_old(i,j,k+1)-u_old(i,j,k-1))*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz/dz*0.25 ! non-conseved term
  
  ! divergence term
  u(i,j,k)=u(i,j,k) +dt*(xnue + xlamda)*(div(i+1,j,k)-div(i-1,j,k))/dx*.5

  ! additional terms by porosity profile
  u(i,j,k)=u(i,j,k)                 &
        +dt*( ( (u_old(i+1,j,k)-u_old(i-1,j,k))/dx*.5+(u_old(i+1,j,k)-u_old(i-1,j,k))/dx*.5) &
                *xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*.5                            &
             +( (u_old(i,j+1,k)-u_old(i,j-1,k))/dy*.5+(v_old(i+1,j,k)-v_old(i-1,j,k))/dx*.5) &
                *xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*.5                            &
             +( (u_old(i,j,k+1)-u_old(i,j,k-1))/dz*.5+(w_old(i+1,j,k)-w_old(i-1,j,k))/dx*.5) &
                *xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*.5                            &
             + div(i,j,k)*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*0.5*xlamda         &
         )/porosity(i,j,k)
  ! force on wall
  if (nonslip) then
    u(i,j,k)=u(i,j,k)- dt*xnue*u_old(i,j,k)/(thickness*dx)**2*alpha*porosity(i,j,k)*(1.-porosity(i,j,k))*(1.-porosity(i,j,k))
  end if
  ! ----------------
  !   velocity v
  ! ----------------
  ! convection_x  (1st upwind scheme)
  v(i,j,k)=v_old(i,j,k)-dt*(          &
        fc *(max(u_old(i,j,k),0.)*(v_old(i,j,k)-v_old(i-1,j,k))/dx   &  ! u>0 1st upwind scheme
            +min(u_old(i,j,k),0.)*(v_old(i+1,j,k)-v_old(i,j,k))/dx)  &  ! u<0 1st upwind scheme
   +(1.-fc)* u_old(i,j,k)*(v_old(i+1,j,k)-v_old(i-1,j,k))/dx*0.5)       ! 2nd central scheme

  ! convection_y
  v(i,j,k)=v(i,j,k)-dt*(          &
        fc *(max(v_old(i,j,k),0.)*(v_old(i,j,k)-v_old(i,j-1,k))/dy   &  ! v>0 
            +min(v_old(i,j,k),0.)*(v_old(i,j+1,k)-v_old(i,j,k))/dy)  &  ! v<0
   +(1.-fc)* v_old(i,j,k)*(v_old(i,j+1,k)-v_old(i,j-1,k))/dy*0.5)       ! 2nd central scheme

  ! convection_z
  v(i,j,k)=v(i,j,k)-dt*(          &
        fc *(max(w_old(i,j,k),0.)*(v_old(i,j,k)-v_old(i,j,k-1))/dz   &  ! v>0 
            +min(w_old(i,j,k),0.)*(v_old(i,j,k+1)-v_old(i,j,k))/dz)  &  ! v<0
   +(1.-fc)* w_old(i,j,k)*(v_old(i,j,k+1)-v_old(i,j,k-1))/dz*0.5)       ! 2nd central scheme
  
  ! diffusion_x
  v(i,j,k)=v(i,j,k) +dt*xnue*(v_old(i+1,j,k)-2.*v_old(i,j,k)+v_old(i-1,j,k))/dx/dx
  !      +dt*xnue/(small+porosity(i,j,k))*(v_old(i+1,j,k)-v_old(i-1,j,k))*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx/dx*0.25 ! non-conseved term
  ! diffusion_y
  v(i,j,k)=v(i,j,k) +dt*xnue*(v_old(i,j+1,k)-2.*v_old(i,j,k)+v_old(i,j-1,k))/dy/dy
  !      +dt*xnue/(small+porosity(i,j,k))*(v_old(i,j+1,k)-v_old(i,j-1,k))*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy/dy*0.25 ! non-conseved term
  ! diffusion_z
  v(i,j,k)=v(i,j,k) +dt*xnue*(v_old(i,j,k+1)-2.*v_old(i,j,k)+v_old(i,j,k-1))/dz/dz
  !      +dt*xnue/(small+porosity(i,j,k))*(v_old(i,j,k+1)-v_old(i,j,k-1))*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz/dz*0.25 ! non-conseved term
  ! divergence term   ! L+(2/3)N = (1/3)N;(2/3) or 0(1/3)
  v(i,j,k)=v(i,j,k) +dt*(xnue + xlamda)*(div(i,j+1,k)-div(i,j-1,k))/dy*.5
  ! additional terms by porosity profile   ! canceled for non-slip condition    ! L+(2/3)N = (1/3)N;(-1/3) or 0:(-2/3)N
  v(i,j,k)=v(i,j,k)               &
        +dt*( ( (v_old(i+1,j,k)-v_old(i-1,j,k))/dx*.5+(u_old(i,j+1,k)-u_old(i,j-1,k))/dy*.5) &
                *xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*.5                            &
             +( (v_old(i,j+1,k)-v_old(i,j-1,k))/dy*.5+(v_old(i,j+1,k)-v_old(i,j-1,k))/dy*.5) &
                *xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*.5                            &
             +( (v_old(i,j,k+1)-v_old(i,j,k-1))/dz*.5+(w_old(i,j+1,k)-w_old(i,j-1,k))/dy*.5) &
                *xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*.5                            &
           + div(i,j,k)*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*0.5*xlamda           &
         )/porosity(i,j,k)
  ! force on wall
  if (nonslip) then
    v(i,j,k)=v(i,j,k)- dt*xnue*v_old(i,j,k)/(thickness*dy)**2*alpha*porosity(i,j,k)*(1.-porosity(i,j,k))*(1.-porosity(i,j,k))
  end if
  ! ----------------
  !   velocity w
  ! ----------------
  ! convection_x  (1st upwind scheme)
  w(i,j,k)=w_old(i,j,k)-dt*(          &
        fc *(max(u_old(i,j,k),0.)*(w_old(i,j,k)-w_old(i-1,j,k))/dx        &  ! w>0 1st upwind scheme
            +min(u_old(i,j,k),0.)*(w_old(i+1,j,k)-w_old(i,j,k))/dx)       &  ! w<0 1st upwind scheme
   +(1.-fc)* u_old(i,j,k)*(w_old(i+1,j,k)-w_old(i-1,j,k))/dx/2.)             ! 2nd central scheme
  
  ! convection_y
  w(i,j,k)=w(i,j,k)-dt*(          &
        fc *(max(v_old(i,j,k),0.)*(w_old(i,j,k)-w_old(i,j-1,k))/dy   &  ! w>0 
            +min(v_old(i,j,k),0.)*(w_old(i,j+1,k)-w_old(i,j,k))/dy)  &  ! w<0
   +(1.-fc)* v_old(i,j,k)*(w_old(i,j+1,k)-w_old(i,j-1,k))/dy/2.)        ! 2nd central scheme

  ! convection_z
  w(i,j,k)=w(i,j,k)-dt*(          &
        fc *(max(w_old(i,j,k),0.)*(w_old(i,j,k)-w_old(i,j,k-1))/dz   &  ! w>0 
            +min(w_old(i,j,k),0.)*(w_old(i,j,k+1)-w_old(i,j,k))/dz)  &  ! w<0
   +(1.-fc)* w_old(i,j,k)*(w_old(i,j,k+1)-w_old(i,j,k-1))/dz/2.   )     ! 2nd central scheme
  
  ! diffusion_x
  w(i,j,k)=w(i,j,k) +dt*xnue*(w_old(i+1,j,k)-2.*w_old(i,j,k)+w_old(i-1,j,k))/dx/dx
  !      +dt*xnue/(small+porosity(i,j,k))*(w_old(i+1,j,k)-w_old(i-1,j,k))*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx/dx*0.25 ! non-conseved term
  ! diffusion_y
  w(i,j,k)=w(i,j,k) +dt*xnue*(w_old(i,j+1,k)-2.*w_old(i,j,k)+w_old(i,j-1,k))/dy/dy
  !      +dt*xnue/(small+porosity(i,j))*(w_old(i,j+1,k)-w_old(i,j-1,k))*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy/dy*0.25 ! non-conseved term
  ! diffusion_z
  w(i,j,k)=w(i,j,k) +dt*xnue*(w_old(i,j,k+1)-2.*w_old(i,j,k)+w_old(i,j,k-1))/dz/dz
  !      +dt*xnue/(small+porosity(i,j))*(w_old(i,j,k+1)-w_old(i,j,k-1))*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz/dz*0.25 ! non-conseved term
  ! divergence term   ! L+(2/3)N = (1/3)N;(2/3) or 0(1/3)
  w(i,j,k)=w(i,j,k) +dt*(xnue + xlamda)*(div(i,j,k+1)-div(i,j,k-1))/dz*.5
  ! additional terms by porosity profile   ! canceled for non-slip condition    ! L+(2/3)N = (1/3)N;(-1/3) or 0:(-2/3)N
  w(i,j,k)=w(i,j,k)               &
        +dt*( ( (w_old(i+1,j,k)-w_old(i-1,j,k))/dx*.5+(u_old(i,j,k+1)-u_old(i,j,k-1))/dz*.5) &
                *xnue*(porosity(i+1,j,k)-porosity(i-1,j,k))/dx*.5                            &
             +( (w_old(i,j+1,k)-w_old(i,j-1,k))/dy*.5+(v_old(i,j,k+1)-v_old(i,j,k-1))/dz*.5) &
                *xnue*(porosity(i,j+1,k)-porosity(i,j-1,k))/dy*.5                            &
             +( (w_old(i,j,k+1)-w_old(i,j,k-1))/dz*.5+(w_old(i,j,k+1)-w_old(i,j,k-1))/dz*.5) &
                *xnue*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*.5                            &
           + div(i,j,k)*(porosity(i,j,k+1)-porosity(i,j,k-1))/dz*0.5*xlamda           &
         )/porosity(i,j,k)
  ! force on wall
  if (nonslip) then
    w(i,j,k)=w(i,j,k)- dt*xnue*w_old(i,j,k)/(thickness*dz)**2*alpha*porosity(i,j,k)*(1.-porosity(i,j,k))*(1.-porosity(i,j,k))
  end if

  end do
  end do
  end do
  
  ! ----------------
  ! matrix solution  !  formulation of porous media
  
  do i = 1, m
  do j = 1, n
  do k = 1, l
  ae(i,j,k)= dt*max(small,(porosity(i+1,j,k)+porosity(i,j,k))*0.5)/dx/dx

  aw(i,j,k)= dt*max(small,(porosity(i,j,k)+porosity(i-1,j,k))*0.5)/dx/dx

  an(i,j,k)= dt*max(small,(porosity(i,j+1,k)+porosity(i,j,k))*0.5)/dy/dy

  as(i,j,k)= dt*max(small,(porosity(i,j,k)+porosity(i,j-1,k))*0.5)/dy/dy

  at(i,j,k)= dt*max(small,(porosity(i,j,k+1)+porosity(i,j,k))*0.5)/dz/dz

  ab(i,j,k)= dt*max(small,(porosity(i,j,k)+porosity(i,j,k-1))*0.5)/dz/dz

  ap(i,j,k)= -ae(i,j,k)-aw(i,j,k)-an(i,j,k)-as(i,j,k)-at(i,j,k)-ab(i,j,k)

  bb(i,j,k)= ((porosity(i+1,j,k)*u(i,j,k)+porosity(i,j,k)*u(i+1,j,k))*0.5             &
             -(porosity(i-1,j,k)*u(i,j,k)+porosity(i,j,k)*u(i-1,j,k))*0.5)*density/dx &
            +((porosity(i,j+1,k)*v(i,j,k)+porosity(i,j,k)*v(i,j+1,k))*0.5             &
             -(porosity(i,j-1,k)*v(i,j,k)+porosity(i,j,k)*v(i,j-1,k))*0.5)*density/dy &
            +((porosity(i,j,k+1)*w(i,j,k)+porosity(i,j,k)*w(i,j,k+1))*0.5             &
             -(porosity(i,j,k-1)*w(i,j,k)+porosity(i,j,k)*w(i,j,k-1))*0.5)*density/dz 

  end do
  end do
  end do
  
  call boundrary_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, l, height, yp)

  ! call solve_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, l)
  call solve_matrix_vec_omp (p, ap, ae, aw, an, as, at, ab, bb, m, n, l)
  ! call solve_matrix_vec_oacc (p, ap, ae, aw, an, as, at, ab, bb, m, n, l)
  ! ----------------
  ! ----------------
  return
  end subroutine solve_p
  !******************
 
  !******************
  ! OpenMP Parallelized
  ! Written only for CPU machine
  ! No efficiency ensured on GPU machine 
  subroutine  solve_matrix_vec_omp (p, ap, ae, aw, an, as, at, ab, bb, m, n, l)
   implicit none
   integer,parameter:: md=300, nd = 300, ld = 300
   real,intent(inout),dimension(0:md,0:nd,0:ld):: p
   real,intent(in),dimension(0:md,0:nd,0:ld):: ap, ae, aw, an, as, at, ab, bb
   integer,intent(in):: m, n, l
  
  ! local variables
  real:: relux_factor, error
  real,dimension(0:md,0:nd,0:ld):: p_old
  integer::i, j, k, iter, iter_max, ii

  !$omp parallel private(iter, i, j, k, ii) &
  !$omp & shared(iter_max, relux_factor, m, n, l) &
  !$omp & shared(error, p_old, p, ap, ae, aw, an, as, at, ab, bb) &
  !$omp & default(none)
  
  ! ----------------
  !   SOR algorithm
  ! ----------------
  !$omp single
  iter_max = 300 ! SOR max interation steps
  relux_factor=1.7 ! SOR reluxation factor
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
  
  !-- EVEN SPACE process
  !$omp do reduction(max: error)
  do ii = 2, m*n*l, 2 ! evenspace
    k = (ii - 1) / (m * n)  + 1
    j = ((ii - 1) / m + 1) - (k - 1) * n
    i = (ii - (j - 1) * m) - (k - 1) * m * n

  !-- IF m is EVEN (Based on Column-Major Order; FORTRAN)
  if(mod(m,2)==0 .and. (mod(j,2)==0 .or. mod(k,2)==0)) i = i - 1
  p(i,j,k) = (bb(i,j,k) &
            - ae(i,j,k)*p_old(i+1,j,k)-aw(i,j,k)*p(i-1,j,k)  &
            - an(i,j,k)*p_old(i,j+1,k)-as(i,j,k)*p(i,j-1,k)  &
            - at(i,j,k)*p_old(i,j,k+1)-ab(i,j,k)*p(i,j,k-1)) &
            / ap(i,j,k)*relux_factor &
            + p_old(i,j,k)*(1.-relux_factor)
  error = max(error, abs(p(i,j,k)-p_old(i,j,k)))
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
  
  !-- ODD SPACE process
  !$omp do reduction(max: error)
  do ii = 1, m*n*l, 2 ! odd space
    k = (ii - 1) / (m * n)  + 1
    j = ((ii - 1) / m + 1) - (k - 1) * n
    i = (ii - (j - 1) * m) - (k - 1) * m * n
  
    !-- IF m is EVEN (Based on Column-Major Order; FORTRAN)
    if(mod(m,2)==0 .and. (mod(j,2)==0 .or. mod(k,2)==0)) i = i + 1
    p(i,j,k) = (bb(i,j,k) &
              - ae(i,j,k)*p_old(i+1,j,k)-aw(i,j,k)*p(i-1,j,k)  &
              - an(i,j,k)*p_old(i,j+1,k)-as(i,j,k)*p(i,j-1,k)  &
              - at(i,j,k)*p_old(i,j,k+1)-ab(i,j,k)*p(i,j,k-1)) &
              / ap(i,j,k)*relux_factor &
              + p_old(i,j,k)*(1.-relux_factor)
    error = max(error, abs(p(i,j,k)-p_old(i,j,k)))
  end do
  !$omp end do

  ! write(*,*)'CHECK iteration no.', iter,'  -- error=', error
  
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
  ! OpenACC Parallelized
  ! Written only for GPU machine
  ! No efficiency ensured on CPU machine 
  subroutine  solve_matrix_vec_oacc (p, ap, ae, aw, an, as, at, ab, bb, m, n, l)
   implicit none
   integer,parameter:: md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(inout),dimension(0:md,0:nd,0:ld):: p
   real,intent(in),dimension(0:md,0:nd,0:ld):: ap, ae, aw, an, as, at, ab, bb
   integer,intent(in):: m, n, l
  
  ! local variables
  real:: relux_factor
  real,dimension(0:md,0:nd,0:ld):: p_old
  integer::i, j, k, iter, iter_max, ii

  !$acc data copy(p_old, p) &
  !$acc & copyin(ap, ae, aw, an, as, at, ab, bb, relux_factor) 
  
  ! ----------------
  !   SOR algorithm
  ! ----------------

  iter_max = 300 ! SOR max interation steps
  relux_factor=1.7 ! SOR reluxation factor

  do iter = 1, iter_max
  ! write(*,*)'CHECK iteration no.'


  ! default periodic condition in yz-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
  !$acc loop independent
    do k = 1, l
      p(i,0,k)  =p(i,n,k)
      p(i,n+1,k)=p(i,1,k)
    end do
  end do
  !$acc end kernels

  ! default periodic condition in xy-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
  !$acc loop independent
    do j = 1, n
      p(i,j,0)  =p(i,j,l)
      p(i,j,l+1)=p(i,j,1)
    end do
  end do
  !$acc end kernels

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

  !-- IF m is EVEN (Based on Column-Major Order; FORTRAN)
  if(mod(m,2)==0 .and. (mod(j,2)==0 .or. mod(k,2)==0)) i = i - 1
  p(i,j,k) = (bb(i,j,k) &
            - ae(i,j,k)*p_old(i+1,j,k)-aw(i,j,k)*p(i-1,j,k)  &
            - an(i,j,k)*p_old(i,j+1,k)-as(i,j,k)*p(i,j-1,k)  &
            - at(i,j,k)*p_old(i,j,k+1)-ab(i,j,k)*p(i,j,k-1)) &
            / ap(i,j,k)*relux_factor &
            + p_old(i,j,k)*(1.-relux_factor)
  
  end do
  !$acc end kernels

  ! default periodic condition in yz-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
  !$acc loop independent
    do k = 1, l
      p(i,0,k)  =p(i,n,k)
      p(i,n+1,k)=p(i,1,k)
    end do
  end do
  !$acc end kernels

  ! default periodic condition in xy-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
  !$acc loop independent
    do j = 1, n
      p(i,j,0)  =p(i,j,l)
      p(i,j,l+1)=p(i,j,1)
    end do
  end do
  !$acc end kernels

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
  
    !-- IF m is EVEN (Based on Column-Major Order; FORTRAN)
    if(mod(m,2)==0 .and. (mod(j,2)==0 .or. mod(k,2)==0)) i = i + 1
    p(i,j,k) = (bb(i,j,k) &
              - ae(i,j,k)*p_old(i+1,j,k)-aw(i,j,k)*p(i-1,j,k)  &
              - an(i,j,k)*p_old(i,j+1,k)-as(i,j,k)*p(i,j-1,k)  &
              - at(i,j,k)*p_old(i,j,k+1)-ab(i,j,k)*p(i,j,k-1)) &
              / ap(i,j,k)*relux_factor &
              + p_old(i,j,k)*(1.-relux_factor)
    
  end do
  !$acc end kernels

  ! write(*,*)'CHECK iteration no.', iter,'  -- error=', error
  
  end do

  ! default periodic condition in yz-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
  !$acc loop independent
    do k = 1, l
      p(i,0,k)  =p(i,n,k)
      p(i,n+1,k)=p(i,1,k)
    end do
  end do
  !$acc end kernels

  ! default periodic condition in xy-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
  !$acc loop independent
    do j = 1, n
      p(i,j,0)  =p(i,j,l)
      p(i,j,l+1)=p(i,j,1)
    end do
  end do
  !$acc end kernels

  !$acc end data 

  ! write(*,*)'SOR iteration no.', iter-1,'  -- error=', error
  
  return
  end subroutine solve_matrix_vec_oacc
  !******************
 
  !******************
  subroutine  solve_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, l)
   implicit none
   integer,parameter:: md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(inout),dimension(0:md,0:nd,0:ld):: p
   real,intent(in),dimension(0:md,0:nd,0:ld):: ap, ae, aw, an, as, at, ab, bb
   integer,intent(in):: m, n, l
  
  ! local variables
  real:: relux_factor, error
  real,dimension(0:md,0:nd,0:ld):: p_old
  integer::i, j, k, iter, iter_max
  
  ! ----------------
  !   SOR algorithm
  ! ----------------
  iter_max = 300 ! SOR max interation steps
  relux_factor=1.7 ! SOR reluxation factor
  
  do iter = 1, iter_max
  ! write(*,*)'CHECK iteration no.'
  error=0.
  
  ! default periodic condition in yz-direction
  do i = 1, m
    do k = 1, l
      p(i,0,k)   = p(i,n,k)
      p(i,n+1,k) = p(i,1,k)
    end do
  end do

  ! default periodic condition in xy-direction
  do i = 1, m
    do j = 1, n
      p(i,j,0)   = p(i,j,l)
      p(i,j,l+1) = p(i,j,1)
    end do
  end do
  
  do i = 0, m+1
    do j = 0, n+1
      do k = 0, l+1
        p_old(i,j,k) = p(i,j,k)
      end do
    end do
  end do
  
  do i = 1, m
  do j = 1, n
  do k = 1, l
  p(i,j,k) = (bb(i,j,k) &
            - ae(i,j,k)*p_old(i+1,j,k)-aw(i,j,k)*p(i-1,j,k)  &
            - an(i,j,k)*p_old(i,j+1,k)-as(i,j,k)*p(i,j-1,k)  &
            - at(i,j,k)*p_old(i,j,k+1)-ab(i,j,k)*p(i,j,k-1)) &
            / ap(i,j,k)*relux_factor &
            + p_old(i,j,k)*(1.-relux_factor)
  
  error = max(error, abs(p(i,j,k)-p_old(i,j,k)))
  end do
  end do
  end do
  
  
  ! write(*,*)'CHECK iteration no.', iter,'  -- error=', error
  
  end do
  
  write(*,*)'SOR iteration no.', iter-1,'  -- error=', error
  if (error > 1e5) then
    write(*,*)'Error value diverges. Terminate the process.'
    call exit(0)
  end if
  
  ! ----------------
  
  return
  end subroutine solve_matrix
  !******************
  
  !******************
  subroutine  boundrary_matrix (p, ap, ae, aw, an, as, at, ab, bb, m, n, l, height, yp)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(in)::height
   real,intent(in),dimension(0:md,0:nd,0:ld)::p
   real,intent(inout),dimension(0:md,0:nd,0:ld)::ap, ae, aw, an, as, at, ab, bb
   real,intent(in),dimension(0:nd)::yp
   integer,intent(in)::m, n, l
  
  ! local variables
  integer::i, j, k
  
  ! ----------------
  ! inlet (dp/x=0 at i=1)
  do j= 1, n
    do k = 1, l
      ae(1,j,k) =ae(1,j,k)+aw(1,j,k)
      aw(1,j,k) =0.
    end do
  end do
  
  ! outlet (p=outlet_pressure at i=m)
  do j= 1, n
    do k = 1, l
      bb(m,j,k) = bb(m,j,k)+ae(m,j,k)*p(m+1,j,k)
      ae(m,j,k) = 0.
      aw(m,j,k) = 0.
      an(m,j,k) = 0.
      as(m,j,k) = 0.
      at(m,j,k) = 0.
      ab(m,j,k) = 0.
    end do
  end do
  ! ----------------
  
  return
  end subroutine  boundrary_matrix
  !******************
  
  !******************
  subroutine  solve_u (p, u, v, w, u_old, v_old, w_old, porosity, xnue, xlamda, density, dx, dy, dz, dt, m, n, l)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(in)::dx, dy, dz, dt
   real,intent(in)::xnue, density, xlamda
   real,intent(inout),dimension(0:md,0:nd,0:ld)::u, v, w, p, u_old, v_old, w_old
   real,intent(in),dimension(0:md,0:nd,0:ld)::porosity
   integer,intent(in)::m, n, l
  
  ! local variables
  integer::i, j, k
  
  ! ----------------
  do i = 1, m
  do j = 1, n
  do k = 1, l
  ! convection_x
  ! (already calculated in solve_p)
  
  ! convection_y
  ! (already calculated in solve_p)

  ! convection_z
  ! (already calculated in solve_p)
  
  ! diffusion_x
  ! (already calculated in solve_p)
  
  ! diffusion_y
  ! (already calculated in solve_p)
  
  ! diffusion_z
  ! (already calculated in solve_p)
  
  ! pressure
  u(i,j,k)=u(i,j,k) -dt/density*(p(i+1,j,k)-p(i-1,j,k))/dx*0.5
  end do
  end do
  end do
  
  ! ----------------
  return
  end subroutine solve_u
  !******************
  
  !******************
  subroutine  solve_v (p, u, v, w, u_old, v_old, w_old, porosity, xnue, xlamda, density, dx, dy, dz, dt, m, n, l)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(in)::dx, dy, dz, dt
   real,intent(in)::xnue, density, xlamda
   real,intent(inout),dimension(0:md,0:nd,0:ld)::u, v, w, p, u_old, v_old, w_old
   real,intent(in),dimension(0:md,0:nd,0:ld)::porosity
   integer,intent(in)::m, n, l
  
  ! local variables
  integer::i, j, k
  
  ! ----------------
  do i = 1, m
  do j = 1, n
  do k = 1, l
  ! convection_x
  ! (already calculated in solve_p)
  
  ! convection_y
  ! (already calculated in solve_p)

  ! convection_z
  ! (already calculated in solve_p)
  
  ! diffusion_x
  ! (already calculated in solve_p)
  
  ! diffusion_y
  ! (already calculated in solve_p)
  
  ! diffusion_z
  ! (already calculated in solve_p)
  
  ! pressure
  v(i,j,k)=v(i,j,k) -dt/density*(p(i,j+1,k)-p(i,j-1,k))/dy*.5
  end do
  end do
  end do
  ! ----------------
  return
  end subroutine solve_v
  !******************

  !******************
  subroutine  solve_w (p, u, v, w, u_old, v_old, w_old, porosity, xnue, xlamda, density, dx, dy, dz, dt, m, n, l)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(in)::dx, dy, dz, dt
   real,intent(in)::xnue, density, xlamda
   real,intent(inout),dimension(0:md,0:nd,0:ld)::u, v, w, p, u_old, v_old, w_old
   real,intent(in),dimension(0:md,0:nd,0:ld)::porosity
   integer,intent(in)::m, n, l
  
  ! local variables
  integer::i, j, k
  
  ! ----------------
  do i = 1, m
  do j = 1, n
  do k = 1, l
  ! convection_x
  ! (already calculated in solve_p)
  
  ! convection_y
  ! (already calculated in solve_p)

  ! convection_z
  ! (already calculated in solve_p)
  
  ! diffusion_x
  ! (already calculated in solve_p)
  
  ! diffusion_y
  ! (already calculated in solve_p)
  
  ! diffusion_z
  ! (already calculated in solve_p)

  ! pressure
  w(i,j,k)=w(i,j,k) -dt/density*(p(i,j,k+1)-p(i,j,k-1))/dz*.5
  end do
  end do
  end do
  ! ----------------
  return
  end subroutine solve_w
  !******************
  
  !  conditions  
  
  !******************
  subroutine  boundary(p, u, v, w, xp, yp, zp, width, height, depth    &
                       , inlet_velocity, outlet_pressure, AoA, porosity, m, n, l)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
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
  
  ! ----------------
  ! inlet (u=inlet_velocity, v=0., dp/dx=0 at i=1)
  do j = 1, n
    do k = 1, l
      u(1,j,k) =inlet_velocity*cos(AoA/1300.*pai)
      v(1,j,k) =inlet_velocity*sin(AoA/1300.*pai)
      w(1,j,k) =0.
      u(0,j,k) =u(1,j,k)    ! dummy
      v(0,j,k) =v(1,j,k)    ! dummy
      w(0,j,k) =w(1,j,k)    ! dummy
      p(0,j,k) =p(2,j,k)
    end do
  end do
  
  ! outlet (du/dx=0., dv/dx=0., p=outlet_pressure at i=m)
  do j = 1, n
    do k = 1, l
      u(m+1,j,k) =u(m-1,j,k)
      v(m+1,j,k) =v(m-1,j,k)
      w(m+1,j,k) =w(m-1,j,k)
      p(m+1,j,k)=outlet_pressure   ! dummy
    end do
  end do

  ! default: periodic condition (xz-direction at j=1 & n)
  do i= 0, m+1
    do k = 0, l+1
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

  ! default: periodic condition (xy-direction at k=1 & l)
  do i= 0, m+1
    do j = 0, n+1
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
  
  return
  end subroutine boundary
  !*****************************
  
  !*****************************
  subroutine physical_conditions(xnue, xlamda, density, width, height, depth, time &
                                , inlet_velocity, outlet_pressure, AoA, m, n, l)
   implicit none
   integer,parameter:: md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(inout):: xnue, xlamda, density
   real,intent(inout):: width, height, depth
   real,intent(inout):: time, inlet_velocity, outlet_pressure, AoA
   integer,intent(in):: m, n, l
  ! local variables
  ! integer:: i, j, k
  
  ! ----------------
  
  ! ----------------
  ! read input file 
  ! by Nobuto Nakamichi 4/7/2023
  namelist /physical/xnue, xlamda, density, width, height, depth, time  &
                   ,inlet_velocity, outlet_pressure, AoA
  
  if (density == 0.)then
  
    open(11,file="config/controlDict.txt",status="old",action="read")
    read(11,nml=physical)
    close(11)
  
  end if
  
  ! ----------------
  
  write(*,*) 
  write(*,*) 'xnue =', xnue
  write(*,*) 'xlamda =', xlamda
  write(*,*) 'density =', density
  write(*,*) 'width =', width
  write(*,*) 'height =', height
  write(*,*) 'depth =', depth
  write(*,*) 'time =', time
  write(*,*) 'inlet_velocity =', inlet_velocity
  write(*,*) 'outlet_pressure =', outlet_pressure
  write(*,*) 'Angle of inlet_velocity (AoA) =', AoA
  
  ! ----------------
  
  return
  end subroutine physical_conditions
  !******************
  
  !******************
  subroutine  grid_conditions (xp, yp, zp, dx, dy, dz, dt, xnue, density, width, height, depth, thickness, time &
                              , inlet_velocity, AoA, porosity, m, n, l, istep_max, iset)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(inout)::dx, dy, dz, dt, AoA, thickness
   real,intent(in)::xnue, density, width, height, depth, time, inlet_velocity
   real,intent(inout),dimension(0:md,0:nd,0:ld)::porosity
   real,intent(inout),dimension(0:md):: xp
   real,intent(inout),dimension(0:nd):: yp
   real,intent(inout),dimension(0:ld):: zp
   integer,intent(inout):: m, n, l, istep_max, iset
   character(len = 50) :: csv_file
   character(len = 50) :: output_folder
  
  ! local variables
  !real,dimension(0:md,0:nd)::	distance
  real:: cfl_no, pecret_no, diffusion_factor
  integer:: i, j, k
  integer:: x, y, z
  real:: poro_val, threshold
  real, parameter:: small=1.e-6, big=1.e6, zero=0.
  ! --- 
  
  ! ----------------
  ! namelist 
  ! by Nobuto Nakamichi 4/7/2023
  namelist /grid_control/istep_max
  namelist /porosity_control/thickness, threshold
  namelist /directory_control/csv_file, output_folder
  open(11,file="config/controlDict.txt",status="old",action="read")
  read(11,nml=grid_control)
  read(11,nml=porosity_control)
  read(11,nml=directory_control)
  close(11)
  !-----------------
  
  ! read pixel data
  open(52,file=csv_file, form='formatted')
  
  read(52,*) m,n,l

  do k = 1, l
    do j = 1, n
      do i = 1, m
        read(52, *) x, y, z, poro_val
        porosity(x, y, z) = max(poro_val, threshold)
      end do
    end do
  end do
  
  close(52) 
  
  dx = width / real(m-1)
  dy = height / real(n-1)
  dz = depth / real(l-1)
  dt = time / real(istep_max)
  
  cfl_no           = inlet_velocity * dt / dx
  pecret_no        = inlet_velocity * dx / xnue
  diffusion_factor = xnue * dt / dy / dy
  
  !----- check print out
  write(*,*)
  write(*,*) 'm, n, l =', m, n, l
  write(*,*) 'istep_max =', istep_max
  write(*,*) 'dx, dy, dz =', dx, dy, dz
  write(*,*) 'dt =', dt
  write(*,*) 'cfl_no =', cfl_no
  write(*,*) 'pecret_no =', pecret_no
  write(*,*) 'diffusion_factor =', diffusion_factor
  write(*,*) 'thickness =', thickness
  write(*,*) 'threshold =', threshold
  
  do i = 0, m+1
    xp(i) = dx * real(i-1) - width*0.5
  end do
  
  do j = 0, n+1
    yp(j) = dy * real(j-1) - height*0.5
  end do

  do k = 0, l+1
    zp(k) = dz * real(k-1) - depth*0.5
  end do
  
  ! default: outlet condtion in x-direction
  do j = 1, n+1
  do k = 1, l+1
  porosity(0,j,k) = porosity(1,j,k)
  porosity(m+1,j,k) = porosity(m,j,k)
  end do
  end do
  
  ! default: periodic condtion in y-direction
  do i = 0, m+1
  do k = 0, l+1
  porosity(i,0,k)   = porosity(i,n,k)
  porosity(i,n+1,k) = porosity(i,1,k)
  end do
  end do
  
  ! ----------------
  return
  end subroutine  grid_conditions
  !******************
  
  !******************
  subroutine  initial_conditions (p, u, v, w, xp, yp, zp, width, height, depth  &
                                 , inlet_velocity, outlet_pressure, AoA, m, n, l)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(in)::width, height, depth, inlet_velocity, outlet_pressure, AoA
   real,intent(out),dimension(0:md,0:nd,0:ld)::u, v, w, p 
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:ld)::zp
   integer,intent(in)::m, n, l
  
  ! local variables
  integer::i, j, k
  real, parameter :: pai=atan(1.)*4.   
  
  ! ----------------
  do j = 1, n
  do i = 1, m
  do k = 1, l
   u(i,j,k) = inlet_velocity*cos(AoA/1300*pai)
   v(i,j,k) = inlet_velocity*sin(AoA/1300*pai)
   w(i,j,k) = 0.
   p(i,j,k) = outlet_pressure
  end do
  end do
  end do
  ! ----------------
  
  return
  end subroutine initial_conditions
  !******************
  
  ! output
  
  !******************
  subroutine  output_solution (p, u, v, w, m, n, l)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md,0:nd,0:ld)::u, v, w, p 
   integer,intent(in)::m, n, l
  
  ! local variables
  integer::i, j, k
  
  ! ----------------
  write(*,*)
  
  write(*,*)'velocity u '
  do k = 0, l+1
  do j = 0, n+1
  write(*,*) (u(i,j,k), i=0,m+1)
  end do
  end do
  write(*,*)
  
  write(*,*)'velocity v '
  do k = 0, l+1
  do j = 0, n+1
  write(*,*) (v(i,j,k), i=0,m+1)
  end do
  end do
  write(*,*)
  
  write(*,*)'velocity w '
  do k = 0, l+1
  do j = 0, n+1
  write(*,*) (w(i,j,k), i=0,m+1)
  end do
  end do
  write(*,*)
  
  write(*,*)'pressure'
  do k = 0, l+1
  do j = 0, n+1
  write(*,*) (p(i,j,k), i=0,m+1)
  end do
  end do
  write(*,*)
  ! ----------------
  
  return
  end subroutine output_solution
  !******************
  
  !******************
  subroutine  output_grid (xp, yp, zp, m, n, l)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:ld)::zp
   integer,intent(in)::m, n, l
  
  ! local variables
  integer::i, j, k
  
  open (63, file='etc/grid.dat', status='replace')
  ! ----------------
  write(63,*)'m, n, l =', m, n, l
  write(63,*)'grid points ='
  write(63,*) (xp(i), i=1,m)
  write(63,*) (yp(j), j=1,n)
  write(63,*) (zp(k), k=1,l)
  ! ----------------
  close (63)
  return
  end subroutine output_grid
  !******************
  
  !******************
  subroutine  output_grid_list (xp, yp, zp, m, n, l, angle_of_attack)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:ld)::zp
   integer,intent(in)::m, n, l
   real,intent(in):: angle_of_attack
  
  ! local variables
  integer::i, j, k
  real::pai=atan(1.)*4.
  real::x, y, z, th
  
  open (67, file='etc/cellcenter.dat', status='replace')
  ! ----------------
  th = angle_of_attack/1300.*pai
  do i=1,m
  do j=1,n
  do k=1,l
  x=xp(i)*cos(th)-yp(j)*sin(th)
  y=xp(i)*sin(th)+yp(j)*cos(th)
  z=zp(k)
  write(67,*) x,y,z
  end do
  end do
  end do
  ! ----------------
  close (67)
  return
  end subroutine output_grid_list
  !******************
  
  !******************
  subroutine  output_solution_post (p, u, v, w, xp, yp, zp, porosity, m, n, l)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md,0:nd,0:ld)::u, v, w, p
   real,intent(in),dimension(0:md,0:nd,0:ld)::porosity
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:ld)::zp
   integer,intent(in)::m, n, l
  
  ! local variables
   real, parameter::small=1.e-6, big=1.e6, zero=0.
   real, parameter::pmin=0.25, pmax=0.75
   integer::i, j, k
   real,dimension(0:md, 0:nd, 0:ld)::u_cnt, v_cnt, w_cnt, p_cnt
  
  open (61, file='etc/solution_uvp.dat', status='replace')
  
  ! ----------------
  ! interpolation at p-center grid
  
  do i = 1, m
  do j = 1, n
  do k = 1, l
   u_cnt(i,j,k)=u(i,j,k)*porosity(i,j,k)
   v_cnt(i,j,k)=v(i,j,k)*porosity(i,j,k)
   w_cnt(i,j,k)=w(i,j,k)*porosity(i,j,k)
   if (porosity(i,j,k) > small) then
    p_cnt(i,j,k)=p(i,j,k)
   else
    p_cnt(i,j,k)=zero
   end if 
  end do
  end do
  end do
  
  do j = 1, n
  do k = 1, l
   u_cnt(0,j,k)=u_cnt(1,j,k)
   v_cnt(0,j,k)=v_cnt(1,j,k)
   w_cnt(0,j,k)=w_cnt(1,j,k)
   p_cnt(0,j,k)=p_cnt(1,j,k)
   u_cnt(m+1,j,k)=u_cnt(m,j,k)
   v_cnt(m+1,j,k)=v_cnt(m,j,k)
   w_cnt(m+1,j,k)=w_cnt(m,j,k)
   p_cnt(m+1,j,k)=p_cnt(m,j,k)
  end do
  end do
  
  do i = 0, m+1
   u_cnt(i,0,k)=u_cnt(i,1,k)
   v_cnt(i,0,k)=v_cnt(i,1,k)
   w_cnt(i,0,k)=w_cnt(i,1,k)
   p_cnt(i,0,k)=p_cnt(i,1,k)
   u_cnt(i,n+1,k)=u_cnt(i,n,k)
   v_cnt(i,n+1,k)=v_cnt(i,n,k)
   w_cnt(i,n+1,k)=w_cnt(i,n,k)
   p_cnt(i,n+1,k)=p_cnt(i,n,k)
  end do
  
  !-----------------
  write(61,*)'m, n, l =', m, n, l
  
  write(61,*)'velocity u_bulk '
  do k = 1, l
  write(61,*) ((u_cnt(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'velocity v_bulk '
  do k = 1, l
  write(61,*) ((v_cnt(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'velocity w_bulk '
  do k = 1, l
  write(61,*) ((w_cnt(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'velocity u_inst '
  do k = 1, l
  write(61,*) ((u(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'velocity v_inst '
  do k = 1, l
  write(61,*) ((v(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'velocity w_inst '
  do k = 1, l
  write(61,*) ((w(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'pressure p_fluid'
  do k = 1, l
  write(61,*) ((p_cnt(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'pressure P_all'
  do k = 1, l
  write(61,*) ((p(i,j,k), i=1,m),j=1,n)
  end do
  
  write(61,*)'porosity'
  do k = 1, l
  write(61,*) ((porosity(i,j,k), i=1,m),j=1,n)
  end do
    
  close (61)
  ! ----------------
  
  ! ----------------
  ! surface profile
  open (62, file='etc/surface_profile.dat', status='replace')
  
  do k=1,l
  do j=1,n
  do i=1,m
  
  if( porosity(i,j,k) < pmax .and. porosity(i,j,k) > pmin )then
   write(62,*) xp(i), yp(j), zp(i), p_cnt(i,j,k), porosity(i,j,k)
  end if
  
  end do
  end do
  end do
  close (62)
  ! ----------------
  
  return
  end subroutine output_solution_post
  !******************
  
  !******************
  subroutine  output_paraview (p, u, v, w, porosity, xp, yp, zp, m, n, l)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:ld)::zp
   real,intent(in),dimension(0:md, 0:nd, 0:ld)::u, v, w, p
   real,intent(in),dimension(0:md,0:nd,0:ld)::porosity
   integer,intent(in)::m, n, l
   integer::i, j, k
  
  ! local variables
   real,dimension(0:md,0:nd,0:ld):: div
  
   character(len=50)::csv_file
   character(len=50)::output_folder
  
   namelist /directory_control/csv_file, output_folder
   open(11,file="config/controlDict.txt",status="old",action="read")
   read(11,nml=directory_control)
   close(11)
  
  open(50,file=trim(output_folder)//'/output_paraview.vtk',status="unknown",form="formatted",position="rewind")
  !open(*,file='solution.vtk',status="replace")
  ! ----------------
  
      write(50,"('# vtk DataFile Version 3.0')")
      write(50,"('3D flow')")
      write(50,"('ASCII ')")
      
      write(50,"('DATASET STRUCTURED_GRID')")
      write(50,"('DIMENSIONS ',3(1x,i4))") m, n, l
      
      write(50,"('POINTS ',i9,' float')") m*n*l
      do k=1,l
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") xp(i), yp(j), zp(k)
      enddo
      enddo
      enddo
      
      write(50,"('POINT_DATA ',i9)") m*n*l
      
  !! velocity vector
      write(50,"('VECTORS velocity float')")
      do k=1,l
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") u(i,j,k), v(i,j,k), w(i,j,k)
      enddo
      enddo
      enddo
  
  !! velocity vector
      write(50,"('VECTORS velocityInFluid float')")
      do k=1,l
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") u(i,j,k)*porosity(i,j,k), v(i,j,k)*porosity(i,j,k), w(i,j,k)*porosity(i,j,k)
      enddo
      enddo
      enddo
        
  !! pressure
      write(50,"('SCALARS pressure float')")
      write(50,"('LOOKUP_TABLE default')")
      do k=1,l
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") p(i,j,k)
      enddo
      enddo
      enddo
  
  do i=1,m
  do j=1,n
  do k=1,l
   div(i,j,k)= (u(i+1,j,k)-u(i-1,j,k))/(xp(i+1)-xp(i-1)) &
              +(v(i,j+1,k)-v(i,j-1,k))/(yp(j+1)-yp(j-1)) &
              +(v(i,j,k+1)-v(i,j,k-1))/(zp(k+1)-zp(k-1))
  end do
  end do
  end do
  
  !! divergent velocity
      write(50,"('SCALARS VelocityDivergent float')")
      write(50,"('LOOKUP_TABLE default')")
      do k=1,l
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") div(i,j,k)
      enddo
      enddo
      enddo
    
  !! porosity
      write(50,"('SCALARS porosity float')")
      write(50,"('LOOKUP_TABLE default')")
      do k=1,l
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") porosity(i,j,k)
      enddo
      enddo
      enddo
  
  ! ----------------
  close(50)
  
  return
  end subroutine  output_paraview
  !******************
  
  !******************
  subroutine  output_divergent (p, u, v, w, porosity, dx, dy, dz, m, n, l)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md,0:nd,0:ld)::u, v, w, p 
   real,intent(in),dimension(0:md,0:nd,0:ld)::porosity
   real,intent(in)::dx, dy, dz
   integer,intent(in)::m, n, l
  
  ! local variables
  integer::i, j, k
  real,dimension(0:md,0:nd,0:ld)::div
  
  open (62, file='etc/divergent.dat', status='replace')
  ! ----------------
  
  do i = 1, m
  do j = 1, n
  do k = 1, l
  div(i,j,k)= ((porosity(i+1,j,k)*u(i,j,k)+porosity(i,j,k)*u(i+1,j,k))/2      &
              -(porosity(i-1,j,k)*u(i,j,k)+porosity(i,j,k)*u(i-1,j,k))/2 )/dx &
             +((porosity(i,j+1,k)*v(i,j,k)+porosity(i,j,k)*v(i,j+1,k))/2      &
              -(porosity(i,j-1,k)*v(i,j,k)+porosity(i,j,k)*v(i,j-1,k))/2 )/dy &
             +((porosity(i,j,k+1)*w(i,j,k)+porosity(i,j,k)*w(i,j,k+1))/2      &
              -(porosity(i,j,k-1)*w(i,j,k)+porosity(i,j,k)*w(i,j,k-1))/2 )/dz
  end do
  end do
  end do
  
  write(62,*)
  write(62,*)'porosity'
  do k = 1, l
  do j = 1, n
  write(62,*) (porosity(i,j,k), i=1,m)
  end do
  end do
  
  write(62,*)
  write(62,*)'divergent velocity'
  do k = 1, l
  do j = 1, n
  write(62,*) (div(i,j,k), i=1,m)
  end do
  end do
  write(62,*)
  
  ! ----------------
  close (62)
  
  end subroutine  output_divergent
  !******************
  
  !******************
  subroutine  output_paraview_temp (p, u, v, w, porosity, xp, yp, zp, m, n, l, istep)
   implicit none
   integer,parameter::md=300, nd = 300, ld = 300     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   real,intent(in),dimension(0:ld)::zp
   real,intent(in),dimension(0:md, 0:nd, 0:ld)::u, v, w, p
   real,intent(in),dimension(0:md, 0:nd, 0:ld)::porosity
   integer,intent(in)::m, n, l, istep
  
  ! -- local variable
   real,dimension(0:md,0:nd,0:ld):: div
   integer::i, j, k
   character(5)::number
   character(len=50)::csv_file
   character(len=50)::output_folder
  ! -- open file
  
   namelist /directory_control/csv_file, output_folder
   open(11,file="config/controlDict.txt",status="old",action="read")
   read(11,nml=directory_control)
   close(11)
  
  write(number,"(I5.5)")istep
  
  open(65,file=trim(output_folder)//"/output_"//number//".vtk",status="unknown",form="formatted",position="rewind")
  !open(*,file='solution.vtk',status="replace")
  ! ----------------
  
  write(65,"('# vtk DataFile Version 3.0')")
  write(65,"('3D flow')")
  write(65,"('ASCII ')")
  
  write(65,"('DATASET STRUCTURED_GRID')")
  write(65,"('DIMENSIONS ',3(1x,i4))") m, n, l
  
  write(65,"('POINTS ',i9,' float')") m*n*l
  do k=1,l
  do j=1,n
  do i=1,m
    write(65,"(3(f16.4,1x))") xp(i), yp(j), zp(k)
  enddo
  enddo
  enddo
  
  write(65,"('POINT_DATA ',i9)") m*n*l
  
!! velocity vector
  write(65,"('VECTORS velocity float')")
  do k=1,l
  do j=1,n
  do i=1,m
    write(65,"(3(f16.4,1x))") u(i,j,k), v(i,j,k), w(i,j,k)
  enddo
  enddo
  enddo

!! velocity vector
  write(65,"('VECTORS velocityInFluid float')")
  do k=1,l
  do j=1,n
  do i=1,m
    write(65,"(3(f16.4,1x))") u(i,j,k)*porosity(i,j,k), v(i,j,k)*porosity(i,j,k), w(i,j,k)*porosity(i,j,k)
  enddo
  enddo
  enddo
    
!! pressure
  write(65,"('SCALARS pressure float')")
  write(65,"('LOOKUP_TABLE default')")
  do k=1,l
  do j=1,n
  do i=1,m
    write(65,"(3(f16.4,1x))") p(i,j,k)
  enddo
  enddo
  enddo

  do i=1,m
  do j=1,n
  do k=1,l
  div(i,j,k)= (u(i+1,j,k)-u(i-1,j,k))/(xp(i+1)-xp(i-1)) &
            +(v(i,j+1,k)-v(i,j-1,k))/(yp(j+1)-yp(j-1)) &
            +(v(i,j,k+1)-v(i,j,k-1))/(zp(k+1)-zp(k-1))
  end do
  end do
  end do

!! divergent velocity
  write(65,"('SCALARS VelocityDivergent float')")
  write(65,"('LOOKUP_TABLE default')")
  do k=1,l
  do j=1,n
  do i=1,m
    write(65,"(3(f16.4,1x))") div(i,j,k)
  enddo
  enddo
  enddo

!! porosity
  write(65,"('SCALARS porosity float')")
  write(65,"('LOOKUP_TABLE default')")
  do k=1,l
  do j=1,n
  do i=1,m
    write(65,"(3(f16.4,1x))") porosity(i,j,k)
  enddo
  enddo
  enddo

  ! ----------------
  close(65)
  
  return
  end subroutine  output_paraview_temp
  !******************
