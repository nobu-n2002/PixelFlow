program main
  !$ use omp_lib
   implicit none
   integer,parameter:: md=4500, nd=4500     ! md, nd > grid size (m,n)
   real:: dx, dy, dt
   real:: xnue, xlamda, density, width, height, time, inlet_velocity, outlet_pressure, AoA, thickness
   real,dimension(0:md,0:nd):: u, v, p, u_old, v_old, p_old
   real,dimension(0:md,0:nd):: porosity
   real,dimension(0:md):: xp
   real,dimension(0:nd):: yp
   real,dimension(0:md,0:nd):: ap, ae, aw, an, as, bb
   real::residual_p, residual_u, residual_v
   integer:: m, n, istep, istep_max, iset, istep_out
   integer:: i, j
   real:: tp, tu, tv, tbc, ttmpout
   character(len=50) :: output_folder
   character(len=50) :: csv_file
  ! ----------------
  ! read input data by using namelist
  ! by Nobuto Nakamichi 4/7/2023
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
  
    call  physical_conditions (xnue, xlamda, density, width, height, time &
                          , inlet_velocity, outlet_pressure, AoA, m, n)
    call  grid_conditions (xp, yp, dx, dy, dt, xnue, xlamda, density, width, height, thickness, time &
                          , inlet_velocity, AoA, porosity, m, n, istep_max, iset)
    ! call  output_grid_list (xp, yp, m, n, angle_of_attack)
    stop
  end if
  
  ! ----------------
  ! calculation start  (if iest=!0)
  ! ----------------
  ! set up condistions
  m=0         ! setup switch for grid conditions
  density=0.  ! setup switch for physical conditions
  
  call  physical_conditions (xnue, xlamda, density, width, height, time &
                         , inlet_velocity, outlet_pressure, AoA, m, n)
  call  grid_conditions (xp, yp, dx, dy, dt, xnue, xlamda, density, width, height,&
               thickness, time, inlet_velocity, AoA, porosity, m, n, istep_max, iset)
  call  output_grid (xp, yp, m, n)
  
  ! write(*,*) "check", (porosity(i,10), i=1,m)
  
  istep = 0
  time = istep * dt
  
  ! ----------------
  
  write(*,*) 'istep_max= ', istep_max,'   istep_out= ', istep_out
  
  call  initial_conditions (p, u, v, xp, yp, width, height &
                         , inlet_velocity, outlet_pressure, AoA, m, n)
  call  boundary (p, u, v, xp, yp, width, height            &
                         , inlet_velocity, outlet_pressure, AoA, porosity, m, n)
  
  ! print initial conditions
  ! call  output_solution (p, u, v, m, n)
  
  ! ----------------
  ! MAC algorithm start
  
  do istep = 1, istep_max
  
    time = istep* dt
    residual_p = 0.0
    residual_u = 0.0
    residual_v = 0.0
    do i = 0, m+1
      do j = 0, n+1
        residual_p = max(residual_p, porosity(i,j) * abs(p(i,j) - p_old(i,j)))
        residual_u = max(residual_u, abs(u(i,j) - u_old(i,j)))
        residual_v = max(residual_v, abs(v(i,j) - v_old(i,j)))
        u_old(i,j) = u(i,j)
        v_old(i,j) = v(i,j)
        p_old(i,j) = p(i,j)
      end do
    end do
    write(*,*)'--- time_steps= ',istep, ' --  time = ',time, 'residual= p:', residual_p, 'u:', residual_u, 'v:', residual_v
  
    call  solve_p (p, u, v, u_old, v_old, porosity, xnue, xlamda, density, height, thickness, yp, dx, dy, dt, m, n)
  
    call  solve_u (p, u, v, u_old, v_old, porosity, xnue, xlamda, density, dx, dy, dt, m, n)
  
    call  solve_v (p, u, v, u_old, v_old, porosity, xnue, xlamda, density, dx, dy, dt, m, n)
  
    call  boundary(p, u, v, xp, yp, width, height    &
                        , inlet_velocity, outlet_pressure, AoA, porosity, m, n)
  
    if(mod(istep,istep_out)==0) call  output_paraview_temp (p, u, v, porosity, xp, yp, m, n, istep)
  
  
  end do
  ! MAC algorithm end
  ! ----------------
  
  ! print conditions (recall)
  call  physical_conditions (xnue, xlamda, density, width, height, time &
                         , inlet_velocity, outlet_pressure, AoA, m, n)
  call  grid_conditions (xp, yp, dx, dy, dt, xnue, xlamda, density, width, height, thickness, time &
  , inlet_velocity, AoA, porosity, m, n, istep_max, iset)
  
  ! print solutions
  call  output_solution_post (p, u, v, xp, yp, porosity, m, n)
  call  output_divergent (p, u, v, porosity, dx, dy, m, n)
  call  output_paraview (p, u, v, porosity, xp, yp, m, n)
  write(*,*) 'program finished'
  
  end program main
  !******************
  
  !  solve variables
  
  !******************
  subroutine  solve_p (p, u, v, u_old, v_old, porosity, xnue, xlamda, density, height,thickness, yp, dx, dy, dt, m, n)
   implicit none
   integer,parameter:: md = 4500, nd = 4500     ! md, nd > grid size (m,n)
   real,intent(in):: dx, dy, dt
   real,intent(in):: xnue, xlamda, density, height, thickness
   real,intent(inout),dimension(0:md,0:nd):: u, v, p, u_old, v_old
   real,intent(in),dimension(0:md,0:nd):: porosity
   real,intent(in),dimension(0:nd):: yp
   integer,intent(in):: m, n
  
  !-----------------
  ! local variables
   real, parameter:: small = 1.e-6, big = 1.e6, zero = 0.
   real, parameter::alpha = 32.0

   real::u_stg, v_stg
   real,dimension(0:md,0:nd):: ap, ae, aw, an, as, bb, div
   integer:: i, j
   real:: fc, poro_grad
  
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
          div(i,j)= (u_old(i+1,j)-u_old(i-1,j))/dx*.5 &
                  + (v_old(i,j+1)-v_old(i,j-1))/dx*.5
      end do
  end do
  
  do j = 1, n
      div(0,j)  = 0.  ! inlet
      div(m+1,j)= 0.  ! outlet
  end do
  
  do i = 1, m
      div(i,0)  = div(i,n)  ! periodic condition
      div(i,n+1)= div(i,1)
  end do
  
  ! ----------------
  fc=0.
  
  do i = 1, m
  do j = 1, n
  ! ----------------
  !   velocity u
  ! ----------------
  ! convection_x  (1st upwind scheme)
  u(i,j)=u_old(i,j)-dt*(                  &
         fc*( max(u_old(i,j),0.)*(u_old(i,j)-u_old(i-1,j))/dx	  &      ! u>0 1st upwind scheme
             +min(u_old(i,j),0.)*(u_old(i+1,j)-u_old(i,j))/dx ) &      ! u<0 1st upwind scheme
   +(1.-fc)* u_old(i,j)*(u_old(i+1,j)-u_old(i-1,j))/dx/2.       )    ! 2nd central scheme
  
  ! convection_y
  u(i,j)=u(i,j)-dt*(						&
         fc*(max(v_old(i,j),0.)*(u_old(i,j)-u_old(i,j-1))/dy	&   ! v>0 1st upwind scheme
            +min(v_old(i,j),0.)*(u_old(i,j+1)-u_old(i,j))/dy) &   ! v<0 1st upwind scheme
   +(1.-fc)* v_old(i,j)*(u_old(i,j+1)-u_old(i,j-1))/dy/2.     ) ! 2nd central scheme
  
  ! diffusion_x
  u(i,j)=u(i,j) +dt*xnue*(u_old(i+1,j)-2.*u_old(i,j)+u_old(i-1,j))/dx/dx
  !      +dt*xnue/(small+porosity(i,j))*(u_old(i+1,j)-u_old(i-1,j))*(porosity(i+1,j)-porosity(i-1,j))/dx/dx*0.25 ! non-conseved term
  ! diffusion_y
  u(i,j)=u(i,j) +dt*xnue*(u_old(i,j+1)-2.*u_old(i,j)+u_old(i,j-1))/dy/dy
  !      +dt*xnue/(small+porosity(i,j))*(u_old(i,j+1)-u_old(i,j-1))*(porosity(i,j+1)-porosity(i,j-1))/dy/dy*0.25 ! non-conseved term
  ! divergence term
  u(i,j)=u(i,j) +dt*(xnue + xlamda)*(div(i+1,j)-div(i-1,j))/dx*.5
  ! additional terms by porosity profile   ! canceled for non-slip condition
  u(i,j)=u(i,j)							&
        +dt*( ( (u_old(i+1,j)-u_old(i-1,j))/dx*.5+(u_old(i+1,j)-u_old(i-1,j))/dx*.5) &
                *xnue*(porosity(i+1,j)-porosity(i-1,j))/dx*.5                        &
             +( (u_old(i,j+1)-u_old(i,j-1))/dy*.5+(v_old(i+1,j)-v_old(i-1,j))/dx*.5) &
                *xnue*(porosity(i,j+1)-porosity(i,j-1))/dy*.5                        &
              + div(i,j)*(porosity(i+1,j)-porosity(i-1,j))/dx*0.5*xlamda             &
             )/porosity(i,j)
  ! force on wall
  if (nonslip) then
    u(i,j)=u(i,j)- dt*xnue*u_old(i,j)/(thickness*dx)**2 *alpha*porosity(i,j)*(1.-porosity(i,j))*(1.-porosity(i,j))
  end if
  ! ----------------
  !   velocity v
  ! ----------------
  ! convection_x  (1st upwind scheme)
  v(i,j)=v_old(i,j)-dt*(						&
        fc *(max(u_old(i,j),0.)*(v_old(i,j)-v_old(i-1,j))/dx	&  ! u>0 1st upwind scheme
            +min(u_old(i,j),0.)*(v_old(i+1,j)-v_old(i,j))/dx) &  ! u<0 1st upwind scheme
   +(1.-fc)* u_old(i,j)*(v_old(i+1,j)-v_old(i-1,j))/dx/2.     ) ! 2nd central scheme
  
  ! convection_y
  v(i,j)=v(i,j)-dt*(							&
        fc *(max(v_old(i,j),0.)*(v_old(i,j)-v_old(i,j-1))/dy	&  ! v>0
            +min(v_old(i,j),0.)*(v_old(i,j+1)-v_old(i,j))/dy) &  ! v<0
   +(1.-fc)* v_old(i,j)*(v_old(i,j+1)-v_old(i,j-1))/dy/2.   ) ! 2nd central scheme
  
  ! diffusion_x
  v(i,j)=v(i,j) +dt*xnue*(v_old(i+1,j)-2.*v_old(i,j)+v_old(i-1,j))/dx/dx
  !      +dt*xnue/(small+porosity(i,j))*(v_old(i+1,j)-v_old(i-1,j))*(porosity(i+1,j)-porosity(i-1,j))/dx/dx*0.25 ! non-conseved term
  ! diffusion_y
  v(i,j)=v(i,j) +dt*xnue*(v_old(i,j+1)-2.*v_old(i,j)+v_old(i,j-1))/dy/dy
  !      +dt*xnue/(small+porosity(i,j))*(v_old(i,j+1)-v_old(i,j-1))*(porosity(i,j+1)-porosity(i,j-1))/dy/dy*0.25 ! non-conseved term
  ! divergence term   ! L+(2/3)N = (1/3)N;(2/3) or 0(1/3)
  v(i,j)=v(i,j) +dt*(xnue + xlamda)*(div(i,j+1)-div(i,j-1))/dy*.5
  ! additional terms by porosity profile   ! canceled for non-slip condition    ! L+(2/3)N = (1/3)N;(-1/3) or 0:(-2/3)N
  v(i,j)=v(i,j)							&
        +dt*( ( (v_old(i+1,j)-v_old(i-1,j))/dx*.5+(u_old(i,j+1)-u_old(i,j-1))/dy*.5) &
                *xnue*(porosity(i+1,j)-porosity(i-1,j))/dx*.5                        &
             +( (v_old(i,j+1)-v_old(i,j-1))/dy*.5+(v_old(i,j+1)-v_old(i,j-1))/dy*.5) &
                *xnue*(porosity(i,j+1)-porosity(i,j-1))/dy*.5                        &
             + div(i,j)*(porosity(i,j+1)-porosity(i,j-1))/dy*0.5*xlamda       &
             )/porosity(i,j)
  ! force on wall
  if (nonslip) then
    v(i,j)=v(i,j)- dt*xnue*v_old(i,j)/(thickness*dx)**2*alpha*porosity(i,j)*(1.-porosity(i,j))*(1.-porosity(i,j))
  end if
  end do
  end do
  
  ! ----------------
  ! matrix solution  !  formulation of porous media
  
  do i = 1, m
  do j = 1, n
  ae(i,j)= dt*max(small,(porosity(i+1,j)+porosity(i,j))*0.5)/dx/dx
  aw(i,j)= dt*max(small,(porosity(i,j)+porosity(i-1,j))*0.5)/dx/dx
  an(i,j)= dt*max(small,(porosity(i,j+1)+porosity(i,j))*0.5)/dy/dy
  as(i,j)= dt*max(small,(porosity(i,j)+porosity(i,j-1))*0.5)/dy/dy
  ap(i,j)= -ae(i,j)-aw(i,j)-an(i,j)-as(i,j)
  
  bb(i,j)= ((porosity(i+1,j)*u(i,j)+porosity(i,j)*u(i+1,j))*0.5             &
           -(porosity(i-1,j)*u(i,j)+porosity(i,j)*u(i-1,j))*0.5)*density/dx &
          +((porosity(i,j+1)*v(i,j)+porosity(i,j)*v(i,j+1))*0.5             &
           -(porosity(i,j-1)*v(i,j)+porosity(i,j)*v(i,j-1))*0.5)*density/dy
  
  !if (porosity(i,j) <small) then   !in solid (dummy solution)
  ! ap(i,j)=-1.
  ! bb(i,j)= 0.
  ! ae(i,j)= 0.25
  ! aw(i,j)= 0.25
  ! an(i,j)= 0.25
  ! as(i,j)= 0.25
  !end if
  
  end do
  end do
  
  call boundrary_matrix (p, ap, ae, aw, an, as, bb, m, n, height, yp)
  
  ! call solve_matrix (p, ap, ae, aw, an, as, bb, m, n)
  call solve_matrix_vec_omp (p, ap, ae, aw, an, as, bb, m, n)
  !call solve_matrix_vec_oacc (p, ap, ae, aw, an, as, bb, m, n)
  ! ----------------
  ! ----------------
  return
  end subroutine solve_p
  !******************
  
  !******************
  ! OpenACC Parallelized
  ! Written only for GPU machine
  ! No efficiency ensured on CPU machine
  subroutine  solve_matrix_vec_oacc (p, ap, ae, aw, an, as, bb, m, n)
   implicit none
   integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(inout),dimension(0:md,0:nd)::	p
   real,intent(in),dimension(0:md,0:nd)::		ap, ae, aw, an, as, bb
   integer,intent(in)::	m, n
  
  ! local variables
  real::		relux_factor, error
  real,dimension(0:md, 0:nd)::	p_old
  integer::	i, j, iter, iter_max, k
  
  ! ----------------
  !   SOR algorithm
  
  iter_max = min(100, max(m, n))		! SOR max interation steps
  relux_factor = 1.7 	! SOR reluxation factor
  
  !$acc data copy(p_old, p, error) &
  !$acc    & copyin(ap, ae, aw, an, as, bb, relux_factor)
  
  do iter = 1, iter_max
    ! write(*,*)'CHECK iteration no. ',iter, ' / iter_max', iter_max
    error=0.
  
    ! default periodic condition in y-direction
    !$acc kernels
    !$acc loop independent
    do i = 1, m
      p(i, 0) = p(i, n)
      p(i, n+1) = p(i, 1)
    end do
    !$acc end kernels
  
    !$acc kernels
    !$acc loop independent
    do i = 0, m+1
    !$acc loop independent
      do j = 0, n+1
        p_old(i, j) = p(i, j)
      end do
    end do
    !$acc end kernels
  
    !-- EVEN SPACE process
    !$acc kernels
    !$acc loop reduction(max:error)
    do k = 2, m*n, 2    ! even space
      j = (k - 1) / m + 1
      i = k - (j - 1) * m
  
      !-- IF m is EVEN (Based on Column-Major Order; FORTRAN)
      if(mod(m,2)==0 .and. mod(j,2)==0) i = i - 1
  
      p(i, j) = ( bb(i, j)					                                      &
                 - ae(i, j) * p_old(i+1, j) - aw(i, j) * p_old(i-1, j)    &
                 - an(i, j) * p_old(i, j+1) - as(i, j) * p_old(i, j-1) )  &
                / ap(i, j) * relux_factor                                 &
               + p_old(i, j) * (1. - relux_factor)
      error = max(error, abs(p(i, j) - p_old(i, j)))
    end do
    !$acc end kernels
  
    ! default periodic condition in y-direction
    !$acc kernels
    !$acc loop independent
    do i = 1, m
      p(i, 0)  = p(i, n)
      p(i, n+1) = p(i, 1)
    end do
    !$acc end kernels
  
    !$acc kernels
    !$acc loop independent
    do i = 0, m+1
    !$acc loop independent
      do j = 0, n+1
        p_old(i, j) = p(i, j)
      end do
    end do
    !$acc end kernels
  
    !-- ODD SPACE process
    !$acc kernels
    !$acc loop reduction(max:error)
    do k = 1, m*n, 2    ! odd space
      j = (k - 1) / m + 1
      i = k - (j - 1) * m
  
      !-- IF m is EVEN (Based on Column-Major Order; FORTRAN)
      if(mod(m,2)==0 .and. mod(j,2)==0) i = i + 1
  
      p(i, j) = ( bb(i, j)					                                      &
                 - ae(i, j) * p_old(i+1, j) - aw(i, j) * p_old(i-1, j)    &
                 - an(i, j) * p_old(i, j+1) - as(i, j) * p_old(i, j-1) )  &
                / ap(i, j) * relux_factor                                 &
               + p_old(i, j) * (1. - relux_factor)
      error = max(error, abs(p(i, j) - p_old(i, j)))
    end do
    !$acc end kernels
  
    !if(mod(iter,10)==0) write(*,*)'CHECK iteration no.', iter,'  -- error=', error
  end do
  
  ! default periodic condition in y-direction
  !$acc kernels
  !$acc loop independent
  do i = 1, m
    p(i, 0)   = p(i, n)
    p(i, n+1) = p(i, 1)
  end do
  !$acc end kernels
  
  !$acc end data
  
   write(*,*)'SOR iteration no.', iter-1,'  -- error=', error
  ! write(*,*)' check P(10,10) in solve_matrix ', p(10,10)
  ! ----------------
  
  return
  end subroutine solve_matrix_vec_oacc
  !******************
  
  !******************
  ! OpenMP Parallelized
  ! Written only for CPU machine
  ! No efficiency ensured on GPU machine
  subroutine  solve_matrix_vec_omp (p, ap, ae, aw, an, as, bb, m, n)
    implicit none
    integer,parameter::md=4500, nd=4500     ! md, nd > grid size (m,n)
    real,intent(inout),dimension(0:md,0:nd)::p
    real,intent(in),dimension(0:md,0:nd)::ap, ae, aw, an, as, bb
    integer,intent(in)::m, n
  
   ! local variables
   real::	relux_factor, error
   real,dimension(0:md, 0:nd)::p_old
   integer::	i, j, iter, iter_max, k
  
   !$omp parallel private(iter, i, j, k) &
   !$omp & shared(iter_max, relux_factor, m, n) &
   !$omp & shared(p_old, p, ap, ae, aw, an, as, bb) &
   !$omp & default(none)
  
   ! ----------------
   !   SOR algorithm
   ! ----------------
   !$omp single
   iter_max = min(300, max(m, n))		! SOR max interation steps
   relux_factor = 1.7 	! SOR reluxation factor
   !$omp end single
  
   do iter = 1, iter_max
     ! write(*,*)'CHECK iteration no. ',iter, ' / iter_max', iter_max
  
     ! default periodic condition in y-direction
     !$omp do
     do i = 1, m
       p(i, 0) = p(i, n)
       p(i, n+1) = p(i, 1)
     end do
     !$omp end do
  
     !$omp do
     do i = 0, m+1
       do j = 0, n+1
         p_old(i, j) = p(i, j)
       end do
     end do
     !$omp end do
  
     !-- EVEN SPACE process
     !$omp do
     do k = 2, m*n, 2    ! even space
       j = (k - 1) / m + 1
       i = k - (j - 1) * m
  
       !-- IF m is EVEN (Based on Column-Major Order; FORTRAN)
       if(mod(m,2)==0 .and. mod(j,2)==0) i = i - 1
  
       p(i, j) = ( bb(i, j)					                                       &
                  - ae(i, j) * p_old(i+1, j) - aw(i, j) * p_old(i-1, j)    &
                  - an(i, j) * p_old(i, j+1) - as(i, j) * p_old(i, j-1) )  &
                 / ap(i, j) * relux_factor                                 &
                + p_old(i, j) * (1. - relux_factor)
     end do
    !$omp end do
  
     ! default periodic condition in y-direction
     !$omp do
     do i = 1, m
       p(i, 0)  = p(i, n)
       p(i, n+1) = p(i, 1)
     end do
     !$omp end do
  
     !$omp do
     do i = 0, m+1
       do j = 0, n+1
         p_old(i, j) = p(i, j)
       end do
     end do
     !$omp end do
  
     !-- ODD SPACE process
    !$omp do
     do k = 1, m*n, 2    ! odd space
       j = (k - 1) / m + 1
       i = k - (j - 1) * m
  
       !-- IF m is EVEN (Based on Column-Major Order; FORTRAN)
       if(mod(m,2)==0 .and. mod(j,2)==0) i = i + 1
  
       p(i, j) = ( bb(i, j)					                                       &
                  - ae(i, j) * p_old(i+1, j) - aw(i, j) * p_old(i-1, j)    &
                  - an(i, j) * p_old(i, j+1) - as(i, j) * p_old(i, j-1) )  &
                 / ap(i, j) * relux_factor                                 &
                + p_old(i, j) * (1. - relux_factor)
     end do
    !$omp end do
  
     !if(mod(iter,10)==0) write(*,*)'CHECK iteration no.', iter,'  -- error=', error
   end do
  
   ! default periodic condition in y-direction
   !$omp do
   do i = 1, m
     p(i, 0)   = p(i, n)
     p(i, n+1) = p(i, 1)
   end do
   !$omp end do
  
   !$omp master
    ! write(*,*)'SOR iteration no.', iter-1,'  -- error=', error
   ! write(*,*)' check P(10,10) in solve_matrix ', p(10,10)
   !$omp end master
   ! ----------------
  
   !$omp end parallel
  
   return
   end subroutine solve_matrix_vec_omp
  !******************
  
  !******************
  subroutine  solve_matrix (p, ap, ae, aw, an, as, bb, m, n)
   implicit none
   integer,parameter:: md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(inout),dimension(0:md,0:nd):: p
   real,intent(in),dimension(0:md,0:nd):: ap, ae, aw, an, as, bb
   integer,intent(in):: m, n
  
  ! local variables
  real:: relux_factor, error
  real,dimension(0:md, 0:nd)::p_old
  integer::i, j, iter, iter_max
  
  ! ----------------
  !   SOR algorithm
  ! ----------------
  ! iter_max = min(100,max(m,n)) ! SOR max interation steps
  iter_max = 50
  relux_factor=1.7 ! SOR reluxation factor
  
  do iter = 1, iter_max
  ! write(*,*)'CHECK iteration no.'
  ! error=0.
  
  ! default periodic condition in y-direction
  do i = 1, m
   p(i,0)  =p(i,n)
   p(i,n+1)=p(i,1)
  end do
  
  ! do i = 0, m+1
  ! do j = 0, n+1
  !  p_old(i,j) = p(i,j)
  ! end do
  ! end do
  
  do i = 1, m
  do j = 1, n
  p(i,j) = (  bb(i,j)					&
            - ae(i,j)*p_old(i+1,j) -aw(i,j)*p(i-1,j)	&
            - an(i,j)*p_old(i,j+1) -as(i,j)*p(i,j-1) )	&
           /ap(i,j)    * relux_factor			&
          + p_old(i,j) * (1.-relux_factor)
  end do
  end do
  
  end do
  
  ! ----------------
  
  return
  end subroutine solve_matrix
  !******************
  
  !******************
  subroutine  boundrary_matrix (p, ap, ae, aw, an, as, bb, m, n, height, yp)
   implicit none
   integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(in)::	height
   real,intent(in),dimension(0:md,0:nd)::		p
   real,intent(inout),dimension(0:md,0:nd)::	ap, ae, aw, an, as, bb
   real,intent(in),dimension(0:nd)::		yp
   integer,intent(in)::	m, n
  
  ! local variables
  integer	i, j
  
  ! ----------------
  ! inlet (dp/x=0 at i=1)
  do j= 1, n
   ae(1,j) =ae(1,j)+aw(1,j)
   aw(1,j) =0.
  end do
  
  ! outlet (p=outlet_pressure at i=m)
  do j= 1, n
   bb(m,j) =bb(m,j)+ae(m,j)*p(m+1,j)
   ae(m,j) = 0.
   aw(m,j) = 0.
   an(m,j) = 0.
   as(m,j) = 0.
  end do
  
  ! default : periodic condition in matrix solver
  
  ! symmetry or wall (dp/dy=0. at j=1)   xp>0
  !do i= 1,m
  ! an(i,1) =an(i,1)+as(i,1)
  ! as(i,1) = 0.
  !end do
  
  ! symmetry or wall  (dp/dy=0. at j=n)  xp>0
  !do i= 1,m
  ! as(i,n) =as(i,n)+an(i,n)
  ! an(i,n) = 0.
  !end do
  ! ----------------
  
  return
  end subroutine  boundrary_matrix
  !******************
  
  !******************
  subroutine  solve_u (p, u, v, u_old, v_old, porosity, xnue, xlamda, density, dx, dy, dt, m, n)
   implicit none
   integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(in)::	dx, dy, dt
   real,intent(in)::	xnue, xlamda, density
   real,intent(inout),dimension(0:md,0:nd)::	u, v, p, u_old, v_old
   real,intent(in),dimension(0:md,0:nd)::	porosity
   integer,intent(in)::	m, n
  
  ! local variables
  integer::	i, j
  
  ! ----------------
  do i = 1, m
  do j = 1, n
  ! convection_x  (1st upwind scheme)
  ! (already calculated in solve_p)
  
  ! convection_y
  ! (already calculated in solve_p)
  
  ! diffusion_x
  ! (already calculated in solve_p)
  
  ! diffusion_y
  ! (already calculated in solve_p)
  
  ! pressure
  u(i,j)=u(i,j) -dt/density*(p(i+1,j)-p(i-1,j))/dx*0.5
  
  end do
  end do
  
  ! ----------------
  return
  end subroutine solve_u
  !******************
  
  !******************
  subroutine  solve_v (p, u, v, u_old, v_old, porosity, xnue, xlamda, density, dx, dy, dt, m, n)
   implicit none
   integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(in)::	dx, dy, dt
   real,intent(in)::	xnue, xlamda, density
   real,intent(inout),dimension(0:md,0:nd)::	u, v, p, u_old, v_old
   real,intent(in),dimension(0:md,0:nd)::	porosity
   integer,intent(in)::	m, n
  
  ! local variables
  integer::	i, j
  
  ! ----------------
  do i = 1, m
  do j = 1, n
  ! convection_x  (1st upwind scheme)
  ! (already calculated in solve_p)
  
  ! convection_y
  ! (already calculated in solve_p)
  
  ! diffusion_x
  ! (already calculated in solve_p)
  
  ! diffusion_y
  ! (already calculated in solve_p)
  
  ! pressure
  v(i,j)=v(i,j) -dt/density*(p(i,j+1)-p(i,j-1))/dy*.5
  
  end do
  end do
  ! ----------------
  return
  end subroutine solve_v
  !******************
  
  !  conditions
  
  !******************
  subroutine  boundary(p, u, v, xp, yp, width, height    &
                       , inlet_velocity, outlet_pressure, AoA, porosity, m, n)
   implicit none
   integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(in)::	width, height, inlet_velocity, outlet_pressure, AoA
   real,intent(inout),dimension(0:md,0:nd)::	u, v, p
   real,intent(in),dimension(0:md,0:nd)::	porosity
   real,intent(in),dimension(0:md)::		xp
   real,intent(in),dimension(0:nd)::		yp
   integer,intent(in)::	m, n
  
  ! local variables
   real, parameter::     small=1.e-6, big=1.e6, zero=0., pai=atan(1.)*4.
   integer	i, j
  
  ! ----------------
  ! inlet (u=inlet_velocity, v=0., dp/dx=0 at i=1)
  do j= 1, n
   u(1,j) =inlet_velocity*cos(AoA/180.*pai)
   v(1,j) =inlet_velocity*sin(AoA/180.*pai)
   u(0,j) =u(1,j)		! dummy
   v(0,j) =v(1,j)  	! dummy
   p(0,j) =p(2,j)
  end do
  
  ! outlet (du/dx=0., dv/dx=0., p=outlet_pressure at i=m)
  do j= 1, n
   u(m+1,j) =u(m-1,j)
   v(m+1,j) =v(m-1,j)
  ! p(m,j) =outlet_pressure
   p(m+1,j)=outlet_pressure   ! dummy
  end do
  
  ! default: periodic condition (y-direction at j=1 & n)
  do i= 0, m+1
   u(i,0)   = u(i,n)
   v(i,0)   = v(i,n)
   p(i,0)   = p(i,n)
   u(i,n+1) = u(i,1)
   v(i,n+1) = v(i,1)
   p(i,n+1) = p(i,1)
  end do
  
  ! option: lower wall (u=0., v=0., dp/dy=0. at j=1)
  !do i= 0, m+1
  ! u(i,1) =0.
  ! v(i,1) =0.
  ! u(i,0) =0.					! dummy
  ! v(i,0) = -v(i,2)		  	! dummy
  ! p(i,0) =p(i,2)
  !end do
  
  ! option: symmetry (du/dy=0., v=0., dp/dy=0. at j=1)  xp>0
  !do i= 1, m
  ! u(i,0) = u(i,2)
  ! v(i,1) =0.
  ! v(i,0) = -v(i,2)		  	! dummy
  ! p(i,0) =p(i,2)
  !end do
  
  ! option: symmetry  (du/dy=0., v=0., dp/dy=0. at j=n)   xp>0
  !do i= 1, m
  ! u(i,n+1) = u(i,n-1)
  ! v(i,n) =0.
  ! v(i,n+1) = -v(i,n-1)		! dummy
  ! p(i,n+1) =p(i,n-1)
  !end do
  ! ----------------
  !do i= 0, m+1
  !do j= 0, n+1
  !if (porosity(i,j) <small) then   !in solid (dummy solution)
  ! u(i,j) = 0.
  ! v(i,j) = 0.
  !end if
  !end do
  !end do
  
  return
  end subroutine boundary
  !*****************************
  
  !*****************************
  subroutine physical_conditions(xnue, xlamda, density, width, height, time &
                                , inlet_velocity, outlet_pressure, AoA, m, n)
   implicit none
   integer,parameter:: md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(inout):: xnue, xlamda, density, width, height, time  &
                         ,inlet_velocity, outlet_pressure, AoA
   integer,intent(in):: m, n
  ! local variables
   real:: reynolds_no, wing_length
   integer:: i, j
   real::depth
  
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
  
  !wing_length =1.    ! (m)
  !reynolds_no=wing_length*inlet_velocity/xnue
  
  write(*,*)
  write(*,*) 'xnue ='	, xnue
  write(*,*) 'xlamda ='	, xlamda
  write(*,*) 'density ='	, density
  write(*,*) 'width ='	, width
  write(*,*) 'height ='	, height
  write(*,*) 'time ='	, time
  write(*,*) 'inlet_velocity ='	, inlet_velocity
  write(*,*) 'outlet_pressure ='	, outlet_pressure
  write(*,*) 'Angle of inlet_velocity (AoA) ='	, AoA
  !write(*,*) 'reynolds_no='	, reynolds_no
  
  ! ----------------
  
  return
  end subroutine physical_conditions
  !******************
  
  !******************
  subroutine  grid_conditions (xp, yp, dx, dy, dt, xnue, xlamda, density, width, height, thickness, time &
                              , inlet_velocity, AoA, porosity, m, n, istep_max, iset)
  implicit none
  integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
  real,intent(inout)::	dx, dy, dt, AoA, thickness
  real,intent(in)::	xnue, xlamda, density, width, height, time, inlet_velocity
  real,intent(inout),dimension(0:md,0:nd):: porosity
  real,intent(inout),dimension(0:md):: xp
  real,intent(inout),dimension(0:nd):: yp
  integer,intent(inout):: m, n, istep_max, iset
  character(len = 50) :: csv_file
  character(len = 50) :: output_folder
  
  ! local variables
  !real,dimension(0:md,0:nd)::	distance
  integer:: x, y, z
  real:: val, threshold
  real:: cfl_no, pecret_no, diffusion_factor, reynolds_no
  real:: pai, distance, center_x, center_y, radius
  integer::	i, j
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
  
  read(52,*) m,n
  
  do i=1,m*n
    read(52,*) x, y, z, val
    porosity(x,y) = max(val, threshold)
  end do
  close(52)
  
  ! thickness = 2.5
  dx = width / real(m-1)
  dy = height / real(n-1)
  dt = time / real(istep_max)
  
  radius = 0.25 * height
  
  cfl_no           = inlet_velocity * dt / dx
  pecret_no        = inlet_velocity * dx / xnue
  diffusion_factor = xnue * dt / dy / dy
  reynolds_no      = radius * inlet_velocity / xnue
  
  !----- check print out
  write(*,*)
  write(*,*) 'm, n =', m, n
  write(*,*) 'istep_max =', istep_max
  write(*,*) 'dx, dy =', dx, dy
  write(*,*) 'dt =', dt
  write(*,*) 'cfl_no =', cfl_no
  write(*,*) 'pecret_no =', pecret_no
  write(*,*) 'diffusion_factor =', diffusion_factor
  write(*,*) 'reynolds_no=' , reynolds_no
  write(*,*) 'thickness =', thickness
  
  do i = 0, m+1
  xp(i) = dx * real(i-1) - width*0.5
  end do
  
  do j = 0, n+1
  yp(j) = dy * real(j-1) - height*0.5
  end do
  
  ! default: outlet condtion in x-direction
  do j = 1, n+1
  porosity(0,j) = porosity(1,j)
  porosity(m+1,j) = porosity(m,j)
  end do
  
  ! default: periodic condtion in y-direction
  do i = 0, m+1
  porosity(i,0)   = porosity(i,n)
  porosity(i,n+1) = porosity(i,1)
  end do
  
  ! ----------------
  return
  end subroutine  grid_conditions
  !******************
  
  !******************
  subroutine  initial_conditions (p, u, v, xp, yp, width, height  &
                                 , inlet_velocity, outlet_pressure, AoA, m, n)
   implicit none
   integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(in)::	width, height, inlet_velocity, outlet_pressure, AoA
   real,intent(out),dimension(0:md,0:nd)::	u, v, p
   real,intent(in),dimension(0:md)::		xp
   real,intent(in),dimension(0:nd)::		yp
   integer,intent(in)::	m, n
  
  ! local variables
  integer::	i, j
  real, parameter :: pai=atan(1.)*4.
  
  ! ----------------
  do j = 1, n
  do i = 1, m
   u(i,j)=inlet_velocity*cos(AoA/180*pai)
   v(i,j)=inlet_velocity*sin(AoA/180*pai)
   p(i,j)=outlet_pressure
  end do
  end do
  ! ----------------
  
  return
  end subroutine initial_conditions
  !******************
  
  ! output
  
  !******************
  subroutine  output_solution (p, u, v, m, n)
   implicit none
   integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md,0:nd)::	u, v, p
   integer,intent(in)::	m, n
  
  ! local variables
  integer::	i, j
  
  ! ----------------
  write(*,*)
  
  write(*,*)'velocity u '
  do j = 0, n+1
  write(*,*) (u(i,j), i=0,m+1)
  end do
  write(*,*)
  
  write(*,*)'velocity v '
  do j = 0, n+1
  write(*,*) (v(i,j), i=0,m+1)
  end do
  write(*,*)
  
  write(*,*)'pressure'
  do j = 0, n+1
  write(*,*) (p(i,j), i=0,m+1)
  end do
  write(*,*)
  ! ----------------
  
  return
  end subroutine output_solution
  !******************
  
  !******************
  subroutine  output_grid (xp, yp, m, n)
   implicit none
   integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md)::	xp
   real,intent(in),dimension(0:nd)::	yp
   integer,intent(in)::	m, n
  
  ! local variables
  integer::	i, j
  
  open (60, file='etc/grid.dat', status='replace')
  ! ----------------
  write(60,*)'m, n =', m, n
  write(60,*)'grid points ='
  write(60,*) (xp(i), i=1,m)
  write(60,*) (yp(j), j=1,n)
  ! ----------------
  close (60)
  return
  end subroutine output_grid
  !******************
  
  !******************
  subroutine  output_grid_list (xp, yp, m, n, angle_of_attack)
   implicit none
   integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md)::	xp
   real,intent(in),dimension(0:nd)::	yp
   integer,intent(in)::	m, n
   real,intent(in):: angle_of_attack
  
  ! local variables
  integer::	i, j
  real::      z=0.0, pai=atan(1.)*4.
  real::      x, y, th
  
  open (60, file='etc/cellcenter.dat', status='replace')
  ! ----------------
  th = angle_of_attack/180.*pai
  do i=1,m
  do j=1,n
   x=xp(i)*cos(th)-yp(j)*sin(th)
   y=xp(i)*sin(th)+yp(j)*cos(th)
   write(60,*) x,y,z
  end do
  end do
  ! ----------------
  close (60)
  return
  end subroutine output_grid_list
  !******************
  
  !******************
  subroutine  output_solution_post (p, u, v, xp, yp, porosity, m, n)
   implicit none
   integer,parameter::md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md,0:nd)::u, v, p
   real,intent(in),dimension(0:md,0:nd)::porosity
   real,intent(in),dimension(0:md)::xp
   real,intent(in),dimension(0:nd)::yp
   integer,intent(in)::m, n
  
  ! local variables
   real, parameter::small=1.e-6, big=1.e6, zero=0.
   real, parameter::pmin=0.25, pmax=0.75
   integer::i, j
   real,dimension(0:md, 0:nd)::u_cnt, v_cnt, p_cnt
  
  open (61, file='etc/solution_uvp.dat', status='replace')
  
  ! ----------------
  ! interpolation at p-center grid
  
  do i = 1, m
  do j = 1, n
   u_cnt(i,j)=u(i,j)*porosity(i,j)
   v_cnt(i,j)=v(i,j)*porosity(i,j)
   if (porosity(i,j) > small)then
    p_cnt(i,j)=p(i,j)
   else
    p_cnt(i,j)=zero
   end if
  end do
  end do
  
  do j = 1, n
   u_cnt(0,j)=u_cnt(1,j)
   v_cnt(0,j)=v_cnt(1,j)
   p_cnt(0,j)=p_cnt(1,j)
   u_cnt(m+1,j)=u_cnt(m,j)
   v_cnt(m+1,j)=v_cnt(m,j)
   p_cnt(m+1,j)=p_cnt(m,j)
  end do
  
  do i = 0, m+1
   u_cnt(i,0)=u_cnt(i,1)
   v_cnt(i,0)=v_cnt(i,1)
   p_cnt(i,0)=p_cnt(i,1)
   u_cnt(i,n+1)=u_cnt(i,n)
   v_cnt(i,n+1)=v_cnt(i,n)
   p_cnt(i,n+1)=p_cnt(i,n)
  end do
  
  !-----------------
  write(61,*)'m, n =', m, n
  
  write(61,*)'velocity u_bulk '
  do j = 1, n
  write(61,*) (u_cnt(i,j), i=1,m)
  end do
  
  write(61,*)'velocity v_bulk '
  do j = 1, n
  write(61,*) (v_cnt(i,j), i=1,m)
  end do
  
  write(61,*)'velocity u_inst '
  do j = 1, n
  write(61,*) (u(i,j), i=1,m)
  end do
  
  write(61,*)'velocity v_inst '
  do j = 1, n
  write(61,*) (v(i,j), i=1,m)
  end do
  
  write(61,*)'pressure p_fluid'
  do j = 1, n
  write(61,*) (p_cnt(i,j), i=1,m)
  end do
  
  write(61,*)'pressure P_all'
  do j = 1, n
  write(61,*) (p(i,j), i=1,m)
  end do
  
  write(61,*)'porosity'
  do j = 1, n
  write(61,*) (porosity(i,j), i=1,m)
  end do
  
  !write(*,*) "check", (porosity(i,10), i=1,m)
  
  close (61)
  ! ----------------
  
  ! ----------------
  ! surface profile
  open (62, file='etc/surface_profile.dat', status='replace')
  
  do j=1,n
  do i=1,m
  
  if( porosity(i,j) < pmax .and. porosity(i,j)>pmin )then
   write(62,*) xp(i), yp(j), p_cnt(i,j), porosity(i,j)
  end if
  
  end do
  end do
  close (62)
  ! ----------------
  
  
  return
  end subroutine output_solution_post
  !******************
  
  !******************
  subroutine  output_paraview (p, u, v, porosity, xp, yp, m, n)
   implicit none
   integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md)::	xp
   real,intent(in),dimension(0:nd)::	yp
   real,intent(in),dimension(0:md, 0:nd)::	u, v, p
   real,intent(in),dimension(0:md,0:nd)::	porosity
   integer,intent(in)::	m, n
   integer::	i, j
  
  ! local variables
   real,dimension(0:md,0:nd):: div
  
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
      write(50,"('2D flow')")
      write(50,"('ASCII ')")
  
      write(50,"('DATASET STRUCTURED_GRID')")
      write(50,"('DIMENSIONS ',3(1x,i4))") m, n, 1
  
      write(50,"('POINTS ',i9,' float')") m*n
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") xp(i), yp(j), 0.0d0
      enddo
      enddo
  
      write(50,"('POINT_DATA ',i9)") m*n
  
  !! velocity vector
      write(50,"('VECTORS velocity float')")
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") u(i,j), v(i,j), 0.0d0
      enddo
      enddo
  
  !! velocity vector
      write(50,"('VECTORS velocityInFluid float')")
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") u(i,j)*porosity(i,j), v(i,j)*porosity(i,j), 0.0d0
      enddo
      enddo
  
  !! pressure
      write(50,"('SCALARS pressure float')")
      write(50,"('LOOKUP_TABLE default')")
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") p(i,j)
      enddo
      enddo
  
  do i = 1, m
  do j = 1, n
   div(i,j)= (u(i+1,j)-u(i-1,j))/(xp(j+1)-xp(j-1))+(v(i,j+1)-v(i,j-1))/(yp(j+1)-yp(j-1))
  end do
  end do
  
  !! divergent velocity
      write(50,"('SCALARS VelocityDivergent float')")
      write(50,"('LOOKUP_TABLE default')")
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") div(i,j)
      enddo
      enddo
  
  !! porosity
      write(50,"('SCALARS porosity float')")
      write(50,"('LOOKUP_TABLE default')")
      do j=1,n
      do i=1,m
        write(50,"(3(f16.4,1x))") porosity(i,j)
      enddo
      enddo
  
  ! ----------------
  close(50)
  
  return
  end subroutine  output_paraview
  !******************
  
  !******************
  subroutine  output_divergent (p, u, v, porosity, dx, dy, m, n)
   implicit none
   integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md,0:nd)::	u, v, p
   real,intent(in),dimension(0:md,0:nd)::	porosity
   real,intent(in)::	dx, dy
   integer,intent(in)::	m, n
  
  ! local variables
  integer::	i, j
  real,dimension(0:md,0:nd)::	div
  
  open (63, file='etc/divergent.dat', status='replace')
  ! ----------------
  
  do i = 1, m
  do j = 1, n
  div(i,j)= ((porosity(i+1,j)*u(i,j)+porosity(i,j)*u(i+1,j))/2     &
           -(porosity(i-1,j)*u(i,j)+porosity(i,j)*u(i-1,j))/2 )/dx &
          +((porosity(i,j+1)*v(i,j)+porosity(i,j)*v(i,j+1))/2      &
           -(porosity(i,j-1)*v(i,j)+porosity(i,j)*v(i,j-1))/2 )/dy
  end do
  end do
  
  write(63,*)
  write(63,*)'porosity'
  do j = 1, n
  write(63,*) (porosity(i,j), i=1,m)
  end do
  
  write(63,*)
  write(63,*)'divergent velocity'
  do j = 1, n
  write(63,*) (div(i,j), i=1,m)
  end do
  write(63,*)
  
  ! ----------------
  close (63)
  
  end subroutine  output_divergent
  !******************
  
  !******************
  subroutine  output_paraview_temp (p, u, v, porosity, xp, yp, m, n, istep)
   implicit none
   integer,parameter::	md=4500, nd=4500     ! md, nd > grid size (m,n)
   real,intent(in),dimension(0:md)::	xp
   real,intent(in),dimension(0:nd)::	yp
   real,intent(in),dimension(0:md, 0:nd)::	u, v, p
   real,intent(in),dimension(0:md, 0:nd)::	porosity
   integer,intent(in)::	m, n, istep
  
  ! -- local variable
   real,dimension(0:md,0:nd):: div
   integer::	i, j
   character(5):: number
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
      write(65,"('2D flow')")
      write(65,"('ASCII ')")
  
      write(65,"('DATASET STRUCTURED_GRID')")
      write(65,"('DIMENSIONS ',3(1x,i4))") m, n, 1
  
      write(65,"('POINTS ',i9,' float')") m*n
      do j=1,n
      do i=1,m
        write(65,"(3(f16.4,1x))") xp(i), yp(j), 0.0d0
      enddo
      enddo
  
      write(65,"('POINT_DATA ',i9)") m*n
  
  !! velocity vector
      write(65,"('VECTORS velocity float')")
      do j=1,n
      do i=1,m
        write(65,"(3(f16.4,1x))") u(i,j), v(i,j), 0.0d0
      enddo
      enddo
  
  !! velocity vector no filter
      write(65,"('VECTORS velocityInFluid float')")
      do j=1,n
      do i=1,m
        write(65,"(3(f16.4,1x))") u(i,j)*porosity(i,j), v(i,j)*porosity(i,j), 0.0d0
      enddo
      enddo
  
  !! pressure
      write(65,"('SCALARS pressure float')")
      write(65,"('LOOKUP_TABLE default')")
      do j=1,n
      do i=1,m
        write(65,"(3(f16.4,1x))") p(i,j)
      enddo
      enddo
  
  do i = 1, m
  do j = 1, n
   div(i,j)= (u(i+1,j)-u(i-1,j))/(xp(j+1)-xp(j-1))+(v(i,j+1)-v(i,j-1))/(yp(j+1)-yp(j-1))
  end do
  end do
  
  !! velocity Ddivergent
      write(65,"('SCALARS VelocityDivergent float')")
      write(65,"('LOOKUP_TABLE default')")
      do j=1,n
      do i=1,m
        write(65,"(3(f16.4,1x))") div(i,j)
      enddo
      enddo
  
  !! porosity
      write(65,"('SCALARS porosity float')")
      write(65,"('LOOKUP_TABLE default')")
      do j=1,n
      do i=1,m
        write(65,"(3(f16.4,1x))") porosity(i,j)
      enddo
      enddo
  
  ! ----------------
  close(65)
  
  return
  end subroutine  output_paraview_temp
  !******************