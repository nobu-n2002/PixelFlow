! ibm_2d_omp_cpu.f90
! link: lib/global.f90 lib/output.f90 lib/utils.f90

program main
  !$ use omp_lib
  use global_2d
  use valiables
  use output_2d
  use utils
  use solver
  implicit none
  integer:: istep
  real,dimension(0:md,0:nd):: u, v, p, u_old, v_old
  real,dimension(0:md,0:nd):: porosity
  real:: dx, dy, dt
  real,dimension(0:md):: xp
  real,dimension(0:nd):: yp
  integer:: m, n
  integer:: i, j
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
  ! -----------------

  call grid_conditions (&
  xp, yp, dx, dy, dt, xnue, xlambda, density, width, height, depth,&
  thickness, threshold, radius,&
  center_x, center_y, time,&
  inlet_velocity, AoA, porosity, m, n, istep_max,&
  csv_file)

  ! ----------------
  ! calculation start  (if iest=!0)
  ! ----------------

  call  output_grid_2d (xp, yp, m, n)

  istep = 0
  time = istep * dt

  ! ----------------

  write(*,*) '# istep_max= ', istep_max,'   istep_out= ', istep_out

  call  initial_conditions (p, u, v, xp, yp, width, height &
                        , inlet_velocity, outlet_pressure, AoA, porosity, m, n)
  call  boundary (p, u, v, xp, yp, width, height            &
                        , inlet_velocity, outlet_pressure, AoA, porosity, m, n)

  ! print initial conditions
  ! call  output_solution (p, u, v, m, n)

  ! ----------------
  ! MAC algorithm start
  call get_now_time()
  write(*,*) '# --- MAC algorithm start'

  do istep = 1, istep_max

  time = istep* dt

  !$omp parallel private(i, j) &
  !$omp & shared(m, n) &
  !$omp & shared(u_old, v_old, u, v) &
  !$omp & default(none)
  !$omp do
  do i = 0, m+1
    do j = 0, n+1
    u_old(i,j) = u(i,j)
    v_old(i,j) = v(i,j)
    end do
  end do
  !$omp end do
  !$omp end parallel
  write(*,*)'--- time_steps= ',istep, ' --  time = ',time

  call  solve_p (p, u, v, u_old, v_old, porosity, xnue, xlambda, density, height, thickness, &
  yp, dx, dy, dt, m, n, nonslip, iter_max, relux_factor)

  !-- solve u, v (fix u, v)
  !$omp parallel private(i, j) &
  !$omp & shared(m, n, dt, dx, dy, density) &
  !$omp & shared(p, u, v) &
  !$omp & default(none)
  !$omp do
  do j = 1, n
    do i = 1, m
      u(i,j) = u(i,j) - dt/density*(p(i+1,j)-p(i-1,j))/dx*0.5
      v(i,j) = v(i,j) - dt/density*(p(i,j+1)-p(i,j-1))/dy*0.5
    end do
  end do
  !$omp end do
  !$omp end parallel

  call  boundary(p, u, v, xp, yp, width, height    &
                    , inlet_velocity, outlet_pressure, AoA, porosity, m, n)

  ! call output_force_log_2d (p, u, v, dx, dy, porosity, m, n, xnue, density, thickness, radius, inlet_velocity)

  if(mod(istep,istep_out)==0) then
    call output_paraview_temp_2d (p, u, v, porosity, xp, yp, m, n, inlet_velocity, istep, output_folder)
  end if
  end do
  call get_now_time()
  ! MAC algorithm end
  ! ----------------

  ! print solutions
  call output_solution_post_2d (p, u, v, xp, yp, porosity, m, n)
  call output_divergent_2d (p, u, v, porosity, dx, dy, m, n)
  call output_paraview_2d (p, u, v, porosity, xp, yp, m, n, inlet_velocity, output_folder)
  write(*,*) 'program finished'
  call get_now_time()
  return
end program main
!******************

!  solve variables
!******************
subroutine  solve_p (p, u, v, u_old, v_old, porosity, xnue, xlambda, density, height, thickness, &
  yp, dx, dy, dt, m, n, nonslip, iter_max, relux_factor)
  use global_2d
  use solver
  implicit none
  real,intent(in):: dx, dy, dt
  real,intent(in):: xnue, xlambda, density, height, thickness
  real,intent(inout),dimension(0:md,0:nd):: u, v, p, u_old, v_old
  real,intent(in),dimension(0:md,0:nd):: porosity
  real,intent(in),dimension(0:nd) :: yp
  integer,intent(in):: m, n
  logical,intent(in):: nonslip
  integer,intent(in):: iter_max
  real,intent(in):: relux_factor
  !-----------------
  ! local variables 
  real, parameter :: small = 1.e-6
  real, parameter :: alpha = 32.0

  real,dimension(0:md,0:nd) :: ap, ae, aw, an, as, bb, div
  integer :: i, j

  !$omp parallel private(i, j) &
  !$omp & shared(dx, dy, dt, xnue, density, height, thickness, m, n) &
  !$omp & shared(p, u, v, u_old, v_old, porosity, yp) &
  !$omp & shared(xlambda, nonslip) &
  !$omp & shared(ap, ae, aw, an, as, bb, div) &
  !$omp & default(none)
  !-----------------
  !  divergence term  div(u)
  !-----------------
  !$omp do
  do i = 1, m
    do j = 1, n
        div(i,j)= (u_old(i+1,j)-u_old(i-1,j))/dx*.5 &
                + (v_old(i,j+1)-v_old(i,j-1))/dx*.5
    end do
  end do
  !$omp end do
  
  !$omp do
  do j = 1, n
    div(0,j)  = 0.  ! inlet
    div(m+1,j)= 0.  ! outlet
  end do
  !$omp end do

  !--- periodic condition
  !$omp do
  do i = 1, m
    div(i,0)  = div(i,n) 
    div(i,n+1)= div(i,1)
  end do
  !$omp end do
  ! ----------------

  ! ----------------
  !   velocity u
  ! ----------------
  !$omp do
  do i = 1, m
  do j = 1, n
  ! --- convection_x  2nd central scheme
  u(i,j)=u_old(i,j)-dt*(u_old(i,j)*(u_old(i+1,j)-u_old(i-1,j))/dx/2.)

  ! --- convection_y  2nd central scheme
  u(i,j)=u(i,j)-dt*(v_old(i,j)*(u_old(i,j+1)-u_old(i,j-1))/dy/2.)
  ! --- diffusion_x
  u(i,j)=u(i,j) +dt*xnue*(u_old(i+1,j)-2.*u_old(i,j)+u_old(i-1,j))/dx/dx
  ! --- diffusion_y
  u(i,j)=u(i,j) +dt*xnue*(u_old(i,j+1)-2.*u_old(i,j)+u_old(i,j-1))/dy/dy
  ! --- divergence term
  u(i,j)=u(i,j) +dt*(xnue + xlambda)*(div(i+1,j)-div(i-1,j))/dx*.5
  ! --- additional terms by porosity profile
  u(i,j)=u(i,j)							&
      +dt*( ( (u_old(i+1,j)-u_old(i-1,j))/dx*.5+(u_old(i+1,j)-u_old(i-1,j))/dx*.5) &
              *xnue*(porosity(i+1,j)-porosity(i-1,j))/dx*.5                        &
            +( (u_old(i,j+1)-u_old(i,j-1))/dy*.5+(v_old(i+1,j)-v_old(i-1,j))/dx*.5) &
              *xnue*(porosity(i,j+1)-porosity(i,j-1))/dy*.5                        &
            + div(i,j)*(porosity(i+1,j)-porosity(i-1,j))/dx*0.5*xlambda             &
            )/porosity(i,j)
  ! --- force on wall
  if (nonslip) then
    u(i,j)=u(i,j)- dt*xnue*u_old(i,j)/(thickness*dx)**2*alpha*porosity(i,j)*(1.-porosity(i,j))*(1.-porosity(i,j))
  end if
  end do
  end do
  !$omp end do
  ! ----------------
  !   velocity v
  ! ----------------
  !$omp do
  do i = 1, m
  do j = 1, n
  ! --- convection_x  2nd central scheme
  v(i,j)=v_old(i,j)-dt*(u_old(i,j)*(v_old(i+1,j)-v_old(i-1,j))/dx/2.) 

  ! --- convection_y 2nd central scheme
  v(i,j)=v(i,j)-dt*(v_old(i,j)*(v_old(i,j+1)-v_old(i,j-1))/dy/2.)

  ! --- diffusion_x
  v(i,j)=v(i,j) +dt*xnue*(v_old(i+1,j)-2.*v_old(i,j)+v_old(i-1,j))/dx/dx
  ! --- diffusion_y
  v(i,j)=v(i,j) +dt*xnue*(v_old(i,j+1)-2.*v_old(i,j)+v_old(i,j-1))/dy/dy
  ! --- divergence term
  v(i,j)=v(i,j) +dt*(xnue + xlambda)*(div(i,j+1)-div(i,j-1))/dy*.5
  ! --- additional terms by porosity profile
  v(i,j)=v(i,j)							&
      +dt*( ( (v_old(i+1,j)-v_old(i-1,j))/dx*.5+(u_old(i,j+1)-u_old(i,j-1))/dy*.5) &
              *xnue*(porosity(i+1,j)-porosity(i-1,j))/dx*.5                        &
            +( (v_old(i,j+1)-v_old(i,j-1))/dy*.5+(v_old(i,j+1)-v_old(i,j-1))/dy*.5) &
              *xnue*(porosity(i,j+1)-porosity(i,j-1))/dy*.5                        &
            + div(i,j)*(porosity(i,j+1)-porosity(i,j-1))/dy*0.5*xlambda       &
            )/porosity(i,j)
  ! --- force on wall
  if (nonslip) then
    v(i,j)=v(i,j)- dt*xnue*v_old(i,j)/(thickness*dx)**2*alpha*porosity(i,j)*(1.-porosity(i,j))*(1.-porosity(i,j))
  end if
  end do
  end do
  !$omp end do

  ! ----------------
  ! matrix solution  !  formulation of porous media
  !$omp do
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

  end do
  end do
  !$omp end do
  !$omp end parallel

  call boundrary_matrix (p, ap, ae, aw, an, as, bb, m, n, height, yp)

  call solve_matrix_x_vec_omp (p, ap, ae, aw, an, as, bb, m, n, relux_factor, iter_max)
  ! ----------------
  return
end subroutine solve_p
!******************
  
!******************
! OpenMP Parallelized
! Written only for CPU machine
! No efficiency ensured on GPU machine
subroutine  solve_matrix_vec_omp (p, ap, ae, aw, an, as, bb, m, n, relux_factor, iter_max)
  use global_2d  
  implicit none
  real,intent(inout),dimension(0:md,0:nd)::p
  real,intent(in),dimension(0:md,0:nd)::ap, ae, aw, an, as, bb
  integer,intent(in)::m, n
  real,intent(in)::relux_factor
  integer,intent(in)::iter_max

  ! local variables
  real::error
  real,dimension(0:md, 0:nd)::p_old
  integer::i, j, iter, k

  !$omp parallel private(iter, i, j, k) &
  !$omp & shared(iter_max, relux_factor, m, n) &
  !$omp & shared(error, p_old, p, ap, ae, aw, an, as, bb) &
  !$omp & default(none)

  ! ----------------
  !   SOR algorithm
  ! ----------------
  !$omp single
  error = 0.0
  !$omp end single

  do iter = 1, iter_max

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
    !$omp do reduction(max: error)
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
      error = max(error, abs(p(i,j)-p_old(i,j)))
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
    !$omp do reduction(max: error)
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
      error = max(error, abs(p(i,j)-p_old(i,j)))
    end do
  !$omp end do

  end do

  ! default periodic condition in y-direction
  !$omp do
  do i = 1, m
    p(i, 0)   = p(i, n)
    p(i, n+1) = p(i, 1)
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
subroutine  boundrary_matrix (p, ap, ae, aw, an, as, bb, m, n, height, yp)
  use global_2d
  implicit none
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

!  conditions

!******************
subroutine  boundary(p, u, v, xp, yp, width, height    &
                  , inlet_velocity, outlet_pressure, AoA, porosity, m, n)
  use global_2d
  implicit none
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

!******************
subroutine  grid_conditions (&
xp, yp, dx, dy, dt, xnue, xlambda, density, width, height, depth,&
thickness, threshold, radius,&
center_x, center_y, time,&
inlet_velocity, AoA, porosity, m, n, istep_max,&
csv_file)

  use global_2d
  implicit none
  real,intent(inout),dimension(0:md):: xp
  real,intent(inout),dimension(0:nd):: yp
  real,intent(inout):: dx, dy, dt
  real,intent(in):: xnue, xlambda, density, width, height, depth, time
  real,intent(in):: inlet_velocity, AoA
  real,intent(in):: thickness, threshold, radius, center_x, center_y
  real,intent(inout),dimension(0:md,0:nd)::porosity
  integer,intent(inout):: m, n
  integer,intent(in):: istep_max
  character(len = 50) :: csv_file
  character(len = 50) :: output_folder

  ! local variables
  !real,dimension(0:md,0:nd)::	distance
  real:: cfl_no, pecret_no, diffusion_factor, reynolds_no
  integer:: i, j
  integer:: x, y, z
  real:: poro_val
  ! --- 

  ! read pixel data
  open(52,file=csv_file, form='formatted')

  read(52,*) m,n

  do j = 1, n
    do i = 1, m
    read(52, *) x, y, z, poro_val
    porosity(x, y) = max(poro_val, threshold)
    end do
  end do

  close(52) 

  !--- calc.
  dx = width / real(m-1)
  dy = height / real(n-1)
  dt = time / real(istep_max)

  cfl_no           = inlet_velocity * dt / dx
  pecret_no        = inlet_velocity * dx / xnue
  diffusion_factor = xnue * dt / dy / dy

  reynolds_no      = inlet_velocity * radius * 2.0 / xnue

  !----- check print out
  write(*,*)
  write(*,*) '# --- Grid conditions'
  write(*,*) '# m, n =', m, n
  write(*,*) '# dx, dy =', dx, dy
  write(*,*) '# dt =', dt
  write(*,*) '# cfl_no =', cfl_no
  write(*,*) '# pecret_no =', pecret_no
  write(*,*) '# diffusion_factor =', diffusion_factor
  write(*,*) '# reynolds_no =', reynolds_no
  write(*,*)

  !$omp parallel private(i, j) &
  !$omp & shared(m, n) &
  !$omp & shared(porosity) &
  !$omp & shared(xp, yp) &
  !$omp & shared(dx, dy) &
  !$omp & shared(width, height, center_x, center_y) &
  !$omp & default(none)

  !$omp do      
  do i = 0, m+1
  xp(i) = dx * real(i-1) - width*center_x
  end do
  !$omp end do

  !$omp do      
  do j = 0, n+1
  yp(j) = dy * real(j-1) - height*center_y
  end do
  !$omp end do

  ! default: outlet condtion in x-direction
  !$omp do      
  do j = 1, n+1
  porosity(0,j) = porosity(1,j)
  porosity(m+1,j) = porosity(m,j)
  end do
  !$omp end do

  ! default: periodic condtion in y-direction
  !$omp do      
  do i = 0, m+1
  porosity(i,0)   = porosity(i,n)
  porosity(i,n+1) = porosity(i,1)
  end do
  !$omp end do
  !$omp end parallel
  ! ----------------
  return
end subroutine  grid_conditions
!******************

!******************
subroutine  initial_conditions (p, u, v, xp, yp, width, height  &
                              , inlet_velocity, outlet_pressure, AoA, porosity, m, n)
  use global_2d
  implicit none
  real,intent(in)::	width, height, inlet_velocity, outlet_pressure, AoA
  real,intent(out),dimension(0:md,0:nd)::	u, v, p
  real,intent(in),dimension(0:md,0:nd)::	porosity
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
