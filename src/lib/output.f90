! output.f90
! link: global.f90

module output_2d
    !$use omp_lib
    implicit none
    contains
    subroutine  output_solution_2d (p, u, v, m, n)
        use global_2d
          implicit none
          real,intent(in),dimension(0:md,0:nd)::u, v, p
          integer,intent(in)::m, n
      
        ! local variables
        integer::i, j
      
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
    end subroutine output_solution_2d
      
    subroutine  output_grid_2d (xp, yp, m, n)
        use global_2d
        implicit none
        real,intent(in),dimension(0:md)::xp
        real,intent(in),dimension(0:nd)::yp
        integer,intent(in)::m, n
      
        ! local variables
        integer::i, j
      
        open (60, file='etc/grid.dat', status='replace')
        ! ----------------
        write(60,*)'m, n =', m, n
        write(60,*)'grid points ='
        write(60,*) (xp(i), i=1,m)
        write(60,*) (yp(j), j=1,n)
        ! ----------------
        close (60)
        return
    end subroutine output_grid_2d
      
    subroutine  output_grid_list_2d (xp, yp, m, n, angle_of_attack)
        use global_2d
        implicit none
        real,intent(in),dimension(0:md)::xp
        real,intent(in),dimension(0:nd)::yp
        integer,intent(in)::m, n
        real,intent(in)::angle_of_attack
      
        ! local variables
        integer::i, j
        real, parameter::z=0.0, pai=atan(1.)*4.
        real::x, y, th
      
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
    end subroutine output_grid_list_2d
      
    subroutine  output_solution_post_2d (p, u, v, xp, yp, porosity, m, n)
        use global_2d
        implicit none
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
    end subroutine output_solution_post_2d
      
    subroutine  output_divergent_2d (p, u, v, porosity, dx, dy, m, n)
        use global_2d
        implicit none
        real,intent(in),dimension(0:md,0:nd)::u, v, p
        real,intent(in),dimension(0:md,0:nd)::porosity
        real,intent(in)::dx, dy
        integer,intent(in)::m, n
      
        ! local variables
        integer::i, j
        real,dimension(0:md,0:nd)::div
      
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
      
    end subroutine  output_divergent_2d
      
    subroutine  output_force_log_2d (p, u, v, dx, dy, porosity, m, n, xnue, density, thickness, radius, inlet_velocity)
        use global_2d
        implicit none
        real,intent(in),dimension(0:md,0:nd)::u, v, p
        real,intent(in),dimension(0:md,0:nd)::porosity
        real,intent(in)::dx, dy
        integer,intent(in)::m, n
        real,intent(in)::xnue, density, thickness, radius, inlet_velocity
      
        ! local variables
        real, parameter::small=1.e-6, big=1.e6, zero=0.
        real, parameter::alpha = 32.0
        integer::i, j
        real::unit_normal_x_tmp, unit_normal_y_tmp
        real::delta_force_px_tmp, delta_force_py_tmp
        real::delta_force_vx_tmp, delta_force_vy_tmp
        real::normal_abs
        real::force_x, force_y, force_px, force_vx, force_py, force_vy
        real::cd, cl
      
        ! ----------------
        
        force_px = 0.0
        force_vx = 0.0
        force_py = 0.0
        force_vy = 0.0
      
        !$omp parallel do private(i,j) reduction(+:force_px, force_py, force_vx, force_vy) &
        !$omp shared(u, v, p, porosity, m, n, dx, dy, xnue, density)
        do i = 1, m
        do j = 1, n
          normal_abs = sqrt(((porosity(i+1,j) - porosity(i-1,j))*0.5)**2+((porosity(i,j+1) - porosity(i,j-1))*0.5)**2)
          unit_normal_x_tmp = (porosity(i+1,j) - porosity(i-1,j))*0.5 / max(normal_abs,small)
          unit_normal_y_tmp = (porosity(i,j+1) - porosity(i,j-1))*0.5 / max(normal_abs,small)
      
          delta_force_px_tmp = -dx*dy*p(i,j)*2*porosity(i,j)*(1.0-porosity(i,j))/(thickness*dx)*unit_normal_x_tmp
          delta_force_py_tmp = -dx*dy*p(i,j)*2*porosity(i,j)*(1.0-porosity(i,j))/(thickness*dy)*unit_normal_y_tmp
      
          delta_force_vx_tmp = +dx*dy*alpha*density*xnue*((porosity(i,j)*(1.0-porosity(i,j)))/(thickness*dx))**2*u(i,j)
          delta_force_vy_tmp = +dx*dy*alpha*density*xnue*((porosity(i,j)*(1.0-porosity(i,j)))/(thickness*dy))**2*v(i,j)
      
          force_px = force_px + delta_force_px_tmp
          force_py = force_py + delta_force_py_tmp
          force_vx = force_vx + delta_force_vx_tmp
          force_vy = force_vy + delta_force_vy_tmp
        end do
        end do
        !$omp end parallel do
      
        force_x = force_px + force_vx
        force_y = force_py + force_vy
      
        cd = force_x / (density * inlet_velocity ** 2 * radius)
        cl = force_y / (density * inlet_velocity ** 2 * radius)
      
        write(*,*)'Fp =',force_px,force_py
        write(*,*)'Fv =',force_vx,force_vy
        write(*,*)'F  =',force_x,force_y
        write(*,*)'Cd =', cd, 'Cl =', cl
      
        return
    end subroutine output_force_log_2d
      
    subroutine  output_paraview_2d (p, u, v, porosity, xp, yp, m, n, inlet_velocity, output_folder)
        use global_2d
        implicit none
        real,intent(in),dimension(0:md)::xp
        real,intent(in),dimension(0:nd)::yp
        real,intent(in),dimension(0:md, 0:nd)::u, v, p
        real,intent(in),dimension(0:md,0:nd)::porosity
        integer,intent(in)::m, n
        real,intent(in)::inlet_velocity
      
        ! local variables
        integer::i, j
        real,dimension(0:md,0:nd):: div
        character(len=50)::output_folder
        real, parameter::small=1.e-6, big=1.e6, zero=0.
      
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
      
        !! velocity vector
        write(50,"('VECTORS dimless_v float')")
        do j=1,n
        do i=1,m
        write(50,"(3(f16.4,1x))") u(i,j)*porosity(i,j)/inlet_velocity, v(i,j)*porosity(i,j)/inlet_velocity, 0.0d0
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
      
        !! pressure
        write(50,"('SCALARS pressure float')")
        write(50,"('LOOKUP_TABLE default')")
        do j=1,n
        do i=1,m
        write(50,"(3(f16.4,1x))") p(i,j)
        enddo
        enddo
      
        !$omp parallel private(i, j) &
        !$omp & shared(div, u, v, xp, yp, m, n) &
        !$omp & default(none)
        !$omp do
        do j = 1, n
        do i = 1, m
        div(i,j)= (u(i+1,j)-u(i-1,j))/(xp(i+1)-xp(i-1))+(v(i,j+1)-v(i,j-1))/(yp(j+1)-yp(j-1))
        end do
        end do
        !$omp end do
        !$omp end parallel
      
        !! divergent velocity
        write(50,"('SCALARS VelocityDivergent float')")
        write(50,"('LOOKUP_TABLE default')")
        do j=1,n
        do i=1,m
        write(50,"(3(f16.4,1x))") div(i,j)
        enddo
        enddo
      
        !! dimless_velocity
        write(50,"('SCALARS abs_dimless_v float')")
        write(50,"('LOOKUP_TABLE default')")
        do j=1,n
        do i=1,m
        write(50,"(3(f16.4,1x))") sqrt((u(i,j)*porosity(i,j)/inlet_velocity)**2+(v(i,j)*porosity(i,j)/inlet_velocity)**2)
        enddo
        enddo
      
        ! ----------------
        close(50)
      
        return
    end subroutine  output_paraview_2d
      
    subroutine  output_paraview_temp_2d (p, u, v, porosity, xp, yp, m, n, inlet_velocity, istep, output_folder)
        use global_2d
        implicit none
        real,intent(in),dimension(0:md)::xp
        real,intent(in),dimension(0:nd)::yp
        real,intent(in),dimension(0:md, 0:nd)::u, v, p
        real,intent(in),dimension(0:md,0:nd)::porosity
        integer,intent(in)::m, n, istep
        real,intent(in)::inlet_velocity
      
        ! local variables
        integer::i, j
        real,dimension(0:md,0:nd):: div
        character(5):: number
        character(len=50)::output_folder
        real, parameter::small=1.e-6, big=1.e6, zero=0.
      
        write(number,"(I5.5)")istep
      
        ! -- open file
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
      
        !! velocity vector
        write(65,"('VECTORS velocityInFluid float')")
        do j=1,n
        do i=1,m
        write(65,"(3(f16.4,1x))") u(i,j)*porosity(i,j), v(i,j)*porosity(i,j), 0.0d0
        enddo
        enddo
      
        !! velocity vector
        write(65,"('VECTORS dimless_v float')")
        do j=1,n
        do i=1,m
        write(65,"(3(f16.4,1x))") u(i,j)*porosity(i,j)/inlet_velocity, v(i,j)*porosity(i,j)/inlet_velocity, 0.0d0
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

        !! pressure
        write(65,"('SCALARS pressure float')")
        write(65,"('LOOKUP_TABLE default')")
        do j=1,n
        do i=1,m
        write(65,"(3(f16.4,1x))") p(i,j)
        enddo
        enddo
      
        !$omp parallel private(i, j) &
        !$omp & shared(div, u, v, xp, yp, m, n) &
        !$omp & default(none)
        !$omp do
        do j = 1, n
        do i = 1, m
        div(i,j)= (u(i+1,j)-u(i-1,j))/(xp(i+1)-xp(i-1))+(v(i,j+1)-v(i,j-1))/(yp(j+1)-yp(j-1))
        end do
        end do
        !$omp end do
        !$omp end parallel
      
        !! divergent velocity
        write(65,"('SCALARS VelocityDivergent float')")
        write(65,"('LOOKUP_TABLE default')")
        do j=1,n
        do i=1,m
        write(65,"(3(f16.4,1x))") div(i,j)
        enddo
        enddo
      
        !! dimless_velocity
        write(65,"('SCALARS abs_dimless_v float')")
        write(65,"('LOOKUP_TABLE default')")
        do j=1,n
        do i=1,m
        write(65,"(3(f16.4,1x))") sqrt((u(i,j)*porosity(i,j)/inlet_velocity)**2+(v(i,j)*porosity(i,j)/inlet_velocity)**2)
        enddo
        enddo
      
        ! ----------------
        close(65)
      
        return
    end subroutine  output_paraview_temp_2d
end module output_2d

module output_3d
    !$use omp_lib
    implicit none
    contains
    subroutine  output_solution_3d (p, u, v, w, m, n, l)

        use global_3d
        implicit none
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
    
    end subroutine output_solution_3d

    subroutine  output_grid_3d (xp, yp, zp, m, n, l)

        use global_3d
        implicit none
        real,intent(in),dimension(0:md)::xp
        real,intent(in),dimension(0:nd)::yp
        real,intent(in),dimension(0:ld)::zp
        integer,intent(in)::m, n, l

        ! local variables
        integer::i, j, k
        
        open (63, file='etc/grid.dat', status='replace')

        write(63,*)'m, n, l =', m, n, l
        write(63,*)'grid points ='
        write(63,*) (xp(i), i=1,m)
        write(63,*) (yp(j), j=1,n)
        write(63,*) (zp(k), k=1,l)

        close (63)

    end subroutine output_grid_3d

    subroutine  output_grid_list_3d (xp, yp, zp, m, n, l, angle_of_attack)

        use global_3d
        implicit none
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

        th = angle_of_attack/1300.*pai
        do k = 1, l
            do j = 1, n
            do i = 1, m
                x = xp(i)*cos(th)-yp(j)*sin(th)
                y = xp(i)*sin(th)+yp(j)*cos(th)
                z = zp(k)
                write(67,*) x, y, z
            end do
            end do
        end do

        close (67)

    end subroutine output_grid_list_3d

    subroutine  output_solution_post_3d (p, u, v, w, xp, yp, zp, porosity, m, n, l)

        use global_3d
        implicit none
        real,intent(in),dimension(0:md,0:nd,0:ld)::u, v, w, p
        real,intent(in),dimension(0:md,0:nd,0:ld)::porosity
        real,intent(in),dimension(0:md)::xp
        real,intent(in),dimension(0:nd)::yp
        real,intent(in),dimension(0:ld)::zp
        integer,intent(in)::m, n, l

        ! local variables
        real, parameter::small=1.e-6, zero=0.
        real, parameter::pmin=0.25, pmax=0.75
        integer::i, j, k
        real,dimension(0:md, 0:nd, 0:ld)::u_cnt, v_cnt, w_cnt, p_cnt
        
        !$omp parallel private(i, j, k) &
        !$omp & shared(m, n, l) &
        !$omp & shared(p, u, v, w, xp, yp, zp, porosity) &
        !$omp & shared(u_cnt, v_cnt, w_cnt, p_cnt) &
        !$omp & default(none)
        
        !$omp do 
        do k = 1, l
            do j = 1, n
            do i = 1, m
                u_cnt(i,j,k) = u(i,j,k)*porosity(i,j,k)
                v_cnt(i,j,k) = v(i,j,k)*porosity(i,j,k)
                w_cnt(i,j,k) = w(i,j,k)*porosity(i,j,k)
                if (porosity(i,j,k) > small) then
                p_cnt(i,j,k) = p(i,j,k)
                else
                p_cnt(i,j,k) = zero
                end if 
            end do
            end do
        end do
        !$omp end do
        
        !$omp do
        do k = 1, l
            do j = 1, n
            u_cnt(0,j,k) = u_cnt(1,j,k)
            v_cnt(0,j,k) = v_cnt(1,j,k)
            w_cnt(0,j,k) = w_cnt(1,j,k)
            p_cnt(0,j,k) = p_cnt(1,j,k)
            u_cnt(m+1,j,k) = u_cnt(m,j,k)
            v_cnt(m+1,j,k) = v_cnt(m,j,k)
            w_cnt(m+1,j,k) = w_cnt(m,j,k)
            p_cnt(m+1,j,k) = p_cnt(m,j,k)
            end do
        end do
        !$omp end do
        
        !$omp do
        do i = 0, m+1
            u_cnt(i,0,k) = u_cnt(i,1,k)
            v_cnt(i,0,k) = v_cnt(i,1,k)
            w_cnt(i,0,k) = w_cnt(i,1,k)
            p_cnt(i,0,k) = p_cnt(i,1,k)
            u_cnt(i,n+1,k) = u_cnt(i,n,k)
            v_cnt(i,n+1,k) = v_cnt(i,n,k)
            w_cnt(i,n+1,k) = w_cnt(i,n,k)
            p_cnt(i,n+1,k) = p_cnt(i,n,k)
        end do
        !$omp end do

        !$omp end parallel

        
        open (61, file='etc/solution_uvp.dat', status='replace')

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
    
    end subroutine output_solution_post_3d

    subroutine  output_paraview_3d (p, u, v, w, porosity, xp, yp, zp, m, n, l)

        use global_3d
        implicit none
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
        do k = 1, l
            do j = 1, n
            do i = 1, m
                write(50,"(3(f16.4,1x))") xp(i), yp(j), zp(k)
            enddo
            enddo
        enddo
        
        write(50,"('POINT_DATA ',i9)") m*n*l
            
        !-- velocity vector
        write(50,"('VECTORS velocity float')")
        do k = 1, l
            do j = 1, n
            do i = 1, m
                write(50,"(3(f16.4,1x))") u(i,j,k), v(i,j,k), w(i,j,k)
            enddo
            enddo
        enddo
        
        !-- velocity vector
        write(50,"('VECTORS velocityInFluid float')")
        do k = 1, l
            do j = 1, n
            do i = 1, m
                write(50,"(3(f16.4,1x))") u(i,j,k)*porosity(i,j,k), v(i,j,k)*porosity(i,j,k), w(i,j,k)*porosity(i,j,k)
            enddo
            enddo
        enddo
            
        !-- pressure
        write(50,"('SCALARS pressure float')")
        write(50,"('LOOKUP_TABLE default')")
        do k = 1, l
            do j = 1, n
            do i = 1, m
                write(50,"(3(f16.4,1x))") p(i,j,k)
            enddo
            enddo
        enddo
        
        !$omp parallel private(i, j, k) &
        !$omp & shared(div, u, v, w, xp, yp, zp, m, n, l) &
        !$omp & default(none)
        !$omp do
        do k = 1, l
            do j = 1, n
            do i = 1, m
                div(i,j,k) = (u(i+1,j,k)-u(i-1,j,k))/(xp(i+1)-xp(i-1)) &
                            +(v(i,j+1,k)-v(i,j-1,k))/(yp(j+1)-yp(j-1)) &
                            +(w(i,j,k+1)-w(i,j,k-1))/(zp(k+1)-zp(k-1))
            end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        
        !-- divergent velocity
        write(50,"('SCALARS VelocityDivergent float')")
        write(50,"('LOOKUP_TABLE default')")
        do k = 1, l
            do j = 1, n
            do i = 1, m
                write(50,"(3(f16.4,1x))") div(i,j,k)
            end do
            end do
        end do
        
        !-- porosity
        write(50,"('SCALARS porosity float')")
        write(50,"('LOOKUP_TABLE default')")
        do k = 1, l
            do j = 1, n
            do i = 1, m
                write(50,"(3(f16.4,1x))") porosity(i,j,k)
            end do
            end do
        end do
        
        ! ----------------
        close(50)

    end subroutine  output_paraview_3d

    subroutine  output_divergent_3d (p, u, v, w, porosity, dx, dy, dz, m, n, l)

        use global_3d
        implicit none
        real,intent(in),dimension(0:md,0:nd,0:ld)::u, v, w, p 
        real,intent(in),dimension(0:md,0:nd,0:ld)::porosity
        real,intent(in)::dx, dy, dz
        integer,intent(in)::m, n, l

        ! local variables
        integer::i, j, k
        real,dimension(0:md,0:nd,0:ld)::div
        
        !$omp parallel private(i, j, k) &
        !$omp & shared(div, u, v, w, porosity, m, n, l, dx, dy, dz) &
        !$omp & default(none)
        !$omp do
        do k = 1, l
            do j = 1, n
            do i = 1, m
            div(i,j,k)= ((porosity(i+1,j,k)*u(i,j,k)+porosity(i,j,k)*u(i+1,j,k))/2      &
                        -(porosity(i-1,j,k)*u(i,j,k)+porosity(i,j,k)*u(i-1,j,k))/2 )/dx &
                        +((porosity(i,j+1,k)*v(i,j,k)+porosity(i,j,k)*v(i,j+1,k))/2      &
                        -(porosity(i,j-1,k)*v(i,j,k)+porosity(i,j,k)*v(i,j-1,k))/2 )/dy &
                        +((porosity(i,j,k+1)*w(i,j,k)+porosity(i,j,k)*w(i,j,k+1))/2      &
                        -(porosity(i,j,k-1)*w(i,j,k)+porosity(i,j,k)*w(i,j,k-1))/2 )/dz
            end do
            end do
        end do
        !$omp end do
        !$omp end parallel
        
        open (62, file='etc/divergent.dat', status='replace')

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
    
    end subroutine  output_divergent_3d

    subroutine  output_paraview_temp_3d (p, u, v, w, porosity, xp, yp, zp, m, n, l, istep)

        use global_3d
        implicit none
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
        do k = 1, l
            do j = 1, n
            do i = 1, m
                write(65,"(3(f16.4,1x))") xp(i), yp(j), zp(k)
            enddo
            enddo
        enddo
        
        write(65,"('POINT_DATA ',i9)") m*n*l
        
        !-- velocity vector
        write(65,"('VECTORS velocity float')")
        do k = 1, l
            do j = 1, n
            do i = 1, m
                write(65,"(3(f16.4,1x))") u(i,j,k), v(i,j,k), w(i,j,k)
            enddo
            enddo
        enddo

        !-- velocity vector
        write(65,"('VECTORS velocityInFluid float')")
        do k = 1, l
            do j = 1, n
            do i = 1, m
                write(65,"(3(f16.4,1x))") u(i,j,k)*porosity(i,j,k), v(i,j,k)*porosity(i,j,k), w(i,j,k)*porosity(i,j,k)
            enddo
            enddo
        enddo

        !-- porosity
        write(65,"('SCALARS porosity float')")
        write(65,"('LOOKUP_TABLE default')")
        do k = 1, l
            do j = 1, n
            do i = 1, m
                write(65,"(3(f16.4,1x))") porosity(i,j,k)
            end do
            end do
        end do

        !-- pressure
        write(65,"('SCALARS pressure float')")
        write(65,"('LOOKUP_TABLE default')")
        do k = 1, l
            do j = 1, n
            do i = 1, m
                write(65,"(3(f16.4,1x))") p(i,j,k)
            enddo
            enddo
        enddo

        !$omp parallel private(i, j, k) &
        !$omp & shared(div, u, v, w, xp, yp, zp, m, n, l) &
        !$omp & default(none)
        !$omp do
        do k = 1, l
            do j = 1, n
            do i = 1, m
                div(i,j,k) = (u(i+1,j,k)-u(i-1,j,k))/(xp(i+1)-xp(i-1)) &
                        + (v(i,j+1,k)-v(i,j-1,k))/(yp(j+1)-yp(j-1)) &
                        + (w(i,j,k+1)-w(i,j,k-1))/(zp(k+1)-zp(k-1))
            end do
            end do
        end do
        !$omp end do
        !$omp end parallel

        !-- divergent velocity
        write(65,"('SCALARS VelocityDivergent float')")
        write(65,"('LOOKUP_TABLE default')")
        do k = 1, l
            do j = 1, n
            do i = 1, m
                write(65,"(3(f16.4,1x))") div(i,j,k)
            end do
            end do
        end do

        ! ----------------
        close(65)
    
    end subroutine  output_paraview_temp_3d

    subroutine  output_force_log_3d (p, u, v, w, dx, dy, dz, porosity, m, n, l,&
         xnue, density, thickness, radius, inlet_velocity)
        use global_3d
        implicit none
        real,intent(in),dimension(0:md,0:nd,0:ld)::u, v, w, p
        real,intent(in),dimension(0:md,0:nd,0:ld)::porosity
        real,intent(in)::dx, dy, dz
        integer,intent(in)::m, n, l
        real,intent(in)::xnue, density, thickness, radius, inlet_velocity
      
        ! local variables
        real, parameter::small=1.e-6, big=1.e6, zero=0.
        real, parameter::alpha = 32.0
        integer::i, j, k
        real::unit_normal_x_tmp, unit_normal_y_tmp, unit_normal_z_tmp
        real::delta_force_px_tmp, delta_force_py_tmp, delta_force_pz_tmp
        real::delta_force_vx_tmp, delta_force_vy_tmp, delta_force_vz_tmp
        real::normal_abs
        real::force_x, force_y, force_z, force_px, force_vx, force_py, force_vy, force_pz, force_vz
        real::cdx, cdz, cl
      
        force_px = 0.0
        force_vx = 0.0
        force_py = 0.0
        force_vy = 0.0
        force_pz = 0.0
        force_vz = 0.0
      
        !$omp parallel do private(i,j) reduction(+:force_px, force_py, force_vx, force_vy) &
        !$omp shared(u, v, p, porosity, m, n, dx, dy, xnue, density)
        do k = 1, l
        do j = 1, n
        do i = 1, m
          normal_abs = sqrt(&
           ((porosity(i+1,j,k) - porosity(i-1,j,k))*0.5)**2&
          +((porosity(i,j+1,k) - porosity(i,j-1,k))*0.5)**2&
          +((porosity(i,j,k+1) - porosity(i,j,k-1))*0.5)**2)

          unit_normal_x_tmp = (porosity(i+1,j,k) - porosity(i-1,j,k))*0.5 / max(normal_abs,small)
          unit_normal_y_tmp = (porosity(i,j+1,k) - porosity(i,j-1,k))*0.5 / max(normal_abs,small)
          unit_normal_z_tmp = (porosity(i,j,k+1) - porosity(i,j,k-1))*0.5 / max(normal_abs,small)
      
          delta_force_px_tmp = -dx*dy*dz*p(i,j,k)*2*porosity(i,j,k)*(1.0-porosity(i,j,k))/(thickness*dx)*unit_normal_x_tmp
          delta_force_py_tmp = -dx*dy*dz*p(i,j,k)*2*porosity(i,j,k)*(1.0-porosity(i,j,k))/(thickness*dy)*unit_normal_y_tmp
          delta_force_pz_tmp = -dx*dy*dz*p(i,j,k)*2*porosity(i,j,k)*(1.0-porosity(i,j,k))/(thickness*dz)*unit_normal_z_tmp
      
          delta_force_vx_tmp = +dx*dy*dz*alpha*density*xnue*((porosity(i,j,k)*(1.0-porosity(i,j,k)))/(thickness*dx))**2*u(i,j,k)
          delta_force_vy_tmp = +dx*dy*dz*alpha*density*xnue*((porosity(i,j,k)*(1.0-porosity(i,j,k)))/(thickness*dy))**2*v(i,j,k)
          delta_force_vz_tmp = +dx*dy*dz*alpha*density*xnue*((porosity(i,j,k)*(1.0-porosity(i,j,k)))/(thickness*dz))**2*w(i,j,k)
      
          force_px = force_px + delta_force_px_tmp
          force_py = force_py + delta_force_py_tmp
          force_pz = force_pz + delta_force_pz_tmp
          force_vx = force_vx + delta_force_vx_tmp
          force_vy = force_vy + delta_force_vy_tmp
          force_vz = force_vz + delta_force_vz_tmp
        end do
        end do
        end do
        !$omp end parallel do
      
        force_x = force_px + force_vx
        force_y = force_py + force_vy
        force_z = force_pz + force_vz
      
        cdx = force_x / (density * inlet_velocity ** 2 * radius)
        cl  = force_y / (density * inlet_velocity ** 2 * radius)
        cdz = force_z / (density * inlet_velocity ** 2 * radius)
      
        write(*,*)'Fp =',force_px,force_py,force_pz
        write(*,*)'Fv =',force_vx,force_vy,force_vz
        write(*,*)'F  =',force_x,force_y,force_z
        write(*,*)'Cd(x) =', cdx, 'Cl =', cl, 'Cd(z) =', cdz
      
        return
    end subroutine output_force_log_3d
end module output_3d