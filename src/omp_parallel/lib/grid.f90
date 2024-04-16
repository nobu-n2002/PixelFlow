module grid_2d
    implicit none

    contains

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
end module grid_2d

module grid_3d
    implicit none
    contains
    subroutine  grid_conditions_wall (&
        xp, yp, zp, dx, dy, dz, dt, xnue, xlambda, density, width, height, depth, &
        thickness, threshold, radius,&
        center_x, center_y, center_z, time,&
        inlet_velocity, AoA, porosity, m, n, l, istep_max, &
        csv_file)
    
        use global_3d
        implicit none
        real,intent(inout),dimension(0:md):: xp
        real,intent(inout),dimension(0:nd):: yp
        real,intent(inout),dimension(0:ld):: zp
        real,intent(inout):: dx, dy, dz, dt
        real,intent(in):: xnue, xlambda, density, width, height, depth, time
        real,intent(in):: inlet_velocity, AoA
        real,intent(in):: thickness, threshold, radius, center_x, center_y, center_z
        real,intent(inout),dimension(0:md,0:nd,0:ld)::porosity
        integer,intent(inout):: m, n, l
        integer,intent(in):: istep_max
        character(len = 50) :: csv_file
        character(len = 50) :: output_folder
    
        ! local variables
        !real,dimension(0:md,0:nd)::	distance
        real:: cfl_no, pecret_no, diffusion_factor, reynolds_no
        integer:: i, j, k
        integer:: x, y, z
        real:: poro_val
        ! --- 
        
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
    
        !--- calc.
        dx = width / real(m-1)
        dy = height / real(n-1)
        dz = depth / real(l-1)
        dt = time / real(istep_max)
        
        cfl_no           = inlet_velocity * dt / dx
        pecret_no        = inlet_velocity * dx / xnue
        diffusion_factor = xnue * dt / dy / dy
    
        reynolds_no      = inlet_velocity * radius * 2.0 / xnue
    
        !----- check print out
        write(*,*)
        write(*,*) '# --- Grid conditions'
        write(*,*) '# m, n, l =', m, n, l
        write(*,*) '# istep_max =', istep_max
        write(*,*) '# dx, dy, dz =', dx, dy, dz
        write(*,*) '# dt =', dt
        write(*,*) '# cfl_no =', cfl_no
        write(*,*) '# pecret_no =', pecret_no
        write(*,*) '# diffusion_factor =', diffusion_factor
        write(*,*) '# reynolds_no =', reynolds_no
        write(*,*) '# thickness =', thickness
        write(*,*) '# threshold =', threshold
        write(*,*)
        
        !$omp parallel private(i, j, k) &
        !$omp & shared(m, n, l) &
        !$omp & shared(porosity) &
        !$omp & shared(xp, yp, zp) &
        !$omp & shared(dx, dy, dz) &
        !$omp & shared(width, height, depth, center_x, center_y, center_z) &
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
    
        !$omp do      
        do k = 0, l+1
        zp(k) = dz * real(k-1) - depth*center_z
        end do
        !$omp end do
    
        ! default: outlet condtion in x-direction
        !$omp do
        do j = 0, n+1
        do k = 0, l+1
            porosity(0,j,k) = porosity(1,j,k)
            porosity(m+1,j,k) = porosity(m,j,k)
        end do
        end do
        !$omp end do
    
        ! default: outlet condtion in y-direction
        !$omp do      
        do i = 0, m+1
        do k = 0, l+1
            porosity(i,0,k) = porosity(i,1,k)
            porosity(i,n+1,k) = porosity(i,n,k)
        end do
        end do
        !$omp end do
    
        ! default: outlet condtion in z-direction
        !$omp do      
        do i = 0, m+1
        do j = 0, n+1
            porosity(i,j,0) = porosity(i,j,1)
            porosity(i,j,l+1) = porosity(i,j,l)
        end do
        end do
        !$omp end do
    
        !$omp end parallel
        ! ----------------
        return
    end subroutine  grid_conditions_wall

    subroutine  grid_conditions_yz_periodic (&
        xp, yp, zp, dx, dy, dz, dt, xnue, xlambda, density, width, height, depth, &
        thickness, threshold, radius,&
        center_x, center_y, center_z, time,&
        inlet_velocity, AoA, porosity, m, n, l, istep_max, &
        csv_file)
      
        use global_3d
        implicit none
        real,intent(inout),dimension(0:md):: xp
        real,intent(inout),dimension(0:nd):: yp
        real,intent(inout),dimension(0:ld):: zp
        real,intent(inout):: dx, dy, dz, dt
        real,intent(in):: xnue, xlambda, density, width, height, depth, time
        real,intent(in):: inlet_velocity, AoA
        real,intent(in):: thickness, threshold, radius, center_x, center_y, center_z
        real,intent(inout),dimension(0:md,0:nd,0:ld)::porosity
        integer,intent(inout):: m, n, l
        integer,intent(in):: istep_max
        character(len = 50) :: csv_file
        character(len = 50) :: output_folder
      
        ! local variables
        !real,dimension(0:md,0:nd)::	distance
        real:: cfl_no, pecret_no, diffusion_factor, reynolds_no
        integer:: i, j, k
        integer:: x, y, z
        real:: poro_val
        ! --- 
        
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
      
        !--- calc.
        dx = width / real(m-1)
        dy = height / real(n-1)
        dz = depth / real(l-1)
        dt = time / real(istep_max)
        
        cfl_no           = inlet_velocity * dt / dx
        pecret_no        = inlet_velocity * dx / xnue
        diffusion_factor = xnue * dt / dy / dy
      
        reynolds_no      = inlet_velocity * radius * 2.0 / xnue
      
        !----- check print out
        write(*,*)
        write(*,*) '# --- Grid conditions'
        write(*,*) '# m, n, l =', m, n, l
        write(*,*) '# istep_max =', istep_max
        write(*,*) '# dx, dy, dz =', dx, dy, dz
        write(*,*) '# dt =', dt
        write(*,*) '# cfl_no =', cfl_no
        write(*,*) '# pecret_no =', pecret_no
        write(*,*) '# diffusion_factor =', diffusion_factor
        write(*,*) '# reynolds_no =', reynolds_no
        write(*,*) '# thickness =', thickness
        write(*,*) '# threshold =', threshold
        write(*,*)
        
        !$omp parallel private(i, j, k) &
        !$omp & shared(m, n, l) &
        !$omp & shared(porosity) &
        !$omp & shared(xp, yp, zp) &
        !$omp & shared(dx, dy, dz) &
        !$omp & shared(width, height, depth, center_x, center_y, center_z) &
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
      
        !$omp do      
        do k = 0, l+1
        zp(k) = dz * real(k-1) - depth*center_z
        end do
        !$omp end do
      
        ! default: outlet condtion in x-direction
        !$omp do      
        do j = 1, n+1
        do k = 1, l+1
        porosity(0,j,k) = porosity(1,j,k)
        porosity(m+1,j,k) = porosity(m,j,k)
        end do
        end do
        !$omp end do
      
        ! default: periodic condtion in y-direction
        !$omp do      
        do i = 0, m+1
        do k = 0, l+1
        porosity(i,0,k)   = porosity(i,n,k)
        porosity(i,n+1,k) = porosity(i,1,k)
        end do
        end do
        !$omp end do
        !$omp end parallel
      
        ! default: periodic condtion in z-direction
        !$omp do      
        do i = 0, m+1
        do j = 0, n+1
        porosity(i,j,0)   = porosity(i,j,l)
        porosity(i,j,l+1) = porosity(i,j,1)
        end do
        end do
        !$omp end do
        !$omp end parallel
        ! ----------------
        return
      end subroutine  grid_conditions_yz_periodic
end module grid_3d