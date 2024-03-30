module scalar_transport
    implicit none
    contains
    !******************
    subroutine solve_T(T, dx, dy, dt, alpha, m, n)
    implicit none
    real, dimension(0:m+1, 0:n+1), intent(inout) :: T
    real, intent(in) :: dx, dy, dt, alpha
    integer, intent(in) :: m, n
    real, dimension(0:m+1, 0:n+1) :: T_new
    integer :: i, j
    
    ! Calculate the solution at the next time step
    do i = 1, m
        do j = 1, n
            T_new(i,j) = T(i,j) + alpha * dt * ((u(i+1,j) - 2*u(i,j) + u(i-1,j)) / dx**2 + &
                                                    (T(i,j+1) - 2*u(i,j) + u(i,j-1)) / dy**2)
        end do
    end do
    
    ! Update the solution for the next time step
    T = T_new
    end subroutine solve_T
    !******************
end module scalar_transport