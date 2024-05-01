! src/omp_parallel/lib/utils.f90

module utils
    implicit none
    contains
    
    subroutine get_now_time()
        implicit none
        character(len=20) :: current_time
        integer ::values(8)
        call date_and_time(values=values)
        write(current_time, '(I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2,":",I2.2)') &
            values(1), values(2), values(3), values(5), values(6), values(7)
        write(*,*) '# --- TIME: ', trim(current_time)
        return
    end subroutine get_now_time
  
end module utils
