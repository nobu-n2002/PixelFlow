! global.f90

module global_2d
    implicit none
    !--- md, nd > grid size (m,n)
    integer, parameter:: md = 1500, nd = 1500
end module global_2d

module global_3d
    implicit none
    !--- md, nd, ld > grid size (m,n,l)
    integer, parameter:: md = 256, nd = 180, ld = 80
end module global_3d

module valiables
    implicit none
    !--- namelist valiables
    real:: xnue, xlambda, density, width, height, depth, time
    real:: inlet_velocity, outlet_pressure, AoA
    integer:: istep_max, istep_out
    real:: thickness, threshold, radius, center_x, center_y, center_z
    logical:: nonslip
    character(len=50) :: output_folder
    character(len=50) :: csv_file
    integer:: iter_max
    real:: relux_factor
    contains
    subroutine read_settings(&
        xnue, xlambda, density, width, height, depth, time,&
        inlet_velocity, outlet_pressure, AoA,&
        istep_max, istep_out,&
        thickness, threshold, radius, center_x, center_y, center_z,&
        nonslip,&
        output_folder,csv_file,&
        iter_max, relux_factor)
      
        real, intent(out):: xnue, xlambda, density, width, height, depth, time
        real, intent(out):: inlet_velocity, outlet_pressure, AoA
        integer, intent(out):: istep_max, istep_out
        real, intent(out):: thickness, threshold, radius, center_x, center_y, center_z
        logical, intent(out):: nonslip
        character(len=50), intent(out) :: output_folder
        character(len=50), intent(out) :: csv_file
        integer, intent(out):: iter_max
        real, intent(out):: relux_factor
      
        namelist /physical/xnue, xlambda, density, width, height, depth, time
        namelist /physical/inlet_velocity, outlet_pressure, AoA
        namelist /file_control/istep_out
        namelist /grid_control/istep_max
        namelist /porosity_control/thickness, threshold, radius, center_x, center_y, center_z
        namelist /calculation_method/nonslip
        namelist /directory_control/output_folder, csv_file
        namelist /solver_control/iter_max, relux_factor
        open(11,file="config/controlDict.txt",status="old",action="read")
        read(11,nml=physical)
        read(11,nml=file_control)
        read(11,nml=grid_control)
        read(11,nml=porosity_control)
        read(11,nml=calculation_method)
        read(11,nml=directory_control)
        read(11,nml=solver_control)
        close(11)
      
        !--- check
        write(*,*) '#'
        write(*,*) '# --- Physical conditions'
        write(*,*) '# xnue =', xnue
        write(*,*) '# xlambda =', xlambda
        write(*,*) '# density =', density
        write(*,*) '# width =', width
        write(*,*) '# height =', height
        write(*,*) '# depth =', depth
        write(*,*) '# time =', time
        write(*,*) '# inlet_velocity =', inlet_velocity
        write(*,*) '# outlet_pressure =', outlet_pressure
        write(*,*) '# Angle of inlet_velocity (AoA) =', AoA
        write(*,*) '#'
        write(*,*) '# --- Porosity information'
        write(*,*) '# thickness =', thickness
        write(*,*) '# threshold =', threshold
        write(*,*) '# radius =', radius
        write(*,*) '#'
        write(*,*) '# --- Directory information'
        write(*,*) '# output_folder =', output_folder  
        write(*,*) '# input_porosity_file =', csv_file  
        write(*,*) '#'
        write(*,*) '# --- Solver information'
        write(*,*) '# SOR max iteration steps =', iter_max
        write(*,*) '# SOR reluxation factor =', relux_factor
        return
    end subroutine read_settings
end module valiables
