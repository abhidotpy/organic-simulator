module state_data_reporter
    use system
    implicit none
    character(len=256) :: file

    contains
    subroutine open_state_data_file( unit, filename )
        implicit none
        integer, intent(in) :: unit
        character(len=*), intent(in) :: filename

        file = filename

        open(unit=unit, file=file, action='write')

        write(unit, "('Step', 5x, 'Potential_Energy', 5x, 'Kinetic_Energy', 5x, 'Total_Energy', 5x, 'Temperature', 5x, 'Pressure', 5x, 'Density')")

    end subroutine open_state_data_file

    subroutine report_state_data( unit, step )
        implicit none
        integer, intent(in) :: unit, step

        write(unit, "(I5, 6(F20.10))") &
        step, potential_energy / dble(N), kinetic_energy / dble(N), &
        (potential_energy + kinetic_energy) / dble(N), temperature, pressure, density
    
    end subroutine report_state_data

    subroutine close_state_data_file( unit )
        implicit none
        integer, intent(in) :: unit

        close(unit=unit)

    end subroutine close_state_data_file

end module state_data_reporter

module trajectory_reporter
    use system
    implicit none
    character(len=256) :: file

    contains
    subroutine open_trajectory_file( unit, filename )
        implicit none
        integer, intent(in) :: unit
        character(len=*), intent(in) :: filename

        file = filename

        open(unit=unit, file=file, action='write')

    end subroutine open_trajectory_file

    subroutine report_trajectory( unit, step )
        implicit none
        integer, intent(in) :: unit, step
        integer :: I

        write(unit, '(A)') "ITEM: TIMESTEP"
        write(unit, *) step
        write(unit, '(A)') "ITEM: NUMBER OF ATOMS"
        write(unit, *) n
        write(unit, "(A)") "ITEM: BOX BOUNDS pp pp pp"
        write(unit, "(1x, 2f20.10)") -box_length / 2.0, box_length / 2.0
        write(unit, "(1x, 2f20.10)") -box_length / 2.0, box_length / 2.0
        write(unit, "(1x, 2f20.10)") -box_length / 2.0, box_length / 2.0
        write(unit, "(A)") "ITEM: ATOMS id radius x y z vx vy vz fx fy fz"

        do I = 1, N
            write(unit, "(1x, i5, 50f20.10)") I, 0.5, RX(I), RY(I), RZ(I), VX(I), VY(I), VZ(I), FX(I), FY(I), FZ(I)
        enddo

    end subroutine report_trajectory

    subroutine close_trajectory_file( unit )
        implicit none
        integer, intent(in) :: unit

        close(unit=unit)

    end subroutine close_trajectory_file

end module trajectory_reporter

module config_reporter
    use system
    implicit none
    character(len=256) :: file

    contains
    subroutine open_confwriter_file( unit, filename )
        implicit none
        integer, intent(in) :: unit
        character(len=*), intent(in) :: filename

        file = filename

        open(unit=unit, file=file, action='write')

    end subroutine open_confwriter_file

    subroutine report_confwriter( unit )
        implicit none
        integer, intent(in) :: unit
        integer :: I

        do I = 1, N
            write(unit, "(1x, 50f20.10)") RX(I), RY(I), RZ(I), VX(I), VY(I), VZ(I)
        enddo

    end subroutine report_confwriter

    subroutine close_confwriter_file( unit )
        implicit none
        integer, intent(in) :: unit

        close(unit=unit)

    end subroutine close_confwriter_file
    
end module config_reporter