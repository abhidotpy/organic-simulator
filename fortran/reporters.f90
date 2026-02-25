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

        write(unit, "(I5, 100F20.10)") &
        step, potential_energy, kinetic_energy, total_energy, temperature, pressure, density
    
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
        write(unit, "(1x, 2f20.10)") -box_length(1) / 2.0, box_length(1) / 2.0
        write(unit, "(1x, 2f20.10)") -box_length(2) / 2.0, box_length(2) / 2.0
        write(unit, "(1x, 2f20.10)") -box_length(3) / 2.0, box_length(3) / 2.0
        write(unit, "(A)") "ITEM: ATOMS id type shapex shapey shapez x y z vx vy vz fx fy fz quatw quati quatj quatk"

        do I = 1, N
            write(unit, "(1x, 2i5, 50f20.10)") I, rtype(I), shape_x(I), shape_y(I), shape_z(I), RX(I), RY(I), RZ(I), VX(I), VY(I), VZ(I), FX(I), FY(I), FZ(I), &
            QW(I), QX(I), QY(I), QZ(I)
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
            write(unit, "(1x, 50f20.10)") RX(I), RY(I), RZ(I), VX(I), VY(I), VZ(I), LX(I), LY(I), LZ(I), QW(I), QX(I), QY(I), QZ(I)
        enddo

    end subroutine report_confwriter

    subroutine close_confwriter_file( unit )
        implicit none
        integer, intent(in) :: unit

        close(unit=unit)

    end subroutine close_confwriter_file
    
end module config_reporter

module checkpoint_reporter
    use system
    implicit none
    character(len=256) :: file
    integer :: config_index = 0

    contains
    subroutine open_checkpoint_file( unit, filename )
        implicit none
        integer, intent(in) :: unit
        character(len=*), intent(in) :: filename

        file = filename

        open(unit=unit, file=file)
        config_index = 0

    end subroutine open_checkpoint_file

    subroutine report_checkpoint( unit )
        implicit none
        integer, intent(in) :: unit
        integer :: I

        write(unit, "('ATOMS ', I0)") N
        write(unit, *)
        write(unit, "('BONDS ', I0)") NC
        write(unit, *)
        write(unit, "('BOX ', 3F10.5)") box_length(1), box_length(2), box_length(3)
        write(unit, *)
        write(unit, "('ATOM: ID TYPE SHAPEX SHAPEY SHAPEZ ESHAPEX ESHAPEY ESHAPEZ RX RY RZ VX VY VZ LX LY LZ QW QX QY QZ')")
        write(unit, *)
        do I = 1, N
            write(unit, "(1x, 2I10, 50F20.10)") I, rtype(I), shape_x(I), shape_y(I), shape_z(I), eshape_x(I), eshape_y(I), eshape_z(I), &
                                          RX(I), RY(I), RZ(I), VX(I), VY(I), VZ(I), LX(I), LY(I), LZ(I), QW(I), QX(I), QY(I), QZ(I)
        enddo
        write(unit, *)
        write(unit, "('BOND: ID I J LEN')")
        write(unit, *)
        do I = 1, NC
            write(unit, "(1x, 3I10, 50F20.10)") I, BBI(I), BBJ(I), BB_len(I)
        enddo
        write(unit, *)
        write(unit, *) config_index

    end subroutine report_checkpoint

    subroutine load_checkpoint( unit )
        implicit none
        integer, intent(in) :: unit
        integer :: I, J, na, nb, ia, ja
        real(8) :: la, bx, by, bz
        character(len=256) :: label

        read(unit, *) label, na
        call set_num_atoms(na)
        read(unit, *)
        read(unit, *) label, nb
        call set_num_constraints(nb)
        read(unit, *)
        read(unit, *) label, bx, by, bz
        call set_box( bx, by, bz )
        read(unit, *)
        read(unit, *) label
        read(unit, *)
        do I = 1, N
            read(unit, *) J, rtype(I), shape_x(I), shape_y(I), shape_z(I), eshape_x(I), eshape_y(I), eshape_z(I), &                  
                                          RX(I), RY(I), RZ(I), VX(I), VY(I), VZ(I), LX(I), LY(I), LZ(I), QW(I), QX(I), QY(I), QZ(I)
        enddo
        read(unit, *)
        read(unit, *) label
        read(unit, *)
        do I = 1, NC
            read(unit, *) J, ia, ja, la
            call add_constraint(J, ia, ja, la)
        enddo
        read(unit, *)
        read(unit, *) config_index

        default_shape_set = .FALSE.
        default_energy_set = .FALSE.
        box_lengths_set = .TRUE.
        
    end subroutine load_checkpoint

    subroutine close_checkpoint_file( unit )
        implicit none
        integer, intent(in) :: unit

        close(unit=unit)

    end subroutine close_checkpoint_file
    
end module checkpoint_reporter

module reporters
    use state_data_reporter
    use trajectory_reporter
    use config_reporter
    use checkpoint_reporter
end module reporters