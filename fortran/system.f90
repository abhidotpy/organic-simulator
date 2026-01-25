module atom
    implicit none

    real(8) :: mass, charge

    real(8), allocatable, dimension(:) :: RX, RY, RZ
    real(8), allocatable, dimension(:) :: VX, VY, VZ
    real(8), allocatable, dimension(:) :: FX, FY, FZ

    real(8), allocatable, dimension(:) :: QW, QX, QY, QZ
    real(8), allocatable, dimension(:) :: QDW, QDX, QDY, QDZ
    real(8), allocatable, dimension(:) :: WX, WY, WZ
    real(8), allocatable, dimension(:) :: LX, LY, LZ
    real(8), allocatable, dimension(:) :: TX, TY, TZ

end module atom

module box
    implicit none

    real(8) :: box_length, volume
    real(8) :: box_xlen, box_ylen, box_zlen
    logical :: is_periodic

end module box

module system
    use atom
    use box
    use xmath
    implicit none

    integer :: N, dof
    
    real(8) :: potential_energy, kinetic_energy, total_energy
    real(8) :: temperature, pressure, density, virial
    
    logical :: default_mass = .TRUE., default_charge=.TRUE.
    logical ::rotaions_enabled = .FALSE.
    logical :: box_lengths_set = .FALSE.
    logical :: num_atoms_set = .FALSE.
    logical :: density_set = .FALSE.
    logical :: periodic_set = .FALSE.
    logical :: fcc_lattice_generated = .FALSE.

    contains

    subroutine set_num_atoms( num_atoms )
        implicit none
        integer, intent(in) :: num_atoms

        N = num_atoms
        dof = 3 * N
        num_atoms_set = .TRUE.

        allocate( RX(N), RY(N), RZ(N) )
        allocate( VX(N), VY(N), VZ(N) )
        allocate( FX(N), FY(N), FZ(N) )

        RX = 0.0; RY = 0.0; RZ = 0.0
        VX = 0.0; VY = 0.0; VZ = 0.0
        FX = 0.0; FY = 0.0; FZ = 0.0

    end subroutine set_num_atoms

    subroutine set_density( target_density )
        implicit none
        
        real(8), intent(in) :: target_density
        density_set = .TRUE.
        density = target_density

    end subroutine set_density

    subroutine set_periodic_boundary
        implicit none

        periodic_set = .TRUE.
    end subroutine set_periodic_boundary

    subroutine add_atom_position(id, pos_x, pos_y, pos_z)
        implicit none

        integer, intent(in) :: id
        real(8), intent(in) :: pos_x, pos_y, pos_z

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "(1x, 'Error: System has no atoms. Set number of atoms first.')")
            stop
        else if ( id < 1 .or. id > N ) then
            write(*, "(1x, 'Error: Atom index ', I10,' does not exist &
            in system of ',I10,' atoms.')") id, N
        else
            RX(id) = pos_x
            RY(id) = pos_y
            RZ(id) = pos_z
        endif

    end subroutine add_atom_position

    subroutine add_atom_velocity(id, vel_x, vel_y, vel_z)
        implicit none

        integer, intent(in) :: id
        real(8), intent(in) :: vel_x, vel_y, vel_z

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "(1x, 'Error: System has no atoms. Set number of atoms first.')")
            stop
        else if ( id < 1 .or. id > N ) then
            write(*, "(1x, 'Error: Atom index ', I10,' does not exist &
            in system of ',I10,' atoms.')") id, N
        else
            VX(id) = vel_x
            VY(id) = vel_y
            VZ(id) = vel_z
        endif

    end subroutine add_atom_velocity

    subroutine set_velocty_to_temperature( rtemp )
        implicit none

        real(8), intent(in) :: rtemp
        real(8) :: target_velocity
        real(8) :: vx_gen, vy_gen, vz_gen
        real(8) :: vx_com, vy_com, vz_com
        integer :: i

        if ( num_atoms_set .eqv. .FALSE. ) then
            write(*, "(1x, 'Error: System has no atoms. Set number of atoms first.')")
            stop
        endif

        vx_gen = 0.0; vy_gen = 0.0; vz_gen = 0.0
        vx_com = 0.0; vy_com = 0.0; vz_com = 0.0

        call random_seed()

        do I = 1, N
            call normal_distribution( vx_gen )
            call normal_distribution( vy_gen )
            call normal_distribution( vz_gen )

            VX(I) = sqrt( rtemp ) * vx_gen
            VY(I) = sqrt( rtemp ) * vy_gen
            VZ(I) = sqrt( rtemp ) * vz_gen
        enddo

        vx_com = sum(vx) / dble(N)
        vy_com = sum(vy) / dble(N)
        vz_com = sum(vz) / dble(N)

        VX = VX - vx_com
        VY = VY - vy_com
        VZ = VZ - vz_com

    end subroutine set_velocty_to_temperature

    subroutine set_rotational_dof()
        implicit none

        dof = dof + 3 * N

        allocate( QW(N), QX(N), QY(N), QZ(N) )
        allocate( QDW(N), QDX(N), QDY(N), QDZ(N) )
        allocate( WX(N), WY(N), WZ(N) )
        allocate( LX(N), LY(N), LZ(N) )
        allocate( TX(N), TY(N), TZ(N) )

        QW = 1.0; QX = 0.0; QY = 0.0; QZ = 0.0
        QDW = 0.0; QDX = 0.0; QDY = 0.0; QDZ = 0.0
        WX = 0.0; WY = 0.0; WZ = 0.0
        LX = 0.0; LY = 0.0; LZ = 0.0
        TX = 0.0; TY = 0.0; TZ = 0.0

    end subroutine set_rotational_dof

    subroutine intiialize_forces
        implicit none

        FX = 0.0; FY = 0.0; FZ = 0.0
        TX = 0.0; TY = 0.0; TZ = 0.0

        potential_energy = 0.0
        kinetic_energy = 0.0
        temperature = 0.0
        pressure = 0.0
        virial = 0.0
    end subroutine intiialize_forces

    subroutine generate_fcc_lattice( unit_cells )
        implicit none
        integer, intent(in) :: unit_cells
        integer :: num_atoms, mtemp, I, J, K, IREF
        real    :: cell, half_cell

        fcc_lattice_generated = .TRUE.
        num_atoms = 4 * unit_cells ** 3
        call set_num_atoms( num_atoms )

        if ( box_lengths_set .eqv. .FALSE. ) then

            if (density_set .eqv. .FALSE.) then
                write(*, "(1x, 'Error: Density not set. Set density before &
                generating FCC lattice.')")
                stop
            endif

            box_lengths_set = .TRUE.
            volume = real(N) / density
            box_length = volume ** (1.0/3.0)
            box_xlen = box_length
            box_ylen = box_length
            box_zlen = box_length
        endif

        cell = box_length / real(unit_cells)
        half_cell = cell / 2.0

        RX(1) = 0.0
        RY(1) = 0.0
        RZ(1) = 0.0

        RX(2) = half_cell
        RY(2) = half_cell
        RZ(2) = 0.0

        RX(3) = 0.0
        RY(3) = half_cell
        RZ(3) = half_cell

        RX(4) = half_cell
        RY(4) = 0.0
        RZ(4) = half_cell

        mtemp = 0
        do I = 1, unit_cells
            do J = 1, unit_cells
                do K = 1, unit_cells
                    do iref = 1, 4

                        RX( IREF + MTEMP ) = RX( IREF ) + cell * ( I - 1 )
                        RY( IREF + MTEMP ) = RY( IREF ) + cell * ( J - 1 )
                        RZ( IREF + MTEMP ) = RZ( IREF ) + cell * ( K - 1 )

                    enddo
                    mtemp = mtemp + 4
                enddo
            enddo
        enddo

        RX = RX - box_length / 2.0
        RY = RY - box_length / 2.0
        RZ = RZ - box_length / 2.0

    end subroutine generate_fcc_lattice

    subroutine initialize_system( )
        implicit none

        write(*, *)
        write(*, "(1x, 'This is SUCROSE v1.0')")

        if (num_atoms_set) then
            write(*, "(1x, 'Number of particles in system set to ', i10)") N
        else
            write(*, "(1x, 'Error: Number of particles not set. Set number of &
            particles before initializing system.')")
            stop
        endif

        if (density_set) then
            write(*, "(1x, 'Density of system set to ', f10.5, ' atoms per cubic unit cell')") density
        else
            write(*, "(1x, 'Error: Density not set. Set density before initializing system.')")
            stop
        endif

        if (box_lengths_set) then
            write(*, "(1x, 'Box length set to ', f10.5)") box_length
        else
            box_lengths_set = .TRUE.
            volume = real(N) / density
            box_length = volume ** (1.0/3.0)
            box_xlen = box_length
            box_ylen = box_length
            box_zlen = box_length
        endif

        if (periodic_set) then
            write(*, "(1x, 'Periodic boundary conditions enabled for current system')")
        else
            write(*, "(1x, 'Periodic boundary conditions disabled for current system')")
        endif
        
        if (rotaions_enabled) then
            write(*, "(1x, 'Rotational degrees of freedom enabled for current system')")
        else
            write(*, "(1x, 'Rotational degrees of freedom disabled for current system')")
        endif

        if (default_mass) then
            mass = 1.0
        endif
        if (default_charge) then
            charge = 0.0
        endif

        ! write(*, "(1x, 'Mass of particles set to ', f10.5)") mass
        ! write(*, "(1x, 'Charge of particles set to ', f10.5)") charge

        call intiialize_forces

        if (fcc_lattice_generated) then
            write(*, "(1x, 'FCC lattice generated with ', i10, ' particles')") N
        endif

        write(*, "(1x, 'System initialization complete.')")
        write(*, *)

    end subroutine initialize_system

    subroutine calculate_kinetic_energy()
        implicit none
        integer :: I

        kinetic_energy = 0.5 * mass * SUM( VX ** 2 + VY ** 2 + VZ ** 2 )
    end subroutine calculate_kinetic_energy

    subroutine calculate_temperature()
        implicit none
    
        temperature = 0.0
        if (kinetic_energy > 0.0) then
            temperature = 2.0 * kinetic_energy / dble(dof)
        endif
    
    end subroutine calculate_temperature

    subroutine calculate_pressure()
        implicit none

        if (box_lengths_set .and. volume > 0.0 .and. temperature > 0.0) then
            pressure = (density * temperature) + (virial / volume)
        endif
    
    end subroutine calculate_pressure

    subroutine calculate_state_variables()
        implicit none

        kinetic_energy = 0.0
        temperature = 0.0
        pressure = 0.0

        kinetic_energy = 0.5 * mass * SUM( VX ** 2 + VY ** 2 + VZ ** 2 )
        total_energy = potential_energy + kinetic_energy

        temperature = 2.0 * kinetic_energy / dble(dof)

        if (volume > 0.0) then
            pressure = (density * temperature) + (virial / (3.0 * volume) )
        endif

        potential_energy = potential_energy / dble(N)
        kinetic_energy = kinetic_energy / dble(N)
        total_energy = total_energy / dble(N)


    end subroutine calculate_state_variables

end module system