module atom
    implicit none

    real(8) :: mass, charge
    real(8), dimension(3) :: rad, erad, exc_rad, MOI

    real(8), allocatable, dimension(:) :: RX, RY, RZ
    real(8), allocatable, dimension(:) :: VX, VY, VZ
    real(8), allocatable, dimension(:) :: FX, FY, FZ

    real(8), allocatable, dimension(:) :: QW, QX, QY, QZ
    real(8), allocatable, dimension(:) :: WX, WY, WZ
    real(8), allocatable, dimension(:) :: LX, LY, LZ
    real(8), allocatable, dimension(:) :: TX, TY, TZ

    real(8), allocatable, dimension(:) :: RX_old, RY_old, RZ_old

end module atom

module box
    implicit none

    real(8) :: box_length(3), volume
    logical :: is_periodic

end module box

module system
    use atom
    use box
    use xmath
    implicit none

    integer :: N, NC, dof
    
    real(8) :: potential_energy, kinetic_energy, total_energy
    real(8) :: temperature, pressure, density, virial
    real(8) :: target_temp, target_pres

    integer, dimension(:), allocatable :: BBI, BBJ
    real(8), dimension(:), allocatable :: BB_len
    
    logical :: default_mass_set     = .TRUE.
    logical :: default_charge_set   = .TRUE.
    logical :: rotations_set        = .FALSE.
    logical :: box_lengths_set      = .FALSE.
    logical :: num_atoms_set        = .FALSE.
    logical :: density_set          = .FALSE.
    logical :: fcc_lattice_set      = .FALSE.
    logical :: constraints_set      = .FALSE.
    logical :: target_temp_set      = .FALSE.
    logical :: target_pres_set      = .FALSE.

    contains

    subroutine set_num_atoms( num_atoms )
        implicit none
        integer, intent(in) :: num_atoms

        N = num_atoms
        dof = 6 * N
        num_atoms_set = .TRUE.

        allocate( RX(N), RY(N), RZ(N) )
        allocate( VX(N), VY(N), VZ(N) )
        allocate( FX(N), FY(N), FZ(N) )
        allocate( LX(N), LY(N), LZ(N) )
        allocate( TX(N), TY(N), TZ(N) )
        allocate( QW(N), QX(N), QY(N), QZ(N) )
        allocate( WX(N), WY(N), WZ(N) )

        RX = 0.0; RY = 0.0; RZ = 0.0
        VX = 0.0; VY = 0.0; VZ = 0.0
        FX = 0.0; FY = 0.0; FZ = 0.0
        LX = 0.0; LY = 0.0; LZ = 0.0
        TX = 0.0; TY = 0.0; TZ = 0.0
        QW = 1.0; QX = 0.0; QY = 0.0; QZ = 0.0
        Wx = 0.0; Wy = 0.0; Wz = 0.0

    end subroutine set_num_atoms

    subroutine set_num_constraints( num_constraints )
        implicit none
        integer, intent(in) :: num_constraints

        constraints_set = .TRUE.
        NC = num_constraints
        dof = dof - NC

        allocate( BBI(NC), BBJ(NC), BB_len(NC) )
        BBI = 0; BBJ = 0; BB_len = 0.0

    end subroutine set_num_constraints

    subroutine set_density( target_density )
        implicit none
        
        real(8), intent(in) :: target_density
        density_set = .TRUE.
        density = target_density

    end subroutine set_density

    subroutine set_box( length, width, height )
        implicit none
        real(8), intent(in) :: length
        real(8), optional, intent(in) :: width, height

        box_lengths_set = .TRUE.

        if ( present(width) .and. present(height) ) then
            box_length(1) = length
            box_length(2) = width
            box_length(3) = height
        else
            box_length(1) = length
            box_length(2) = length
            box_length(3) = length
        endif

    end subroutine set_box

    subroutine set_periodic_boundary
        implicit none

        is_periodic = .TRUE.
    end subroutine set_periodic_boundary

    subroutine set_target_temperature( temp )
        implicit none
        real(8), intent(in) :: temp

        target_temp = temp
        target_temp_set = .TRUE.

    end subroutine set_target_temperature

    subroutine set_target_pressure( pres )
        implicit none
        real(8), intent(in) :: pres

        target_pres = pres
        target_pres_set = .TRUE.

    end subroutine set_target_pressure

    subroutine set_atom_shape( length, width, height )
        implicit none
        real(8), intent(in) :: length, width, height

        rad(1) = length
        rad(2) = width
        rad(3) = height

        MOI(1) = (rad(2)**2 + rad(3)**2) / 20.0
        MOI(2) = (rad(1)**2 + rad(3)**2) / 20.0
        MOI(3) = (rad(1)**2 + rad(2)**2) / 20.0

        exc_rad = max( length, width, height )
    
    end subroutine set_atom_shape

    subroutine set_atom_energy( eng_x, eng_y, eng_z )
        implicit none
        real(8), intent(in) :: eng_x, eng_y, eng_z

        erad(1) = eng_x
        erad(2) = eng_y
        erad(3) = eng_z
        
    end subroutine set_atom_energy

    subroutine set_atom_position(id, pos_x, pos_y, pos_z)
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

    end subroutine set_atom_position

    subroutine set_atom_velocity(id, vel_x, vel_y, vel_z)
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

    end subroutine set_atom_velocity

    subroutine add_constraint(id, I, J, length)
        implicit none
        integer, intent(in) :: id
        integer, intent(in) :: I, J
        real(8), intent(in) :: length

        if (constraints_set .eqv. .FALSE.) then
            write(*, "(1x, 'Error: System has no constraints. Set number of constraints first.')")
            stop
        else if ( id < 1 .or. id > NC ) then
            write(*, "(1x, 'Error: Constraint index ', I10,' does not exist &
            in system of ',I10,' constraints.')") id, N
        else
            BBI(id) = I
            BBJ(id) = J
            BB_len(id) = length
        endif

    end subroutine add_constraint

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

        dof = dof - 3

    end subroutine set_velocty_to_temperature

    subroutine initialize_forces
        implicit none

        FX = 0.0; FY = 0.0; FZ = 0.0
        TX = 0.0; TY = 0.0; TZ = 0.0

        potential_energy = 0.0
        kinetic_energy = 0.0
        temperature = 0.0
        pressure = 0.0
        virial = 0.0
    end subroutine initialize_forces

    subroutine generate_fcc_lattice( unit_cells )
        implicit none
        integer, intent(in) :: unit_cells
        integer :: num_atoms, mtemp, I, J, K, IREF
        real(8)    :: cell, half_cell, rroot3

        fcc_lattice_set = .TRUE.
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
        endif

        cell = box_length(3) / real(unit_cells)
        half_cell = cell / 2.0
        rroot3 = 1.0 / sqrt(3.0)

        RX(1) = 0.0
        RY(1) = 0.0
        RZ(1) = 0.0
        QW(1) = sqrt( ( 1.0 + rroot3 ) / 2.0 )
        QX(1) = sqrt( ( 1.0 - rroot3 ) / 2.0 ) * (  rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QY(1) = sqrt( ( 1.0 - rroot3 ) / 2.0 ) * ( -rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QZ(1) = 0.0

        RX(2) = half_cell
        RY(2) = half_cell
        RZ(2) = 0.0
        QW(2) = sqrt( ( 1.0 - rroot3 ) / 2.0 )
        QX(2) = sqrt( ( 1.0 + rroot3 ) / 2.0 ) * ( -rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QY(2) = sqrt( ( 1.0 + rroot3 ) / 2.0 ) * ( -rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QZ(2) = 0.0

        RX(3) = 0.0
        RY(3) = half_cell
        RZ(3) = half_cell
        QW(3) = sqrt( ( 1.0 - rroot3 ) / 2.0 )
        QX(3) = sqrt( ( 1.0 + rroot3 ) / 2.0 ) * (  rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QY(3) = sqrt( ( 1.0 + rroot3 ) / 2.0 ) * (  rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QZ(3) = 0.0

        RX(4) = half_cell
        RY(4) = 0.0
        RZ(4) = half_cell
        QW(4) = sqrt( ( 1.0 + rroot3 ) / 2.0 )
        QX(4) = sqrt( ( 1.0 - rroot3 ) / 2.0 ) * ( -rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QY(4) = sqrt( ( 1.0 - rroot3 ) / 2.0 ) * (  rroot3 / sqrt( 1.0 - rroot3 ** 2.0 ) )
        QZ(4) = 0.0

        mtemp = 0
        do I = 1, unit_cells
            do J = 1, unit_cells
                do K = 1, unit_cells
                    do iref = 1, 4

                        RX( IREF + MTEMP ) = RX( IREF ) + cell * ( I - 1 )
                        RY( IREF + MTEMP ) = RY( IREF ) + cell * ( J - 1 )
                        RZ( IREF + MTEMP ) = RZ( IREF ) + cell * ( K - 1 )

                        QW( IREF + MTEMP ) = QW( IREF )
                        QX( IREF + MTEMP ) = QX( IREF )
                        QY( IREF + MTEMP ) = QY( IREF )
                        QZ( IREF + MTEMP ) = QZ( IREF )

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

        if (box_lengths_set) then
            write(*, "(1x, 'Box length set to ', f10.5)") box_length
        else if (density_set) then
            write(*, "(1x, 'Density of system set to ', f10.5, ' atoms per cubic unit cell')") density

            box_lengths_set = .TRUE.
            volume = real(N) / density
            box_length = volume ** (1.0/3.0)
        else
            write(*, "(1x, 'Error: Density not set. Set density before initializing system.')")
            stop
        endif

        if (is_periodic) then
            write(*, "(1x, 'Periodic boundary conditions enabled for current system')")
        else
            write(*, "(1x, 'Periodic boundary conditions disabled for current system')")
        endif

        if (default_mass_set) then
            mass = 1.0
        endif
        if (default_charge_set) then
            charge = 0.0
        endif

        ! write(*, "(1x, 'Mass of particles set to ', f10.5)") mass
        ! write(*, "(1x, 'Charge of particles set to ', f10.5)") charge

        call initialize_forces

        if (fcc_lattice_set) then
            write(*, "(1x, 'FCC lattice generated with ', i10, ' particles')") N
        endif

        write(*, "(1x, 'System initialization complete.')")
        write(*, *)

    end subroutine initialize_system

    subroutine calculate_state_variables()
        implicit none

        kinetic_energy = kinetic_energy + 0.5 * mass * SUM( VX ** 2 + VY ** 2 + VZ ** 2 )
        total_energy = potential_energy + kinetic_energy

        temperature = 2.0 * kinetic_energy / dble(dof)

        virial = virial + SUM( RX * FX + RY * FY + RZ * FZ )

        if (volume > 0.0) then
            pressure = (density * temperature) + (virial / (3.0 * volume) )
        endif

        potential_energy = potential_energy / dble(N)
        kinetic_energy = kinetic_energy / dble(N)
        total_energy = total_energy / dble(N)

    end subroutine calculate_state_variables

end module system