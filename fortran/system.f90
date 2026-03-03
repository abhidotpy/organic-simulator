module atom
    implicit none

    real(8), allocatable, dimension(:) :: mass, charge
    real(8), allocatable, dimension(:) :: shape_x, shape_y, shape_z, exc_rad
    real(8), allocatable, dimension(:) :: eshape_x, eshape_y, eshape_z
    real(8), allocatable, dimension(:) :: Ixx, Iyy, Izz
    integer, allocatable, dimension(:) :: rtype

    real(8), allocatable, dimension(:) :: RX, RY, RZ
    real(8), allocatable, dimension(:) :: VX, VY, VZ
    real(8), allocatable, dimension(:) :: FX, FY, FZ

    real(8), allocatable, dimension(:) :: QW, QX, QY, QZ
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
    logical :: default_shape_set    = .TRUE.
    logical :: default_energy_set   = .TRUE.
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

        allocate( mass(N), charge(N), rtype(N) )
        allocate( shape_x(N), shape_y(N), shape_z(N), exc_rad(N) )
        allocate( eshape_x(N), eshape_y(N), eshape_z(N) )
        allocate( Ixx(N), Iyy(N), Izz(N) )

        allocate( RX(N), RY(N), RZ(N) )
        allocate( VX(N), VY(N), VZ(N) )
        allocate( FX(N), FY(N), FZ(N) )
        allocate( LX(N), LY(N), LZ(N) )
        allocate( TX(N), TY(N), TZ(N) )
        allocate( QW(N), QX(N), QY(N), QZ(N) )

        mass = 1.0; charge = 0.0; rtype = 1;
        shape_x = 1.0; shape_y = 1.0; shape_z = 1.0; exc_rad = 1.0
        eshape_x = 1.0; eshape_y = 1.0; eshape_z = 1.0
        Ixx = 0.1; Iyy = 0.1; Izz = 0.1

        RX = 0.0; RY = 0.0; RZ = 0.0
        VX = 0.0; VY = 0.0; VZ = 0.0
        FX = 0.0; FY = 0.0; FZ = 0.0
        LX = 0.0; LY = 0.0; LZ = 0.0
        TX = 0.0; TY = 0.0; TZ = 0.0
        QW = 1.0; QX = 0.0; QY = 0.0; QZ = 0.0

    end subroutine set_num_atoms

    subroutine set_num_constraints( num_constraints )
        implicit none
        integer, intent(in) :: num_constraints

        constraints_set = .TRUE.
        NC = num_constraints
        dof = dof - NC

        allocate( RX_old(N), RY_old(N), RZ_old(N) )
        allocate( BBI(NC), BBJ(NC), BB_len(NC) )
        BBI = 0; BBJ = 0; BB_len = 0.0
        RX_old = 0.0; RY_old = 0.0; RZ_old = 0.0

    end subroutine set_num_constraints

    subroutine set_density( dens )
        implicit none
        real(8), intent(in) :: dens

        density_set = .TRUE.
        density = dens
    end subroutine set_density

    subroutine set_atom_type( id, atom_type )
        implicit none
        integer, intent(in) :: id
        integer, intent(in) :: atom_type

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "(1x, 'Error: System has no atoms. Set number of atoms first.')")
            stop
        else if ( id < 1 .or. id > N ) then
            write(*, "(1x, 'Error: Atom index ', I0,' does not exist &
            in system of ',I0,' atoms.')") id, N
        else
            rtype(id) = atom_type
        endif

    end subroutine set_atom_type

    subroutine set_atom_mass( id, atom_mass )
        implicit none
        integer, intent(in) :: id
        real(8), intent(in) :: atom_mass

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "(1x, 'Error: System has no atoms. Set number of atoms first.')")
            stop
        else if ( id < 1 .or. id > N ) then
            write(*, "(1x, 'Error: Atom index ', I0,' does not exist &
            in system of ',I0,' atoms.')") id, N
        else
            default_mass_set = .FALSE.
            mass(id) = atom_mass
        endif

    end subroutine set_atom_mass

    subroutine set_atom_charge( id, atom_charge )
        implicit none
        integer, intent(in) :: id
        real(8), intent(in) :: atom_charge

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "(1x, 'Error: System has no atoms. Set number of atoms first.')")
            stop
        else if ( id < 1 .or. id > N ) then
            write(*, "(1x, 'Error: Atom index ', I0,' does not exist &
            in system of ',I0,' atoms.')") id, N
        else
            default_charge_set = .FALSE.
            charge(id) = atom_charge
        endif

    end subroutine set_atom_charge

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

        volume = box_length(1) * box_length(2) * box_length(3)

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

    subroutine set_atom_shape( id, length, width, height )
        implicit none
        integer, intent(in) :: id
        real(8), intent(in) :: length, width, height

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "(1x, 'Error: System has no atoms. Set number of atoms first.')")
            stop
        else if ( id < 1 .or. id > N ) then
            write(*, "(1x, 'Error: Atom index ', I0,' does not exist &
            in system of ',I0,' atoms.')") id, N
        else
            default_shape_set = .FALSE.

            shape_x(id) = length
            shape_y(id) = width
            shape_z(id) = height

            Ixx(id) = ( shape_y(id) ** 2.0 + shape_z(id) ** 2.0 ) / 20.0
            Iyy(id) = ( shape_x(id) ** 2.0 + shape_z(id) ** 2.0 ) / 20.0
            Izz(id) = ( shape_x(id) ** 2.0 + shape_y(id) ** 2.0 ) / 20.0

            exc_rad(id) = max( length, width, height )
        endif

    
    end subroutine set_atom_shape

    subroutine set_atom_energy( id, eng_x, eng_y, eng_z )
        implicit none
        integer, intent(in) :: id
        real(8), intent(in) :: eng_x, eng_y, eng_z

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "(1x, 'Error: System has no atoms. Set number of atoms first.')")
            stop
        else if ( id < 1 .or. id > N ) then
            write(*, "(1x, 'Error: Atom index ', I0,' does not exist &
            in system of ',I0,' atoms.')") id, N
        else
            default_energy_set = .FALSE.

            eshape_x(id) = eng_x
            eshape_y(id) = eng_y
            eshape_z(id) = eng_z
        endif
        
    end subroutine set_atom_energy

    subroutine set_atom_position(id, pos_x, pos_y, pos_z)
        implicit none

        integer, intent(in) :: id
        real(8), intent(in) :: pos_x, pos_y, pos_z

        if (num_atoms_set .eqv. .FALSE.) then
            write(*, "(1x, 'Error: System has no atoms. Set number of atoms first.')")
            stop
        else if ( id < 1 .or. id > N ) then
            write(*, "(1x, 'Error: Atom index ', I0,' does not exist &
            in system of ',I0,' atoms.')") id, N
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
            write(*, "(1x, 'Error: Atom index ', I0,' does not exist &
            in system of ',I0,' atoms.')") id, N
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
            write(*, "(1x, 'Error: Constraint index ', I0,' does not exist &
            in system of ',I0,' constraints.')") id, NC
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

    subroutine generate_fcc_lattice( unit_cells, dens )
        implicit none
        integer, intent(in) :: unit_cells
        real(8), intent(in) :: dens
        integer             :: num_atoms, mtemp, I, J, K, IREF
        real(8)             :: cell, half_cell, rroot3

        fcc_lattice_set = .TRUE.
        num_atoms = 4 * unit_cells ** 3
        call set_num_atoms( num_atoms )

        if ( box_lengths_set .eqv. .FALSE. ) then

            box_lengths_set = .TRUE.
            
            if (dens == 0.0) then
                write(*, "(1x, 'Error: Density cannot be zero.')")
            else
                density_set = .TRUE.
                density = dens
                volume = real(N) / density
                box_length = volume ** (1.0/3.0)
            endif

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

        RX = RX - box_length(1) / 2.0
        RY = RY - box_length(2) / 2.0
        RZ = RZ - box_length(3) / 2.0

    end subroutine generate_fcc_lattice

    subroutine initialize_system( m )
        implicit none
        logical, optional :: m
        logical :: message

        if (not(present(m))) then
            message = .TRUE.
        else
            message = m
        endif

        if(message) write(*, *)
        if(message) write(*, "(1x, 'This is SUCROSE v1.0')")

        if (num_atoms_set) then
            if(message) write(*, "(1x, 'Number of particles in system set to ', I0)") N
        else
            if(message) write(*, "(1x, 'Error: Number of particles not set. Set number of &
            particles before initializing system.')")
            stop
        endif

        if (box_lengths_set) then
            if(message) write(*, "(1x, 'Box dimensions set to ', F0.2, 3X, F0.2, 3X, F0.2)") box_length

            if (density_set) then
                if(message) write(*, "(1x, 'Density of system set to ', F0.5)") density
            else
                density = dble(N) / volume
                density_set = .TRUE.
                if(message) write(*, "(1x, 'Density of system set to ', F0.5)") density
            endif

        else
            if(message) write(*, "(1x, 'Error: Box dimensions not set. Set box dimensions &
            before initializing system.')")
            stop
        endif
            
        if (is_periodic) then
            if(message) write(*, "(1x, 'Periodic boundary conditions enabled for current system')")
        else
            if(message) write(*, "(1x, 'Periodic boundary conditions disabled for current system')")
        endif

        if (default_mass_set) then
            mass = 1.0
            if(message) write(*, "(1x, 'Default mass for atoms set to ', G0.2)") mass(1)
        endif
        if (default_charge_set) then
            charge = 0.0
            if(message) write(*, "(1x, 'Default charge for atoms set to ', G0.2)") charge(1)
        endif
        if (default_shape_set) then
            shape_x = 1.0; shape_y = 1.0; shape_z = 1.0

            Ixx = ( shape_y(1) ** 2.0 + shape_z(1) ** 2.0 ) / 20.0
            Iyy = ( shape_x(1) ** 2.0 + shape_z(1) ** 2.0 ) / 20.0
            Izz = ( shape_x(1) ** 2.0 + shape_y(1) ** 2.0 ) / 20.0

            exc_rad = 1.0
            if(message) write(*, "(1x, 'Default shape for atoms set to ', G0.2)") shape_x(1)
        endif
        if (default_energy_set) then
            eshape_x = 1.0; eshape_y = 1.0; eshape_z = 1.0
            if(message) write(*, "(1x, 'Default energy for atoms set to ', G0.2)") eshape_x(1)
        endif

        if (fcc_lattice_set) then
            if(message) write(*, "(1x, 'FCC lattice generated with ', I0, ' particles')") N
        endif

        Ixx = Ixx * mass
        Iyy = Iyy * mass
        Izz = Izz * mass
        
        if(message) write(*, "(1x, 'System initialization complete.')")

        call initialize_forces

        if(message) write(*, *)

    end subroutine initialize_system

    subroutine calculate_state_variables()
        implicit none

        kinetic_energy = kinetic_energy + 0.5 * SUM( mass * ( VX ** 2 + VY ** 2 + VZ ** 2 ) )
        total_energy = potential_energy + kinetic_energy

        temperature = 2.0 * kinetic_energy / dble(dof)

        virial = virial + SUM( RX * FX + RY * FY + RZ * FZ )

        if (volume > 0.0) then
            pressure = (density * temperature) + (virial / (3.0 * volume) )
        endif

        density = dble(N) / volume

        potential_energy = potential_energy / dble(N)
        kinetic_energy = kinetic_energy / dble(N)
        total_energy = total_energy / dble(N)

    end subroutine calculate_state_variables

    subroutine copy_transform( from, to, xoffset, yoffset, zoffset )
        implicit none
        integer, intent(in) :: from, to
        real(8), intent(in) :: xoffset, yoffset, zoffset
        real(8)             :: rot(3, 3), offset(3), rot_offset(3)

        offset = (/ xoffset, yoffset, zoffset /)

        rot(1, 1) = QW(from)**2 + QX(from)**2 - QY(from)**2 - QZ(from)**2
        rot(1, 2) = 2.0 * ( QX(from) * QY(from) + QW(from) * QZ(from) )
        rot(1, 3) = 2.0 * ( QX(from) * QZ(from) - QW(from) * QY(from) )
        rot(2, 1) = 2.0 * ( QX(from) * QY(from) - QW(from) * QZ(from) )
        rot(2, 2) = QW(from)**2 - QX(from)**2 + QY(from)**2 - QZ(from)**2
        rot(2, 3) = 2.0 * ( QY(from) * QZ(from) + QW(from) * QX(from) )
        rot(3, 1) = 2.0 * ( QX(from) * QZ(from) + QW(from) * QY(from) )
        rot(3, 2) = 2.0 * ( QY(from) * QZ(from) - QW(from) * QX(from) )
        rot(3, 3) = QW(from)**2 - QX(from)**2 - QY(from)**2 + QZ(from)**2

        rot_offset(1) = rot(1, 1) * offset(1) + rot(1, 2) * offset(2) + rot(1, 3) * offset(3)
        rot_offset(2) = rot(2, 1) * offset(1) + rot(2, 2) * offset(2) + rot(2, 3) * offset(3)
        rot_offset(3) = rot(3, 1) * offset(1) + rot(3, 2) * offset(2) + rot(3, 3) * offset(3)

        RX(to) = RX(from) + rot_offset(1)
        RY(to) = RY(from) + rot_offset(2)
        RZ(to) = RZ(from) + rot_offset(3)

    end subroutine copy_transform

end module system