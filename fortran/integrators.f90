
module verlet_integrator
    use system
    implicit none

    contains
    subroutine vv_initial_step( DT )
        implicit none
        real(8), intent(in) :: DT

        VX = VX + 0.5 * DT * FX
        VY = VY + 0.5 * DT * FY
        VZ = VZ + 0.5 * DT * FZ

        RX = RX + DT * VX
        RY = RY + DT * VY
        RZ = RZ + DT * VZ

        if (periodic_set) then
            RX = RX - ANINT( RX / box_xlen ) * box_xlen
            RY = RY - ANINT( RY / box_ylen ) * box_ylen
            RZ = RZ - ANINT( RZ / box_zlen ) * box_zlen
        endif

    end subroutine vv_initial_step

    subroutine vv_final_step( DT )
        implicit none
        real(8), intent(in) :: DT

        VX = VX + 0.5 * DT * FX
        VY = VY + 0.5 * DT * FY
        VZ = VZ + 0.5 * DT * FZ

    end subroutine vv_final_step

end module verlet_integrator
    