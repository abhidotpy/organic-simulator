module lennard_jones
    use system
    implicit none

    real(8) :: epsilon = 1.0
    real(8) :: sigma = 1.0

    contains
    subroutine lj_calculate_forces( I, J )
        implicit none
        integer, intent(in) :: I, J
        real(8) :: rxij, ryij, rzij, fxij, fyij, fzij
        real(8) :: rij_sq, sr2_lj, vir

        RXIJ = RX(I) - RX(J)
        RYIJ = RY(I) - RY(J)
        RZIJ = RZ(I) - RZ(J)

        if (periodic_set) then
            RXIJ = RXIJ - ANINT( RXIJ / box_xlen ) * box_xlen
            RYIJ = RYIJ - ANINT( RYIJ / box_ylen ) * box_ylen
            RZIJ = RZIJ - ANINT( RZIJ / box_zlen ) * box_zlen
        endif

        rij_sq = RXIJ**2 + RYIJ**2 + RZIJ**2
        sr2_lj = ( sigma ** 2.0 ) / rij_sq
        potential_energy = potential_energy + 4.0 * epsilon * ( sr2_lj**6 - sr2_lj**3 )

        vir = 24.0 * epsilon * (2.0 * sr2_lj**6 - sr2_lj**3)
        FXIJ = vir * RXIJ / rij_sq
        FYIJ = vir * RYIJ / rij_sq
        FZIJ = vir * RZIJ / rij_sq

        virial = virial + vir

        FX(I) = FX(I) + FXIJ
        FY(I) = FY(I) + FYIJ
        FZ(I) = FZ(I) + FZIJ

        FX(J) = FX(J) - FXIJ
        FY(J) = FY(J) - FYIJ
        FZ(J) = FZ(J) - FZIJ

    end subroutine lj_calculate_forces

end module lennard_jones

module lennard_jones12
    use system
    implicit none

    real(8) :: epsilon = 1.0
    real(8) :: sigma = 1.0

    contains
    subroutine lj12_calculate_forces( I, J )
        implicit none
        integer, intent(in) :: I, J
        real(8) :: RXIJ, RYIJ, RZIJ, FXIJ, FYIJ, FZIJ
        real(8) :: rij_sq, sr2_lj, vir

        RXIJ = RX(I) - RX(J)
        RYIJ = RY(I) - RY(J)
        RZIJ = RZ(I) - RZ(J)

        if (periodic_set) then
            RXIJ = RXIJ - ANINT( RXIJ / box_xlen ) * box_xlen
            RYIJ = RYIJ - ANINT( RYIJ / box_ylen ) * box_ylen
            RZIJ = RZIJ - ANINT( RZIJ / box_zlen ) * box_zlen
        endif

        rij_sq = RXIJ ** 2 + RYIJ ** 2 + RZIJ ** 2
        sr2_lj = ( sigma ** 2.0 ) / rij_sq

        potential_energy = potential_energy + epsilon * ( sr2_lj**6.0 - 2.0*sr2_lj**3.0 )
        vir = 12.0 * epsilon * (sr2_lj**6 - sr2_lj**3)

        FXIJ = vir * RXIJ / rij_sq
        FYIJ = vir * RYIJ / rij_sq
        FZIJ = vir * RZIJ / rij_sq

        virial = virial + vir

        FX(I) = FX(I) + FXIJ
        FY(I) = FY(I) + FYIJ
        FZ(I) = FZ(I) + FZIJ

        FX(J) = FX(J) - FXIJ
        FY(J) = FY(J) - FYIJ
        FZ(J) = FZ(J) - FZIJ

    end subroutine lj12_calculate_forces

end module lennard_jones12