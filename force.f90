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

module morse
    use system
    implicit none

    real(8) :: DE, alpha, r0, rp

    contains
    subroutine morse_calculate_forces( I, J )
        implicit none
        integer, intent(in) :: I, J
        real(8) :: RXIJ, RYIJ, RZIJ, FXIJ, FYIJ, FZIJ
        real(8) :: rij, exp_term, vir

        RXIJ = RX(I) - RX(J)
        RYIJ = RY(I) - RY(J)
        RZIJ = RZ(I) - RZ(J)

        if (periodic_set) then
            RXIJ = RXIJ - ANINT( RXIJ / box_xlen ) * box_xlen
            RYIJ = RYIJ - ANINT( RYIJ / box_ylen ) * box_ylen
            RZIJ = RZIJ - ANINT( RZIJ / box_zlen ) * box_zlen
        endif

        rij = SQRT( RXIJ**2 + RYIJ**2 + RZIJ**2 )
        exp_term = EXP( -alpha * ( rij - r0 ) )

        potential_energy = potential_energy + DE * ( (1.0 - exp_term)**2 - 1.0 )
        vir = -2.0 * alpha * DE * (1.0 - exp_term) * exp_term / rij

        FXIJ = vir * RXIJ
        FYIJ = vir * RYIJ
        FZIJ = vir * RZIJ

        virial = virial + vir * rij

        FX(I) = FX(I) + FXIJ
        FY(I) = FY(I) + FYIJ
        FZ(I) = FZ(I) + FZIJ

        FX(J) = FX(J) - FXIJ
        FY(J) = FY(J) - FYIJ
        FZ(J) = FZ(J) - FZIJ

    end subroutine morse_calculate_forces

end module morse

module gay_berne
    use system
    implicit none
    
    real(8), dimension(3), public :: U1, U2, RIJ, FIJ, TI, TJ, TORQ1, TORQ2
    public :: gb_calculate_forces
    
    private
    real(8) :: meu = 1.0, neu = 2.0
    real(8) :: ru1, ru2, uu, chi, xhi, rij_sq, rij_mag
    real(8) :: eps1, eps2, sr, sigma, sigma_0 = 1.0


    contains
    function g_func( chie ) result (res)
        implicit none
        real(8), intent(in) :: chie
        real(8)             :: term1, term2, res

        term1 = (ru1 + ru2)**2 / (1 + chie * uu)
        term2 = (ru1 - ru2)**2 / (1 - chie * uu)
        
        res = 1 - (chie / 2.0 / rij_sq) * (term1 + term2)

    end function g_func

    function dG_dr( chie ) result (res)
        implicit none
        real(8), intent(in)     :: chie
        real(8), dimension(3)   :: term1, term2, res

        term1 = ((ru1 + ru2) / (1 + chie * uu)) * (U1 + U2)
        term2 = ((ru1 - ru2) / (1 - chie * uu)) * (U1 - U2)
        
        res = (-chie / RIJ_SQ ) * (term1 + term2) + (2 * RIJ / RIJ_SQ) * (1 - g_func( chie ))

    end function dG_dr

    function dG_du1( chie ) result(res)
        implicit none
        real(8), intent(in)    :: chie
        real(8)                :: term1, term2, res(3)
        
        term1 = ((ru1 + ru2) / (1 + chie * uu))
        term2 = ((ru1 - ru2) / (1 - chie * uu))
        
        res = (-chie * RIJ / rij_sq ) * (term1 + term2) + ((chie**2 * U2) / (2 * rij_sq)) * (term1 ** 2 - term2 ** 2)
        
    end function dG_du1

    function dG_du2( chie ) result(res)
        implicit none
        real(8), intent(in)    :: chie
        real(8)                :: term1, term2, res(3)
        
        term1 = ((ru1 + ru2) / (1 + chie * uu))
        term2 = ((ru1 - ru2) / (1 - chie * uu))
        
        res = (-chie * RIJ / RIJ_SQ ) * (term1 - term2) + ((chie**2 * U1) / (2 * rij_sq)) * (term1 ** 2 - term2 ** 2)
        
    end function dG_du2

    subroutine gb_calculate_forces( I, J )
        implicit none
        integer, intent(in)     :: I, J
        real(8), dimension(3)   :: dR_dr, de2_dr, dR_du1, dR_du2
        real(8), dimension(3)   :: de1_du1, de2_du1, de1_du2, de2_du2
        real(8)                 :: pot, vir

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (periodic_set) then
            RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
        endif

        RIJ_SQ = SUM( RIJ**2 )
        RIj_MAG = SQRT( RIJ_SQ )

        U1(1) = 2.0 *  ( QX(I) * QZ(I) + QW(I) * QY(I) )
        U1(2) = 2.0 *  ( QY(I) * QZ(I) - QW(I) * QX(I) )
        U1(3) = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

        U2(1) = 2.0 *  ( QX(J) * QZ(J) + QW(J) * QY(J) )
        U2(2) = 2.0 *  ( QY(J) * QZ(J) - QW(J) * QX(J) )
        U2(3) = QW(J)**2 - QX(J)**2 - QY(J)**2 + QZ(J)**2

        ru1 = dot( RIJ, U1 )
        ru2 = dot( RIJ, U2 )
        uu  = dot( U1, U2 )
        chi = ( rad(3) ** 2 - rad(1) ** 2 ) / ( rad(3) ** 2 + rad(1) ** 2 )
        xhi = ( erad(1) ** (1.0/meu) - erad(3) ** (1.0/meu) ) / ( erad(1) ** (1.0/meu) + erad(3) ** (1.0/meu) )

        eps1 = ( 1 - ( chi * uu ) ** 2) ** (-1.0 / 2.0)
        eps2 = g_func( xhi )
        sigma = rad(1) * g_func( chi ) ** (-1.0 / 2.0)
        sr = sigma_0 / ( rij_mag - sigma + sigma_0 )

        dR_dr   = (1.0 / sigma_0) * (( RIJ / rij_mag ) + 0.5 * rad(1) * dG_dr( chi ) * g_func( chi ) ** (-3.0 / 2.0))
        de2_dr  = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_dr( xhi )

        de1_du1 = (eps2 ** meu) * neu * eps1 ** (neu + 2) * (chi ** 2 * uu * U2)
        de2_du1 = (eps1 ** neu) * meu * eps2 ** (meu - 1) * dG_du1( xhi )
        dR_du1  = (12.0 * sr**7 - 12.0 * sr**13) * 0.5 * (rad(1) / sigma_0) * dG_du1( chi ) * g_func( chi ) ** (-3.0 / 2.0)

        de1_du2 = (eps2 ** meu) * neu * eps1 ** (neu + 2) * (chi ** 2 * uu * U1)
        de2_du2 = (eps1 ** neu) * meu * eps2 ** (meu - 1) * dG_du2( xhi )
        dR_du2  = (12.0 * sr**7 - 12.0 * sr**13) * 0.5 * (rad(1) / sigma_0) * dG_du2( chi ) * g_func( chi ) ** (-3.0 / 2.0)
    
        pot = (eps1 ** neu) * (eps2 ** meu) * (sr ** 12.0 - 2.0 * sr ** 6.0)
        FIJ = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * (12.0 * sr**7 - 12.0 * sr**13) * dR_dr  + de2_dr * (sr**12 - 2.0*sr**6))
        TI  = (-1.0) * ((sr**12 - 2.0*sr**6) * (de1_du1 + de2_du1) + ((eps1 ** neu) * (eps2 ** meu) * dR_du1))
        TJ  = (-1.0) * ((sr**12 - 2.0*sr**6) * (de1_du2 + de2_du2) + ((eps1 ** neu) * (eps2 ** meu) * dR_du2))

        potential_energy = potential_energy + pot
        virial = virial + dot( FIJ, RIJ )
        TORQ1 = cross( U1, TI )
        TORQ2 = cross( U2, TJ )

        FX(I) = FX(I) + FIJ(1)
        FY(I) = FY(I) + FIJ(2)
        FZ(I) = FZ(I) + FIJ(3)

        FX(J) = FX(J) - FIJ(1)
        FY(J) = FY(J) - FIJ(2)
        FZ(J) = FZ(J) - FIJ(3)

        TX(I) = TX(I) + TORQ1(1)
        TY(I) = TY(I) + TORQ1(2)
        TZ(I) = TZ(I) + TORQ1(3)

        TX(J) = TX(J) + TORQ2(1)
        TY(J) = TY(J) + TORQ2(2)
        TZ(J) = TZ(J) + TORQ2(3)

    end subroutine gb_calculate_forces

end module gay_berne