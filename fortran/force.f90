module lennard_jones
    use system
    implicit none

    real(8) :: epsilon = 1.0
    real(8) :: sigma = 1.0

    contains
    subroutine lj_calculate_forces( I, J )
        implicit none
        integer, intent(in) :: I, J
        real(8)             :: RXIJ, RYIJ, RZIJ, FXIJ, FYIJ, FZIJ
        real(8)             :: rij_sq, sr2_lj, vir

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

module ecp
    use system
    implicit none

    real(8), dimension(3, 3) :: U1, U2, S1_sq, S2_sq, E1, E2, EM
    real(8), dimension(3, 3) :: AE, BE, GE, AR, BR, GR, MM1, MM2
    real(8), dimension(3) :: KE, KR
    real(8) :: lambda_E, lambda_R
    
    real(8) :: meu = 1.0, neu = 2.0
    real(8), dimension(3) :: RIJ, RIJ_hat, FIJ, TORQ1, TORQ2
    real(8) :: RIJ_SQ, RIJ_mag
    real(8) :: eps1, eps2, sr, phi, sigma, sigma_0 = 1.0

    contains
    real(8) function brent(func) result (res)
        implicit none
        integer, parameter :: max_iter = 10000
        real(8), parameter :: eps = 1e-8, psi = 0.5 * ( 3.0 - sqrt(5.0) )
        real(8)            :: a, b, c, x, w, v, u
        real(8)            :: fa, fb, fc, ft, fw, fv, fu
        real(8)            :: deltax = 0.0, atol = 1e-11, rtol = 1e-8
        real(8)            :: tol1, tol2, xmid, tmp1, tmp2, rat, p, dx_temp
        integer            :: I, iter = 0

        interface
            real(8) function func(xx)
                real(8), intent(in) :: xx
            end function func
        end interface

        a = 0.0; b = 1.0 
        c = (1 - psi) * a + psi * b
        fa = func(a); fb = func(b); fc = func(c)
        v = c; w = v; x = w
        fv = fc; fw = fv; ft = fw

        do while (iter < max_iter)
            tol1 = rtol * abs(x) + atol
            tol2 = 2.0 * tol1
            xmid = 0.5 * (a + b)

            if (abs(x - xmid) < (tol2 - 0.5 * (b - a))) exit

            if (abs(deltax) <= tol1) then
                if (x >= xmid) then
                    deltax = a - x
                else
                    deltax = b - x
                endif
                rat = psi * deltax

            else
                tmp1 = (x - w) * (ft - fv)
                tmp2 = (x - v) * (ft - fw)
                p = (x - v) * tmp2 - (x - w) * tmp1
                tmp2 = 2.0 * (tmp2 - tmp1)

                if (tmp2 > 0.0) p = -p

                tmp2 = abs(tmp2)
                dx_temp = deltax
                deltax = rat
                
                if ((p > tmp2 * (a - x)) .and. (p < tmp2 * (b - x)) .and. (abs(p) < abs(0.5 * tmp2 * dx_temp))) then
                    rat = p * 1.0 / tmp2
                    u = x + rat
                    if ((u - a) < tol2 .or. (b - u) < tol2) then
                        if (xmid - x >= 0) then
                            rat = tol1
                        else
                            rat = -tol1
                        endif
                    endif
                else
                    if (x >= xmid) then
                        deltax = a - x 
                    else
                        deltax = b - x
                    endif
                    rat = psi * deltax
                endif
            endif

            if (abs(rat) < tol1) then
                if (rat >= 0) then
                    u = x + tol1
                else
                    u = x - tol1
                endif
            else
                u = x + rat

            endif

            fu = func(u)

            if (fu > ft) then

                if (u < x) then
                    a = u
                else
                    b = u
                endif

                if ((fu <= fw) .or. (w == x)) then
                    v = w
                    w = u
                    fv = fw
                    fw = fu

                else if ((fu <= fv) .or. (v == x) .or. (v == w)) then
                    v = u
                    fv = fu
                endif

            else

                if (u >= x) then
                    a = x
                else
                    b = x
                endif

                v = w
                w = x
                x = u
                fv = fw
                fw = ft
                ft = fu

            endif
            iter = iter + 1
    
        enddo

        res = x
    end function brent

    real(8) function optim_eps( L )
        implicit none
        real(8), intent(in) :: L
        real(8), dimension(3, 3) :: GM
        real(8), dimension(3) :: KM

        GM = ( 1 - L ) * AE + L * BE
        KM = inverse(GM) .x. RIJ

        optim_eps = -L * ( 1 - L ) * dot( RIJ, KM )
        
    end function optim_eps

    real(8) function optim_dist( L )
        implicit none
        real(8), intent(in) :: L
        real(8), dimension(3, 3) :: GM_R
        real(8), dimension(3) :: KM_R

        GM_R = ( 1 - L ) * AR + L * BR
        KM_R = inverse(GM_R) .x. RIJ

        optim_dist = -L * ( 1 - L ) * dot( RIJ, KM_R )
        
    end function optim_dist

    pure function dF_dr( lam, k_vec ) result(res)
        implicit none
        real(8), intent(in) :: k_vec(3), lam
        real(8)             :: coeff, res(3)

        coeff = 2.0 * lam * ( 1 - lam ) / RIJ_SQ
        res = coeff * ( k_vec - dot_product( RIJ_hat, k_vec ) * RIJ_hat )
    end function dF_dr

    pure function dF_du( lam, K_vec, M_mat ) result(res)
        implicit none
        real(8), intent(in) :: M_mat(3, 3), K_vec(3), lam
        real(8)             :: coeff, res(3)

        coeff = 2.0 * lam * ( 1 - lam ) / RIJ_SQ
        res = coeff * cross( K_vec, matmul( K_vec, M_mat ) )
    end function dF_du

    subroutine ecp_calculate_forces( I, J )
        implicit none
        integer, intent(in)     :: I, J
        real(8), dimension(3)   :: dR_dr, de2_dr, dR_du1, dR_du2
        real(8), dimension(3)   :: de1_du1, de2_du1, de1_du2, de2_du2
        real(8)                 :: e0, sigma_t1, sigma_t2, sigma_t
        real(8)                 :: pot, vir
        integer                 :: IX
        
        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (periodic_set) then
            RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
        endif

        RIJ_SQ = SUM( RIJ**2 )
        RIj_MAG = SQRT( RIJ_SQ )
        RIJ_hat = RIJ / RIj_MAG

        U1(1, 1) = QW(I)**2 + QX(I)**2 - QY(I)**2 - QZ(I)**2
        U1(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
        U1(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
        U1(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
        U1(2, 2) = QW(I)**2 - QX(I)**2 + QY(I)**2 - QZ(I)**2
        U1(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
        U1(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
        U1(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
        U1(3, 3) = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

        U2(1, 1) = QW(J)**2 + QX(J)**2 - QY(J)**2 - QZ(J)**2
        U2(1, 2) = 2.0 * ( QX(J) * QY(J) + QW(J) * QZ(J) )
        U2(1, 3) = 2.0 * ( QX(J) * QZ(J) - QW(J) * QY(J) )
        U2(2, 1) = 2.0 * ( QX(J) * QY(J) - QW(J) * QZ(J) )
        U2(2, 2) = QW(J)**2 - QX(J)**2 + QY(J)**2 - QZ(J)**2
        U2(2, 3) = 2.0 * ( QY(J) * QZ(J) + QW(J) * QX(J) )
        U2(3, 1) = 2.0 * ( QX(J) * QZ(J) + QW(J) * QY(J) )
        U2(3, 2) = 2.0 * ( QY(J) * QZ(J) - QW(J) * QX(J) )
        U2(3, 3) = QW(J)**2 - QX(J)**2 - QY(J)**2 + QZ(J)**2

        e0 = maxval( erad )
        E1 = 0.0_8; E2 = 0.0_8; S1_sq = 0.0_8; S2_sq = 0.0_8

        E1(1, 1) = ( (e0 / erad(1)) ** (1/meu) ) / 4.0
        E1(2, 2) = ( (e0 / erad(2)) ** (1/meu) ) / 4.0
        E1(3, 3) = ( (e0 / erad(3)) ** (1/meu) ) / 4.0

        E2(1, 1) = ( (e0 / erad(1)) ** (1/meu) ) / 4.0
        E2(2, 2) = ( (e0 / erad(2)) ** (1/meu) ) / 4.0
        E2(3, 3) = ( (e0 / erad(3)) ** (1/meu) ) / 4.0

        S1_sq(1, 1) = ( (rad(1) / 2.0) ** 2.0 )
        S1_sq(2, 2) = ( (rad(2) / 2.0) ** 2.0 )
        S1_sq(3, 3) = ( (rad(3) / 2.0) ** 2.0 )

        S2_sq(1, 1) = ( (rad(1) / 2.0) ** 2.0 )
        S2_sq(2, 2) = ( (rad(2) / 2.0) ** 2.0 )
        S2_sq(3, 3) = ( (rad(3) / 2.0) ** 2.0 )

        AR = transpose(U1) .x. S1_sq .x. U1
        BR = transpose(U2) .x. S2_sq .x. U2
        AE = transpose(U1) .x. E1 .x. U1
        BE = transpose(U2) .x. E2 .x. U2

        sigma_t1 = (rad(1) * rad(2) + rad(3) ** 2) * sqrt(2.0 * rad(1) * rad(2)) / 8.0
        sigma_t2 = (rad(1) * rad(2) + rad(3) ** 2) * sqrt(2.0 * rad(1) * rad(2)) / 8.0
        sigma_t = sqrt( sigma_t1 * sigma_t2 )
        
        EM = inverse( AR + BR )
        eps1 = sigma_t * sqrt( determinant( EM ) )
        
        lambda_E = brent( optim_eps )
        if ( lambda_E < 0.0 .or. lambda_E > 1.0 ) then
            write(*, *) "LE = ", lambda_E, " not within [0, 1]"
            stop
        endif
        
        GE = ( 1 - lambda_E ) * AE + lambda_E * BE
        KE = inverse( GE ) .x. RIJ
        eps2 = lambda_E * (1.0 - lambda_E) * dot( RIJ_hat, KE ) / RIJ_mag
        
        lambda_R = brent( optim_dist )
        if ( lambda_R < 0.0 .or. lambda_R > 1.0 ) then
            write(*, *) "LR = ", lambda_R, " not within [0, 1]"
            stop
        endif
        
        GR = ( 1 - lambda_R ) * AR + lambda_R * BR
        KR = inverse( GR ) .x. RIJ
        phi = lambda_R * (1.0 - lambda_R) * dot( RIJ_hat, KR ) / RIJ_mag
        sigma = phi ** (-1.0/2.0)
        sr = sigma_0 / ( RIJ_mag - sigma + sigma_0 )
        

        dR_dr = (1.0 / sigma_0) * (( RIJ_hat ) + 0.5 * dF_dr( lambda_R, KR ) * phi ** (-3.0 / 2.0))
        de2_dr = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dF_dr( lambda_E, KE )
        
        de1_du1 = 0.0;
        MM1 = transpose(U1) .x. sqrt(S1_sq)
        do IX = 1, 3
            de1_du1 = de1_du1 + cross( MM1(:, IX), EM .x. MM1(:, IX) )
        enddo
        
        de1_du1 = -(eps2 ** meu) * (neu * eps1 ** neu) * de1_du1
        de2_du1 = (eps1 ** neu) * meu * eps2 ** (meu - 1) * (1 - lambda_E) * dF_du( lambda_E, KE, AE )
        dR_du1 =  (0.5 / sigma_0) * (1 - lambda_R) * dF_du( lambda_R, KR, AR ) * phi ** (-3.0 / 2.0)
        
        de1_du2 = 0.0
        MM2 = transpose(U2) .x. sqrt(S2_sq)
        do IX = 1, 3
            de1_du2 = de1_du2 + cross( MM2(:, IX), EM .x. MM2(:, IX) )
        enddo
        
        de1_du2 = -(eps2 ** meu) * (neu * eps1 ** neu) * de1_du2
        de2_du2 = (eps1 ** neu) * meu * eps2 ** (meu - 1) * lambda_E * dF_du( lambda_E, KE, BE )
        dR_du2 =  (0.5 / sigma_0) * lambda_R * dF_du( lambda_R, KR, BR ) * phi ** (-3.0 / 2.0)
        
        pot   = (eps1 ** neu) * (eps2 ** meu) * (sr ** 12.0 - 2.0 * sr ** 6.0)
        FIJ   = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * (12.0 * sr**7 - 12.0 * sr**13) * dR_dr  + (sr**12 - 2.0*sr**6) * de2_dr)
        TORQ1 = (-1.0) * ((sr**12 - 2.0*sr**6) * (de1_du1 + de2_du1) + ((eps1 ** neu) * (eps2 ** meu) * (12.0 * sr**7 - 12.0 * sr**13) * dR_du1))
        TORQ2 = (-1.0) * ((sr**12 - 2.0*sr**6) * (de1_du2 + de2_du2) + ((eps1 ** neu) * (eps2 ** meu) * (12.0 * sr**7 - 12.0 * sr**13) * dR_du2))
        
        potential_energy = potential_energy + pot
        virial = virial + dot( FIJ, RIJ )
        
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
        
    end subroutine ecp_calculate_forces
end module ecp