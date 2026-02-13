module lennard_jones
    use system
    implicit none

    real(8) :: epsilon = 1.0
    real(8) :: sigma = 1.0

    contains
    subroutine lj_calculate_forces( I, J, cutoff, eng )
        implicit none
        integer, intent(in)           :: I, J
        real(8), optional, intent(in) :: eng
        real(8), intent(in)           :: cutoff
        real(8), dimension(3)         :: RIJ, FIJ
        real(8)                       :: rij_sq, rcut_sq, sr2_lj, coeff
        real(8)                       :: pot, pot_cut, sr2_cut

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (is_periodic) then
            RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
        endif

        if (present(eng)) then
            epsilon = eng
        else
            epsilon = 1.0
        endif

        rij_sq = SUM( RIJ ** 2.0 )
        rcut_sq = cutoff ** 2.0

        if (rij_sq < rcut_sq) then
            sigma = ( exc_rad(I) + exc_rad(J) ) / 2.0
            sr2_lj = ( sigma ** 2.0 ) / rij_sq
            sr2_cut = ( sigma ** 2.0 ) / rcut_sq

            pot = 4.0 * epsilon * ( sr2_lj**6 - sr2_lj**3 )
            pot_cut = 4.0 * epsilon * ( sr2_cut**6 - sr2_cut**3 )
            coeff = 24.0 * epsilon * (2.0 * sr2_lj**6 - sr2_lj**3)

            potential_energy = potential_energy + pot - pot_cut
            FIJ = coeff * RIJ / rij_sq

            FX(I) = FX(I) + FIJ(1)
            FY(I) = FY(I) + FIJ(2)
            FZ(I) = FZ(I) + FIJ(3)

            FX(J) = FX(J) - FIJ(1)
            FY(J) = FY(J) - FIJ(2)
            FZ(J) = FZ(J) - FIJ(3)
        endif

    end subroutine lj_calculate_forces

end module lennard_jones

module lennard_jones12
    use system
    implicit none

    real(8) :: epsilon = 1.0
    real(8) :: sigma = 1.0

    contains
    subroutine lj12_calculate_forces( I, J, cutoff, eng )
        implicit none
        integer, intent(in)           :: I, J
        real(8), optional, intent(in) :: eng 
        real(8), intent(in)           :: cutoff
        real(8), dimension(3)         :: RIJ, FIJ
        real(8)                       :: rij_sq, rcut_sq, sr2_lj, coeff
        real(8)                       :: pot, pot_cut, sr2_cut


        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (is_periodic) then
            RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
        endif

        if (present(eng)) then
            epsilon = eng
        else
            epsilon = 1.0
        endif

        rij_sq = SUM( RIJ ** 2.0 )
        rcut_sq = cutoff ** 2.0

        if (rij_sq < rcut_sq) then
            sigma = ( exc_rad(I) + exc_rad(J) ) / 2.0
            sr2_lj = ( sigma ** 2.0 ) / rij_sq
            sr2_cut = ( sigma ** 2.0 ) / rcut_sq

            pot = epsilon * ( sr2_lj**6.0 - 2.0*sr2_lj**3.0 )
            pot_cut = epsilon * ( sr2_cut**6.0 - 2.0*sr2_cut**3.0 )
            coeff = 12.0 * epsilon * (sr2_lj**6 - sr2_lj**3)

            potential_energy = potential_energy + pot - pot_cut
            FIJ = coeff * RIJ / rij_sq

            FX(I) = FX(I) + FIJ(1)
            FY(I) = FY(I) + FIJ(2)
            FZ(I) = FZ(I) + FIJ(3)

            FX(J) = FX(J) - FIJ(1)
            FY(J) = FY(J) - FIJ(2)
            FZ(J) = FZ(J) - FIJ(3)
        endif

    end subroutine lj12_calculate_forces

end module lennard_jones12

module morse
    use system
    implicit none

    contains
    subroutine morse_calculate_forces( I, J, dissoc, width, rmin, cutoff )
        implicit none
        integer, intent(in)   :: I, J
        real(8), dimension(3) :: RIJ, FIJ
        real(8), intent(in)   :: dissoc, width, rmin, cutoff
        real(8)               :: rij_mag, exp_term
        real(8)               :: exp_cut, pot, pot_cut

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (is_periodic) then
            RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
        endif

        rij_mag = SQRT( SUM( RIJ**2.0 ) )

        if (rij_mag < cutoff) then
            exp_term = EXP( -width * ( rij_mag - rmin ) )
            exp_cut  = EXP( -width * ( cutoff - rmin ) )

            pot = dissoc * ( (1.0 - exp_term)**2.0 - 1.0 )
            pot_cut = dissoc * ( (1.0 - exp_cut)**2.0 - 1.0 )

            potential_energy = potential_energy + pot - pot_cut
            FIJ = -2.0 * width * dissoc * (1.0 - exp_term) * exp_term * RIJ / rij_mag

            FX(I) = FX(I) + FIJ(1)
            FY(I) = FY(I) + FIJ(2)
            FZ(I) = FZ(I) + FIJ(3)

            FX(J) = FX(J) - FIJ(1)
            FY(J) = FY(J) - FIJ(2)
            FZ(J) = FZ(J) - FIJ(3)
        endif

    end subroutine morse_calculate_forces

end module morse

module gay_berne
    use system
    implicit none
    
    real(8), dimension(3), public :: U1, U2, RIJ, RIJ_cut
    real(8), dimension(3), public :: FIJ, TI, TJ, TORQ1, TORQ2
    real(8), dimension(3), public :: FIJ_cut, TI_cut, TJ_cut, TORQ1_cut, TORQ2_cut
    public :: gb_calculate_forces
    
    private
    real(8) :: meu = 1.0, neu = 2.0
    real(8) :: ru1, ru2, uu, chi, xhi, rij_sq, rij_mag
    real(8) :: eps1, eps2, sr, sigma, sigma_0 = 1.0
    real(8) :: sr_cut


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

    subroutine gb_calculate_forces( I, J, cutoff, wf )
        implicit none
        integer, intent(in)             :: I, J
        real(8), optional               :: wf
        real(8), dimension(3)           :: dR_dr, de2_dr, dR_du1, dR_du2
        real(8), dimension(3)           :: de1_du1, de2_du1, de1_du2, de2_du2
        real(8)                         :: pot, pot_cut, cutoff

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (is_periodic) then
            RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
        endif

        RIJ_SQ = SUM( RIJ**2 )
        RIJ_MAG = SQRT( RIJ_SQ )

        U1(1) = 2.0 *  ( QX(I) * QZ(I) + QW(I) * QY(I) )
        U1(2) = 2.0 *  ( QY(I) * QZ(I) - QW(I) * QX(I) )
        U1(3) = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

        U2(1) = 2.0 *  ( QX(J) * QZ(J) + QW(J) * QY(J) )
        U2(2) = 2.0 *  ( QY(J) * QZ(J) - QW(J) * QX(J) )
        U2(3) = QW(J)**2 - QX(J)**2 - QY(J)**2 + QZ(J)**2

        ! CALCULATE GAY BERNE PARAMETERS

        ru1 = dot( RIJ, U1 )
        ru2 = dot( RIJ, U2 )
        uu  = dot( U1, U2 )
        chi = ( shape_z(I) ** 2 - shape_x(I) ** 2 ) / ( shape_z(I) ** 2 + shape_x(I) ** 2 )
        xhi = ( eshape_x(I) ** (1.0/meu) - eshape_z(I) ** (1.0/meu) ) / ( eshape_x(I) ** (1.0/meu) + eshape_z(I) ** (1.0/meu) )

        eps1 = ( 1 - ( chi * uu ) ** 2) ** (-1.0 / 2.0)
        eps2 = g_func( xhi )
        sigma = shape_x(I) * g_func( chi ) ** (-1.0 / 2.0)

        if (present(wf)) then
            if ( rij_mag < sigma ) then
                sr = sigma_0 / ( rij_mag - sigma + sigma_0 )
            else if (rij_mag .ge. sigma .and. rij_mag .le. (sigma + wf)) then
                sr = sigma_0 / ( sigma - sigma + sigma_0 )
            else
                sr = sigma_0 / ( rij_mag - sigma + sigma_0 - wf )
            endif
        else
            sr = sigma_0 / ( rij_mag - sigma + sigma_0 )
        endif

        dR_dr   = (1.0 / sigma_0) * (( RIJ / rij_mag ) + 0.5 * shape_x(I) * dG_dr( chi ) * g_func( chi ) ** (-3.0 / 2.0))
        de2_dr  = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_dr( xhi )

        de1_du1 = (eps2 ** meu) * neu * eps1 ** (neu + 2) * (chi ** 2 * uu * U2)
        de2_du1 = (eps1 ** neu) * meu * eps2 ** (meu - 1) * dG_du1( xhi )
        dR_du1  = (12.0 * sr**7 - 12.0 * sr**13) * 0.5 * (shape_x(I) / sigma_0) * dG_du1( chi ) * g_func( chi ) ** (-3.0 / 2.0)

        de1_du2 = (eps2 ** meu) * neu * eps1 ** (neu + 2) * (chi ** 2 * uu * U1)
        de2_du2 = (eps1 ** neu) * meu * eps2 ** (meu - 1) * dG_du2( xhi )
        dR_du2  = (12.0 * sr**7 - 12.0 * sr**13) * 0.5 * (shape_x(I) / sigma_0) * dG_du2( chi ) * g_func( chi ) ** (-3.0 / 2.0)
    
        pot = (eps1 ** neu) * (eps2 ** meu) * (sr ** 12.0 - 2.0 * sr ** 6.0)
        FIJ = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * (12.0 * sr**7 - 12.0 * sr**13) * dR_dr  + de2_dr * (sr ** 12.0 - 2.0 * sr ** 6.0))
        TI  = (-1.0) * ((sr ** 12.0 - 2.0 * sr ** 6.0) * (de1_du1 + de2_du1) + ((eps1 ** neu) * (eps2 ** meu) * dR_du1))
        TJ  = (-1.0) * ((sr ** 12.0 - 2.0 * sr ** 6.0) * (de1_du2 + de2_du2) + ((eps1 ** neu) * (eps2 ** meu) * dR_du2))
        TORQ1 = cross( U1, TI )
        TORQ2 = cross( U2, TJ )

        ! CALCULATE CUTOFF PARAMETERS

        RIJ_cut = cutoff * (RIJ / rij_mag)
        RIJ = RIJ_cut

        if (is_periodic) then
            RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
        endif

        RIJ_SQ = cutoff ** 2.0
        RIJ_MAG = cutoff
        
        ru1 = dot( RIJ, U1 )
        ru2 = dot( RIJ, U2 )
        eps2 = g_func( xhi )
        sigma = shape_x(I) * g_func( chi ) ** (-1.0 / 2.0)

        if (present(wf)) then
            if ( rij_mag < sigma ) then
                sr = sigma_0 / ( rij_mag - sigma + sigma_0 )
            else if (rij_mag .ge. sigma .and. rij_mag .le. (sigma + wf)) then
                sr = sigma_0 / ( sigma - sigma + sigma_0 )
            else
                sr = sigma_0 / ( rij_mag - sigma + sigma_0 - wf )
            endif
        else
            sr = sigma_0 / ( rij_mag - sigma + sigma_0 )
        endif

        dR_dr   = (1.0 / sigma_0) * (( RIJ / rij_mag ) + 0.5 * shape_x(I) * dG_dr( chi ) * g_func( chi ) ** (-3.0 / 2.0))
        de2_dr  = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_dr( xhi )

        de1_du1 = (eps2 ** meu) * neu * eps1 ** (neu + 2) * (chi ** 2 * uu * U2)
        de2_du1 = (eps1 ** neu) * meu * eps2 ** (meu - 1) * dG_du1( xhi )
        dR_du1  = (12.0 * sr**7 - 12.0 * sr**13) * 0.5 * (shape_x(I) / sigma_0) * dG_du1( chi ) * g_func( chi ) ** (-3.0 / 2.0)

        de1_du2 = (eps2 ** meu) * neu * eps1 ** (neu + 2) * (chi ** 2 * uu * U1)
        de2_du2 = (eps1 ** neu) * meu * eps2 ** (meu - 1) * dG_du2( xhi )
        dR_du2  = (12.0 * sr**7 - 12.0 * sr**13) * 0.5 * (shape_x(I) / sigma_0) * dG_du2( chi ) * g_func( chi ) ** (-3.0 / 2.0)
    
        pot_cut   = (eps1 ** neu) * (eps2 ** meu) * (sr ** 12.0 - 2.0 * sr ** 6.0)
        FIJ_cut   = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * (12.0 * sr**7 - 12.0 * sr**13) * dR_dr  + de2_dr * (sr ** 12.0 - 2.0 * sr ** 6.0))
        TI_cut    = (-1.0) * ((sr ** 12.0 - 2.0 * sr ** 6.0) * (de1_du1 + de2_du1) + ((eps1 ** neu) * (eps2 ** meu) * dR_du1))
        TJ_cut    = (-1.0) * ((sr ** 12.0 - 2.0 * sr ** 6.0) * (de1_du2 + de2_du2) + ((eps1 ** neu) * (eps2 ** meu) * dR_du2))
        TORQ1_cut = cross( U1, TI_cut )
        TORQ2_cut = cross( U2, TJ_cut )

        potential_energy = potential_energy + pot - pot_cut

        FX(I) = FX(I) + ( FIJ(1) - FIJ_cut(1) )
        FY(I) = FY(I) + ( FIJ(2) - FIJ_cut(2) )
        FZ(I) = FZ(I) + ( FIJ(3) - FIJ_cut(3) )

        FX(J) = FX(J) - ( FIJ(1) - FIJ_cut(1) )
        FY(J) = FY(J) - ( FIJ(2) - FIJ_cut(2) )
        FZ(J) = FZ(J) - ( FIJ(3) - FIJ_cut(3) )

        TX(I) = TX(I) + ( TORQ1(1) - TORQ1_cut(1) )
        TY(I) = TY(I) + ( TORQ1(2) - TORQ1_cut(2) )
        TZ(I) = TZ(I) + ( TORQ1(3) - TORQ1_cut(3) )

        TX(J) = TX(J) + ( TORQ2(1) - TORQ2_cut(1) )
        TY(J) = TY(J) + ( TORQ2(2) - TORQ2_cut(2) )
        TZ(J) = TZ(J) + ( TORQ2(3) - TORQ2_cut(3) )

    end subroutine gb_calculate_forces

end module gay_berne

module ecp
    use system
    implicit none

    real(8), dimension(3, 3) :: U1, U2, S1, S2, E1, E2, EM
    real(8), dimension(3, 3) :: AE, BE, GE, AR, BR, GR, MM1, MM2
    real(8), dimension(3)    :: KE, KR
    real(8)                  :: lambda_E, lambda_R
    
    real(8)                  :: meu = 1.0, neu = 2.0
    real(8), dimension(3)    :: RIJ, RIJ_hat, FIJ, TORQ1, TORQ2
    real(8), dimension(3)    :: RIJ_cut, FIJ_cut, TORQ1_cut, TORQ2_cut
    real(8)                  :: RIJ_SQ, RIJ_mag
    real(8)                  :: eps1, eps2, sr, phi, sigma, sigma_0 = 1.0

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
        real(8), intent(in)         :: L
        real(8), dimension(3, 3)    :: GM
        real(8), dimension(3)       :: KM

        GM = ( 1 - L ) * AE + L * BE
        KM = inverse(GM) .x. RIJ

        optim_eps = -L * ( 1 - L ) * dot( RIJ, KM )
        
    end function optim_eps

    real(8) function optim_dist( L )
        implicit none
        real(8), intent(in)         :: L
        real(8), dimension(3, 3)    :: GM_R
        real(8), dimension(3)       :: KM_R

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

    subroutine ecp_calculate_forces( I, J, cutoff, wf )
        implicit none
        integer, intent(in)             :: I, J
        real(8), optional, intent(in)   :: wf
        real(8), dimension(3)           :: dR_dr, de2_dr, dR_du1, dR_du2
        real(8), dimension(3)           :: de1_du1, de2_du1, de1_du2, de2_du2
        real(8)                         :: e10, e20, sigma_t1, sigma_t2, sigma_t
        real(8)                         :: pot, pot_cut, cutoff
        integer                         :: IX
        
        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (is_periodic) then
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

        e10 = max( eshape_x(I), eshape_y(I), eshape_z(I) )
        e20 = max( eshape_x(J), eshape_y(J), eshape_z(J) )
        E1 = 0.0_8; E2 = 0.0_8; S1 = 0.0_8; S2 = 0.0_8

        E1(1, 1) = ( (e10 / eshape_x(I)) ** (1/meu) ) / 4.0
        E1(2, 2) = ( (e10 / eshape_y(I)) ** (1/meu) ) / 4.0
        E1(3, 3) = ( (e10 / eshape_z(I)) ** (1/meu) ) / 4.0

        E2(1, 1) = ( (e20 / eshape_x(J)) ** (1/meu) ) / 4.0
        E2(2, 2) = ( (e20 / eshape_y(J)) ** (1/meu) ) / 4.0
        E2(3, 3) = ( (e20 / eshape_z(J)) ** (1/meu) ) / 4.0

        S1(1, 1) = ( shape_x(I) / 2.0 )
        S1(2, 2) = ( shape_y(I) / 2.0 )
        S1(3, 3) = ( shape_z(I) / 2.0 )

        S2(1, 1) = ( shape_x(J) / 2.0 )
        S2(2, 2) = ( shape_y(J) / 2.0 )
        S2(3, 3) = ( shape_z(J) / 2.0 )

        AR = transpose(U1) .x. ( S1 ** 2.0 ) .x. U1
        BR = transpose(U2) .x. ( S2 ** 2.0 ) .x. U2
        AE = transpose(U1) .x. E1 .x. U1
        BE = transpose(U2) .x. E2 .x. U2

        sigma_t1 = (shape_x(I) * shape_y(I) + shape_z(I) ** 2) * sqrt( 2.0 * shape_x(I) * shape_y(I) ) / 8.0
        sigma_t2 = (shape_x(J) * shape_y(J) + shape_z(J) ** 2) * sqrt( 2.0 * shape_x(J) * shape_y(J) ) / 8.0
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
        
        if (present(wf)) then
            if ( rij_mag < sigma ) then
                sr = sigma_0 / ( rij_mag - sigma + sigma_0 )
            else if (rij_mag .ge. sigma .and. rij_mag .le. (sigma + wf)) then
                sr = sigma_0 / ( sigma - sigma + sigma_0 )
            else
                sr = sigma_0 / ( rij_mag - sigma + sigma_0 - wf )
            endif
        else
            sr = sigma_0 / ( rij_mag - sigma + sigma_0 )
        endif
        
        dR_dr = (1.0 / sigma_0) * (( RIJ_hat ) + 0.5 * dF_dr( lambda_R, KR ) * phi ** (-3.0 / 2.0))
        de2_dr = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dF_dr( lambda_E, KE )
        
        de1_du1 = 0.0;
        MM1 = transpose(U1) .x. S1
        do IX = 1, 3
            de1_du1 = de1_du1 + cross( MM1(:, IX), EM .x. MM1(:, IX) )
        enddo
        
        de1_du1 = -(eps2 ** meu) * (neu * eps1 ** neu) * de1_du1
        de2_du1 = (eps1 ** neu) * meu * eps2 ** (meu - 1) * (1 - lambda_E) * dF_du( lambda_E, KE, AE )
        dR_du1 =  (0.5 / sigma_0) * (1 - lambda_R) * dF_du( lambda_R, KR, AR ) * phi ** (-3.0 / 2.0)
        
        de1_du2 = 0.0
        MM2 = transpose(U2) .x. S2
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

module bonded_interactions
    use system
    implicit none

    contains
    subroutine bend_calculate_forces( I, J, K, coeff )
        implicit none
        integer, intent(in) :: I, J, K
        real(8), intent(in) :: coeff
        real(8) :: VI(3), VJ(3), FI(3), FJ(3), FK(3)
        real(8) :: CC11, CC12, CC22
        real(8) :: prefac, fac, fac1, fac2, pot

        VI(1) = RX(J) - RX(I)
        VI(2) = RY(J) - RY(I)
        VI(3) = RZ(J) - RZ(I)

        VJ(1) = RX(K) - RX(J)
        VJ(2) = RY(K) - RY(J)
        VJ(3) = RZ(K) - RZ(J)

        if (is_periodic) then
            VI = VI - ANINT( VI / box_length ) * box_length
            VJ = VJ - ANINT( VJ / box_length ) * box_length
        endif

        CC11 = dot( VI, VI )
        CC12 = dot( VI, VJ )
        CC22 = dot( VJ, VJ )

        prefac = 1.0 / SQRT( CC11 * CC22 )
        fac    = CC12
        fac1   = fac / CC22
        fac2   = fac / CC11

        pot = -coeff * prefac * fac

        FK = -coeff * prefac * ( fac1 * VJ - VI )
        FJ =  coeff * prefac * ( fac1 * VJ - fac2 * VI + VJ - VI )
        FI =  coeff * prefac * ( fac2 * VI - VJ )

        potential_energy = potential_energy + pot

        FX(I) = FX(I) + FI(1)
        FY(I) = FY(I) + FI(2)
        FZ(I) = FZ(I) + FI(3)

        FX(J) = FX(J) + FJ(1)
        FY(J) = FY(J) + FJ(2)
        FZ(J) = FZ(J) + FJ(3)

        FX(K) = FX(K) + FK(1)
        FY(K) = FY(K) + FK(2)
        FZ(K) = FZ(K) + FK(3)
        
    end subroutine bend_calculate_forces

    subroutine dih_calculate_forces( I, J, K, L, coeff )
        implicit none
        integer, intent(in) :: I, J, K, L
        real(8), intent(in) :: coeff
        real(8) :: VI(3), VJ(3), VK(3), FI(3), FJ(3), FK(3), FL(3)
        real(8) :: CC11, CC12, CC13, CC22, CC23, CC33
        real(8) :: DD11, DD12, DD13, DD22, DD23, DD33
        real(8) :: prefac, fac, fac1, fac2, fac3
        real(8) :: pot

        VI(1) = RX(J) - RX(I)
        VI(2) = RY(J) - RY(I)
        VI(3) = RZ(J) - RZ(I)

        VJ(1) = RX(K) - RX(J)
        VJ(2) = RY(K) - RY(J)
        VJ(3) = RZ(K) - RZ(J)

        VK(1) = RX(L) - RX(K)
        VK(2) = RY(L) - RY(K)
        VK(3) = RZ(L) - RZ(K)

        if (is_periodic) then
            VI = VI - ANINT( VI / box_length ) * box_length
            VJ = VJ - ANINT( VJ / box_length ) * box_length
            VK = VK - ANINT( VK / box_length ) * box_length
        endif

        CC11 = dot( VI, VI )
        CC12 = dot( VI, VJ )
        CC13 = dot( VI, VK )
        CC22 = dot( VJ, VJ )
        CC23 = dot( VJ, VK )
        CC33 = dot( VK, VK )

        DD11 = CC11 * CC11 - CC11 ** 2.0
        DD12 = CC11 * CC22 - CC12 ** 2.0
        DD13 = CC11 * CC33 - CC13 ** 2.0
        DD22 = CC22 * CC22 - CC22 ** 2.0
        DD23 = CC22 * CC33 - CC23 ** 2.0
        DD33 = CC33 * CC33 - CC33 ** 2.0

        prefac = 1.0 / SQRT( DD23 * DD12 )
        fac = CC23 * CC12 - CC13 * CC22
        fac1 = fac / DD23
        fac2 = fac / DD12

        pot = coeff * prefac * fac

        FL = -coeff * prefac * ( CC12 * VJ - CC22 * VI - fac1 * ( CC22 * VK- CC23 * VJ ) )

        FK = -coeff * prefac * ( CC12 * VK - CC12 * VJ + CC23 * VI + CC22 * VI - 2.0 * CC13 * VJ &
                    - fac2  * ( CC11 * VJ - CC12 * VI ) - fac1 * ( CC33 * VJ - CC22 * VK - CC23 * VK + CC23 * VJ ) )

        FJ = -coeff * prefac * ( -CC12 * VK + CC23 * VJ - CC23 * VI - CC22 * VK + 2.0 * CC13 * VJ &
                    - fac2  * (  CC22 * VI - CC11 * VJ - CC12 * VJ + CC12 * VI ) - fac1 * ( -CC33 * VJ + CC23 * VK ) )

        FI = -coeff * prefac * ( -CC23 * VJ + CC22 * VK - fac2 * ( -CC22 * VI + CC12 * VJ ) )

        potential_energy = potential_energy + pot

        FX(I) = FX(I) + FI(1)
        FY(I) = FY(I) + FI(2)
        FZ(I) = FZ(I) + FI(3)

        FX(J) = FX(J) + FJ(1)
        FY(J) = FY(J) + FJ(2)
        FZ(J) = FZ(J) + FJ(3)    

        FX(K) = FX(K) + FK(1)
        FY(K) = FY(K) + FK(2)
        FZ(K) = FZ(K) + FK(3)

        FX(L) = FX(L) + FL(1)
        FY(L) = FY(L) + FL(2)
        FZ(L) = FZ(L) + FL(3)

    end subroutine dih_calculate_forces

    subroutine angle_bend_calculate_forces( I, J, K, coeff, ang )
        implicit none
        integer, intent(in) :: I, J, K
        real(8), intent(in) :: coeff, ang
        real(8) :: VI(3), VJ(3), FI(3), FJ(3), FK(3)
        real(8) :: CC11, CC12, CC22
        real(8) :: prefac, fac, fac1, fac2, pot

        VI(1) = RX(J) - RX(I)
        VI(2) = RY(J) - RY(I)
        VI(3) = RZ(J) - RZ(I)

        VJ(1) = RX(K) - RX(J)
        VJ(2) = RY(K) - RY(J)
        VJ(3) = RZ(K) - RZ(J)

        CC11 = dot( VI, VI )
        CC12 = dot( VI, VJ )
        CC22 = dot( VJ, VJ )

        prefac = 1.0 / SQRT( CC11 * CC22 )
        fac    = CC12
        fac1   = fac / CC22
        fac2   = fac / CC11

        pot = -coeff * prefac * fac

        FK = -coeff * prefac * ( fac1 * VJ - VI )
        FJ =  coeff * prefac * ( fac1 * VJ - fac2 * VI + VJ - VI )
        FI =  coeff * prefac * ( fac2 * VI - VJ )

        potential_energy = potential_energy + pot

        FX(I) = FX(I) + FI(1)
        FY(I) = FY(I) + FI(2)
        FZ(I) = FZ(I) + FI(3)

        FX(J) = FX(J) + FJ(1)
        FY(J) = FY(J) + FJ(2)
        FZ(J) = FZ(J) + FJ(3)

        FX(K) = FX(K) + FK(1)
        FY(K) = FY(K) + FK(2)
        FZ(K) = FZ(K) + FK(3)
        
    end subroutine angle_bend_calculate_forces

end module bonded_interactions

module constraints
    use system
    implicit none

    contains
    subroutine apply_constraints_a( DT )
        implicit none
        real(8), intent(in) :: DT
        real(8) :: RIJ(3), RIJ_old(3), DR(3)   
        real(8) :: L_ij, chi, vdot, r_tol, len
        integer :: I, J, IX, rattle, max_rattle = 100
        logical :: moved(N)

        L_ij = 0.0; rattle = 0; 
        r_tol = 1.0e-5; max_rattle = 100
        moved = .TRUE.

        do while (any(moved) .and. rattle < max_rattle)

            do IX = 1, NC
                I   = BBI(IX)
                J   = BBJ(IX)
                len = BB_len(IX)
                
                if ( moved(I) .or. moved(J) ) then

                    RIJ(1) = RX(I) - RX(J)
                    RIJ(2) = RY(I) - RY(J)
                    RIJ(3) = RZ(I) - RZ(J)

                    if (is_periodic) then
                        RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
                    endif

                    chi = len ** 2 - SUM( RIJ ** 2)

                    if (abs(chi) > 2 * r_tol * len**2.0) then 

                        RIJ_old(1) = RX_old(I) - RX_old(J)
                        RIJ_old(2) = RY_old(I) - RY_old(J)
                        RIJ_old(3) = RZ_old(I) - RZ_old(J)

                        if (is_periodic) then
                            RIJ_old = RIJ_old - ANINT( RIJ_old / box_length ) * box_length
                        endif

                        vdot = dot( RIJ, RIJ_old )

                        ! write(*, *) chi, vdot, len, sum(RIJ**2), sum(RIJ_old**2)

                        if (vdot < r_tol * len ** 2.0) then
                            stop "Constraint failure. SB position multipliers exploded !!"
                        endif

                        L_ij = chi / ( 4.0 * vdot )
                        DR = L_ij * rij_old

                        RX(I) = RX(I) + DR(1)
                        RY(I) = RY(I) + DR(2)
                        RZ(I) = RZ(I) + DR(3)

                        RX(J) = RX(J) - DR(1)
                        RY(J) = RY(J) - DR(2)
                        RZ(J) = RZ(J) - DR(3)

                        VX(I) = VX(I) + DR(1) / DT
                        VY(I) = VY(I) + DR(2) / DT
                        VZ(I) = VZ(I) + DR(3) / DT

                        VX(J) = VX(J) - DR(1) / DT
                        VY(J) = VY(J) - DR(2) / DT
                        VZ(J) = VZ(J) - DR(3) / DT

                        moved(i) = .TRUE.
                        moved(j) = .TRUE.
                    else
                        moved(i) = .FALSE.
                        moved(j) = .FALSE.
                    endif

                endif

            enddo

            rattle = rattle + 1

        enddo

    end subroutine apply_constraints_a

    subroutine apply_constraints_b( DT )
        implicit none
        real(8), intent(in) :: DT
        real(8) :: RIJ(3), VIJ(3), DV(3)   
        real(8) :: L_ij, chi, vdot, r_tol, len
        integer :: I, J, IX, rattle, max_rattle = 100
        logical :: moved(N)

        L_ij = 0.0; rattle = 0; 
        r_tol = 1.0e-5; max_rattle = 100
        moved = .TRUE.

        do while (any(moved) .and. rattle < max_rattle)

            do IX = 1, NC
                I   = BBI(IX)
                J   = BBJ(IX)
                len = BB_len(IX)

                if( moved(I) .or. moved(J) ) then

                    RIJ(1) = RX(I) - RX(J)
                    RIJ(2) = RY(I) - RY(J)
                    RIJ(3) = RZ(I) - RZ(J)

                    if (is_periodic) then
                        RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
                    endif

                    VIJ(1) = VX(I) - VX(J)
                    VIJ(2) = VY(I) - VY(J)
                    VIJ(3) = VZ(I) - VZ(J)

                    vdot = dot( RIJ, VIJ )
                    L_ij = -vdot / ( 2.0 * len**2 )
                    
                    if (abs(L_ij) > r_tol) then

                        virial = virial + L_ij * len**2
                        DV = L_ij * RIJ

                        VX(I) = VX(I) + DV(1)
                        VY(I) = VY(I) + DV(2)
                        VZ(I) = VZ(I) + DV(3)

                        VX(J) = VX(J) - DV(1)
                        VY(J) = VY(J) - DV(2)
                        VZ(J) = VZ(J) - DV(3)

                        moved(i) = .TRUE.
                        moved(j) = .TRUE.

                    else

                        moved(i) = .FALSE.
                        moved(j) = .FALSE.

                    endif

                endif
                
            enddo

            rattle = rattle + 1

        enddo

    end subroutine apply_constraints_b
    
end module constraints