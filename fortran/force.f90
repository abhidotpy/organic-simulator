module lennard_jones
    use system
    implicit none

    contains
    pure subroutine lj_calculate_forces( I, J, cutoff, PE, FXIJ)
        implicit none
        integer, intent(in)         :: I, J
        real(real64), intent(in)    :: cutoff
        real(real64), intent(out)   :: PE, FXIJ(3)

        real(real64)                :: RIJ(3), epsilon, sigma
        real(real64)                :: rij_sq, rcut_sq, sr2_lj, coeff
        real(real64)                :: pot, pot_cut, sr2_cut

        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (is_periodic) then
            RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
        endif

        epsilon = sqrt( eshape_z(I) * eshape_z(J) )
        rij_sq = SUM( RIJ ** 2.0 )
        rcut_sq = cutoff ** 2.0

        if (rij_sq < rcut_sq) then
            sigma = ( shape_z(I) + shape_z(J) ) / 2.0
            sr2_lj = ( sigma ** 2.0 ) / rij_sq
            sr2_cut = ( sigma ** 2.0 ) / rcut_sq

            pot = 4.0 * epsilon * ( sr2_lj**6 - sr2_lj**3 )
            pot_cut = 4.0 * epsilon * ( sr2_cut**6 - sr2_cut**3 )
            coeff = 24.0 * epsilon * (2.0 * sr2_lj**6 - sr2_lj**3)

            PE = pot - pot_cut
            FXIJ = coeff * RIJ / rij_sq
        
        else
            PE = 0.0
            FXIJ = 0.0
        endif

    end subroutine lj_calculate_forces

end module lennard_jones

module lennard_jones12
    use system
    implicit none

    contains
    pure subroutine lj12_calculate_forces( I, J, cutoff, PE, FXIJ )
        implicit none
        integer, intent(in)        :: I, J
        real(real64), intent(in)   :: cutoff
        real(real64), intent(out)  :: PE, FXIJ(3)

        real(real64)               :: RIJ(3), epsilon, sigma
        real(real64)               :: rij_sq, rcut_sq, sr2_lj, coeff
        real(real64)               :: pot, pot_cut, sr2_cut


        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (is_periodic) then
            RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
        endif

        epsilon = sqrt( eshape_z(I) * eshape_z(J) )
        rij_sq = SUM( RIJ ** 2.0 )
        rcut_sq = cutoff ** 2.0

        if (rij_sq < rcut_sq) then
            sigma = ( shape_z(I) + shape_z(J) ) / 2.0
            sr2_lj = ( sigma ** 2.0 ) / rij_sq
            sr2_cut = ( sigma ** 2.0 ) / rcut_sq

            pot = epsilon * ( sr2_lj**6.0 - 2.0*sr2_lj**3.0 )
            pot_cut = epsilon * ( sr2_cut**6.0 - 2.0*sr2_cut**3.0 )
            coeff = 12.0 * epsilon * (sr2_lj**6 - sr2_lj**3)

            PE = pot - pot_cut
            FXIJ = coeff * RIJ / rij_sq

        else
            PE = 0.0
            FXIJ = 0.0
        endif

    end subroutine lj12_calculate_forces

end module lennard_jones12

module wca
    use system
    implicit none

    contains
    pure subroutine wca_calculate_forces( I, J, cutoff, PE, FXIJ )
        implicit none
        integer, intent(in)             :: I, J
        real(real64), intent(in)        :: cutoff
        real(real64), intent(out)       :: PE, FXIJ(3)
        real(real64), dimension(3)      :: RIJ, FIJ
        real(real64)                    :: rij_sq, rcut_sq, sigma_sq, sr2_lj
        real(real64)                    :: coeff, pot


        RIJ(1) = RX(I) - RX(J)
        RIJ(2) = RY(I) - RY(J)
        RIJ(3) = RZ(I) - RZ(J)

        if (is_periodic) then
            RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
        endif

        rij_sq = SUM( RIJ ** 2.0 )
        rcut_sq = cutoff ** 2.0

        if (rij_sq < rcut_sq) then
            sigma_sq = ( (shape_z(I) + shape_z(J)) / 2.0 ) ** 2.0
            sr2_lj = sigma_sq / rij_sq

            if ( rij_sq <= sigma_sq ) then
                PE = sigma_sq * (( sr2_lj**6.0 - 2.0*sr2_lj**3.0 ) + 1.0)
                coeff = 12.0 * sigma_sq * (sr2_lj**6 - sr2_lj**3)
                FXIJ = coeff * RIJ / rij_sq
            else
                PE = 0.0
                FXIJ = 0.0
            endif
        endif

    end subroutine wca_calculate_forces

end module wca

module morse
    use system
    implicit none

    contains
    pure subroutine morse_calculate_forces( I, J, dissoc, width, rmin, cutoff, PE, FXIJ )
        implicit none
        integer, intent(in)        :: I, J
        real(real64), dimension(3) :: RIJ, FIJ
        real(real64), intent(in)   :: dissoc, width, rmin, cutoff
        real(real64), intent(out)  :: PE, FXIJ(3)
        real(real64)               :: rij_mag, exp_term
        real(real64)               :: exp_cut, pot, pot_cut

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

            PE = pot - pot_cut
            FXIJ = -2.0 * width * dissoc * (1.0 - exp_term) * exp_term * RIJ / rij_mag

        else
            PE = 0.0
            FXIJ = 0.0
        endif

    end subroutine morse_calculate_forces

end module morse

module gay_berne
    use system
    implicit none
    
    type :: gb_parameters
        real(real64) :: ru1, ru2, uu, chi, xhi
    end type
    
    type :: gb_pair
        real(real64), dimension(3)   :: U1, U2, RIJ
        real(real64)                 :: rij_sq, rij_mag
    end type

    contains
    pure function g_func( chie, param, pair ) result (res)
        implicit none
        real(real64), intent(in)             :: chie
        type(gb_parameters), intent(in)      :: param
        type(gb_pair), intent(in)            :: pair
        real(real64)                         :: term1, term2, res

        term1 = (param % ru1 + param % ru2)**2 / (1 + chie * param % uu)
        term2 = (param % ru1 - param % ru2)**2 / (1 - chie * param % uu)
        
        res = 1 - (chie / 2.0 / pair % rij_sq) * (term1 + term2)

    end function g_func

    pure function dG_dr( chie, param, pair ) result (res)
        implicit none
        real(real64), intent(in)             :: chie
        type(gb_parameters), intent(in) :: param
        type(gb_pair), intent(in)       :: pair
        real(real64), dimension(3)           :: term1, term2, res

        term1 = ((param % ru1 + param % ru2) / (1 + chie * param % uu)) * (pair % U1 + pair % U2)
        term2 = ((param % ru1 - param % ru2) / (1 - chie * param % uu)) * (pair % U1 - pair % U2)
        
        res = (-chie / pair % RIJ_SQ ) * (term1 + term2) + (2 * pair % RIJ / pair % RIJ_SQ) * (1 - g_func( chie, param, pair ))

    end function dG_dr

    pure function dG_du1( chie, param, pair ) result(res)
        implicit none
        real(real64), intent(in)             :: chie
        type(gb_parameters), intent(in) :: param
        type(gb_pair), intent(in)       :: pair
        real(real64)                         :: term1, term2, res(3)
        
        term1 = ((param % ru1 + param % ru2) / (1 + chie * param % uu))
        term2 = ((param % ru1 - param % ru2) / (1 - chie * param % uu))
        
        res = (-chie * pair % RIJ / pair % rij_sq ) * (term1 + term2) + ((chie**2 * pair % U2) / (2 * pair % rij_sq)) * (term1 ** 2 - term2 ** 2)
        
    end function dG_du1

    pure function dG_du2( chie, param, pair ) result(res)
        implicit none
        real(real64), intent(in)             :: chie
        type(gb_parameters), intent(in) :: param
        type(gb_pair), intent(in)       :: pair
        real(real64)                         :: term1, term2, res(3)
        
        term1 = ((param % ru1 + param % ru2) / (1 + chie * param % uu))
        term2 = ((param % ru1 - param % ru2) / (1 - chie * param % uu))
        
        res = (-chie * pair % RIJ / pair % RIJ_SQ ) * (term1 - term2) + ((chie**2 * pair % U1) / (2 * pair % rij_sq)) * (term1 ** 2 - term2 ** 2)
        
    end function dG_du2

    pure subroutine gb_calculate_forces( I, J, cutoff, PE, FXIJ, TXI, TXJ, width )
        implicit none
        integer, intent(in)             :: I, J
        real(real64), intent(in)             :: cutoff
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        real(real64), optional, intent(in)   :: width
        
        type(gb_parameters) :: param
        type(gb_pair)       :: pair
        
        real(real64), dimension(3)           :: dR_dr, de2_dr, dR_du1, dR_du2
        real(real64), dimension(3)           :: de1_du1, de2_du1, de1_du2, de2_du2
        real(real64), dimension(3)           :: FIJ, TI, TJ, TORQ1, TORQ2

        real(real64)                         :: meu, neu
        real(real64)                         :: eps1, eps2, sr, sr_cut, sigma, sigma_0
        real(real64)                         :: pot, wf, srx_fac, srd_fac
        real(real64)                         :: switch, xx, sx, sdx, factor
        logical                         :: imp_wat

        meu = 1.0; neu = 2.0; sigma_0 = 1.0; imp_wat = .FALSE.

        if (present(width)) then
            wf = width
            imp_wat = .TRUE.
        else
            wf = 0.0
            imp_wat = .FALSE.
        endif

        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (is_periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box_length ) * box_length
        endif

        pair % RIJ_SQ = SUM( pair % RIJ**2 )
        pair % RIJ_MAG = SQRT( pair % RIJ_SQ )

        pair % U1(1) = 2.0 *  ( QX(I) * QZ(I) + QW(I) * QY(I) )
        pair % U1(2) = 2.0 *  ( QY(I) * QZ(I) - QW(I) * QX(I) )
        pair % U1(3) = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

        pair % U2(1) = 2.0 *  ( QX(J) * QZ(J) + QW(J) * QY(J) )
        pair % U2(2) = 2.0 *  ( QY(J) * QZ(J) - QW(J) * QX(J) )
        pair % U2(3) = QW(J)**2 - QX(J)**2 - QY(J)**2 + QZ(J)**2

        ! CALCULATE GAY BERNE PARAMETERS

        param % ru1 = dot( pair % RIJ, pair % U1 )
        param % ru2 = dot( pair % RIJ, pair % U2 )
        param % uu  = dot( pair % U1, pair % U2 )
        param % chi = ( shape_z(I) ** 2 - shape_x(I) ** 2 ) / ( shape_z(I) ** 2 + shape_x(I) ** 2 )
        param % xhi = ( eshape_x(I) ** (1.0/meu) - eshape_z(I) ** (1.0/meu) ) / ( eshape_x(I) ** (1.0/meu) + eshape_z(I) ** (1.0/meu) )

        eps1 = ( 1 - ( param % chi * param % uu ) ** 2) ** (-1.0 / 2.0)
        eps2 = g_func( param % xhi, param, pair )
        sigma = shape_x(I) * g_func( param % chi, param, pair ) ** (-1.0 / 2.0)

        switch = ( cutoff + wf ) - 2.0 * sigma_0    ! Apply switch function 2sigma_0 inside cutoff
        factor = 1.0 / ( ( cutoff + wf ) - switch )
        xx = ( pair % RIJ_MAG - switch ) * factor
        sx = 1.0 - 6.0 * xx ** 5.0 + 15.0 * xx ** 4.0 - 10.0 * xx ** 3.0
        sdx = -30.0 * xx ** 4.0 + 60.0 * xx ** 3.0 - 30.0 * xx ** 2.0

        if ( pair % RIJ_MAG < ( sigma + cutoff + wf ) ) then

            if (imp_wat) then
                if ( pair % RIJ_MAG < sigma ) then
                    sr = sigma_0 / ( pair % RIJ_MAG - sigma + sigma_0 )
                else if ((pair % RIJ_MAG .ge. sigma) .and. (pair % RIJ_MAG .le. (sigma + wf))) then
                    sr = sigma_0 / ( sigma - sigma + sigma_0 )
                else
                    sr = sigma_0 / ( pair % RIJ_MAG - sigma + sigma_0 - wf )
                endif
            else
                sr = sigma_0 / ( pair % RIJ_MAG - sigma + sigma_0 )
            endif

            srx_fac = (sr ** 12.0 - 2.0 * sr ** 6.0)
            srd_fac = (12.0 * sr**7 - 12.0 * sr**13)

            dR_dr   = (1.0 / sigma_0) * (( pair % RIJ / pair % RIJ_MAG ) + 0.5 * shape_x(I) * dG_dr( param % chi, param, pair ) * g_func( param % chi, param, pair ) ** (-3.0 / 2.0))
            de2_dr  = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_dr( param % xhi, param, pair )

            de1_du1 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (param % chi ** 2 * param % uu * pair % U2)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du1( param % xhi, param, pair )
            dR_du1  = srd_fac * 0.5 * (shape_x(I) / sigma_0) * dG_du1( param % chi, param, pair ) * g_func( param % chi, param, pair ) ** (-3.0 / 2.0)

            de1_du2 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (param % chi ** 2 * param % uu * pair % U1)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du2( param % xhi, param, pair )
            dR_du2  = srd_fac * 0.5 * (shape_x(I) / sigma_0) * dG_du2( param % chi, param, pair ) * g_func( param % chi, param, pair ) ** (-3.0 / 2.0)
        
            pot = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * srd_fac * dR_dr  + de2_dr * srx_fac)
            TI  = (-1.0) * (srx_fac * (de1_du1 + de2_du1) + ((eps1 ** neu) * (eps2 ** meu) * dR_du1))
            TJ  = (-1.0) * (srx_fac * (de1_du2 + de2_du2) + ((eps1 ** neu) * (eps2 ** meu) * dR_du2))
            TORQ1 = cross( pair % U1, TI )
            TORQ2 = cross( pair % U2, TJ )

            if ((pair % rij_mag .ge. switch) .and. (pair % rij_mag .le. (cutoff + wf))) then
                pot   = sx * pot
                FIJ   = sx * FIJ   - sdx * pot * factor
                TORQ1 = sx * TORQ1 - sdx * pot * factor
                TORQ2 = sx * TORQ2 - sdx * pot * factor
            endif

        else
            pot = 0.0
            FIJ = 0.0
            TORQ1 = 0.0
            TORQ2 = 0.0
        endif

        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2

    end subroutine gb_calculate_forces

end module gay_berne

module ecp
    use system
    implicit none

    type :: ecp_pair
        real(real64), dimension(3, 3) :: AE, BE, AR, BR
        real(real64), dimension(3)    :: RIJ, RIJ_hat
        real(real64)                  :: RIJ_SQ, RIJ_MAG
    end type

    contains
    pure real(real64) function brent(func, pair) result (res)
        implicit none
        type(ecp_pair), intent(in)  :: pair
        integer, parameter          :: max_iter = 10000
        real(real64), parameter     :: eps = 1e-8, psi = 0.5 * ( 3.0 - sqrt(5.0) )
        real(real64)                :: a, b, c, x, w, v, u
        real(real64)                :: fa, fb, fc, ft, fw, fv, fu
        real(real64)                :: deltax, atol, rtol
        real(real64)                :: tol1, tol2, xmid, tmp1, tmp2, rat, p, dx_temp
        integer                     :: I, iter

        interface
            pure real(real64) function func(xx, pp)
                import ecp_pair, real64
                real(real64), intent(in)        :: xx
                type(ecp_pair), intent(in) :: pp
            end function func
        end interface

        deltax = 0.0; atol = 1e-11; rtol = 1e-8; iter = 0
        a = 0.0; b = 1.0 
        c = (1 - psi) * a + psi * b
        fa = func(a, pair); fb = func(b, pair); fc = func(c, pair)
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

            fu = func(u, pair)

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

    pure real(real64) function optim_eps( L, pair )
        implicit none
        real(real64), intent(in)         :: L
        type(ecp_pair), intent(in)  :: pair
        real(real64), dimension(3, 3)    :: GM
        real(real64), dimension(3)       :: KM

        GM = ( 1 - L ) * pair % AE + L * pair % BE
        KM = inverse(GM) .x. pair % RIJ

        optim_eps = -L * ( 1 - L ) * dot( pair % RIJ, KM )
        
    end function optim_eps

    pure real(real64) function optim_dist( L, pair )
        implicit none
        real(real64), intent(in)         :: L
        type(ecp_pair), intent(in)  :: pair
        real(real64), dimension(3, 3)    :: GM_R
        real(real64), dimension(3)       :: KM_R

        GM_R = ( 1 - L ) * pair % AR + L * pair % BR
        KM_R = inverse(GM_R) .x. pair % RIJ

        optim_dist = -L * ( 1 - L ) * dot( pair % RIJ, KM_R )
        
    end function optim_dist

    pure subroutine ecp_calculate_forces( I, J, cutoff, PE, FXIJ, TXI, TXJ, width )
        implicit none
        integer, intent(in)             :: I, J
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        real(real64), intent(in)             :: cutoff
        real(real64), optional, intent(in)   :: width

        type(ecp_pair) :: pair

        real(real64), dimension(3, 3) :: U1, U2, S1, S2, E1, E2, EM
        real(real64), dimension(3, 3) :: GE, GR, MM1, MM2
        real(real64), dimension(3)    :: KE, KR
        real(real64)                  :: lambda_E, lambda_R
        real(real64)                  :: meu, neu

        real(real64), dimension(3)    :: dR_dr, de2_dr
        real(real64), dimension(3)    :: de1_du1x, de1_du1y, de1_du1z, de1_du1, de2_du1, dR_du1
        real(real64), dimension(3)    :: de1_du2x, de1_du2y, de1_du2z, de1_du2, de2_du2, dR_du2
        real(real64), dimension(3)    :: FIJ, TORQ1, TORQ2

        real(real64)                  :: eps1, eps2, sr, phi, sigma, sigma_0
        real(real64)                  :: e10, e20, sigma_t1, sigma_t2, sigma_t
        real(real64)                  :: pot, wf, efac, rfac, srx_fac, srd_fac
        real(real64)                  :: switch, xx, sx, sdx, factor
        logical                  :: imp_wat, ierr_e, ierr_r

        meu = 1.0; neu = 2.0; sigma_0 = 1.0; 
        imp_wat = .FALSE.; ierr_e = .FALSE.; ierr_r = .FALSE.
        if (present(width)) then
            wf = width
            imp_wat = .TRUE.
        else
            wf = 0.0
            imp_wat = .FALSE.
        endif
        
        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (is_periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box_length ) * box_length
        endif

        pair % RIJ_SQ = SUM( pair % RIJ**2 )
        pair % RIj_MAG = SQRT( pair % RIJ_SQ )
        pair % RIJ_hat = pair % RIJ / pair % RIj_MAG

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

        MM1 = transpose(U1) .x. S1
        MM2 = transpose(U2) .x. S2
        pair % AR = transpose(U1) .x. ( S1 ** 2.0 ) .x. U1
        pair % BR = transpose(U2) .x. ( S2 ** 2.0 ) .x. U2
        pair % AE = transpose(U1) .x. E1 .x. U1
        pair % BE = transpose(U2) .x. E2 .x. U2

        sigma_t1 = (shape_x(I) * shape_y(I) + shape_z(I) ** 2) * sqrt( 2.0 * shape_x(I) * shape_y(I) ) / 8.0
        sigma_t2 = (shape_x(J) * shape_y(J) + shape_z(J) ** 2) * sqrt( 2.0 * shape_x(J) * shape_y(J) ) / 8.0
        sigma_t = sqrt( sigma_t1 * sigma_t2 )
        
        EM = inverse( pair % AR + pair % BR )
        eps1 = sigma_t * sqrt( determinant( EM ) )
        
        lambda_E = brent( optim_eps, pair )
        if ( lambda_E < 0.0 .or. lambda_E > 1.0 ) then
            ierr_e = .TRUE.
        endif
        
        GE = ( 1 - lambda_E ) * pair % AE + lambda_E * pair % BE
        KE = inverse( GE ) .x. pair % RIJ
        eps2 = lambda_E * (1.0 - lambda_E) * dot( pair % RIJ_hat, KE ) / pair % RIJ_mag
        
        lambda_R = brent( optim_dist, pair )
        if ( lambda_R < 0.0 .or. lambda_R > 1.0 ) then
            ierr_r = .TRUE.
        endif
        
        GR = ( 1 - lambda_R ) * pair % AR + lambda_R * pair % BR
        KR = inverse( GR ) .x. pair % RIJ
        phi = lambda_R * (1.0 - lambda_R) * dot( pair % RIJ_hat, KR ) / pair % RIJ_mag
        sigma = phi ** (-1.0/2.0)
        efac = 2.0 * lambda_E * ( 1.0 - lambda_E ) / pair % RIJ_SQ
        rfac = 2.0 * lambda_R * ( 1.0 - lambda_R ) / pair % RIJ_SQ

        switch = ( cutoff + wf ) - 2.0 * sigma_0    ! Apply switch function 2sigma_0 inside cutoff
        factor = 1.0 / ( ( cutoff + wf ) - switch )
        xx = ( pair % RIJ_MAG - switch ) * factor
        sx = 1.0 - 6.0 * xx ** 5.0 + 15.0 * xx ** 4.0 - 10.0 * xx ** 3.0
        sdx = -30.0 * xx ** 4.0 + 60.0 * xx ** 3.0 - 30.0 * xx ** 2.0

        if ( pair % rij_mag < ( sigma + cutoff + wf ) ) then
        
            if (imp_wat) then
                if ( pair % rij_mag < sigma ) then
                    sr = sigma_0 / ( pair % rij_mag - sigma + sigma_0 )
                else if ((pair % rij_mag .ge. sigma) .and. (pair % rij_mag .le. (sigma + wf))) then
                    sr = sigma_0 / ( sigma - sigma + sigma_0 )
                else
                    sr = sigma_0 / ( pair % rij_mag - sigma + sigma_0 - wf )
                endif
            else
                sr = sigma_0 / ( pair % rij_mag - sigma + sigma_0 )
            endif
            
            srx_fac = (sr ** 12.0 - 2.0 * sr ** 6.0)
            srd_fac = (12.0 * sr**7 - 12.0 * sr**13)

            dR_dr = (1.0 / sigma_0) * (  ( pair % RIJ_hat )     + 0.5 * rfac * ( KR - dot_product( pair % RIJ_hat, KR ) * pair % RIJ_hat ) * phi ** (-3.0 / 2.0)  )
            de2_dr = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * efac * ( KE - dot_product( pair % RIJ_hat, KE ) * pair % RIJ_hat )

            de1_du1x = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM1(:, 1), ( EM .x. MM1(:, 1) ) )
            de1_du1y = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM1(:, 2), ( EM .x. MM1(:, 2) ) )
            de1_du1z = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM1(:, 3), ( EM .x. MM1(:, 3) ) )

            de1_du1 = (de1_du1x + de1_du1y + de1_du1z)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * srx_fac * (1 - lambda_E) * efac * cross( KE .x. pair % AE, KE )
            dR_du1 = (eps1 ** neu) * (eps2 ** meu) * srd_fac * (0.5 / sigma_0) * (1 - lambda_R) * rfac * cross( KR .x. pair % AR, KR ) * phi ** (-3.0 / 2.0)

            de1_du2x = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM2(:, 1), ( EM .x. MM2(:, 1) ) )
            de1_du2y = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM2(:, 2), ( EM .x. MM2(:, 2) ) )
            de1_du2z = (eps2 ** meu) * (neu * eps1 ** neu) * srx_fac * cross( MM2(:, 3), ( EM .x. MM2(:, 3) ) )

            de1_du2 = (de1_du2x + de1_du2y + de1_du2z)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * srx_fac * lambda_E * efac * cross( KE .x. pair % BE, KE )
            dR_du2 = (eps1 ** neu) * (eps2 ** meu) * srd_fac * (0.5 / sigma_0) * lambda_R * rfac * cross( KR .x. pair % BR, KR ) * phi ** (-3.0 / 2.0)
            
            pot   = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ   = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * srd_fac * dR_dr  + srx_fac * de2_dr)
            TORQ1 = (de1_du1 + de2_du1 + dR_du1)
            TORQ2 = (de1_du2 + de2_du2 + dR_du2)

            if ((pair % rij_mag .ge. switch) .and. (pair % rij_mag .le. (cutoff + wf))) then
                pot   = sx * pot
                FIJ   = sx * FIJ   - sdx * pot * factor
                TORQ1 = sx * TORQ1 - sdx * pot * factor
                TORQ2 = sx * TORQ2 - sdx * pot * factor
            endif

        else
            pot = 0.0
            FIJ = 0.0
            TORQ1 = 0.0
            TORQ2 = 0.0
        endif
        
        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2
        
    end subroutine ecp_calculate_forces
end module ecp

module gay_berne_chiral
    use system
    implicit none
    
    type :: gb_parameters
        real(real64) :: ru1, ru2, uu, chi, xhi
    end type
    
    type :: gb_pair
        real(real64), dimension(3)   :: U1, U2, RIJ, RIJ_hat
        real(real64)                 :: rij_sq, rij_mag
    end type

    contains
    pure function g_func( chie, param, pair ) result (res)
        implicit none
        real(real64), intent(in)            :: chie
        type(gb_parameters), intent(in)     :: param
        type(gb_pair), intent(in)           :: pair
        real(real64)                        :: term1, term2, res

        term1 = (param % ru1 + param % ru2)**2 / (1 + chie * param % uu)
        term2 = (param % ru1 - param % ru2)**2 / (1 - chie * param % uu)
        
        res = 1 - (chie / 2.0 / pair % rij_sq) * (term1 + term2)

    end function g_func

    pure function dG_dr( chie, param, pair ) result (res)
        implicit none
        real(real64), intent(in)             :: chie
        type(gb_parameters), intent(in) :: param
        type(gb_pair), intent(in)       :: pair
        real(real64), dimension(3)           :: term1, term2, res

        term1 = ((param % ru1 + param % ru2) / (1 + chie * param % uu)) * (pair % U1 + pair % U2)
        term2 = ((param % ru1 - param % ru2) / (1 - chie * param % uu)) * (pair % U1 - pair % U2)
        
        res = (-chie / pair % RIJ_SQ ) * (term1 + term2) + (2 * pair % RIJ / pair % RIJ_SQ) * (1 - g_func( chie, param, pair ))

    end function dG_dr

    pure function dG_du1( chie, param, pair ) result(res)
        implicit none
        real(real64), intent(in)             :: chie
        type(gb_parameters), intent(in) :: param
        type(gb_pair), intent(in)       :: pair
        real(real64)                         :: term1, term2, res(3)
        
        term1 = ((param % ru1 + param % ru2) / (1 + chie * param % uu))
        term2 = ((param % ru1 - param % ru2) / (1 - chie * param % uu))
        
        res = (-chie * pair % RIJ / pair % rij_sq ) * (term1 + term2) + ((chie**2 * pair % U2) / (2 * pair % rij_sq)) * (term1 ** 2 - term2 ** 2)
        
    end function dG_du1

    pure function dG_du2( chie, param, pair ) result(res)
        implicit none
        real(real64), intent(in)             :: chie
        type(gb_parameters), intent(in) :: param
        type(gb_pair), intent(in)       :: pair
        real(real64)                         :: term1, term2, res(3)
        
        term1 = ((param % ru1 + param % ru2) / (1 + chie * param % uu))
        term2 = ((param % ru1 - param % ru2) / (1 - chie * param % uu))
        
        res = (-chie * pair % RIJ / pair % RIJ_SQ ) * (term1 - term2) + ((chie**2 * pair % U1) / (2 * pair % rij_sq)) * (term1 ** 2 - term2 ** 2)
        
    end function dG_du2

    pure subroutine gbc_calculate_forces( I, J, cutoff, PE, FXIJ, TXI, TXJ, width, chirality )
        implicit none
        integer, intent(in)             :: I, J
        real(real64), intent(in)             :: cutoff
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        real(real64), optional, intent(in)   :: width, chirality
        
        type(gb_parameters)             :: param
        type(gb_pair)                   :: pair
        
        real(real64), dimension(3)           :: dR_dr, de2_dr, dR_du1, dR_du2, U1XU2
        real(real64), dimension(3)           :: de1_du1, de2_du1, de1_du2, de2_du2
        real(real64), dimension(3)           :: duur_dr, duur_du1, duur_du2, duu_du1, duu_du2
        real(real64), dimension(3)           :: FIJ, TI, TJ, TORQ1, TORQ2
        real(real64), dimension(3)           :: FIJ_chiral, TI_chiral, TJ_chiral, TORQ1_chiral, TORQ2_chiral

        real(real64)                         :: meu, neu, uur, uu
        real(real64)                         :: eps1, eps2, sr, sr_cut, sigma, sigma_0
        real(real64)                         :: pot, pot_chiral, wf, ch_str
        real(real64)                         :: switch, xx, sx, sdx, factor, srx_fac, srd_fac
        logical                         :: imp_wat

        meu = 1.0; neu = 2.0; sigma_0 = 1.0; imp_wat = .FALSE.

        if (present(width)) then
            wf = width
            imp_wat = .TRUE.
        else
            wf = 0.0
            imp_wat = .FALSE.
        endif

        if (present(chirality)) then
            ch_str = chirality
        else
            ch_str = 0.0
        endif

        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (is_periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box_length ) * box_length
        endif

        pair % RIJ_SQ = SUM( pair % RIJ**2 )
        pair % RIJ_MAG = SQRT( pair % RIJ_SQ )
        pair % RIJ_hat = pair % RIJ / pair % RIJ_MAG

        pair % U1(1) = 2.0 *  ( QX(I) * QZ(I) + QW(I) * QY(I) )
        pair % U1(2) = 2.0 *  ( QY(I) * QZ(I) - QW(I) * QX(I) )
        pair % U1(3) = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

        pair % U2(1) = 2.0 *  ( QX(J) * QZ(J) + QW(J) * QY(J) )
        pair % U2(2) = 2.0 *  ( QY(J) * QZ(J) - QW(J) * QX(J) )
        pair % U2(3) = QW(J)**2 - QX(J)**2 - QY(J)**2 + QZ(J)**2

        ! CALCULATE GAY BERNE PARAMETERS

        param % ru1 = dot( pair % RIJ, pair % U1 )
        param % ru2 = dot( pair % RIJ, pair % U2 )
        param % uu  = dot( pair % U1, pair % U2 )
        param % chi = ( shape_z(I) ** 2 - shape_x(I) ** 2 ) / ( shape_z(I) ** 2 + shape_x(I) ** 2 )
        param % xhi = ( eshape_x(I) ** (1.0/meu) - eshape_z(I) ** (1.0/meu) ) / ( eshape_x(I) ** (1.0/meu) + eshape_z(I) ** (1.0/meu) )

        eps1 = ( 1 - ( param % chi * param % uu ) ** 2) ** (-1.0 / 2.0)
        eps2 = g_func( param % xhi, param, pair )
        sigma = shape_x(I) * g_func( param % chi, param, pair ) ** (-1.0 / 2.0)

        switch = ( cutoff + wf ) - 2.0 * sigma_0    ! Apply switch function 2sigma_0 inside cutoff
        factor = 1.0 / ( ( cutoff + wf ) - switch )
        xx = ( pair % RIJ_MAG - switch ) * factor
        sx = 1.0 - 6.0 * xx ** 5.0 + 15.0 * xx ** 4.0 - 10.0 * xx ** 3.0
        sdx = -30.0 * xx ** 4.0 + 60.0 * xx ** 3.0 - 30.0 * xx ** 2.0

        if ( pair % RIJ_MAG < ( sigma + cutoff + wf ) ) then

            if (imp_wat) then
                if ( pair % RIJ_MAG < sigma ) then
                    sr = sigma_0 / ( pair % RIJ_MAG - sigma + sigma_0 )
                else if ((pair % RIJ_MAG .ge. sigma) .and. (pair % RIJ_MAG .le. (sigma + wf))) then
                    sr = sigma_0 / ( sigma - sigma + sigma_0 )
                else
                    sr = sigma_0 / ( pair % RIJ_MAG - sigma + sigma_0 - wf )
                endif
            else
                sr = sigma_0 / ( pair % RIJ_MAG - sigma + sigma_0 )
            endif

            ! EVALUATE NORMAL GB POTENTIAL

            srx_fac = (sr ** 12.0 - 2.0 * sr ** 6.0)
            srd_fac = (12.0 * sr**7 - 12.0 * sr**13)

            dR_dr   = (eps1 ** neu) * (eps2 ** meu) * (1.0 / sigma_0) * (( pair % RIJ / pair % RIJ_MAG ) + 0.5 * shape_x(I) * dG_dr( param % chi, param, pair ) * g_func( param % chi, param, pair ) ** (-3.0 / 2.0))
            de2_dr  = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_dr( param % xhi, param, pair )

            de1_du1 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (param % chi ** 2 * param % uu * pair % U2)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du1( param % xhi, param, pair )
            dR_du1  = (eps1 ** neu) * (eps2 ** meu) * 0.5 * (shape_x(I) / sigma_0) * dG_du1( param % chi, param, pair ) * g_func( param % chi, param, pair ) ** (-3.0 / 2.0)

            de1_du2 = (eps2 ** meu) * (neu * eps1 ** (neu + 2)) * (param % chi ** 2 * param % uu * pair % U1)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * dG_du2( param % xhi, param, pair )
            dR_du2  = (eps1 ** neu) * (eps2 ** meu) * 0.5 * (shape_x(I) / sigma_0) * dG_du2( param % chi, param, pair ) * g_func( param % chi, param, pair ) ** (-3.0 / 2.0)
        
            pot = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ = (-1.0) * ( srx_fac *  de2_dr             + srd_fac * dR_dr  )
            TI  = (-1.0) * ( srx_fac * (de1_du1 + de2_du1) + srd_fac * dR_du1 )
            TJ  = (-1.0) * ( srx_fac * (de1_du2 + de2_du2) + srd_fac * dR_du2 )
            TORQ1 = cross( pair % U1, TI )
            TORQ2 = cross( pair % U2, TJ )

            ! EVALUATE CHIRAL TERMS

            ! srx_fac = (sr ** 7.0)
            ! srd_fac = (-7.0 * sr ** 8.0)
            u1xu2   = cross( pair % u1, pair % u2 )
            uur     = dot_product( U1XU2, pair % RIJ_hat )
            uu      = param % uu

            dR_dr   = dR_dr  * srd_fac * uur * uu
            de2_dr  = de2_dr * srx_fac * uur * uu
            duur_dr = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * ( U1XU2 - dot( U1XU2, pair % RIJ_hat ) * pair % RIJ_hat ) / pair % rij_mag

            de1_du1  = de1_du1 * srx_fac * uur * uu
            de2_du1  = de2_du1 * srx_fac * uur * uu
            dR_du1   = dR_du1  * srd_fac * uur * uu
            duur_du1 = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * cross( pair % U2, pair % RIJ_hat )
            duu_du1  = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * pair % U2

            de1_du2  = de1_du2 * srx_fac * uur * uu
            de2_du2  = de2_du2 * srx_fac * uur * uu
            dR_du2   = dR_du2  * srd_fac * uur * uu
            duur_du2 = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * cross( pair % RIJ_hat, pair % U1 )
            duu_du2  = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * pair % U1

            pot_chiral   = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * uu
            FIJ_chiral   = (-1.0) * (          de2_dr  + dR_dr  + duur_dr           )
            TI           = (-1.0) * (de1_du1 + de2_du1 + dR_du1 + duur_du1 + duu_du1)
            TJ           = (-1.0) * (de1_du2 + de2_du2 + dR_du2 + duur_du2 + duu_du2)
            TORQ1_chiral = cross( pair % U1, TI )
            TORQ2_chiral = cross( pair % U2, TJ )

            pot   = pot   + ch_str * pot_chiral
            FIJ   = FIJ   + ch_str * FIJ_chiral
            TORQ1 = TORQ1 + ch_str * TORQ1_chiral
            TORQ2 = TORQ2 + ch_str * TORQ2_chiral

            if ((pair % rij_mag .ge. switch) .and. (pair % rij_mag .le. (cutoff + wf))) then
                pot   = sx * pot
                FIJ   = sx * FIJ   - sdx * pot * factor
                TORQ1 = sx * TORQ1 - sdx * pot * factor
                TORQ2 = sx * TORQ2 - sdx * pot * factor
            endif

        else
            pot = 0.0
            FIJ = 0.0
            TORQ1 = 0.0
            TORQ2 = 0.0
        endif

        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2

    end subroutine gbc_calculate_forces

end module gay_berne_chiral

module ecp_chiral
    use system
    implicit none

    type :: ecp_pair
        real(real64), dimension(3, 3) :: AE, BE, AR, BR
        real(real64), dimension(3)    :: RIJ, RIJ_hat
        real(real64)                  :: RIJ_SQ, RIJ_MAG
    end type

    contains
    pure real(real64) function brent(func, pair) result (res)
        implicit none
        type(ecp_pair), intent(in)  :: pair
        integer, parameter          :: max_iter = 10000
        real(real64), parameter     :: eps = 1e-8, psi = 0.5 * ( 3.0 - sqrt(5.0) )
        real(real64)                :: a, b, c, x, w, v, u
        real(real64)                :: fa, fb, fc, ft, fw, fv, fu
        real(real64)                :: deltax, atol, rtol
        real(real64)                :: tol1, tol2, xmid, tmp1, tmp2, rat, p, dx_temp
        integer                     :: I, iter

        interface
            pure real(real64) function func(xx, pp)
                import ecp_pair, real64
                real(real64), intent(in)        :: xx
                type(ecp_pair), intent(in) :: pp
            end function func
        end interface

        deltax = 0.0; atol = 1e-11; rtol = 1e-8; iter = 0
        a = 0.0; b = 1.0 
        c = (1 - psi) * a + psi * b
        fa = func(a, pair); fb = func(b, pair); fc = func(c, pair)
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

            fu = func(u, pair)

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

    pure real(real64) function optim_eps( L, pair )
        implicit none
        real(real64), intent(in)         :: L
        type(ecp_pair), intent(in)  :: pair
        real(real64), dimension(3, 3)    :: GM
        real(real64), dimension(3)       :: KM

        GM = ( 1 - L ) * pair % AE + L * pair % BE
        KM = inverse(GM) .x. pair % RIJ

        optim_eps = -L * ( 1 - L ) * dot( pair % RIJ, KM )
        
    end function optim_eps

    pure real(real64) function optim_dist( L, pair )
        implicit none
        real(real64), intent(in)         :: L
        type(ecp_pair), intent(in)  :: pair
        real(real64), dimension(3, 3)    :: GM_R
        real(real64), dimension(3)       :: KM_R

        GM_R = ( 1 - L ) * pair % AR + L * pair % BR
        KM_R = inverse(GM_R) .x. pair % RIJ

        optim_dist = -L * ( 1 - L ) * dot( pair % RIJ, KM_R )
        
    end function optim_dist

    pure subroutine ecpc_calculate_forces( I, J, cutoff, PE, FXIJ, TXI, TXJ, width, chirality )
        implicit none
        integer, intent(in)             :: I, J
        real(real64), intent(out)            :: PE, FXIJ(3), TXI(3), TXJ(3)
        real(real64), intent(in)             :: cutoff
        real(real64), optional, intent(in)   :: width, chirality

        type(ecp_pair) :: pair

        real(real64), dimension(3, 3) :: U1, U2, S1, S2, E1, E2, EM
        real(real64), dimension(3, 3) :: GE, GR, MM1, MM2
        real(real64), dimension(3)    :: KE, KR, U1XU2
        real(real64)                  :: lambda_E, lambda_R
        real(real64)                  :: meu, neu, uur, uu

        real(real64), dimension(3)    :: dR_dr, de2_dr
        real(real64), dimension(3)    :: de1_du1x, de1_du1y, de1_du1z, de1_du1, de2_du1, dR_du1
        real(real64), dimension(3)    :: de1_du2x, de1_du2y, de1_du2z, de1_du2, de2_du2, dR_du2
        real(real64), dimension(3)    :: duur_dr, duur_du1z, duur_du2z, duur_du1, duur_du2
        real(real64), dimension(3)    :: duu_du1z, duu_du2z, duu_du1, duu_du2
        real(real64), dimension(3)    :: FIJ, TORQ1, TORQ2
        real(real64), dimension(3)    :: FIJ_chiral, TORQ1_chiral, TORQ2_chiral

        real(real64)                  :: eps1, eps2, sr, phi, sigma, sigma_0, ch_str
        real(real64)                  :: e10, e20, sigma_t1, sigma_t2, sigma_t
        real(real64)                  :: pot, wf, efac, rfac, srx_fac, srd_fac
        real(real64)                  :: switch, xx, sx, sdx, factor, pot_chiral
        logical                       :: imp_wat, ierr_e, ierr_r

        meu = 1.0; neu = 2.0; sigma_0 = 1.0; 
        imp_wat = .FALSE.; ierr_e = .FALSE.; ierr_r = .FALSE.
        if (present(width)) then
            wf = width
            imp_wat = .TRUE.
        else
            wf = 0.0
            imp_wat = .FALSE.
        endif

        if (present(chirality)) then
            ch_str = chirality
        else
            ch_str = 0.0
        endif
        
        pair % RIJ(1) = RX(I) - RX(J)
        pair % RIJ(2) = RY(I) - RY(J)
        pair % RIJ(3) = RZ(I) - RZ(J)

        if (is_periodic) then
            pair % RIJ = pair % RIJ - ANINT( pair % RIJ / box_length ) * box_length
        endif

        pair % RIJ_SQ = SUM( pair % RIJ**2 )
        pair % RIj_MAG = SQRT( pair % RIJ_SQ )
        pair % RIJ_hat = pair % RIJ / pair % RIj_MAG

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

        MM1 = transpose(U1) .x. S1
        MM2 = transpose(U2) .x. S2
        pair % AR = transpose(U1) .x. ( S1 ** 2.0 ) .x. U1
        pair % BR = transpose(U2) .x. ( S2 ** 2.0 ) .x. U2
        pair % AE = transpose(U1) .x. E1 .x. U1
        pair % BE = transpose(U2) .x. E2 .x. U2

        sigma_t1 = (shape_x(I) * shape_y(I) + shape_z(I) ** 2) * sqrt( 2.0 * shape_x(I) * shape_y(I) ) / 8.0
        sigma_t2 = (shape_x(J) * shape_y(J) + shape_z(J) ** 2) * sqrt( 2.0 * shape_x(J) * shape_y(J) ) / 8.0
        sigma_t = sqrt( sigma_t1 * sigma_t2 )
        
        EM = inverse( pair % AR + pair % BR )
        eps1 = sigma_t * sqrt( determinant( EM ) )
        
        lambda_E = brent( optim_eps, pair )
        if ( lambda_E < 0.0 .or. lambda_E > 1.0 ) then
            ierr_e = .TRUE.
        endif
        
        GE = ( 1 - lambda_E ) * pair % AE + lambda_E * pair % BE
        KE = inverse( GE ) .x. pair % RIJ
        eps2 = lambda_E * (1.0 - lambda_E) * dot( pair % RIJ_hat, KE ) / pair % RIJ_mag
        
        lambda_R = brent( optim_dist, pair )
        if ( lambda_R < 0.0 .or. lambda_R > 1.0 ) then
            ierr_r = .TRUE.
        endif
        
        GR = ( 1 - lambda_R ) * pair % AR + lambda_R * pair % BR
        KR = inverse( GR ) .x. pair % RIJ
        phi = lambda_R * (1.0 - lambda_R) * dot( pair % RIJ_hat, KR ) / pair % RIJ_mag
        sigma = phi ** (-1.0/2.0)
        efac = 2.0 * lambda_E * ( 1.0 - lambda_E ) / pair % RIJ_SQ
        rfac = 2.0 * lambda_R * ( 1.0 - lambda_R ) / pair % RIJ_SQ

        switch = ( cutoff + wf ) - 2.0 * sigma_0    ! Apply switch function 2sigma_0 inside cutoff
        factor = 1.0 / ( ( cutoff + wf ) - switch )
        xx = ( pair % RIJ_MAG - switch ) * factor
        sx = 1.0 - 6.0 * xx ** 5.0 + 15.0 * xx ** 4.0 - 10.0 * xx ** 3.0
        sdx = -30.0 * xx ** 4.0 + 60.0 * xx ** 3.0 - 30.0 * xx ** 2.0

        if ( pair % rij_mag < ( sigma + cutoff + wf ) ) then
        
            if (imp_wat) then
                if ( pair % rij_mag < sigma ) then
                    sr = sigma_0 / ( pair % rij_mag - sigma + sigma_0 )
                else if ((pair % rij_mag .ge. sigma) .and. (pair % rij_mag .le. (sigma + wf))) then
                    sr = sigma_0 / ( sigma - sigma + sigma_0 )
                else
                    sr = sigma_0 / ( pair % rij_mag - sigma + sigma_0 - wf )
                endif
            else
                sr = sigma_0 / ( pair % rij_mag - sigma + sigma_0 )
            endif

            ! EVALUATE NORMAL ECP POTENTIAL
            
            srx_fac = (sr ** 12.0 - 2.0 * sr ** 6.0)
            srd_fac = (12.0 * sr**7 - 12.0 * sr**13)

            dR_dr = (1.0 / sigma_0) * (  ( pair % RIJ_hat )     + 0.5 * rfac * ( KR - dot_product( pair % RIJ_hat, KR ) * pair % RIJ_hat ) * phi ** (-3.0 / 2.0)  )
            de2_dr = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * efac * ( KE - dot_product( pair % RIJ_hat, KE ) * pair % RIJ_hat )

            de1_du1x = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM1(:, 1), ( EM .x. MM1(:, 1) ) )
            de1_du1y = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM1(:, 2), ( EM .x. MM1(:, 2) ) )
            de1_du1z = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM1(:, 3), ( EM .x. MM1(:, 3) ) )

            de1_du1 = (de1_du1x + de1_du1y + de1_du1z)
            de2_du1 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * (1 - lambda_E) * efac * cross( KE .x. pair % AE, KE )
            dR_du1 = (eps1 ** neu) * (eps2 ** meu) * (0.5 / sigma_0) * (1 - lambda_R) * rfac * cross( KR .x. pair % AR, KR ) * phi ** (-3.0 / 2.0)

            de1_du2x = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM2(:, 1), ( EM .x. MM2(:, 1) ) )
            de1_du2y = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM2(:, 2), ( EM .x. MM2(:, 2) ) )
            de1_du2z = (eps2 ** meu) * (neu * eps1 ** neu) * cross( MM2(:, 3), ( EM .x. MM2(:, 3) ) )

            de1_du2 = (de1_du2x + de1_du2y + de1_du2z)
            de2_du2 = (eps1 ** neu) * (meu * eps2 ** (meu - 1)) * lambda_E * efac * cross( KE .x. pair % BE, KE )
            dR_du2 = (eps1 ** neu) * (eps2 ** meu) * (0.5 / sigma_0) * lambda_R * rfac * cross( KR .x. pair % BR, KR ) * phi ** (-3.0 / 2.0)

            pot   = (eps1 ** neu) * (eps2 ** meu) * srx_fac
            FIJ   = (-1.0) * ((eps1 ** neu) * (eps2 ** meu) * srd_fac * dR_dr  + srx_fac * de2_dr)
            
            TORQ1 = (srx_fac * (de1_du1 + de2_du1) + srd_fac * dR_du1)
            TORQ2 = (srx_fac * (de1_du2 + de2_du2) + srd_fac * dR_du2)

            ! EVALUATE CHIRAL TERMS

            ! srx_fac = (sr ** 7.0)
            ! srd_fac = (-7.0 * sr ** 8.0)
            u1xu2   = cross( U1(3, :), U2(3, :) )
            uur     = dot_product( u1xu2, pair % RIJ_hat )
            uu      = dot( U1(3, :), U2(3, :) )

            dR_dr   = (eps1 ** neu) * (eps2 ** meu) * dR_dr  * srd_fac * uur * uu
            de2_dr  = de2_dr * srx_fac * uur * uu
            duur_dr = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * ( U1XU2 - dot( U1XU2, pair % RIJ_hat ) * pair % RIJ_hat ) / pair % rij_mag

            de1_du1   = de1_du1 * srx_fac * uur * uu
            de2_du1   = de2_du1 * srx_fac * uur * uu
            dR_du1    = dR_du1  * srd_fac * uur * uu
            duur_du1z = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * cross( U2(3, :), pair % RIJ_hat )
            duu_du1z  = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * U2(3, :)

            de1_du2   = de1_du2 * srx_fac * uur * uu
            de2_du2   = de2_du2 * srx_fac * uur * uu
            dR_du2    = dR_du2  * srd_fac * uur * uu
            duur_du2z = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uu * cross( pair % RIJ_hat, U1(3, :) )
            duu_du2z  = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * U1(3, :)

            duur_du1  = (cross( U1(3, :), duur_du1z ))
            duu_du1   = (cross( U1(3, :), duu_du1z  ))

            duur_du2  = (cross( U2(3, :), duur_du2z ))
            duu_du2   = (cross( U2(3, :), duu_du2z  ))

            pot_chiral   = (eps1 ** neu) * (eps2 ** meu) * srx_fac * uur * uu
            FIJ_chiral   = (-1.0) * (de2_dr  + dR_dr  + duur_dr)
            TORQ1_chiral = (de1_du1 + de2_du1 + dR_du1 - duur_du1 - duu_du1)
            TORQ2_chiral = (de1_du2 + de2_du2 + dR_du2 - duur_du2 - duu_du2)

            pot   = pot   + ch_str * pot_chiral
            FIJ   = FIJ   + ch_str * FIJ_chiral
            TORQ1 = TORQ1 + ch_str * TORQ1_chiral
            TORQ2 = TORQ2 + ch_str * TORQ2_chiral

            if ((pair % rij_mag .ge. switch) .and. (pair % rij_mag .le. (cutoff + wf))) then
                pot   = sx * pot
                FIJ   = sx * FIJ   - sdx * pot * factor
                TORQ1 = sx * TORQ1 - sdx * pot * factor
                TORQ2 = sx * TORQ2 - sdx * pot * factor
            endif

        else
            pot = 0.0
            FIJ = 0.0
            TORQ1 = 0.0
            TORQ2 = 0.0
        endif
        
        PE = pot
        FXIJ = FIJ
        TXI = TORQ1
        TXJ = TORQ2
        
    end subroutine ecpc_calculate_forces
end module ecp_chiral

module bonded_interactions
    use system
    implicit none

    contains
    subroutine bend_calculate_forces( I, J, K, coeff )
        implicit none
        integer, intent(in)         :: I, J, K
        real(real64), intent(in)    :: coeff
        real(real64)                :: VI(3), VJ(3), FI(3), FJ(3), FK(3)
        real(real64)                :: CC11, CC12, CC22
        real(real64)                :: prefac, fac, fac1, fac2, pot

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
        integer, intent(in)         :: I, J, K, L
        real(real64), intent(in)    :: coeff
        real(real64)                :: VI(3), VJ(3), VK(3), FI(3), FJ(3), FK(3), FL(3)
        real(real64)                :: CC11, CC12, CC13, CC22, CC23, CC33
        real(real64)                :: DD11, DD12, DD13, DD22, DD23, DD33
        real(real64)                :: prefac, fac, fac1, fac2, fac3
        real(real64)                :: pot

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
        integer, intent(in)         :: I, J, K
        real(real64), intent(in)    :: coeff, ang
        real(real64)                :: VI(3), VJ(3), FI(3), FJ(3), FK(3)
        real(real64)                :: CC11, CC12, CC22
        real(real64)                :: prefac, fac, fac1, fac2, pot

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
        real(real64), intent(in)    :: DT
        real(real64)                :: RIJ(3), RIJ_old(3), DR(3)   
        real(real64)                :: L_ij, chi, vdot, r_tol, len
        real(real64)                :: inv_mass_a, inv_mass_b
        integer                     :: I, J, IX, rattle, max_rattle = 100
        logical                     :: moved(N)

        L_ij = 0.0; rattle = 0; 
        r_tol = 1.0e-5; max_rattle = 100
        moved = .TRUE.

        do while (any(moved) .and. rattle < max_rattle)

            do IX = 1, NC
                I   = BBI(IX)
                J   = BBJ(IX)
                len = BB_len(IX)

                inv_mass_a = 1.0 / mass(I)
                inv_mass_b = 1.0 / mass(J)
                
                if ( moved(I) .or. moved(J) ) then

                    RIJ(1) = RX(I) - RX(J)
                    RIJ(2) = RY(I) - RY(J)
                    RIJ(3) = RZ(I) - RZ(J)

                    if (is_periodic) then
                        RIJ = RIJ - ANINT( RIJ / box_length ) * box_length
                    endif

                    chi = len ** 2 - SUM( RIJ ** 2 )

                    if (abs(chi) > 2 * r_tol * len**2.0) then 

                        RIJ_old(1) = RX_old(I) - RX_old(J)
                        RIJ_old(2) = RY_old(I) - RY_old(J)
                        RIJ_old(3) = RZ_old(I) - RZ_old(J)

                        if (is_periodic) then
                            RIJ_old = RIJ_old - ANINT( RIJ_old / box_length ) * box_length
                        endif

                        vdot = dot( RIJ, RIJ_old )

                        if (vdot < r_tol * len ** 2.0) then
                            write(*, "(5F10.5)") chi, vdot, len, sum(RIJ**2), sum(RIJ_old**2)
                            write(*, "('Constraint failure between atoms ', I0, ' and ', I0)") I, J
                            stop
                        endif

                        L_ij = chi / ( 2.0 * vdot * ( inv_mass_a + inv_mass_b ) )

                        virial = virial + ( L_ij * len ** 2 ) / DT ** 2

                        RX(I) = RX(I) + ( L_ij * inv_mass_a ) * RIJ_old(1)
                        RY(I) = RY(I) + ( L_ij * inv_mass_a ) * RIJ_old(2)
                        RZ(I) = RZ(I) + ( L_ij * inv_mass_a ) * RIJ_old(3)

                        RX(J) = RX(J) - ( L_ij * inv_mass_b ) * RIJ_old(1)
                        RY(J) = RY(J) - ( L_ij * inv_mass_b ) * RIJ_old(2)
                        RZ(J) = RZ(J) - ( L_ij * inv_mass_b ) * RIJ_old(3)

                        VX(I) = VX(I) + ( L_ij * inv_mass_a ) * RIJ_old(1) / DT
                        VY(I) = VY(I) + ( L_ij * inv_mass_a ) * RIJ_old(2) / DT
                        VZ(I) = VZ(I) + ( L_ij * inv_mass_a ) * RIJ_old(3) / DT

                        VX(J) = VX(J) - ( L_ij * inv_mass_b ) * RIJ_old(1) / DT
                        VY(J) = VY(J) - ( L_ij * inv_mass_b ) * RIJ_old(2) / DT
                        VZ(J) = VZ(J) - ( L_ij * inv_mass_b ) * RIJ_old(3) / DT

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
        real(real64), intent(in)    :: DT
        real(real64)                :: RIJ(3), VIJ(3)
        real(real64)                :: L_ij, chi, vdot, r_tol, len
        real(real64)                :: inv_mass_a, inv_mass_b
        integer                     :: I, J, IX, rattle, max_rattle = 100
        logical                     :: moved(N)

        L_ij = 0.0; rattle = 0; 
        r_tol = 1.0e-5; max_rattle = 100
        moved = .TRUE.

        do while (any(moved) .and. rattle < max_rattle)

            do IX = 1, NC
                I   = BBI(IX)
                J   = BBJ(IX)
                len = BB_len(IX)

                inv_mass_a = 1.0 / mass(I)
                inv_mass_b = 1.0 / mass(J)

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
                    L_ij = -vdot / ( ( inv_mass_a + inv_mass_b ) * len**2 )
                    
                    if (abs(L_ij) > r_tol) then

                        virial = virial + ( L_ij * len**2 ) / DT

                        VX(I) = VX(I) + ( L_ij * inv_mass_a ) * RIJ(1)
                        VY(I) = VY(I) + ( L_ij * inv_mass_a ) * RIJ(2)
                        VZ(I) = VZ(I) + ( L_ij * inv_mass_a ) * RIJ(3)

                        VX(J) = VX(J) - ( L_ij * inv_mass_b ) * RIJ(1)
                        VY(J) = VY(J) - ( L_ij * inv_mass_b ) * RIJ(2)
                        VZ(J) = VZ(J) - ( L_ij * inv_mass_b ) * RIJ(3)

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