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

        if (is_periodic) then
            RX = RX - ANINT( RX / box_length(1) ) * box_length(1)
            RY = RY - ANINT( RY / box_length(2) ) * box_length(2)
            RZ = RZ - ANINT( RZ / box_length(3) ) * box_length(3)
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

module quaternion_integrator
    use xmath
    use system
    implicit none

    real(8), dimension(3, 3), private :: RM
    real(8), private :: QDW, QDX, QDY, QDZ
    real(8), private :: LX_body, LY_body, LZ_body
    real(8), private :: WX_body, WY_body, WZ_body

    contains
    subroutine qq_initial_step( DT )
        implicit none
        real(8), intent(in) :: DT
        real(8) :: QW_Old, QX_Old, QY_Old, QZ_Old
        real(8) :: QW_New, QX_New, QY_New, QZ_New
        real(8) :: Q_mag
        integer :: I, J

        do I = 1, N
            QW_New = QW(I); QX_New = QX(I)
            QY_New = QY(I); QZ_New = QZ(I)

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            do J = 1, 2
                QW_Old = QW_New; QX_Old = QX_New
                QY_Old = QY_New; QZ_Old = QZ_New
                
                RM(1, 1) = QW_old**2 + QX_old**2 - QY_old**2 - QZ_old**2
                RM(1, 2) = 2.0 * ( QX_old * QY_old + QW_old * QZ_old )
                RM(1, 3) = 2.0 * ( QX_old * QZ_old - QW_old * QY_old )
                RM(2, 1) = 2.0 * ( QX_old * QY_old - QW_old * QZ_old )
                RM(2, 2) = QW_old**2 - QX_old**2 + QY_old**2 - QZ_old**2
                RM(2, 3) = 2.0 * ( QY_old * QZ_old + QW_old * QX_old )
                RM(3, 1) = 2.0 * ( QX_old * QZ_old + QW_old * QY_old )
                RM(3, 2) = 2.0 * ( QY_old * QZ_old - QW_old * QX_old )
                RM(3, 3) = QW_old**2 - QX_old**2 - QY_old**2 + QZ_old**2

                WX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) ) / MOI(1)
                WY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) ) / MOI(2)
                WZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) ) / MOI(3)

                QDW = 0.5 * ( - QX_Old * WX_body - QY_Old * WY_body - QZ_Old * WZ_body )
                QDX = 0.5 * (   QW_Old * WX_body + QY_Old * WZ_body - QZ_Old * WY_body )
                QDY = 0.5 * (   QW_Old * WY_body - QX_Old * WZ_body + QZ_Old * WX_body )
                QDZ = 0.5 * (   QW_Old * WZ_body + QX_Old * WY_body - QY_Old * WX_body )

                QW_New = QW(I) + 0.5 * DT * QDW
                QX_New = QX(I) + 0.5 * DT * QDX
                QY_New = QY(I) + 0.5 * DT * QDY
                QZ_New = QZ(I) + 0.5 * DT * QDZ
            
            enddo

            QW(I) = QW(I) + DT * QDW
            QX(I) = QX(I) + DT * QDX
            QY(I) = QY(I) + DT * QDY
            QZ(I) = QZ(I) + DT * QDZ

            Q_mag = SQRT( QW(I)**2 + QX(I)**2 + QY(I)**2 + QZ(I)**2 )
            QW(I) = QW(I) / Q_mag
            QX(I) = QX(I) / Q_mag
            QY(I) = QY(I) / Q_mag
            QZ(I) = QZ(I) / Q_mag

        enddo

    end subroutine qq_initial_step

    subroutine qq_final_step( DT )
        implicit none
        real(8), intent(in) :: DT
        integer :: I

        do I = 1, N

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            RM(1, 1) = QW(I)**2 + QX(I)**2 - QY(I)**2 - QZ(I)**2
            RM(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
            RM(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
            RM(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
            RM(2, 2) = QW(I)**2 - QX(I)**2 + QY(I)**2 - QZ(I)**2
            RM(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
            RM(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
            RM(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
            RM(3, 3) = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

            LX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
            LY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
            LZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

            kinetic_energy = kinetic_energy + 0.5 * ( ( LX_body ** 2.0) / MOI(1) + &
                                                      ( LY_body ** 2.0) / MOI(2) + &
                                                      ( LZ_body ** 2.0) / MOI(3) )

        enddo

    end subroutine qq_final_step

end module quaternion_integrator

module nose_hoover
    use system
    implicit none

    integer, parameter       :: N_th = 2
    real(8)                  :: tau = 5.0
    real(8), dimension(N_th) :: HB = 1.0, eta = 0.0, p_eta = 0.0

    contains
    ! subroutine nht_initialize( chain_length )
    !     implicit none

    !     integer, optional, intent(in) :: chain_length

    !     if (present(chain_length)) then
    !         N_th = chain_length
    !     else
    !         N_th = 2
    !     endif

    !     allocate( HB(N_th), eta(N_th), p_eta(N_th) )
    !     HB = 0.0; eta = 0.0; p_eta = 0.0

    !     if (not(target_temp_set)) then
    !         write(*,*) "ERROR: Target temperature has not been set"
    !         stop
    !     endif

    !     HB = target_temp * tau ** 2.0
    !     HB(1) = dof * target_temp * tau ** 2.0

    ! end subroutine nht_initialize

    SUBROUTINE u4_propagator ( t, j_start, j_stop, j_stride )
        implicit none
        real(8), intent(in) :: t               
        integer, intent(in) :: j_start, j_stop

        integer :: j, j_stride
        REAL(8) :: gj, x, c

        do j = j_start, j_stop, j_stride

            if ( j == 1 ) then
                gj = SUM(VX**2 + VY**2 + VZ**2) - (dof-3*N) * target_temp
            else
                gj = ( p_eta(j-1)**2 / HB(j-1) ) - target_temp
            endif

            if ( j == N_th ) then
                p_eta(j)  = p_eta(j) + t * gj
            else
                x = t * p_eta(j+1)/HB(j+1)
                c = exprel(-x) ! (1-exp(-x))/x, preserving accuracy for small x

                p_eta(j) = p_eta(j)*EXP(-x) + t * gj * c
            endif

        enddo

    end subroutine u4_propagator

    subroutine nht_initial_step( DT )
        implicit none
        real(8), intent(in) :: DT

        call u4_propagator( DT / 4.0, N_th, 1, -1)

        VX = VX * EXP( -0.5 * DT * p_eta(1) / HB(1) )
        VY = VY * EXP( -0.5 * DT * p_eta(1) / HB(1) )
        VZ = VZ * EXP( -0.5 * DT * p_eta(1) / HB(1) )

        eta = eta + 0.5 * DT * p_eta / HB

        CALL u4_propagator( DT / 4.0, 1, N_th, 1)

        VX = VX + 0.5 * DT * FX
        VY = VY + 0.5 * DT * FY
        VZ = VZ + 0.5 * DT * FZ

        RX = RX + DT * VX
        RY = RY + DT * VY
        RZ = RZ + DT * VZ

        if (is_periodic) then
            RX = RX - ANINT( RX / box_length(1) ) * box_length(1)
            RY = RY - ANINT( RY / box_length(2) ) * box_length(2)
            RZ = RZ - ANINT( RZ / box_length(3) ) * box_length(3)
        endif

    end subroutine nht_initial_step

    subroutine nht_final_step( DT )
        implicit none
        real(8), intent(in) :: DT

        VX = VX + 0.5 * DT * FX
        VY = VY + 0.5 * DT * FY
        VZ = VZ + 0.5 * DT * FZ

        call u4_propagator( DT / 4.0, N_th, 1, -1)

        VX = VX * EXP( -0.5 * DT * p_eta(1) / HB(1) )
        VY = VY * EXP( -0.5 * DT * p_eta(1) / HB(1) )
        VZ = VZ * EXP( -0.5 * DT * p_eta(1) / HB(1) )

        eta = eta + 0.5 * DT * p_eta / HB

        CALL u4_propagator( DT / 4.0, 1, N_th, 1)

    end subroutine nht_final_step

end module nose_hoover

module langevin_integrator
    use system
    implicit none
    real(8), dimension(3, 3), private   :: RM
    real(8), dimension(3)               :: noise1, noise2
    real(8), private                    :: QDW, QDX, QDY, QDZ
    real(8), private                    :: LX_body, LY_body, LZ_body
    real(8), private                    :: WX_body, WY_body, WZ_body
    real(8)                             :: kk1, kk2, damp = 1.0
    
    contains
    subroutine lgv_initial_step( DT )
        implicit none
        real(8), intent(in)  :: DT
        real(8)              :: QW_Old, QX_Old, QY_Old, QZ_Old
        real(8)              :: QW_New, QX_New, QY_New, QZ_New
        real(8)              :: Q_mag
        integer              :: I, J

        kk1 = 1.0 - (damp * DT) / 2.0
        kk2 = 1.0 / (1.0 + (damp * DT) / 2.0)

        do I = 1, N

            VX(I) = VX(I) + 0.5 * DT * FX(I)
            VY(I) = VY(I) + 0.5 * DT * FY(I)
            VZ(I) = VZ(I) + 0.5 * DT * FZ(I)

            RX(I) = RX(I) + 0.5 * DT * VX(I)
            RY(I) = RY(I) + 0.5 * DT * VY(I)
            RZ(I) = RZ(I) + 0.5 * DT * VZ(I)

            call normal_sequences( 3, noise1, noise2 )

            VX(I) = EXP( -damp * DT ) * VX(I) + sqrt( target_temp * ( 1 - EXP( -2.0 * damp * DT ))) * noise1(1)
            VY(I) = EXP( -damp * DT ) * VY(I) + sqrt( target_temp * ( 1 - EXP( -2.0 * damp * DT ))) * noise1(2)
            VZ(I) = EXP( -damp * DT ) * VZ(I) + sqrt( target_temp * ( 1 - EXP( -2.0 * damp * DT ))) * noise1(3)

            RX(I) = RX(I) + 0.5 * DT * VX(I)
            RY(I) = RY(I) + 0.5 * DT * VY(I)
            RZ(I) = RZ(I) + 0.5 * DT * VZ(I)

            if (is_periodic) then
                RX = RX - ANINT( RX / box_length(1) ) * box_length(1)
                RY = RY - ANINT( RY / box_length(2) ) * box_length(2)
                RZ = RZ - ANINT( RZ / box_length(3) ) * box_length(3)
            endif

            QW_New = QW(I); QX_New = QX(I)
            QY_New = QY(I); QZ_New = QZ(I)

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            do J = 1, 2
                QW_Old = QW_New; QX_Old = QX_New
                QY_Old = QY_New; QZ_Old = QZ_New
                
                RM(1, 1) = QW_old**2 + QX_old**2 - QY_old**2 - QZ_old**2
                RM(1, 2) = 2.0 * ( QX_old * QY_old + QW_old * QZ_old )
                RM(1, 3) = 2.0 * ( QX_old * QZ_old - QW_old * QY_old )
                RM(2, 1) = 2.0 * ( QX_old * QY_old - QW_old * QZ_old )
                RM(2, 2) = QW_old**2 - QX_old**2 + QY_old**2 - QZ_old**2
                RM(2, 3) = 2.0 * ( QY_old * QZ_old + QW_old * QX_old )
                RM(3, 1) = 2.0 * ( QX_old * QZ_old + QW_old * QY_old )
                RM(3, 2) = 2.0 * ( QY_old * QZ_old - QW_old * QX_old )
                RM(3, 3) = QW_old**2 - QX_old**2 - QY_old**2 + QZ_old**2

                LX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
                LY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
                LZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

                WX_body = LX_body / MOI(1)
                WY_body = LY_body / MOI(2)
                WZ_body = LZ_body / MOI(3)

                QDW = 0.5 * ( - QX_Old * WX_body - QY_Old * WY_body - QZ_Old * WZ_body )
                QDX = 0.5 * (   QW_Old * WX_body + QY_Old * WZ_body - QZ_Old * WY_body )
                QDY = 0.5 * (   QW_Old * WY_body - QX_Old * WZ_body + QZ_Old * WX_body )
                QDZ = 0.5 * (   QW_Old * WZ_body + QX_Old * WY_body - QY_Old * WX_body )

                QW_New = QW(I) + 0.5 * DT * QDW
                QX_New = QX(I) + 0.5 * DT * QDX
                QY_New = QY(I) + 0.5 * DT * QDY
                QZ_New = QZ(I) + 0.5 * DT * QDZ
            
            enddo

            LX_body = kk2 * kk1 * LX_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * MOI(1) ) * noise2(1)
            LY_body = kk2 * kk1 * LY_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * MOI(2) ) * noise2(2)
            LZ_body = kk2 * kk1 * LZ_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * MOI(3) ) * noise2(3)

            LX(I) = ( RM(1,1) * LX_body + RM(2,1) * LY_body + RM(3,1) * LZ_body )
            LY(I) = ( RM(1,2) * LX_body + RM(2,2) * LY_body + RM(3,2) * LZ_body )
            LZ(I) = ( RM(1,3) * LX_body + RM(2,3) * LY_body + RM(3,3) * LZ_body )

            QW(I) = QW(I) + DT * QDW
            QX(I) = QX(I) + DT * QDX
            QY(I) = QY(I) + DT * QDY
            QZ(I) = QZ(I) + DT * QDZ

            Q_mag = SQRT( QW(I)**2 + QX(I)**2 + QY(I)**2 + QZ(I)**2 )
            QW(I) = QW(I) / Q_mag
            QX(I) = QX(I) / Q_mag
            QY(I) = QY(I) / Q_mag
            QZ(I) = QZ(I) / Q_mag

        enddo

    end subroutine lgv_initial_step

    subroutine lgv_final_step( DT )
        implicit none
        real(8), intent(in) :: DT
        integer :: I

        do I = 1, N

            VX(I) = VX(I) + 0.5 * DT * FX(I)
            VY(I) = VY(I) + 0.5 * DT * FY(I)
            VZ(I) = VZ(I) + 0.5 * DT * FZ(I)

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            RM(1, 1) = QW(I)**2 + QX(I)**2 - QY(I)**2 - QZ(I)**2
            RM(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
            RM(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
            RM(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
            RM(2, 2) = QW(I)**2 - QX(I)**2 + QY(I)**2 - QZ(I)**2
            RM(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
            RM(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
            RM(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
            RM(3, 3) = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

            LX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
            LY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
            LZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

            kinetic_energy = kinetic_energy + 0.5 * ( ( LX_body ** 2.0) / MOI(1) + &
                                                      ( LY_body ** 2.0) / MOI(2) + &
                                                      ( LZ_body ** 2.0) / MOI(3) )

        enddo

    end subroutine lgv_final_step

end module langevin_integrator

module qlangevin_integrator
    use system
    implicit none

    real(8), dimension(3, 3), private   :: RM
    real(8), private                    :: QDW, QDX, QDY, QDZ
    real(8), private                    :: LX_body, LY_body, LZ_body
    real(8), private                    :: WX_body, WY_body, WZ_body
    real(8)                             :: kk1, kk2, damp = 1.0
    real(8), dimension(3)               :: noise1, noise2

    contains
    subroutine qlgv_initial_step( DT )
        implicit none
        real(8), intent(in) :: DT
        real(8) :: QW_Old, QX_Old, QY_Old, QZ_Old
        real(8) :: QW_New, QX_New, QY_New, QZ_New
        real(8) :: Q_mag
        integer :: I, J

        call normal_sequences( 3, noise1, noise2 )

        kk1 = 1.0 - (damp * DT) / 2.0
        kk2 = 1.0 / (1.0 + (damp * DT) / 2.0)

        do I = 1, N
            QW_New = QW(I); QX_New = QX(I)
            QY_New = QY(I); QZ_New = QZ(I)

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            do J = 1, 2
                QW_Old = QW_New; QX_Old = QX_New
                QY_Old = QY_New; QZ_Old = QZ_New
                
                RM(1, 1) = QW_old**2 + QX_old**2 - QY_old**2 - QZ_old**2
                RM(1, 2) = 2.0 * ( QX_old * QY_old + QW_old * QZ_old )
                RM(1, 3) = 2.0 * ( QX_old * QZ_old - QW_old * QY_old )
                RM(2, 1) = 2.0 * ( QX_old * QY_old - QW_old * QZ_old )
                RM(2, 2) = QW_old**2 - QX_old**2 + QY_old**2 - QZ_old**2
                RM(2, 3) = 2.0 * ( QY_old * QZ_old + QW_old * QX_old )
                RM(3, 1) = 2.0 * ( QX_old * QZ_old + QW_old * QY_old )
                RM(3, 2) = 2.0 * ( QY_old * QZ_old - QW_old * QX_old )
                RM(3, 3) = QW_old**2 - QX_old**2 - QY_old**2 + QZ_old**2

                LX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
                LY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
                LZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

                WX_body = LX_body / MOI(1)
                WY_body = LY_body / MOI(2)
                WZ_body = LZ_body / MOI(3)

                QDW = 0.5 * ( - QX_Old * WX_body - QY_Old * WY_body - QZ_Old * WZ_body )
                QDX = 0.5 * (   QW_Old * WX_body + QY_Old * WZ_body - QZ_Old * WY_body )
                QDY = 0.5 * (   QW_Old * WY_body - QX_Old * WZ_body + QZ_Old * WX_body )
                QDZ = 0.5 * (   QW_Old * WZ_body + QX_Old * WY_body - QY_Old * WX_body )

                QW_New = QW(I) + 0.5 * DT * QDW
                QX_New = QX(I) + 0.5 * DT * QDX
                QY_New = QY(I) + 0.5 * DT * QDY
                QZ_New = QZ(I) + 0.5 * DT * QDZ
            
            enddo

            LX_body = kk2 * kk1 * LX_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * MOI(1) ) * noise2(1)
            LY_body = kk2 * kk1 * LY_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * MOI(2) ) * noise2(2)
            LZ_body = kk2 * kk1 * LZ_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * MOI(3) ) * noise2(3)

            LX(I) = ( RM(1,1) * LX_body + RM(2,1) * LY_body + RM(3,1) * LZ_body )
            LY(I) = ( RM(1,2) * LX_body + RM(2,2) * LY_body + RM(3,2) * LZ_body )
            LZ(I) = ( RM(1,3) * LX_body + RM(2,3) * LY_body + RM(3,3) * LZ_body )

            QW(I) = QW(I) + DT * QDW
            QX(I) = QX(I) + DT * QDX
            QY(I) = QY(I) + DT * QDY
            QZ(I) = QZ(I) + DT * QDZ

            Q_mag = SQRT( QW(I)**2 + QX(I)**2 + QY(I)**2 + QZ(I)**2 )
            QW(I) = QW(I) / Q_mag
            QX(I) = QX(I) / Q_mag
            QY(I) = QY(I) / Q_mag
            QZ(I) = QZ(I) / Q_mag

        enddo

    end subroutine qlgv_initial_step

    subroutine qlgv_final_step( DT )
        implicit none
        real(8), intent(in) :: DT
        integer :: I

        do I = 1, N

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            RM(1, 1) = QW(I)**2 + QX(I)**2 - QY(I)**2 - QZ(I)**2
            RM(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
            RM(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
            RM(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
            RM(2, 2) = QW(I)**2 - QX(I)**2 + QY(I)**2 - QZ(I)**2
            RM(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
            RM(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
            RM(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
            RM(3, 3) = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

            LX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
            LY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
            LZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

            kinetic_energy = kinetic_energy + 0.5 * ( ( LX_body ** 2.0) / MOI(1) + &
                                                      ( LY_body ** 2.0) / MOI(2) + &
                                                      ( LZ_body ** 2.0) / MOI(3) )

        enddo

    end subroutine qlgv_final_step

end module qlangevin_integrator

module verlet_constraint_integrator
    use system
    use constraints
    implicit none

    contains
    subroutine vvc_initial_step( DT )
        implicit none
        real(8), intent(in) :: DT

        VX = VX + 0.5 * DT * FX
        VY = VY + 0.5 * DT * FY
        VZ = VZ + 0.5 * DT * FZ

        RX_old = RX
        RY_old = RY
        RZ_old = RZ

        RX = RX + DT * VX
        RY = RY + DT * VY
        RZ = RZ + DT * VZ

        if (is_periodic) then
            RX = RX - ANINT( RX / box_length(1) ) * box_length(1)
            RY = RY - ANINT( RY / box_length(2) ) * box_length(2)
            RZ = RZ - ANINT( RZ / box_length(3) ) * box_length(3)
        endif

        call apply_constraints_a( DT )

    end subroutine vvc_initial_step

    subroutine vvc_final_step( DT )
        implicit none
        real(8), intent(in) :: DT

        VX = VX + 0.5 * DT * FX
        VY = VY + 0.5 * DT * FY
        VZ = VZ + 0.5 * DT * FZ

        call apply_constraints_b( DT )

    end subroutine vvc_final_step

end module verlet_constraint_integrator

module langevin_constraint_integrator
    use system
    use constraints
    implicit none
    real(8), dimension(3, 3), private   :: RM
    real(8), dimension(3)               :: noise1, noise2
    real(8), private                    :: QDW, QDX, QDY, QDZ
    real(8), private                    :: LX_body, LY_body, LZ_body
    real(8), private                    :: WX_body, WY_body, WZ_body
    real(8)                             :: kk1, kk2, damp = 1.0
    
    contains
    subroutine lgvc_initial_step( DT )
        implicit none
        real(8), intent(in)  :: DT
        real(8)              :: QW_Old, QX_Old, QY_Old, QZ_Old
        real(8)              :: QW_New, QX_New, QY_New, QZ_New
        real(8)              :: Q_mag
        integer              :: I, J

        kk1 = 1.0 - (damp * DT) / 2.0
        kk2 = 1.0 / (1.0 + (damp * DT) / 2.0)

        do I = 1, N

            VX(I) = VX(I) + 0.5 * DT * FX(I)
            VY(I) = VY(I) + 0.5 * DT * FY(I)
            VZ(I) = VZ(I) + 0.5 * DT * FZ(I)

            RX_old(I) = RX(I)
            RY_old(I) = RY(I)
            RZ_old(I) = RZ(I)

            RX(I) = RX(I) + 0.5 * DT * VX(I)
            RY(I) = RY(I) + 0.5 * DT * VY(I)
            RZ(I) = RZ(I) + 0.5 * DT * VZ(I)

            call normal_sequences( 3, noise1, noise2 )

            VX(I) = EXP( -damp * DT ) * VX(I) + sqrt( target_temp * ( 1 - EXP( -2.0 * damp * DT ))) * noise1(1)
            VY(I) = EXP( -damp * DT ) * VY(I) + sqrt( target_temp * ( 1 - EXP( -2.0 * damp * DT ))) * noise1(2)
            VZ(I) = EXP( -damp * DT ) * VZ(I) + sqrt( target_temp * ( 1 - EXP( -2.0 * damp * DT ))) * noise1(3)

            RX(I) = RX(I) + 0.5 * DT * VX(I)
            RY(I) = RY(I) + 0.5 * DT * VY(I)
            RZ(I) = RZ(I) + 0.5 * DT * VZ(I)

            if (is_periodic) then
                RX = RX - ANINT( RX / box_length(1) ) * box_length(1)
                RY = RY - ANINT( RY / box_length(2) ) * box_length(2)
                RZ = RZ - ANINT( RZ / box_length(3) ) * box_length(3)
            endif

            QW_New = QW(I); QX_New = QX(I)
            QY_New = QY(I); QZ_New = QZ(I)

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            do J = 1, 2
                QW_Old = QW_New; QX_Old = QX_New
                QY_Old = QY_New; QZ_Old = QZ_New
                
                RM(1, 1) = QW_old**2 + QX_old**2 - QY_old**2 - QZ_old**2
                RM(1, 2) = 2.0 * ( QX_old * QY_old + QW_old * QZ_old )
                RM(1, 3) = 2.0 * ( QX_old * QZ_old - QW_old * QY_old )
                RM(2, 1) = 2.0 * ( QX_old * QY_old - QW_old * QZ_old )
                RM(2, 2) = QW_old**2 - QX_old**2 + QY_old**2 - QZ_old**2
                RM(2, 3) = 2.0 * ( QY_old * QZ_old + QW_old * QX_old )
                RM(3, 1) = 2.0 * ( QX_old * QZ_old + QW_old * QY_old )
                RM(3, 2) = 2.0 * ( QY_old * QZ_old - QW_old * QX_old )
                RM(3, 3) = QW_old**2 - QX_old**2 - QY_old**2 + QZ_old**2

                LX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
                LY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
                LZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

                WX_body = LX_body / MOI(1)
                WY_body = LY_body / MOI(2)
                WZ_body = LZ_body / MOI(3)

                QDW = 0.5 * ( - QX_Old * WX_body - QY_Old * WY_body - QZ_Old * WZ_body )
                QDX = 0.5 * (   QW_Old * WX_body + QY_Old * WZ_body - QZ_Old * WY_body )
                QDY = 0.5 * (   QW_Old * WY_body - QX_Old * WZ_body + QZ_Old * WX_body )
                QDZ = 0.5 * (   QW_Old * WZ_body + QX_Old * WY_body - QY_Old * WX_body )

                QW_New = QW(I) + 0.5 * DT * QDW
                QX_New = QX(I) + 0.5 * DT * QDX
                QY_New = QY(I) + 0.5 * DT * QDY
                QZ_New = QZ(I) + 0.5 * DT * QDZ
            
            enddo

            LX_body = kk2 * kk1 * LX_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * MOI(1) ) * noise2(1)
            LY_body = kk2 * kk1 * LY_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * MOI(2) ) * noise2(2)
            LZ_body = kk2 * kk1 * LZ_body + kk2 * sqrt( 2.0 * damp * target_temp * DT * MOI(3) ) * noise2(3)

            LX(I) = ( RM(1,1) * LX_body + RM(2,1) * LY_body + RM(3,1) * LZ_body )
            LY(I) = ( RM(1,2) * LX_body + RM(2,2) * LY_body + RM(3,2) * LZ_body )
            LZ(I) = ( RM(1,3) * LX_body + RM(2,3) * LY_body + RM(3,3) * LZ_body )

            QW(I) = QW(I) + DT * QDW
            QX(I) = QX(I) + DT * QDX
            QY(I) = QY(I) + DT * QDY
            QZ(I) = QZ(I) + DT * QDZ

            Q_mag = SQRT( QW(I)**2 + QX(I)**2 + QY(I)**2 + QZ(I)**2 )
            QW(I) = QW(I) / Q_mag
            QX(I) = QX(I) / Q_mag
            QY(I) = QY(I) / Q_mag
            QZ(I) = QZ(I) / Q_mag

        enddo

        call apply_constraints_a( DT )

    end subroutine lgvc_initial_step

    subroutine lgvc_final_step( DT )
        implicit none
        real(8), intent(in) :: DT
        integer :: I

        do I = 1, N

            VX(I) = VX(I) + 0.5 * DT * FX(I)
            VY(I) = VY(I) + 0.5 * DT * FY(I)
            VZ(I) = VZ(I) + 0.5 * DT * FZ(I)

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            RM(1, 1) = QW(I)**2 + QX(I)**2 - QY(I)**2 - QZ(I)**2
            RM(1, 2) = 2.0 * ( QX(I) * QY(I) + QW(I) * QZ(I) )
            RM(1, 3) = 2.0 * ( QX(I) * QZ(I) - QW(I) * QY(I) )
            RM(2, 1) = 2.0 * ( QX(I) * QY(I) - QW(I) * QZ(I) )
            RM(2, 2) = QW(I)**2 - QX(I)**2 + QY(I)**2 - QZ(I)**2
            RM(2, 3) = 2.0 * ( QY(I) * QZ(I) + QW(I) * QX(I) )
            RM(3, 1) = 2.0 * ( QX(I) * QZ(I) + QW(I) * QY(I) )
            RM(3, 2) = 2.0 * ( QY(I) * QZ(I) - QW(I) * QX(I) )
            RM(3, 3) = QW(I)**2 - QX(I)**2 - QY(I)**2 + QZ(I)**2

            LX_body = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
            LY_body = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
            LZ_body = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

            kinetic_energy = kinetic_energy + 0.5 * ( ( LX_body ** 2.0) / MOI(1) + &
                                                      ( LY_body ** 2.0) / MOI(2) + &
                                                      ( LZ_body ** 2.0) / MOI(3) )

        enddo

        call apply_constraints_b( DT )

    end subroutine lgvc_final_step

end module langevin_constraint_integrator
    