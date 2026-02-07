
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
    