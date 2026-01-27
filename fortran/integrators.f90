
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

    contains
    subroutine qq_initial_step( DT )
        implicit none
        real(8), intent(in) :: DT
        real(8) :: QW_Old, QX_Old, QY_Old, QZ_Old, QW_New, QX_New, QY_New, QZ_New
        real(8) :: QDW, QDX, QDY, QDZ
        real(8), dimension(3) :: L_lab, L_body, W_body
        integer :: I, J

        do I = 1, N
            QW_New = QW(I); QX_New = QX(I)
            QY_New = QY(I); QZ_New = QZ(I)

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            L_lab = (/ LX(I), LY(I), LZ(I) /)

            do J = 1, 3
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

                WX(I) = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) ) / MOI(1)
                WY(I) = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) ) / MOI(2)
                WZ(I) = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) ) / MOI(3)

                QDW = 0.5 * ( - QX_Old * WX(I) - QY_Old * WY(I) - QZ_Old * WZ(I) )
                QDX = 0.5 * (   QW_Old * WX(I) + QY_Old * WZ(I) - QZ_Old * WY(I) )
                QDY = 0.5 * (   QW_Old * WY(I) - QX_Old * WZ(I) + QZ_Old * WX(I) )
                QDZ = 0.5 * (   QW_Old * WZ(I) + QX_Old * WY(I) - QY_Old * WX(I) )
                
                QW_New = QW(I) + 0.5 * DT * QDW
                QX_New = QX(I) + 0.5 * DT * QDX
                QY_New = QY(I) + 0.5 * DT * QDY
                QZ_New = QZ(I) + 0.5 * DT * QDZ
            
            enddo

            QW(I) = QW(I) + DT * QDW
            QX(I) = QX(I) + DT * QDX
            QY(I) = QY(I) + DT * QDY
            QZ(I) = QZ(I) + DT * QDZ

        enddo

    end subroutine qq_initial_step

    subroutine qq_final_step( DT )
        implicit none
        real(8), intent(in) :: DT
        real(8) :: LXB, LYB, LZB
        integer :: I

        do I = 1, N

            LX(I) = LX(I) + 0.5 * DT * TX(I)
            LY(I) = LY(I) + 0.5 * DT * TY(I)
            LZ(I) = LZ(I) + 0.5 * DT * TZ(I)

            LXB = ( RM(1,1) * LX(I) + RM(1,2) * LY(I) + RM(1,3) * LZ(I) )
            LYB = ( RM(2,1) * LX(I) + RM(2,2) * LY(I) + RM(2,3) * LZ(I) )
            LZB = ( RM(3,1) * LX(I) + RM(3,2) * LY(I) + RM(3,3) * LZ(I) )

            kinetic_energy = kinetic_energy + 0.5 * (LXB ** 2.0 / MOI(1) + LYB ** 2.0 / MOI(2) + LZB ** 2.0 / MOI(3) )

        enddo


    end subroutine qq_final_step

end module quaternion_integrator
    