module lcg_random
    implicit none
    integer, parameter :: i64 = selected_int_kind(18)
    integer(8), parameter, private :: a = 1103515245
    integer(8), parameter, private :: b = 214013
    integer(8), parameter, private :: c = 2147483648_i64

    real, parameter, private :: A1 = 3.949846138
    real, parameter, private :: A3 = 0.252408784
    real, parameter, private :: A5 = 0.076542912
    real, parameter, private :: A7 = 0.008355968
    real, parameter, private :: A9 = 0.029899776

    integer, save :: seed

    contains
    subroutine init_rand()
        implicit none

        real :: sead
        call random_seed()
        call random_number(sead)
        seed = int(sead * 1E+5)

    end subroutine init_rand

    subroutine rand_lcg(p)
        implicit none
        real, intent(out):: p

        p = mod(a * seed + b, c)
        seed = p
        p = p / c
    end subroutine rand_lcg

    subroutine rand_array(arr, size)
        implicit none
        integer, intent(in) :: size
        real, dimension(size), intent(out) :: arr
        integer :: loop

        do loop = 1, size
            call rand_lcg(arr(loop))
        enddo
    end subroutine rand_array

    subroutine randrange_lcg(out, lower, upper)
        implicit none
        real, intent(in) :: lower, upper
        real, intent(out) :: out
        real :: yy

        call rand_lcg(yy)
        out = yy * (upper - lower) + lower
    end subroutine randrange_lcg

    subroutine normal_lcg(out)
        implicit none
        real, intent(out) :: out
        real :: sum_randnormal, dummy_randnormal, r_randnormal, r2_randnormal
        integer :: i_rnc

        sum_randnormal = 0.0

        do i_rnc = 1, 12
            call random_number(dummy_randnormal)
            sum_randnormal = sum_randnormal + dummy_randnormal
        enddo

        r_randnormal  = ( sum_randnormal - 6.0 ) / 4.0
        r2_randnormal = r_randnormal * r_randnormal

        out = (((( A9 * r2_randnormal + A7 ) * r2_randnormal + A5 ) * r2_randnormal + A3 ) * r2_randnormal + A1 ) * r_randnormal

    end subroutine normal_lcg

    function normal_distribution( rand ) result(out)
        implicit none
        real, intent(in) :: rand
        real :: out, saead
        real :: sum_randnormal, dummy_randnormal, r_randnormal, r2_randnormal
        integer :: i_rnc

        sum_randnormal = 0.0

        do i_rnc = 1, 12
            call random_number(dummy_randnormal)
            sum_randnormal = sum_randnormal + dummy_randnormal
        enddo

        r_randnormal  = ( sum_randnormal - 6.0 ) / 4.0
        r2_randnormal = r_randnormal * r_randnormal

        out = (((( A9 * r2_randnormal + A7 ) * r2_randnormal + A5 ) * r2_randnormal + A3 ) * r2_randnormal + A1 ) * r_randnormal


    end function normal_distribution
end module lcg_random

module xmath
    use iso_fortran_env, only: real64
    implicit none
    real, parameter :: pi = 3.1415926535897932384626433832795

    interface operator(.x.)
        module procedure vmdot
        module procedure mvdot
        module procedure mmdot
    end interface

    contains

    pure function mag(v) result (out)
        implicit none
        real(real64), intent(in), dimension(3) :: v
        real(real64) :: out

        out = sqrt(sum(v**2))

    end function mag

    pure function mag2(v) result (out)
        implicit none
        real(real64), intent(in), dimension(3) :: v
        real(real64) :: out

        out = sum(v**2)

    end function mag2

    pure function normalize(v) result (out)
        implicit none
        real(real64), intent(in), dimension(3) :: v
        real(real64), dimension(size(v)) :: out

        out = v / mag(v)

    end function normalize

    pure function dot(v1, v2) result (out)
        implicit none
        real(real64), intent(in), dimension(3) :: v1, v2
        real(real64) :: out

        out = sum(v1 * v2)

    end function

    pure function cross(a, b) result (out)
        implicit none
        real(real64), intent(in), dimension(3) :: a, b
        real(real64), dimension(3) :: out

        out(1) = a(2)*b(3) - a(3)*b(2)
        out(2) = a(3)*b(1) - a(1)*b(3)
        out(3) = a(1)*b(2) - a(2)*b(1)

    end function cross

    pure function determinant(M) result (out)
        implicit none
        real(real64), intent(in), dimension(3,3) :: M
        real(real64) :: out

        out = M(1, 1) * ( M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3) ) &
            - M(1, 2) * ( M(2, 1) * M(3, 3) - M(3, 1) * M(2, 3) ) &
            + M(1, 3) * ( M(2, 1) * M(3, 2) - M(3, 1) * M(2, 2) )

    end function determinant

    pure function inverse(M) result (out)
        implicit none
        real(real64), intent(in), dimension(3,3) :: M
        real(real64), dimension(3,3) :: out
        real(real64) :: det

        out(1, 1) = M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3)
        out(1, 2) = M(1, 3) * M(3, 2) - M(1, 2) * M(3, 3)
        out(1, 3) = M(1, 2) * M(2, 3) - M(1, 3) * M(2, 2)

        out(2, 1) = M(2, 3) * M(3, 1) - M(2, 1) * M(3, 3)
        out(2, 2) = M(1, 1) * M(3, 3) - M(1, 3) * M(3, 1)    
        out(2, 3) = M(1, 3) * M(2, 1) - M(1, 1) * M(2, 3)

        out(3, 1) = M(2, 1) * M(3, 2) - M(2, 2) * M(3, 1)
        out(3, 2) = M(1, 2) * M(3, 1) - M(1, 1) * M(3, 2)
        out(3, 3) = M(1, 1) * M(2, 2) - M(1, 2) * M(2, 1)

        det = determinant(M)
        out = out / det

    end function inverse

    pure function as_matrix( w, x, y, z ) result (out)
        implicit none
        real(real64), intent(in)     :: w, x, y, z
        real(real64), dimension(3,3) :: out

        out(1, 1) = 1 - 2 * y ** 2 - 2 * z ** 2
        out(1, 2) = 2 * x * y - 2 * w * z
        out(1, 3) = 2 * w * y + 2 * x * z

        out(2, 1) = 2 * x * y + 2 * w * z
        out(2, 2) = 1 - 2 * x ** 2 - 2 * z ** 2
        out(2, 3) = 2 * y * z - 2 * w * x

        out(3, 1) = 2 * w * z - 2 * x * y
        out(3, 2) = 2 * y * z + 2 * w * x
        out(3, 3) = 1 - 2 * x ** 2 - 2 * y ** 2

    end function as_matrix

    pure function rotation_matrix( w, x, y, z ) result (rot_matrix)
        real(real64), intent(in)   :: w, x, y, z
        real, dimension(3, 3) :: rot_matrix

        rot_matrix(1, 1) = w ** 2 + x ** 2 - y ** 2 - z ** 2
        rot_matrix(1, 2) = 2.0 * ( x * y + w * z )
        rot_matrix(1, 3) = 2.0 * ( x * z - w * y )
        rot_matrix(2, 1) = 2.0 * ( x * y - w * z )
        rot_matrix(2, 2) = w ** 2 - x ** 2 + y ** 2 - z ** 2
        rot_matrix(2, 3) = 2.0 * ( y * z + w * x )
        rot_matrix(3, 1) = 2.0 * ( x * z + w * y )
        rot_matrix(3, 2) = 2.0 * ( y * z - w * x )
        rot_matrix(3, 3) = w ** 2 - x ** 2 - y ** 2 + z ** 2

    end function rotation_matrix

    pure function vvdot(a, b) result (out)
        implicit none
        real(real64), intent(in), dimension(3) :: a, b
        real(real64) :: out

        out = sum(a * b)

    end function vvdot

    pure function vmdot(a, b) result (out)
        implicit none
        real(real64), intent(in), dimension(3) :: a
        real(real64), intent(in), dimension(3, 3) :: b
        real(real64), dimension(3) :: out

        out(1) = sum(a(:) * b(:, 1))
        out(2) = sum(a(:) * b(:, 2))
        out(3) = sum(a(:) * b(:, 3))

    end function vmdot

    pure function mvdot(a, b) result (out)
        implicit none
        real(real64), intent(in), dimension(3, 3) :: a
        real(real64), intent(in), dimension(3) :: b
        real(real64), dimension(3) :: out

        out(1) = sum(a(1, :) * b(:))
        out(2) = sum(a(2, :) * b(:))
        out(3) = sum(a(3, :) * b(:))
    end function mvdot

    pure function mmdot(a, b) result (out)
        implicit none
        real(real64), intent(in), dimension(3, 3) :: a, b
        real(real64), dimension(3, 3) :: out

        out(1, 1) = sum(a(1, :) * b(:, 1))
        out(1, 2) = sum(a(1, :) * b(:, 2))
        out(1, 3) = sum(a(1, :) * b(:, 3))

        out(2, 1) = sum(a(2, :) * b(:, 1))
        out(2, 2) = sum(a(2, :) * b(:, 2))
        out(2, 3) = sum(a(2, :) * b(:, 3))

        out(3, 1) = sum(a(3, :) * b(:, 1))
        out(3, 2) = sum(a(3, :) * b(:, 2))
        out(3, 3) = sum(a(3, :) * b(:, 3))

    end function mmdot

    pure function quat_product(q1, q2) result (out)
        implicit none
        real(real64), dimension(4), intent(in) :: q1, q2
        real(real64), dimension(4)             :: out

        out(1) = q1(1) * q2(1) - q1(2) * q2(2) - q1(3) * q2(3) - q1(4) * q2(4)
        out(2) = q1(1) * q2(2) + q1(2) * q2(1) + q1(3) * q2(4) - q1(4) * q2(3)
        out(3) = q1(1) * q2(3) - q1(2) * q2(4) + q1(3) * q2(1) + q1(4) * q2(2)
        out(4) = q1(1) * q2(4) + q1(2) * q2(3) - q1(3) * q2(2) + q1(4) * q2(1)

    end function quat_product

    pure function degrees_to_radians(angle) result (out)
        implicit none
        real(real64), intent(in) :: angle
        real(real64) :: out

        out = angle * pi / 180.0

    end function degrees_to_radians

    pure function radians_to_degrees(angle) result (out)
        implicit none
        real(real64), intent(in) :: angle
        real(real64) :: out

        out = angle * 180.0 / pi

    end function radians_to_degrees

    subroutine normal_distribution( rand )
        implicit none
        real(real64), intent(out) :: rand
        real(real64) :: sum_randnormal, dummy_randnormal, r_randnormal, r2_randnormal
        integer :: i_rnc

        real(real64), parameter :: A1 = 3.949846138
        real(real64), parameter :: A3 = 0.252408784
        real(real64), parameter :: A5 = 0.076542912
        real(real64), parameter :: A7 = 0.008355968
        real(real64), parameter :: A9 = 0.029899776

        sum_randnormal = 0.0

        do i_rnc = 1, 12
            call random_number(dummy_randnormal)
            sum_randnormal = sum_randnormal + dummy_randnormal
        enddo

        r_randnormal  = ( sum_randnormal - 6.0 ) / 4.0
        r2_randnormal = r_randnormal * r_randnormal

        rand = (((( A9 * r2_randnormal + A7 ) * r2_randnormal + A5 ) * r2_randnormal + A3 ) * r2_randnormal + A1 ) * r_randnormal

    end subroutine normal_distribution

    subroutine normal_sequences( num, noise1, noise2 )
        implicit none
        integer, intent(in) :: num
        real(real64), dimension(num), intent(out) :: noise1, noise2
        real(real64), dimension(num) :: unf_noise1, unf_noise2

        call random_number(unf_noise1)
        call random_number(unf_noise2)

        noise1 = sqrt(-2.0 * log(unf_noise1)) * cos(2.0 * pi * unf_noise2)
        noise2 = sqrt(-2.0 * log(unf_noise1)) * sin(2.0 * pi * unf_noise2)

    end subroutine normal_sequences

    PURE FUNCTION polyval ( x, c ) RESULT ( f )
        IMPLICIT NONE
        real(real64)                            :: f ! Returns polynomial in ...
        real(real64),                INTENT(in) :: x ! argument
        real(real64), DIMENSION(0:), INTENT(in) :: c ! given coefficients (ascending powers of x)

        ! Uses Horner's rule
        
        INTEGER :: i, upper

        upper = UBOUND(c,1)
        f = c(upper)
        DO i = upper - 1, 0, -1
        f = f * x + c(i)
        END DO
    END FUNCTION polyval

    PURE FUNCTION exprel ( x ) RESULT ( f )
        IMPLICIT NONE
        real(real64), INTENT(in) :: x ! Argument
        real(real64)             :: f ! Returns value of (exp(x)-1)/x

        ! At small x, we must guard against the ratio of imprecise small values.
        ! There are various ways of doing this.
        ! We follow some others and use the identity: (exp(x)-1)/x = exp(x/2)*[sinh(x/2)/(x/2)].
        ! For small x, sinh(x)/x = g0 + g1*x**2 + g2*x**4 + ...
        ! where the coefficient of x**(2n) is gn = 1/(2*n+1)!
        ! Alternatively, the exprel function is available in some math and scientific libraries.

        real(real64), DIMENSION(0:4), PARAMETER :: g = 1.0 / [1,6,120,5040,362880]
        real(real64),                 PARAMETER :: tol = 0.01

        IF ( ABS(x) > tol ) THEN
        f = ( EXP(x) - 1.0 ) / x
        ELSE
        f = EXP(x/2) * polyval ( (x/2)**2, g )
        END IF

    END FUNCTION exprel


end module