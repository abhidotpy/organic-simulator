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
    implicit none
    real, parameter :: pi = 3.1415926535897932384626433832795
    real, parameter :: e = 2.7182818284590452353602874713527

    type, public :: Quaternion
        real(8), dimension(4) :: q
    end type

    interface operator(.x.)
        module procedure vmdot
        module procedure mvdot
        module procedure mmdot
    end interface

    interface operator(*)
        module procedure quat_product
        module procedure quat_vec_product
    end interface

    contains

    pure function mag(v) result (out)
        implicit none
        real(8), intent(in), dimension(3) :: v
        real(8) :: out

        out = sqrt(sum(v**2))

    end function

    pure function mag2(v) result (out)
        implicit none
        real(8), intent(in), dimension(3) :: v
        real(8) :: out

        out = sum(v**2)

    end function

    pure function normalize(v) result (out)
        implicit none
        real(8), intent(in), dimension(3) :: v
        real(8), dimension(size(v)) :: out

        out = v / mag(v)

    end function normalize

    pure function dot(v1, v2) result (out)
        implicit none
        real(8), intent(in), dimension(3) :: v1, v2
        real(8) :: out

        out = sum(v1 * v2)

    end function

    pure function cross(a, b) result (out)
        implicit none
        real(8), intent(in), dimension(3) :: a, b
        real(8), dimension(3) :: out

        out(1) = a(2)*b(3) - a(3)*b(2)
        out(2) = a(3)*b(1) - a(1)*b(3)
        out(3) = a(1)*b(2) - a(2)*b(1)

    end function cross

    pure function determinant(M) result (out)
        implicit none
        real(8), intent(in), dimension(3,3) :: M
        real(8) :: out

        out = M(1, 1) * ( M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3) ) &
            - M(1, 2) * ( M(2, 1) * M(3, 3) - M(3, 1) * M(2, 3) ) &
            + M(1, 3) * ( M(2, 1) * M(3, 2) - M(3, 1) * M(2, 2) )

    end function determinant

    pure function inverse(M) result (out)
        implicit none
        real(8), intent(in), dimension(3,3) :: M
        real(8), dimension(3,3) :: out
        real(8) :: det

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

    pure function quat_conjugate(quat) result (out)
        implicit none
        type(Quaternion), intent(in) :: quat
        type(Quaternion) :: out

        out % q(1) =  quat % q(1)
        out % q(2) = -quat % q(2)
        out % q(3) = -quat % q(3)
        out % q(4) = -quat % q(4)

    end function quat_conjugate

    pure function as_matrix(quat) result (out)
        implicit none
        type(Quaternion), intent(in) :: quat
        real(8), dimension(3,3) :: out

        out(1, 1) = 1 - 2 * quat % q(3)**2 - 2 * quat % q(4)**2
        out(1, 2) = 2 * quat % q(2) * quat % q(3) - 2 * quat % q(1) * quat % q(4)
        out(1, 3) = 2 * quat % q(1) * quat % q(3) + 2 * quat % q(2) * quat % q(4)

        out(2, 1) = 2 * quat % q(2) * quat % q(3) + 2 * quat % q(1) * quat % q(4)
        out(2, 2) = 1 - 2 * quat % q(2)**2 - 2 * quat % q(4)**2
        out(2, 3) = 2 * quat % q(3) * quat % q(4) - 2 * quat % q(1) * quat % q(2)

        out(3, 1) = 2 * quat % q(1) * quat % q(4) - 2 * quat % q(2) * quat % q(3)
        out(3, 2) = 2 * quat % q(3) * quat % q(4) + 2 * quat % q(1) * quat % q(2)
        out(3, 3) = 1 - 2 * quat % q(2)**2 - 2 * quat % q(3)**2

    end function as_matrix

    pure function vvdot(a, b) result (out)
        implicit none
        real(8), intent(in), dimension(3) :: a, b
        real(8) :: out

        out = sum(a * b)

    end function

    pure function vmdot(a, b) result (out)
        implicit none
        real(8), intent(in), dimension(3) :: a
        real(8), intent(in), dimension(3, 3) :: b
        real(8), dimension(3) :: out

        out(1) = sum(a(:) * b(:, 1))
        out(2) = sum(a(:) * b(:, 2))
        out(3) = sum(a(:) * b(:, 3))

    end function vmdot

    pure function mvdot(a, b) result (out)
        implicit none
        real(8), intent(in), dimension(3, 3) :: a
        real(8), intent(in), dimension(3) :: b
        real(8), dimension(3) :: out

        out(1) = sum(a(1, :) * b(:))
        out(2) = sum(a(2, :) * b(:))
        out(3) = sum(a(3, :) * b(:))
    end function mvdot

    pure function mmdot(a, b) result (out)
        implicit none
        real(8), intent(in), dimension(3, 3) :: a, b
        real(8), dimension(3, 3) :: out

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

    pure function quat_product(a, b) result (out)
        implicit none
        type(Quaternion), intent(in) :: a, b
        type(Quaternion) :: out

        out % q(1) = a % q(1) * b % q(1) - a % q(2) * b % q(2) - a % q(3) * b % q(3) - a % q(4) * b % q(4)
        out % q(2) = a % q(1) * b % q(2) + a % q(2) * b % q(1) + a % q(3) * b % q(4) - a % q(4) * b % q(3)
        out % q(3) = a % q(1) * b % q(3) - a % q(2) * b % q(4) + a % q(3) * b % q(1) + a % q(4) * b % q(2)
        out % q(4) = a % q(1) * b % q(4) + a % q(2) * b % q(3) - a % q(3) * b % q(2) + a % q(4) * b % q(1)

    end function quat_product

    pure function quat_vec_product(a, b) result (out)
        implicit none
        type(Quaternion), intent(in) :: a
        real(8), intent(in), dimension(3) :: b
        real(8), dimension(3) :: out

        out(1) = a % q(1) * b(1) + a % q(2) * b(2) + a % q(3) * b(3)
        out(2) = a % q(2) * b(1) - a % q(1) * b(2) - a % q(4) * b(3)
        out(3) = a % q(3) * b(1) + a % q(4) * b(2) - a % q(1) * b(3)

    end function quat_vec_product

    pure function degrees_to_radians(angle) result (out)
        implicit none
        real(8), intent(in) :: angle
        real(8) :: out

        out = angle * pi / 180.0

    end function degrees_to_radians

    pure function radians_to_degrees(angle) result (out)
        implicit none
        real(8), intent(in) :: angle
        real(8) :: out

        out = angle * 180.0 / pi

    end function radians_to_degrees

    pure function golden_section_search(func, args) result(res)
        implicit none
        integer, parameter :: max_iter = 1000
        real(8) :: eps, phi, a, b, x1, x2, fx1, fx2
        real(8) :: res
        real(8), dimension(:), intent(in) :: args
        integer :: count

        interface
            pure function func(x, argv) result (out)
                implicit none
                real(8), intent(in) :: x
                real(8), dimension(:), intent(in) :: argv
                real(8) :: out
            end function func
        end interface

        eps = 1e-8; phi = 0.5 * ( 3.0 - sqrt(5.0) )
        a = 0.0; b = 1.0;
        x1 = (1 - phi) * a + phi * b
        x2 = phi * a + (1 - phi) * b
        fx1 = func(x1, args)
        fx2 = func(x2, args)
        count = 0;

        do while ( abs(x2 - x1) > eps .and. count < max_iter )
            count = count + 1

            if ( fx1 < fx2 ) then
                b = x2
                x2 = x1;  fx2 = fx1
                x1 = (1 - phi) * a + phi * b
                fx1 = func(x1, args)
            
            else
                a = x1
                x1 = x2; fx1 = fx2
                x2 = phi * a + (1 - phi) * b
                fx2 = func(x2, args)
            endif
        enddo

        res = (x1 + x2) / 2.0

    end function golden_section_search

    ! real function brent(func, args) result (res)
    !     implicit none
    !     integer, parameter :: max_iter = 10000
    !     real(8), parameter    :: eps = 1e-8, phi = 0.5 * ( 3.0 - sqrt(5.0) )
    !     real(8)               :: a, b, c, x, w, v, u
    !     real(8)               :: fa, fb, fc, fx, fw, fv, fu
    !     real(8)               :: deltax = 0.0, atol = 1e-11, rtol = 1e-8
    !     real(8)               :: tol1, tol2, xmid, tmp1, tmp2, rat, p, dx_temp
    !     integer            :: I, iter = 0

    !     interface
    !         real function func(xx, argv)
    !             real(8), intent(in) :: xx
    !             real(8), intent(in) :: argv(:)
    !         end function func
    !     end interface

    !     a = 0.0; b = 1.0 
    !     c = (1 - phi) * a + phi * b
    !     fa = func(a); fb = func(b); fc = func(c)
    !     v = c; w = v; x = w
    !     fv = fc; fw = fv; fx = fw

    !     do while (iter < max_iter)
    !         tol1 = rtol * abs(x) + atol
    !         tol2 = 2.0 * tol1
    !         xmid = 0.5 * (a + b)

    !         if (abs(x - xmid) < (tol2 - 0.5 * (b - a))) exit

    !         if (abs(deltax) <= tol1) then
    !             if (x >= xmid) then
    !                 deltax = a - x
    !             else
    !                 deltax = b - x
    !             endif
    !             rat = phi * deltax

    !         else
    !             tmp1 = (x - w) * (fx - fv)
    !             tmp2 = (x - v) * (fx - fw)
    !             p = (x - v) * tmp2 - (x - w) * tmp1
    !             tmp2 = 2.0 * (tmp2 - tmp1)

    !             if (tmp2 > 0.0) p = -p

    !             tmp2 = abs(tmp2)
    !             dx_temp = deltax
    !             deltax = rat
                
    !             if ((p > tmp2 * (a - x)) .and. (p < tmp2 * (b - x)) .and. (abs(p) < abs(0.5 * tmp2 * dx_temp))) then
    !                 rat = p * 1.0 / tmp2
    !                 u = x + rat
    !                 if ((u - a) < tol2 .or. (b - u) < tol2) then
    !                     if (xmid - x >= 0) then
    !                         rat = tol1
    !                     else
    !                         rat = -tol1
    !                     endif
    !                 endif
    !             else
    !                 if (x >= xmid) then
    !                     deltax = a - x 
    !                 else
    !                     deltax = b - x
    !                 endif
    !                 rat = phi * deltax
    !             endif
    !         endif

    !         if (abs(rat) < tol1) then
    !             if (rat >= 0) then
    !                 u = x + tol1
    !             else
    !                 u = x - tol1
    !             endif
    !         else
    !             u = x + rat

    !         endif

    !         fu = func(u)

    !         if (fu > fx) then

    !             if (u < x) then
    !                 a = u
    !             else
    !                 b = u
    !             endif

    !             if ((fu <= fw) .or. (w == x)) then
    !                 v = w
    !                 w = u
    !                 fv = fw
    !                 fw = fu

    !             else if ((fu <= fv) .or. (v == x) .or. (v == w)) then
    !                 v = u
    !                 fv = fu
    !             endif

    !         else

    !             if (u >= x) then
    !                 a = x
    !             else
    !                 b = x
    !             endif

    !             v = w
    !             w = x
    !             x = u
    !             fv = fw
    !             fw = fx
    !             fx = fu

    !         endif
    !         iter = iter + 1
    
    !     enddo

    !     res = x
    ! end function brent

    subroutine normal_distribution( rand )
        implicit none
        real(8), intent(out) :: rand
        real(8) :: sum_randnormal, dummy_randnormal, r_randnormal, r2_randnormal
        integer :: i_rnc

        real(8), parameter :: A1 = 3.949846138
        real(8), parameter :: A3 = 0.252408784
        real(8), parameter :: A5 = 0.076542912
        real(8), parameter :: A7 = 0.008355968
        real(8), parameter :: A9 = 0.029899776

        sum_randnormal = 0.0

        do i_rnc = 1, 12
            call random_number(dummy_randnormal)
            sum_randnormal = sum_randnormal + dummy_randnormal
        enddo

        r_randnormal  = ( sum_randnormal - 6.0 ) / 4.0
        r2_randnormal = r_randnormal * r_randnormal

        rand = (((( A9 * r2_randnormal + A7 ) * r2_randnormal + A5 ) * r2_randnormal + A3 ) * r2_randnormal + A1 ) * r_randnormal

    end subroutine normal_distribution

    subroutine normal_sequences(num, noise1, noise2 )
        implicit none
        integer, intent(in) :: num
        real, dimension(num), intent(out) :: noise1, noise2
        real, dimension(num) :: unf_noise1, unf_noise2

        call random_number(unf_noise1)
        call random_number(unf_noise2)

        noise1 = sqrt(-2.0 * log(unf_noise1)) * cos(2.0 * pi * unf_noise2)
        noise2 = sqrt(-2.0 * log(unf_noise1)) * sin(2.0 * pi * unf_noise2)

    end subroutine

end module