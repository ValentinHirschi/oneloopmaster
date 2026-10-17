! Independent value driver for the unchanged original OneLOop tensor routines.
! Build against an original double-precision OneLOop module/object, e.g.
!   gfortran -O2 -fno-fast-math -J "$D" -c /path/to/avh_olo.f90 -o "$D/avh_olo.o"
!   gfortran -O2 -I "$D" tests/support/tensor_oracle.f90 "$D/avh_olo.o" -o tensor_oracle
! Input: one row per line (blank lines and lines starting with # are skipped):
!   5 mu_squared m_re m_im                          olo_an, rank 4: A0, A00, A0000
!   6 mu_squared p_squared m0_re m0_im m1_re m1_im  olo_b11: B0, B1, B11
!   7 mu_squared p_squared m0_re m0_im m1_re m1_im  olo_bn, rank 3: B0, B1, B11, B111
! The scale is supplied squared and passed to OneLOop unsquared (rmu).
! Output per row: VALUE kind, then real/imaginary parts of the finite, simple-
! pole and double-pole coefficients of each tensor coefficient in order.
program tensor_oracle
    use avh_olo
    implicit none
    integer, parameter :: dp = kind(1d0)
    integer :: family, status, j, k
    real(dp) :: values(7), mu
    complex(dp) :: p, m0, m1, an(0:2, 0:2), bn(0:2, 0:5)
    complex(dp) :: b11(0:2), b00(0:2), b1(0:2), b0(0:2)
    character(4096) :: line

    do
        read(*, '(a)', iostat=status) line
        if (status /= 0) exit
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#') cycle
        read(line, *) family
        select case (family)
        case (5)
            read(line, *) values(1:4)
            mu = sqrt(values(2))
            m0 = cmplx(values(3), values(4), dp)
            call olo_an(an, 4, m0, mu)
            write(*, '(a,1x,i0,18(1x,es26.17))') 'VALUE', family, &
                ((real(an(j, k), dp), aimag(an(j, k)), j = 0, 2), k = 0, 2)
        case (6)
            read(line, *) values(1:7)
            mu = sqrt(values(2))
            p = cmplx(values(3), 0.0_dp, dp)
            m0 = cmplx(values(4), values(5), dp)
            m1 = cmplx(values(6), values(7), dp)
            call olo_b11(b11, b00, b1, b0, p, m0, m1, mu)
            write(*, '(a,1x,i0,18(1x,es26.17))') 'VALUE', family, &
                (real(b0(j), dp), aimag(b0(j)), j = 0, 2), &
                (real(b1(j), dp), aimag(b1(j)), j = 0, 2), &
                (real(b11(j), dp), aimag(b11(j)), j = 0, 2)
        case (7)
            read(line, *) values(1:7)
            mu = sqrt(values(2))
            p = cmplx(values(3), 0.0_dp, dp)
            m0 = cmplx(values(4), values(5), dp)
            m1 = cmplx(values(6), values(7), dp)
            call olo_bn(bn, 3, p, m0, m1, mu)
            write(*, '(a,1x,i0,24(1x,es26.17))') 'VALUE', family, &
                (real(bn(j, 0), dp), aimag(bn(j, 0)), j = 0, 2), &
                (real(bn(j, 1), dp), aimag(bn(j, 1)), j = 0, 2), &
                (real(bn(j, 3), dp), aimag(bn(j, 3)), j = 0, 2), &
                (real(bn(j, 5), dp), aimag(bn(j, 5)), j = 0, 2)
        case default
            error stop 'Unsupported tensor oracle kind'
        end select
    end do
end program tensor_oracle
