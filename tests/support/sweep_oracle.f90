! Batch driver for permutation/regression sweeps against unchanged OneLOop.
! Build (module and archive paths refer to an existing original build):
! gfortran -O2 -I /path/to/OneLOop tests/support/sweep_oracle.f90 \
!   /path/to/OneLOop/libavh_olo.a -o /tmp/oneloop-sweep-oracle
! Input: one row per line in the scalar fixture order without references:
!   kind mu_squared momenta... (squared_mass_re squared_mass_im)...
! Kinds 1,2,3,4 are A0,B0,C0,D0; -2 is dB0/dp^2.
! Output per row: "dp" and "qp" lines with the three Laurent coefficients
! [finite, 1/eps, 1/eps^2] as real/imaginary pairs. The quadruple-precision
! library evaluates the identical binary64 inputs; disagreement between both
! lines flags points where the binary64 reference itself is unstable.
program oneloop_sweep_oracle
    use avh_olo, only: olo, olo_db0, olo_scale, olo_unit
    use iso_fortran_env, only: real64
    implicit none
    integer, parameter :: qp = 16
    integer :: kind, np, nm, nfields, i, status
    real(real64) :: values(22), mu2
    complex(real64) :: p(6), m(4), dp_result(0:2)
    complex(qp) :: pq(6), mq(4), qp_result(0:2)
    character(4096) :: line

    call olo_unit(-1)
    do
        read(*, '(a)', iostat=status) line
        if (status /= 0) exit
        if (len_trim(line) == 0) cycle
        read(line, *, iostat=status) kind
        if (status /= 0) error stop 'Missing row kind'
        select case (kind)
        case (1)
            np = 0; nm = 1
        case (2, -2)
            np = 1; nm = 2
        case (3)
            np = 3; nm = 3
        case (4)
            np = 6; nm = 4
        case default
            error stop 'Unsupported sweep kind'
        end select
        nfields = 2 + np + 2*nm
        read(line, *, iostat=status) values(1:nfields)
        if (status /= 0) error stop 'Malformed sweep row'
        mu2 = values(2)
        p = (0.0_real64, 0.0_real64)
        m = (0.0_real64, 0.0_real64)
        do i = 1, np
            p(i) = cmplx(values(2 + i), 0.0_real64, real64)
        end do
        do i = 1, nm
            m(i) = cmplx(values(2 + np + 2*i - 1), values(2 + np + 2*i), real64)
        end do
        pq = cmplx(real(p, real64), aimag(p), qp)
        mq = cmplx(real(m, real64), aimag(m), qp)
        call olo_scale(sqrt(mu2))
        select case (kind)
        case (1)
            call olo(dp_result, m(1))
            call olo(qp_result, mq(1))
        case (2)
            call olo(dp_result, p(1), m(1), m(2))
            call olo(qp_result, pq(1), mq(1), mq(2))
        case (-2)
            call olo_db0(dp_result, p(1), m(1), m(2))
            call olo_db0(qp_result, pq(1), mq(1), mq(2))
        case (3)
            call olo(dp_result, p(1), p(2), p(3), m(1), m(2), m(3))
            call olo(qp_result, pq(1), pq(2), pq(3), mq(1), mq(2), mq(3))
        case (4)
            call olo(dp_result, p(1), p(2), p(3), p(4), p(5), p(6), m(1), m(2), m(3), m(4))
            call olo(qp_result, pq(1), pq(2), pq(3), pq(4), pq(5), pq(6), &
                mq(1), mq(2), mq(3), mq(4))
        end select
        write(*, '(a,6es26.17e3)') 'dp ', (real(dp_result(i)), aimag(dp_result(i)), i = 0, 2)
        write(*, '(a,6es26.17e3)') 'qp ', &
            (real(real(qp_result(i), qp), real64), real(aimag(qp_result(i)), real64), i = 0, 2)
        flush(6)
    end do
end program oneloop_sweep_oracle
