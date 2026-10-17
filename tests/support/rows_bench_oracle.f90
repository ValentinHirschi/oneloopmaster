! Warmed throughput driver for unchanged original OneLOop over mixed master rows,
! including the tensor routines olo_an, olo_b11 and olo_bn.
! Build: gfortran -O2 -I "$D" tests/support/rows_bench_oracle.f90 "$D/avh_olo.o" -o rows_bench
! Input (stdin): OLO_ROWS_BENCH_V1, then "nrows warmup calls repetitions real",
! then nrows rows "kind mu_squared arguments...":
!   1 A0:   m_re m_im                 5 An (rank 4):  m_re m_im
!   2 B0:   p m0 m1                   6 olo_b11:      p m0 m1
!   3 C0:   p1 p2 p3 m0 m1 m2         7 olo_bn rank 3: p m0 m1
!   4 D0:   p1 p2 p3 p4 s12 s23 m0 m1 m2 m3
! Masses are re/im pairs and the scale is squared. With real=1 and all masses
! of a row real, the real-argument interfaces (Ninja's *_rm wrappers) are used.
! Workloads: SAME repeats each row, HETERO cycles the rows of one kind and ALL
! cycles all rows in input order. Every call computes all requested outputs.
! A timed block has at least `calls` calls, raised after an untimed calibration
! block so that it lasts at least 25 ms.
program rows_bench_oracle
    use avh_olo
    use iso_fortran_env, only: int64, real64
    implicit none
    type row_type
        integer :: kind
        logical :: real_masses
        real(real64) :: mu
        real(real64) :: p(6), mr(4)
        complex(real64) :: m(4)
    end type row_type
    type(row_type), allocatable :: rows(:)
    integer, allocatable :: selection(:)
    integer :: nrows, repetitions, real_flag, row, kind, nselected, np, nm, j, status
    integer(int64) :: warmup, calls, clock_rate, clock_max
    real(real64) :: values(24)
    character(4096) :: line

    read(*, '(a)', iostat=status) line
    if (status /= 0 .or. trim(line) /= 'OLO_ROWS_BENCH_V1') error stop 'Expected OLO_ROWS_BENCH_V1'
    read(*, *) nrows, warmup, calls, repetitions, real_flag
    allocate(rows(nrows), selection(nrows))
    do row = 1, nrows
        read(*, '(a)') line
        read(line, *) rows(row)%kind
        call dimensions(rows(row)%kind, np, nm)
        read(line, *) values(1:2 + np + 2*nm)
        rows(row)%mu = sqrt(values(2))
        rows(row)%p = 0
        rows(row)%p(1:np) = values(3:2 + np)
        do j = 1, nm
            rows(row)%m(j) = cmplx(values(1 + np + 2*j), values(2 + np + 2*j), real64)
            rows(row)%mr(j) = values(1 + np + 2*j)
        end do
        rows(row)%real_masses = real_flag /= 0 .and. all(aimag(rows(row)%m(1:nm)) == 0)
    end do
    do row = 1, nrows
        call print_values(row, rows(row))
    end do
    call system_clock(count_rate=clock_rate, count_max=clock_max)
    write(*, '(a,1x,i0,1x,i0)') 'CLOCK', clock_rate, clock_max

    do row = 1, nrows
        selection(1) = row
        call workload('SAME', rows(row)%kind, row, selection(1:1))
    end do
    do kind = 1, 7
        nselected = 0
        do row = 1, nrows
            if (rows(row)%kind /= kind) cycle
            nselected = nselected + 1
            selection(nselected) = row
        end do
        if (nselected > 0) call workload('HETERO', kind, 0, selection(1:nselected))
    end do
    selection = [(row, row = 1, nrows)]
    call workload('ALL', 0, 0, selection)

contains

    subroutine dimensions(kind, momenta, masses)
        integer, intent(in) :: kind
        integer, intent(out) :: momenta, masses
        select case (kind)
        case (1, 5)
            momenta = 0; masses = 1
        case (2, 6, 7)
            momenta = 1; masses = 2
        case (3)
            momenta = 3; masses = 3
        case (4)
            momenta = 6; masses = 4
        case default
            error stop 'Unsupported row kind'
        end select
    end subroutine dimensions

    ! VALUE row kind, then re/im of [finite, 1/eps, 1/eps^2] per coefficient.
    subroutine print_values(index, r)
        integer, intent(in) :: index
        type(row_type), intent(in) :: r
        complex(real64) :: out(0:2), an(0:2, 0:2), bn(0:2, 0:5), b11(0:2), b00(0:2), b1(0:2), b0(0:2)
        complex(real64) :: flat(12)
        integer :: n, k
        select case (r%kind)
        case (1, 2, 3, 4)
            select case (r%kind)
            case (1)
                call olo_a0(out, r%m(1), r%mu)
            case (2)
                call olo_b0(out, r%p(1), r%m(1), r%m(2), r%mu)
            case (3)
                call olo_c0(out, r%p(1), r%p(2), r%p(3), r%m(1), r%m(2), r%m(3), r%mu)
            case (4)
                call olo_d0(out, r%p(1), r%p(2), r%p(3), r%p(4), r%p(5), r%p(6), &
                    r%m(1), r%m(2), r%m(3), r%m(4), r%mu)
            end select
            flat(1:3) = out
            n = 3
        case (5)
            call olo_an(an, 4, r%m(1), r%mu)
            flat(1:9) = [an(:, 0), an(:, 1), an(:, 2)]
            n = 9
        case (6)
            call olo_b11(b11, b00, b1, b0, cmplx(r%p(1), 0.0_real64, real64), r%m(1), r%m(2), r%mu)
            flat(1:9) = [b0, b1, b11]
            n = 9
        case (7)
            call olo_bn(bn, 3, cmplx(r%p(1), 0.0_real64, real64), r%m(1), r%m(2), r%mu)
            flat(1:12) = [bn(:, 0), bn(:, 1), bn(:, 3), bn(:, 5)]
            n = 12
        end select
        write(*, '(a,1x,i0,1x,i0,24(1x,es26.17))') 'VALUE', index, r%kind, &
            (real(flat(k), real64), aimag(flat(k)), k = 1, n)
    end subroutine print_values

    ! Returns the finite part of the first coefficient, for the checksum.
    function evaluate(r) result(first)
        type(row_type), intent(in) :: r
        complex(real64) :: first
        complex(real64) :: out(0:2), an(0:2, 0:2), bn(0:2, 0:5), b11(0:2), b00(0:2), b1(0:2), b0(0:2)
        complex(real64) :: pc
        select case (r%kind)
        case (1)
            if (r%real_masses) then
                call olo_a0(out, r%mr(1), r%mu)
            else
                call olo_a0(out, r%m(1), r%mu)
            end if
            first = out(0)
        case (2)
            if (r%real_masses) then
                call olo_b0(out, r%p(1), r%mr(1), r%mr(2), r%mu)
            else
                call olo_b0(out, r%p(1), r%m(1), r%m(2), r%mu)
            end if
            first = out(0)
        case (3)
            if (r%real_masses) then
                call olo_c0(out, r%p(1), r%p(2), r%p(3), r%mr(1), r%mr(2), r%mr(3), r%mu)
            else
                call olo_c0(out, r%p(1), r%p(2), r%p(3), r%m(1), r%m(2), r%m(3), r%mu)
            end if
            first = out(0)
        case (4)
            if (r%real_masses) then
                call olo_d0(out, r%p(1), r%p(2), r%p(3), r%p(4), r%p(5), r%p(6), &
                    r%mr(1), r%mr(2), r%mr(3), r%mr(4), r%mu)
            else
                call olo_d0(out, r%p(1), r%p(2), r%p(3), r%p(4), r%p(5), r%p(6), &
                    r%m(1), r%m(2), r%m(3), r%m(4), r%mu)
            end if
            first = out(0)
        case (5)
            if (r%real_masses) then
                call olo_an(an, 4, r%mr(1), r%mu)
            else
                call olo_an(an, 4, r%m(1), r%mu)
            end if
            first = an(0, 0)
        case (6)
            if (r%real_masses) then
                call olo_b11(b11, b00, b1, b0, r%p(1), r%mr(1), r%mr(2), r%mu)
            else
                pc = cmplx(r%p(1), 0.0_real64, real64)
                call olo_b11(b11, b00, b1, b0, pc, r%m(1), r%m(2), r%mu)
            end if
            first = b0(0)
        case (7)
            if (r%real_masses) then
                call olo_bn(bn, 3, r%p(1), r%mr(1), r%mr(2), r%mu)
            else
                pc = cmplx(r%p(1), 0.0_real64, real64)
                call olo_bn(bn, 3, pc, r%m(1), r%m(2), r%mu)
            end if
            first = bn(0, 0)
        end select
    end function evaluate

    subroutine workload(mode, kind, fixture_id, indices)
        character(*), intent(in) :: mode
        integer, intent(in) :: kind, fixture_id, indices(:)
        integer :: cursor, repetition
        integer(int64) :: iteration, start_tick, end_tick, block
        real(real64) :: wall_seconds
        complex(real64) :: checksum

        cursor = 1
        checksum = 0
        do iteration = 1, warmup
            checksum = checksum + evaluate(rows(indices(cursor)))
            cursor = cursor + 1
            if (cursor > size(indices)) cursor = 1
        end do
        block = calls
        do repetition = 0, repetitions
            cursor = 1
            checksum = 0
            call system_clock(start_tick)
            do iteration = 1, block
                checksum = checksum + evaluate(rows(indices(cursor)))
                cursor = cursor + 1
                if (cursor > size(indices)) cursor = 1
            end do
            call system_clock(end_tick)
            wall_seconds = real(end_tick - start_tick, real64)/real(clock_rate, real64)
            if (repetition == 0) then
                ! Calibration only: size the timed blocks to at least 25 ms.
                if (wall_seconds < 0.025_real64) block = &
                    int(real(block, real64)*0.025_real64/max(wall_seconds, 1e-6_real64), int64) + 1
                cycle
            end if
            write(*, '(a,1x,a,4(1x,i0),1x,i0,3(1x,es26.17))') 'SAMPLE', mode, kind, fixture_id, &
                repetition, block, size(indices), wall_seconds, real(checksum, real64), aimag(checksum)
        end do
    end subroutine workload
end program rows_bench_oracle
