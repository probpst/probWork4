!=======================================================================
! stop_criteria.F90 -- stop criteria and Reason_for_stop.md (probCFD
!   style): NaN anywhere, CFL above Stop_CFL_Max, |v| above Stop_Vmax;
!   every stop (criterion, input error, normal end) writes
!   Reason_for_stop.md with the reason, the location and the history
!   of the last steps.
!=======================================================================
module stop_criteria
  use global_vars
  use mpi_vars
  implicit none
  private
  public :: check_nan, check_limits, record_step, stop_report
  integer, parameter :: MXHIST = 40
  integer :: nhist = 0
  integer  :: h_step(MXHIST)
  real(dp) :: h_t(MXHIST), h_dt(MXHIST), h_cfl(MXHIST), h_ek(MXHIST), h_div(MXHIST), h_vmax(MXHIST)

contains

  subroutine record_step
    integer :: m
    if (nhist < MXHIST) then
       nhist = nhist + 1;  m = nhist
    else
       h_step(1:MXHIST-1) = h_step(2:MXHIST);  h_t(1:MXHIST-1) = h_t(2:MXHIST);  h_dt(1:MXHIST-1) = h_dt(2:MXHIST)
       h_cfl(1:MXHIST-1) = h_cfl(2:MXHIST);  h_ek(1:MXHIST-1) = h_ek(2:MXHIST);  h_div(1:MXHIST-1) = h_div(2:MXHIST)
       h_vmax(1:MXHIST-1) = h_vmax(2:MXHIST)
       m = MXHIST
    end if
    h_step(m) = istep;  h_t(m) = t*t0;  h_dt(m) = dt*t0;  h_cfl(m) = cfl_now;  h_ek(m) = e_kin;  h_div(m) = divb_max;  h_vmax(m) = vmax
  end subroutine record_step

  subroutine check_nan(where)
    character(len=*), intent(in) :: where
    integer :: i, j, k, q
    logical :: bad, gbad
    character(len=300) :: msg
    bad = .false.;  msg = ''
    outer: do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             do q = 1, NVAR
                if (u(q,i,j,k) /= u(q,i,j,k)) then
                   bad = .true.
                   write(msg, '(A,I0,A,3I6,A,3ES12.4)') 'NaN in variable ', q, ' at global cell ', ioff+i, joff+j, koff+k, &
                        '  x [m] = ', xc(:,i,j,k)*Lch
                   exit outer
                end if
             end do
          end do
       end do
    end do outer
    gbad = gany(bad)
    if (gbad) then
       if (bad) write(*,'(A)') ' STOP: '//trim(msg)//'  ('//trim(where)//')'
       flush(6)
       call stop_report('NaN detected '//trim(where), msg)
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
       call MPI_Finalize(ierr)
       stop
    end if
  end subroutine check_nan

  subroutine check_limits
    character(len=200) :: msg
    if (cfl_now > stop_cfl_max) then
       write(msg, '(A,ES12.4,A,ES12.4)') 'maximum CFL ', cfl_now, ' exceeds Stop_CFL_Max ', stop_cfl_max
       if (myrank == 0) write(*,'(A)') ' STOP: '//trim(msg)
       call stop_report('CFL limit', msg)
       call MPI_Barrier(MPI_COMM_WORLD, ierr);  call MPI_Finalize(ierr);  stop
    end if
    if (vmax > stop_vmax) then
       write(msg, '(A,ES12.4,A,ES12.4)') 'maximum |v|* ', vmax, ' exceeds Stop_Vmax ', stop_vmax
       if (myrank == 0) write(*,'(A)') ' STOP: '//trim(msg)
       call stop_report('velocity limit', msg)
       call MPI_Barrier(MPI_COMM_WORLD, ierr);  call MPI_Finalize(ierr);  stop
    end if
  end subroutine check_limits

  subroutine stop_report(reason, detail)
    character(len=*), intent(in) :: reason, detail
    integer :: iu, m
    if (myrank /= 0) return
    open(newunit=iu, file='Reason_for_stop.md', status='replace', action='write')
    write(iu,'(A)') '# Reason for stop'
    write(iu,'(A)') ''
    write(iu,'(A,A)') '**Case:** ', trim(casename)
    write(iu,'(A,A)') '**Reason:** ', trim(reason)
    if (len_trim(detail) > 0) write(iu,'(A,A)') '**Detail:** ', trim(detail)
    write(iu,'(A,I0,A,ES12.4,A)') '**Stopped at step** ', istep, ' (t = ', t*t0, ' s)'
    write(iu,'(A)') ''
    write(iu,'(A)') '## History of the last steps'
    write(iu,'(A)') ''
    write(iu,'(A)') '| step | t [s] | dt [s] | CFL | E_kin* | max |divB| | max |v|* |'
    write(iu,'(A)') '|---|---|---|---|---|---|---|'
    do m = 1, nhist
       write(iu,'(A,I0,A,ES11.4,A,ES11.4,A,F8.4,A,ES11.4,A,ES11.4,A,ES11.4,A)') '| ', h_step(m), ' | ', h_t(m), ' | ', h_dt(m), &
            ' | ', h_cfl(m), ' | ', h_ek(m), ' | ', h_div(m), ' | ', h_vmax(m), ' |'
    end do
    write(iu,'(A)') ''
    select case (trim(reason(1:3)))
    case ('NaN')
       write(iu,'(A)') '## Likely causes'
       write(iu,'(A)') '- negative pressure / density from a too large time step: lower CFL (0.2-0.3) or use Time_scheme = 3'
       write(iu,'(A)') '- an unresolved current sheet: raise the resistivity, use Cnvt_scheme = 4 (ENO) instead of 2 / 3'
       write(iu,'(A)') '- a boundary condition inconsistent with the equilibrium (check BC_* and the fixed state)'
    case ('CFL')
       write(iu,'(A)') '## Likely causes'
       write(iu,'(A)') '- Del_T fixed and too large for the fast magnetosonic speed: set Del_T = 0 (automatic) or lower it'
    case ('vel')
       write(iu,'(A)') '## Likely causes'
       write(iu,'(A)') '- a runaway instability or a floor hit: check the equilibrium (grid / Solov_ev parameters) and the floors'
    case default
       write(iu,'(A)') 'Normal end of the run.'
    end select
    close(iu)
  end subroutine stop_report

end module stop_criteria
