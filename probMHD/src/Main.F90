!=======================================================================
! Main.F90 -- probMHD: full (compressible, visco-resistive) MHD solver
!   on structured curvilinear grids, in the style of probCFD_부산.
!
!   Start-up : name list, units, grid, metrics, initial state / equilibrium
!   Loop     : time step (CFL / VN), SSP Runge-Kutta step with constrained
!              transport, diagnostics, monitor, log, output, restart,
!              stop criteria
!=======================================================================
program probMHD_main
  use global_vars
  use mpi_vars
  use input
  use geometry
  use equilibrium
  use boundary
  use convection
  use induction
  use diffusion
  use time_integration
  use output
  use stop_criteria
  implicit none
  integer :: step0, nsteps, k
  real(dp) :: t_end_c, wall, e_tot0, dt_last
  logical :: done

  call mpi_setup
  call print_banner
  call read_namelist
  call set_units
  if (SW_Turb_MHD == 2) call mpi_abort_msg('SW_Turb_MHD = 2 (dynamic Smagorinsky) is reserved for a later update; use 0, 1 or 3')
  if (SW_Hall == 1) call mpi_abort_msg('SW_Hall = 1 is reserved (not implemented in this version)')
  call read_grid
  xn = xn/Lch                     ! code units
  call build_metrics
  call print_input_summary
  call initialize_fields
  call refresh_state(u, bf, eint)
  step0 = 0
  if (run_mode == 1) then
     call read_restart(step0)
     call refresh_state(u, bf, eint)
  end if
  call write_grid_plot3d
  call monitor_open(run_mode == 1)

  ! initial diagnostics and output
  call compute_rhs(u, bf, eint, t)
  call divergence_B(bf, w, divb_max);  divb_max = gmax(divb_max)
  call global_diagnostics
  e_tot0 = e_tot;  e_kin0 = e_kin
  istep = step0
  call compute_dt
  if (myrank == 0) then
     write(*,'(A)') ' ---------------- initial state ----------------'
     write(*,'(A,4ES14.6)') '  E_kin, E_mag, E_th, E_tot* : ', e_kin, e_mag, e_th, e_tot
     write(*,'(A,ES12.4,A,ES12.4)') '  max |v|*  : ', vmax, '   max |div B| (normalised): ', divb_max
     write(*,'(A,ES12.4,A,ES12.4,A)') '  first dt   : ', dt*t0, ' s  (', dt, ' t_A)'
     write(*,'(A)') ' -----------------------------------------------'
     write(*,'(A)') '    step       t [s]        dt [s]     CFL     E_kin*        E_mag*        E_th*       dE_tot/E0    max|v|*    max|divB|   min rho*    wall [s]'
     flush(6)
  end if
  if (run_mode == 0) then
     call write_solution(0)
     call monitor_write
  end if

  t_end_c = End_T/t0
  done = .false.
  do while (.not. done)
     istep = istep + 1
     if (t + dt >= t_end_c) then
        dt = t_end_c - t;  done = .true.
     end if
     call rk_step
     t = t + dt
     call check_nan('after the time step')
     ! RHS of the new state (also the signal speeds for the next dt)
     call compute_rhs(u, bf, eint, t)
     dt_last = dt
     call compute_dt
     if (mod(istep, monitor_interval) == 0 .or. done) then
        call divergence_B(bf, w, divb_max);  divb_max = gmax(divb_max)
        call global_diagnostics
        call monitor_write
        call record_step
     end if
     if (mod(istep, log_interval) == 0 .or. done) then
        if (mod(istep, monitor_interval) /= 0 .and. .not. done) then
           call divergence_B(bf, w, divb_max);  divb_max = gmax(divb_max)
           call global_diagnostics
        end if
        wall = MPI_Wtime() - wall_t0
        if (myrank == 0) then
           write(*,'(I8,2ES13.4,F8.4,3ES14.6,ES13.3,3ES12.3,F10.1)') istep, t*t0, dt_last*t0, cfl_now, e_kin, e_mag, e_th, &
                (e_tot - e_tot0)/max(abs(e_tot0), 1.0e-30_dp), vmax, divb_max, rho_min, wall
           flush(6)
        end if
        call check_limits
     end if
     if (out_interval > 0) then
        if (mod(istep, out_interval) == 0 .or. done) call write_solution(istep)
     else if (done) then
        call write_solution(istep)
     end if
     if (restart_interval > 0) then
        if (mod(istep, restart_interval) == 0 .or. done) call write_restart(istep)
     end if
     if (istep - step0 >= nsteps_max) done = .true.
  end do

  wall = MPI_Wtime() - wall_t0
  if (myrank == 0) then
     write(*,'(A)') ' -----------------------------------------------'
     write(*,'(A,I0,A,ES12.4,A,F10.1,A)') '  normal end after ', istep - step0, ' steps, t = ', t*t0, ' s, wall time ', wall, ' s'
     write(*,'(A,ES12.4,A,ES12.4)') '  dE_tot/E0 = ', (e_tot - e_tot0)/max(abs(e_tot0), 1.0e-30_dp), '   max |div B| = ', divb_max
  end if
  call stop_report('normal end (End_T reached)', '')
  call MPI_Finalize(ierr)

contains

  subroutine print_banner
    if (myrank /= 0) return
    write(*,'(A)') ' ==============================================================================='
    write(*,'(A)') ' ||                                                                           ||'
    write(*,'(A)') ' ||     pppp   r  r   ooo     b       M   M  H   H  DDDD                      ||'
    write(*,'(A)') ' ||     p   p  r r   o   o    b       MM MM  H   H  D   D     probMHD         ||'
    write(*,'(A)') ' ||     p  p   rr   o     o   b       M M M  HHHHH  D    D    full MHD        ||'
    write(*,'(A)') ' ||     ppp    r    o     o   bbbb    M   M  H   H  D    D    curvilinear     ||'
    write(*,'(A)') ' ||     p      r    o     o   b   b   M   M  H   H  D   D     ENO / HLLD / CT ||'
    write(*,'(A)') ' ||     p      r      ooo     bbbb    M   M  H   H  DDDD      MPI             ||'
    write(*,'(A)') ' ||                                                                           ||'
    write(*,'(A)') ' ==============================================================================='
    write(*,'(A)') ' '
    write(*,'(A)') ' probMHD - tokamak full-MHD solver (probCFD_부산 family)'
    write(*,'(A)') ' '
  end subroutine print_banner

end program probMHD_main
