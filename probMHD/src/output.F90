!=======================================================================
! output.F90 -- PLOT3D grid / solution files (Fortran unformatted,
!   single block, cell centres), the .nam variable-name file, the
!   monitor file, the k-mid-plane ASCII slice and restart files.
!
!   Volume_solutions_grids/<case>_cells.grd         cell-centre grid [m]
!   Volume_solutions_grids/<case>_NNNNNNNN.p3d      function file, nvar variables
!   Volume_solutions_grids/<case>.nam               names and SI units
!   Output_data/<case>_monitor.dat                  global history
!   Output_data/<case>_slice_NNNNNNNN.dat           k = mid plane (Output_slice = 1)
!   <case>_NNNNNNNN.rsto                            restart (conserved cells + face fluxes)
!=======================================================================
module output
  use global_vars
  use mpi_vars
  implicit none
  private
  public :: write_grid_plot3d, write_solution, monitor_open, monitor_write, write_restart, read_restart, global_diagnostics
  integer, parameter :: NOUT = 17
  integer :: iu_mon = -1

contains

  subroutine make_dirs
    if (myrank == 0) then
       call execute_command_line('mkdir -p Volume_solutions_grids Output_data', wait=.true.)
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  end subroutine make_dirs

  subroutine write_grid_plot3d
    real(dp), allocatable :: g(:,:,:,:)
    integer :: iu
    call make_dirs
    if (myrank == 0) allocate(g(3, ngx, ngy, ngz))
    if (myrank /= 0) allocate(g(1,1,1,1))
    call gather_cells(xc, 3, g)
    if (myrank == 0) then
       g = g*Lch
       open(newunit=iu, file='Volume_solutions_grids/'//trim(casename)//'_cells.grd', form='unformatted', status='replace')
       write(iu) ngx, ngy, ngz
       write(iu) g(1,:,:,:), g(2,:,:,:), g(3,:,:,:)
       close(iu)
       open(newunit=iu, file='Volume_solutions_grids/'//trim(casename)//'.nam', status='replace')
       write(iu,'(A)') 'rho [kg/m3]';  write(iu,'(A)') 'u [m/s]';  write(iu,'(A)') 'v [m/s]';  write(iu,'(A)') 'w [m/s]'
       write(iu,'(A)') 'p [Pa]';  write(iu,'(A)') 'T [K]';  write(iu,'(A)') 'Bx [T]';  write(iu,'(A)') 'By [T]';  write(iu,'(A)') 'Bz [T]'
       write(iu,'(A)') 'Bmag [T]';  write(iu,'(A)') 'Jx [A/m2]';  write(iu,'(A)') 'Jy [A/m2]';  write(iu,'(A)') 'Jz [A/m2]'
       write(iu,'(A)') 'divB_norm [-]';  write(iu,'(A)') 'nu_t [m2/s]';  write(iu,'(A)') 'eta_t [Ohm m]';  write(iu,'(A)') 'chi_perp [m2/s]'
       close(iu)
    end if
    deallocate(g)
  end subroutine write_grid_plot3d

  subroutine write_solution(step)
    integer, intent(in) :: step
    real(dp), allocatable :: f(:,:,:,:), g(:,:,:,:)
    integer :: i, j, k, iu, q, km
    character(len=8) :: tag
    real(dp) :: wc(NVAR), jscale
    allocate(f(NOUT, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    f = 0.0_dp
    jscale = Bch/(MU0*Lch)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             call cons2prim(u(:,i,j,k), wc)
             f(1,i,j,k) = wc(1)*Rhoch
             f(2:4,i,j,k) = wc(2:4)*vA0
             f(5,i,j,k) = wc(5)*p0
             f(6,i,j,k) = wc(5)/wc(1)*T0_K
             f(7:9,i,j,k) = wc(6:8)*Bch
             f(10,i,j,k) = norm2(wc(6:8))*Bch
             f(11:13,i,j,k) = jc(:,i,j,k)*jscale
             f(14,i,j,k) = divb(i,j,k)
             f(15,i,j,k) = nu_sgs(i,j,k)*Lch*vA0
             f(16,i,j,k) = eta_sgs(i,j,k)*MU0*Lch*vA0
             f(17,i,j,k) = chi_perp_cell(i,j,k)*Lch*vA0
          end do
       end do
    end do
    if (myrank == 0) allocate(g(NOUT, ngx, ngy, ngz))
    if (myrank /= 0) allocate(g(1,1,1,1))
    call gather_cells(f, NOUT, g)
    if (myrank == 0) then
       write(tag, '(I8.8)') step
       open(newunit=iu, file='Volume_solutions_grids/'//trim(casename)//'_'//tag//'.p3d', form='unformatted', status='replace')
       write(iu) ngx, ngy, ngz, NOUT
       write(iu) (g(q,:,:,:), q = 1, NOUT)
       close(iu)
       if (out_slice == 1) then
          km = (ngz + 1)/2
          open(newunit=iu, file='Output_data/'//trim(casename)//'_slice_'//tag//'.dat', status='replace')
          write(iu,'(A,ES16.8,A,I0,A,I0)') '# t [s] = ', t*t0, '   ni = ', ngx, '   nj = ', ngy
          write(iu,'(A)') '# i j rho u v w p T Bx By Bz Jz divB'
          do j = 1, ngy
             do i = 1, ngx
                write(iu,'(2I6,12ES16.8)') i, j, g(1,i,j,km), g(2,i,j,km), g(3,i,j,km), g(4,i,j,km), g(5,i,j,km), g(6,i,j,km), &
                     g(7,i,j,km), g(8,i,j,km), g(9,i,j,km), g(13,i,j,km), g(14,i,j,km)
             end do
          end do
          close(iu)
       end if
    end if
    deallocate(f, g)
  end subroutine write_solution

  !--------------------------------------------------------------------
  subroutine global_diagnostics
    integer :: i, j, k
    real(dp) :: ek, em, et, ms, vm, ep, wc(NVAR), we(NVAR), rmn, pmn
    ek = 0.0_dp;  em = 0.0_dp;  et = 0.0_dp;  ms = 0.0_dp;  vm = 0.0_dp;  ep = 0.0_dp;  rmn = huge(1.0_dp);  pmn = huge(1.0_dp)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             call cons2prim(u(:,i,j,k), wc)
             call cons2prim(ueq(:,i,j,k), we)
             ek = ek + 0.5_dp*wc(1)*sum(wc(2:4)**2)*vol(i,j,k)
             em = em + 0.5_dp*sum(wc(6:8)**2)*vol(i,j,k)
             et = et + wc(5)/(gam - 1.0_dp)*vol(i,j,k)
             ms = ms + wc(1)*vol(i,j,k)
             vm = max(vm, norm2(wc(2:4)))
             ep = ep + 0.5_dp*sum((wc(6:8) - we(6:8))**2)*vol(i,j,k)
             rmn = min(rmn, wc(1));  pmn = min(pmn, wc(5))
          end do
       end do
    end do
    e_kin = gsum(ek);  e_mag = gsum(em);  e_th = gsum(et);  mass_tot = gsum(ms);  vmax = gmax(vm);  e_pert = gsum(ep)
    rho_min = gmin(rmn);  p_min = gmin(pmn)
    e_tot = e_kin + e_mag + e_th
  end subroutine global_diagnostics

  subroutine monitor_open(append)
    logical, intent(in) :: append
    call make_dirs
    if (myrank /= 0) return
    if (append) then
       open(newunit=iu_mon, file='Output_data/'//trim(casename)//'_monitor.dat', status='old', position='append', action='write')
    else
       open(newunit=iu_mon, file='Output_data/'//trim(casename)//'_monitor.dat', status='replace', action='write')
       write(iu_mon,'(A)') '# step  t[s]  dt[s]  CFL  E_kin*  E_mag*  E_th*  E_tot*  mass*  max|v|*  max|divB|  E_magpert*  t/tA  rho_min*  p_min*'
    end if
  end subroutine monitor_open

  subroutine monitor_write
    if (myrank /= 0) return
    write(iu_mon,'(I8,14ES17.9)') istep, t*t0, dt*t0, cfl_now, e_kin, e_mag, e_th, e_tot, mass_tot, vmax, divb_max, e_pert, t, rho_min, p_min
    flush(iu_mon)
  end subroutine monitor_write

  !--------------------------------------------------------------------
  subroutine write_restart(step)
    integer, intent(in) :: step
    real(dp), allocatable :: g(:,:,:,:), gb(:,:,:,:)
    integer :: iu
    character(len=8) :: tag
    if (myrank == 0) then
       allocate(g(NVAR, ngx, ngy, ngz), gb(3, ngx, ngy, ngz))
    else
       allocate(g(1,1,1,1), gb(1,1,1,1))
    end if
    call gather_cells(u, NVAR, g)
    call gather_cells(bf, 3, gb)
    if (myrank == 0) then
       write(tag, '(I8.8)') step
       open(newunit=iu, file=trim(casename)//'_'//tag//'.rsto', form='unformatted', status='replace')
       write(iu) step, t, ngx, ngy, ngz
       write(iu) g
       write(iu) gb
    end if
    deallocate(g);  if (myrank == 0) allocate(g(1, ngx, ngy, ngz));  if (myrank /= 0) allocate(g(1,1,1,1))
    call gather_cells(eint, 1, g)
    if (myrank == 0) then
       write(iu) g
       close(iu)
    end if
    deallocate(g, gb)
  end subroutine write_restart

  ! restart from the latest <case>_NNNNNNNN.rsto in the working folder
  subroutine read_restart(step0)
    integer, intent(out) :: step0
    real(dp), allocatable :: g(:,:,:,:), gb(:,:,:,:)
    integer :: iu, n1, n2, n3, io, s, best
    character(len=8) :: tag
    logical :: ex
    real(dp) :: tr
    best = -1
    if (myrank == 0) then
       do s = 0, 99999999, 1
          write(tag, '(I8.8)') s
          inquire(file=trim(casename)//'_'//tag//'.rsto', exist=ex)
          if (ex) best = s
          if (s > 1000000) exit
       end do
       ! faster scan: use ls through a shell
       call execute_command_line('ls -1 '//trim(casename)//'_????????.rsto 2>/dev/null | sort | tail -1 > .probmhd_last_rsto', wait=.true.)
       open(newunit=iu, file='.probmhd_last_rsto', status='old', action='read', iostat=io)
       if (io == 0) then
          read(iu, '(A)', iostat=io) tag
          close(iu)
       end if
       inquire(file='.probmhd_last_rsto', exist=ex)
    end if
    call MPI_Bcast(best, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    if (best < 0) call mpi_abort_msg('Mode = 1 but no restart file <case>_NNNNNNNN.rsto found')
    if (myrank == 0) then
       allocate(g(NVAR, ngx, ngy, ngz), gb(3, ngx, ngy, ngz))
       write(tag, '(I8.8)') best
       open(newunit=iu, file=trim(casename)//'_'//tag//'.rsto', form='unformatted', status='old', action='read')
       read(iu) step0, tr, n1, n2, n3
       if (n1 /= ngx .or. n2 /= ngy .or. n3 /= ngz) call mpi_abort_msg('restart file grid size differs from the grid')
       read(iu) g
       read(iu) gb
       write(*,'(A,I0,A,ES12.4,A)') '  restart from step ', step0, ' (t = ', tr*t0, ' s)'
    else
       allocate(g(1,1,1,1), gb(1,1,1,1))
    end if
    call MPI_Bcast(step0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(tr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    t = tr
    call scatter_cells(g, NVAR, u)
    call scatter_cells(gb, 3, bf)
    deallocate(g)
    if (myrank == 0) then
       allocate(g(1, ngx, ngy, ngz));  read(iu) g;  close(iu)
    else
       allocate(g(1,1,1,1))
    end if
    call scatter_cells(g, 1, eint)
    deallocate(g, gb)
  end subroutine read_restart

end module output
