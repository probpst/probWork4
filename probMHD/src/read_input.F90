!=======================================================================
! read_input.F90 -- probMHD.nml (grouped name list, probCFD format) and
!   the grid (PLOT3D ASCII / Fortran-unformatted / generated box).
!
!   $Group_name
!    Variable = value(s)        ! comment
!   /end
!
!   Variable names are case-insensitive; unknown names abort.  Every
!   variable has a default (mod_global.F90); a group may be omitted.
!=======================================================================
module input
  use global_vars
  use mpi_vars
  implicit none
  private
  public :: read_namelist, read_grid, print_input_summary, lower

contains

  function lower(s) result(r)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: r
    integer :: i, ic
    r = s
    do i = 1, len(s)
       ic = iachar(s(i:i))
       if (ic >= 65 .and. ic <= 90) r(i:i) = achar(ic + 32)
    end do
  end function lower

  subroutine read_namelist
    integer :: iu, ios, lnum, pe, nv, i, nchar
    character(len=512) :: line, work, name, group
    character(len=256) :: vals(32)
    character(len=256), allocatable :: lines(:)
    logical :: ex
    integer :: units = 0

    if (myrank == 0) then
       inquire(file='probMHD.nml', exist=ex)
       if (.not. ex) call mpi_abort_msg('probMHD.nml not found in the working folder')
       open(newunit=iu, file='probMHD.nml', status='old', action='read')
       nchar = 0
       do
          read(iu, '(A)', iostat=ios) line
          if (ios /= 0) exit
          nchar = nchar + 1
       end do
       close(iu)
    end if
    call MPI_Bcast(nchar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    allocate(lines(nchar))
    if (myrank == 0) then
       open(newunit=iu, file='probMHD.nml', status='old', action='read')
       do i = 1, nchar
          read(iu, '(A)') lines(i)
       end do
       close(iu)
    end if
    call MPI_Bcast(lines, 256*nchar, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

    group = ''
    do lnum = 1, nchar
       line = lines(lnum)
       pe = index(line, '!')
       if (pe > 0) line = line(1:pe-1)
       work = adjustl(line)
       if (len_trim(work) == 0) cycle
       if (work(1:1) == '$') then
          group = lower(trim(work(2:)));  cycle
       end if
       if (lower(work(1:4)) == '/end') then
          group = '';  cycle
       end if
       pe = index(work, '=')
       if (pe == 0) call mpi_abort_msg('probMHD.nml: cannot parse line "'//trim(line)//'"')
       name = lower(trim(adjustl(work(1:pe-1))))
       call split_values(work(pe+1:), vals, nv)
       call assign(name, vals, nv, units)
    end do
    nondim_inputs = (units == 1)
  end subroutine read_namelist

  subroutine split_values(s, vals, nv)
    character(len=*), intent(in)  :: s
    character(len=256), intent(out) :: vals(32)
    integer, intent(out) :: nv
    integer :: i, n
    character(len=512) :: w
    logical :: inq
    w = s;  nv = 0;  n = len_trim(w)
    ! commas to blanks (outside quotes)
    inq = .false.
    do i = 1, n
       if (w(i:i) == '"' .or. w(i:i) == "'") inq = .not. inq
       if (.not. inq .and. w(i:i) == ',') w(i:i) = ' '
    end do
    i = 1
    do while (i <= n)
       if (w(i:i) == ' ') then
          i = i + 1;  cycle
       end if
       nv = nv + 1
       if (w(i:i) == '"' .or. w(i:i) == "'") then
          vals(nv) = ''
          i = i + 1
          do while (i <= n .and. w(i:i) /= '"' .and. w(i:i) /= "'")
             vals(nv) = trim(vals(nv))//w(i:i);  i = i + 1
          end do
          i = i + 1
       else
          vals(nv) = ''
          do while (i <= n .and. w(i:i) /= ' ')
             vals(nv) = trim(vals(nv))//w(i:i);  i = i + 1
          end do
       end if
    end do
  end subroutine split_values

  subroutine assign(name, vals, nv, units)
    character(len=*), intent(in) :: name
    character(len=256), intent(in) :: vals(32)
    integer, intent(in) :: nv
    integer, intent(inout) :: units
    integer :: i
    select case (trim(name))
    ! $Case_title
    case ('mode');              read(vals(1), *) run_mode
    case ('case_name');         casename = vals(1)
    ! $Control
    case ('units');             read(vals(1), *) units
    case ('char_l');            read(vals(1), *) Lch
    case ('char_b');            read(vals(1), *) Bch
    case ('char_rho');          read(vals(1), *) Rhoch
    case ('ion_mass_amu');      read(vals(1), *) Ion_mass_amu
    case ('del_t');             read(vals(1), *) dt_dim
    case ('end_t');             read(vals(1), *) End_T
    case ('cfl');               read(vals(1), *) cfl_max
    case ('vn');                read(vals(1), *) vn_max
    case ('gamma');             read(vals(1), *) gam
    case ('max_steps');         read(vals(1), *) nsteps_max
    ! $Grid
    case ('grid_file');         gridfile = vals(1)
    case ('grid_format');       read(vals(1), *) grid_format
    case ('grid_multiblock');   read(vals(1), *) grid_multiblock
    case ('ncell');             read(vals(1:3), *) ncell_box
    case ('box_min');           read(vals(1:3), *) box_min
    case ('box_max');           read(vals(1:3), *) box_max
    case ('periodic');          read(vals(1:3), *) periodic
    case ('periodic_dup');      read(vals(1:3), *) periodic_dup
    case ('periodic_rot');      read(vals(1:3), *) periodic_rot
    case ('bc_imin');           read(vals(1), *) bc_type(1,1)
    case ('bc_imax');           read(vals(1), *) bc_type(2,1)
    case ('bc_jmin');           read(vals(1), *) bc_type(1,2)
    case ('bc_jmax');           read(vals(1), *) bc_type(2,2)
    case ('bc_kmin');           read(vals(1), *) bc_type(1,3)
    case ('bc_kmax');           read(vals(1), *) bc_type(2,3)
    case ('nproc_dir');         read(vals(1:3), *) nproc_dir
    ! $Initial_Condition
    case ('ic_type');           read(vals(1), *) ic_type
    case ('rho_init');          read(vals(1), *) rho_init
    case ('p_init');            read(vals(1), *) p_init
    case ('v_init');            read(vals(1:3), *) v_init
    case ('b_init');            read(vals(1:3), *) B_init
    case ('ic_param');          do i = 1, min(nv, 10); read(vals(i), *) ic_param(i); end do
    case ('ic_dir');            read(vals(1), *) ic_dir
    case ('perturb_amp');       read(vals(1), *) perturb_amp
    case ('perturb_m');         read(vals(1), *) perturb_m
    case ('perturb_n');         read(vals(1), *) perturb_n
    case ('perturb_seed');      read(vals(1), *) perturb_seed
    ! $Tokamak
    case ('r0');                read(vals(1), *) tok_R0
    case ('a');                 read(vals(1), *) tok_a
    case ('kappa');             read(vals(1), *) tok_kappa
    case ('delta');             read(vals(1), *) tok_delta
    case ('b0');                read(vals(1), *) tok_B0
    case ('a_solovev');         read(vals(1), *) tok_Asol
    case ('psi0');              read(vals(1), *) tok_psi0
    case ('p_floor_ratio');     read(vals(1), *) tok_p_floor
    case ('n_axis');            read(vals(1), *) tok_n_axis
    case ('n_exponent');        read(vals(1), *) tok_n_exp
    case ('n_floor_ratio');     read(vals(1), *) tok_n_floor
    case ('phi_extent');        read(vals(1), *) tok_phi_extent
    ! $Numerical_Scheme
    case ('cnvt_scheme');       read(vals(1), *) Cnvt_scheme
    case ('riemann_solver');    read(vals(1), *) Riemann_solver
    case ('time_scheme');       read(vals(1), *) Time_scheme
    case ('recon_vars');        read(vals(1), *) Recon_vars
    case ('ct_scheme');         read(vals(1), *) CT_scheme
    case ('rho_floor');         read(vals(1), *) rho_floor
    case ('p_floor');           read(vals(1), *) p_floor
    case ('sw_eq_subtract');    read(vals(1), *) SW_Eq_Subtract
    case ('sw_energy');         read(vals(1), *) SW_Energy
    case ('dual_eps');          read(vals(1), *) dual_eps
    ! $Physics
    case ('eta');               read(vals(1), *) eta_dim
    case ('spitzer');           read(vals(1), *) Spitzer
    case ('t_ref_ev');          read(vals(1), *) T_ref_eV
    case ('eta_max_factor');    read(vals(1), *) eta_max_factor
    case ('nu');                read(vals(1), *) nu_dim
    case ('chi_par');           read(vals(1), *) chi_par_dim
    case ('chi_perp');          read(vals(1), *) chi_perp_dim
    case ('chi_ratio_max');     read(vals(1), *) chi_ratio_max
    case ('d_perp');            read(vals(1), *) D_perp_dim
    case ('p_heat');            read(vals(1), *) P_heat_dim
    case ('heat_width');        read(vals(1), *) heat_width
    case ('heat_start');        read(vals(1), *) heat_start
    case ('sw_radiation');      read(vals(1), *) SW_radiation
    case ('z_eff');             read(vals(1), *) Z_eff
    case ('sw_hall');           read(vals(1), *) SW_Hall
    case ('gravity');           read(vals(1:3), *) grav
    ! $Turbulence
    case ('sw_turb_mhd');       read(vals(1), *) SW_Turb_MHD
    case ('cs');                read(vals(1), *) Cs_sgs
    case ('c_eta');             read(vals(1), *) Ceta_sgs
    case ('pr_sgs');            read(vals(1), *) Pr_sgs
    case ('anom_model');        read(vals(1), *) Anom_model
    case ('chi_anom');          read(vals(1), *) chi_anom_dim
    case ('d_anom');            read(vals(1), *) D_anom_dim
    case ('nu_anom');           read(vals(1), *) nu_anom_dim
    case ('bohm_coef');         read(vals(1), *) Bohm_coef
    case ('anom_max');          read(vals(1), *) anom_max_dim
    ! $Export_Interval
    case ('output_intvl');      read(vals(1), *) out_interval
    case ('log_intvl');         read(vals(1), *) log_interval
    case ('monitor_intvl');     read(vals(1), *) monitor_interval
    case ('restart_intvl');     read(vals(1), *) restart_interval
    case ('output_slice');      read(vals(1), *) out_slice
    ! $Stop
    case ('stop_cfl_max');      read(vals(1), *) stop_cfl_max
    case ('stop_vmax');         read(vals(1), *) stop_vmax
    case ('debug_balance');     read(vals(1), *) i;  debug_balance = (i == 1)
    case default
       call mpi_abort_msg('probMHD.nml: unknown variable "'//trim(name)//'"')
    end select
  end subroutine assign

  !--------------------------------------------------------------------
  ! Grid: rank 0 reads (or generates) the global node array, broadcasts
  ! it, every rank extracts its sub-domain with ghost nodes.
  !   Periodic directions: the file holds ng nodes without the closing
  !   node (Periodic_dup = 0, KSTAR file) or ng+1 nodes with the closing
  !   node repeated (Periodic_dup = 1, probCFD grids).  The ghost images
  !   are the wrapped nodes plus the translation (closing - first node),
  !   which is zero for a physically closed grid (torus / cylinder).
  !--------------------------------------------------------------------
  subroutine read_grid
    integer :: nn(3), iu, nb, i, j, k, c, ii, jj, kk, io
    real(dp), allocatable :: g(:,:,:,:)
    real(dp) :: shift(3)
    integer :: ibuf(3), gi(3), ng3(3), im(3), nwrap, ngc, ngl(3)
    character(len=64) :: form

    if (myrank == 0) then
       select case (grid_format)
       case (1)
          open(newunit=iu, file=trim(gridfile), status='old', action='read', iostat=io)
          if (io /= 0) call mpi_abort_msg('grid file not found: '//trim(gridfile))
          if (grid_multiblock == 1) then
             read(iu, *) nb
             if (nb /= 1) call mpi_abort_msg('probMHD handles a single grid block')
          end if
          read(iu, *) nn
          allocate(g(3, nn(1), nn(2), nn(3)))
          do c = 1, 3
             read(iu, *) (((g(c,i,j,k), i=1,nn(1)), j=1,nn(2)), k=1,nn(3))
          end do
          close(iu)
       case (2)
          open(newunit=iu, file=trim(gridfile), status='old', action='read', form='unformatted', access='sequential', iostat=io)
          if (io /= 0) call mpi_abort_msg('grid file not found: '//trim(gridfile))
          if (grid_multiblock == 1) then
             read(iu) nb
             if (nb /= 1) call mpi_abort_msg('probMHD handles a single grid block')
          end if
          read(iu) nn
          allocate(g(3, nn(1), nn(2), nn(3)))
          read(iu) g(1,:,:,:), g(2,:,:,:), g(3,:,:,:)
          close(iu)
       case (3)
          nn = ncell_box + 1
          allocate(g(3, nn(1), nn(2), nn(3)))
          do k = 1, nn(3); do j = 1, nn(2); do i = 1, nn(1)
             g(1,i,j,k) = box_min(1) + (box_max(1) - box_min(1))*real(i-1,dp)/real(ncell_box(1),dp)
             g(2,i,j,k) = box_min(2) + (box_max(2) - box_min(2))*real(j-1,dp)/real(ncell_box(2),dp)
             g(3,i,j,k) = box_min(3) + (box_max(3) - box_min(3))*real(k-1,dp)/real(ncell_box(3),dp)
          end do; end do; end do
          periodic_dup = 1
       case default
          call mpi_abort_msg('Grid_format must be 1, 2 or 3')
       end select
       ibuf = nn
    end if
    call MPI_Bcast(ibuf, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(periodic_dup, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    nn = ibuf
    if (myrank /= 0) allocate(g(3, nn(1), nn(2), nn(3)))
    call MPI_Bcast(g, 3*product(nn), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! global cell counts
    do c = 1, 3
       if (periodic(c) == 1 .and. periodic_dup(c) == 0) then
          ng3(c) = nn(c)
       else
          ng3(c) = nn(c) - 1
       end if
       if (ng3(c) < 1) call mpi_abort_msg('grid has fewer than 2 node planes in a non-periodic direction')
    end do
    ngx = ng3(1);  ngy = ng3(2);  ngz = ng3(3)

    shift = 0.0_dp
    call mpi_decompose
    allocate(xn(3, 1-ngi:nx+ngi+1, 1-ngj:ny+ngj+1, 1-ngk:nz+ngk+1))
    do k = 1-ngk, nz+ngk+1
       do j = 1-ngj, ny+ngj+1
          do i = 1-ngi, nx+ngi+1
             gi = [ioff + i, joff + j, koff + k]
             call node_image(gi, nn, ng3, im, shift)
             xn(:, i, j, k) = g(:, im(1), im(2), im(3))
             ! translation / rotation of periodic images
             do c = 1, 3
                if (periodic(c) == 1 .and. periodic_dup(c) == 1) then
                   nwrap = 0
                   if (gi(c) > ng3(c)) nwrap =  (gi(c) - 1)/ng3(c)
                   if (gi(c) < 1)      nwrap = -((ng3(c) - gi(c))/ng3(c))
                   if (nwrap /= 0) then
                      if (periodic_rot(c) /= 0.0_dp) then
                         xn(:, i, j, k) = rotz(xn(:, i, j, k), real(nwrap, dp)*periodic_rot(c)*PI/180.0_dp)
                      else
                         xn(:, i, j, k) = xn(:, i, j, k) + shift_vec(g, nn, c)*real(nwrap, dp)
                      end if
                   end if
                end if
             end do
          end do
       end do
    end do
    ! non-periodic boundaries: ghost nodes by linear extrapolation of the two boundary planes
    ngl = [ngi, ngj, ngk]
    do c = 1, 3
       if (periodic(c) == 1) cycle
       ngc = ngl(c)
       if (coords(c) == 0) then
          do ii = 0, 1 - ngc, -1
             select case (c)
             case (1); xn(:, ii, :, :) = 2.0_dp*xn(:, ii+1, :, :) - xn(:, ii+2, :, :)
             case (2); xn(:, :, ii, :) = 2.0_dp*xn(:, :, ii+1, :) - xn(:, :, ii+2, :)
             case (3); xn(:, :, :, ii) = 2.0_dp*xn(:, :, :, ii+1) - xn(:, :, :, ii+2)
             end select
          end do
       end if
       if (coords(c) == dims(c) - 1) then
          select case (c)
          case (1); do ii = nx+2, nx+ngi+1; xn(:, ii, :, :) = 2.0_dp*xn(:, ii-1, :, :) - xn(:, ii-2, :, :); end do
          case (2); do ii = ny+2, ny+ngj+1; xn(:, :, ii, :) = 2.0_dp*xn(:, :, ii-1, :) - xn(:, :, ii-2, :); end do
          case (3); do ii = nz+2, nz+ngk+1; xn(:, :, :, ii) = 2.0_dp*xn(:, :, :, ii-1) - xn(:, :, :, ii-2); end do
          end select
       end if
    end do
    deallocate(g)
  end subroutine read_grid

  pure function rotz(x, ang) result(y)
    real(dp), intent(in) :: x(3), ang
    real(dp) :: y(3)
    y(1) = cos(ang)*x(1) - sin(ang)*x(2)
    y(2) = sin(ang)*x(1) + cos(ang)*x(2)
    y(3) = x(3)
  end function rotz

  function shift_vec(g, nn, c) result(s)
    real(dp), intent(in) :: g(:,:,:,:)
    integer,  intent(in) :: nn(3), c
    real(dp) :: s(3)
    select case (c)
    case (1); s = g(:, nn(1), 1, 1) - g(:, 1, 1, 1)
    case (2); s = g(:, 1, nn(2), 1) - g(:, 1, 1, 1)
    case default; s = g(:, 1, 1, nn(3)) - g(:, 1, 1, 1)
    end select
  end function shift_vec

  ! global node index gi (may lie outside 1..ng+1) -> node index in the file array
  subroutine node_image(gi, nn, ng3, im, shift)
    integer,  intent(in)  :: gi(3), nn(3), ng3(3)
    integer,  intent(out) :: im(3)
    real(dp), intent(in)  :: shift(3)
    integer :: c
    do c = 1, 3
       if (periodic(c) == 1) then
          im(c) = modulo(gi(c) - 1, ng3(c)) + 1      ! 1..ng (the closing node = node 1)
       else
          im(c) = min(max(gi(c), 1), nn(c))          ! clamped; overwritten by extrapolation
       end if
    end do
  end subroutine node_image

  !--------------------------------------------------------------------
  subroutine print_input_summary
    character(len=64) :: s
    if (myrank /= 0) return
    write(*,'(A)') ' ---------------- input summary (probMHD.nml) ----------------'
    write(*,'(A,A)')          '  Case_name        : ', trim(casename)
    write(*,'(A,I0,A)')       '  Mode             : ', run_mode, merge(' (restart) ', ' (standard)', run_mode == 1)
    if (nondim_inputs) then
       write(*,'(A)')         '  Units            : 1 (code units; Char_* = 1)'
    else
       write(*,'(A)')         '  Units            : 0 (SI inputs / outputs)'
    end if
    write(*,'(A,ES12.4)')     '  Char_L   [m]     : ', Lch
    write(*,'(A,ES12.4)')     '  Char_B   [T]     : ', Bch
    write(*,'(A,ES12.4)')     '  Char_Rho [kg/m3] : ', Rhoch
    write(*,'(A,ES12.4)')     '  v_A0     [m/s]   : ', vA0
    write(*,'(A,ES12.4)')     '  t_A0     [s]     : ', t0
    write(*,'(A,ES12.4)')     '  p0       [Pa]    : ', p0
    write(*,'(A,ES12.4)')     '  End_T    [s]     : ', End_T
    if (dt_dim > 0.0_dp) then
       write(*,'(A,ES12.4)')  '  Del_T    [s]     : ', dt_dim
    else
       write(*,'(A,F6.3,A,F6.3)') '  Del_T            : automatic, CFL ', cfl_max, ', VN ', vn_max
    end if
    write(*,'(A,F8.4)')       '  gamma            : ', gam
    write(*,'(A,3I6,A,I0)')   '  cells            : ', ngx, ngy, ngz, '  total ', ngx*ngy*ngz
    write(*,'(A,3I6)')        '  MPI decomposition: ', dims
    write(*,'(A,3I2)')        '  periodic (i,j,k) : ', periodic
    if (any(periodic_rot /= 0.0_dp)) write(*,'(A,3F8.3)') '  rotational period [deg]: ', periodic_rot
    write(*,'(A,6I3)')        '  BC (imin,imax,jmin,jmax,kmin,kmax): ', bc_type(1,1), bc_type(2,1), bc_type(1,2), bc_type(2,2), bc_type(1,3), bc_type(2,3)
    select case (Cnvt_scheme)
    case (1); s = '1 (1st-order upwind)'
    case (2); s = '2 (2nd-order upwind)'
    case (3); s = '3 (QUICK)'
    case (4); s = '4 (2nd-order ENO, minmod - probCFD)'
    case (5); s = '5 (3rd-order ENO)'
    case default; s = '?'
    end select
    write(*,'(A,A)')          '  Cnvt_scheme      : ', trim(s)
    select case (Riemann_solver)
    case (1); s = '1 (Rusanov)'
    case (2); s = '2 (HLL)'
    case default; s = '3 (HLLD)'
    end select
    write(*,'(A,A)')          '  Riemann_solver   : ', trim(s)
    select case (Time_scheme)
    case (1); s = '1 (Euler)'
    case (2); s = '2 (SSP RK2)'
    case default; s = '3 (SSP RK3)'
    end select
    write(*,'(A,A)')          '  Time_scheme      : ', trim(s)
    write(*,'(A,I0,A)')       '  CT_scheme        : ', CT_scheme, merge(' (constrained transport, face fluxes)', ' (cell-centred B, no CT)             ', CT_scheme >= 1)
    write(*,'(A,I0)')         '  SW_Eq_Subtract   : ', SW_Eq_Subtract
    select case (SW_Energy)
    case (2); s = '2 (internal-energy equation: p from e_int, low-beta safe)'
    case (3); s = '3 (dual: total energy where the thermal fraction exceeds Dual_eps)'
    case default; s = '1 (total energy, conservative)'
    end select
    write(*,'(A,A)')          '  SW_Energy        : ', trim(s)
    if (eta_c > 0.0_dp) then
       write(*,'(A,ES12.4,A,ES12.4,A)') '  resistivity      : eta* = ', eta_c, '  (Lundquist number S = ', 1.0_dp/eta_c, ')'
    else
       write(*,'(A)') '  resistivity      : 0 (ideal MHD)'
    end if
    if (Spitzer == 1) write(*,'(A,ES10.3,A)') '                     Spitzer T^-3/2, T_ref = ', T_ref_eV, ' eV'
    write(*,'(A,ES12.4)')     '  viscosity nu*    : ', nu_c
    write(*,'(A,2ES12.4)')    '  chi_par*, chi_perp*: ', chi_par_c, chi_perp_c
    write(*,'(A,ES12.4)')     '  D_perp*          : ', D_perp_c
    write(*,'(A,ES12.4,A,ES10.3,A)') '  heating P*       : ', P_heat_c, '  (width ', heat_width, ' a)'
    write(*,'(A,I0)')         '  SW_radiation     : ', SW_radiation
    select case (SW_Turb_MHD)
    case (0); s = '0 (none: ENO dissipation only / implicit LES)'
    case (1); s = '1 (MHD Smagorinsky: nu_sgs, eta_sgs)'
    case (3); s = '3 (anomalous perpendicular transport)'
    case default; s = '?'
    end select
    write(*,'(A,A)')          '  SW_Turb_MHD      : ', trim(s)
    if (SW_Turb_MHD == 1) write(*,'(A,3F8.3)') '                     Cs, C_eta, Pr_sgs = ', Cs_sgs, Ceta_sgs, Pr_sgs
    if (SW_Turb_MHD == 3) then
       if (Anom_model == 1) then
          write(*,'(A,3ES11.3)') '                     constant chi_anom*, D_anom*, nu_anom* = ', chi_anom_c, D_anom_c, nu_anom_c
       else
          write(*,'(A,F7.3,A,ES10.3)') '                     Bohm, coefficient ', Bohm_coef, ', cap* ', anom_max_c
       end if
    end if
    write(*,'(A,I0)')         '  IC_type          : ', ic_type
    write(*,'(A,I0,A,I0,A,I0)') '  Output_Intvl     : ', out_interval, '  Log_Intvl ', log_interval, '  Restart_Intvl ', restart_interval
    write(*,'(A)') ' -------------------------------------------------------------'
    flush(6)
  end subroutine print_input_summary

end module input
