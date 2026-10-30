!=======================================================================
! mod_global.F90 -- global data of probMHD
!
!   probMHD: full (compressible, visco-resistive, single-fluid) MHD
!   solver on structured curvilinear grids, written in the style of
!   probCFD_부산 (grouped name list, PLOT3D in/out, MPI, stop criteria).
!
!   Conserved cell variables  U(1:8) = rho, rho u, rho v, rho w, E,
!                                      Bx, By, Bz        (Cartesian)
!   with E = p/(gam-1) + rho |v|^2/2 + |B|^2/2   (code units: mu0 = 1)
!   Face variables            bf_c  = magnetic flux through the face
!                                     (B . S_f), constrained transport
!
!   NON-DIMENSIONALISATION (inputs and outputs are dimensional, SI):
!     length   L0 = Char_L [m]
!     field    B0 = Char_B [T]
!     density  rho0 = Char_Rho [kg/m3]
!     velocity vA0 = B0 / sqrt(mu0 rho0)      time  t0 = L0 / vA0
!     pressure p0 = B0^2 / mu0                temperature T* = p*/rho*
!       (T [K] = m_i p0 / (2 k_B rho0) * T* for a quasi-neutral
!        hydrogenic plasma with T_i = T_e, ion mass Ion_mass_amu)
!     resistivity  eta* = eta / (mu0 L0 vA0) = 1/S (Lundquist number)
!     diffusivity  nu*, chi* = nu, chi / (L0 vA0)
!     power density  P* = P L0 / (p0 vA0)
!=======================================================================
module global_vars
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: NVAR = 8
  integer, parameter :: NG   = 4          ! ghost layers (3rd-order ENO on the ghost faces)
  integer :: ngi = NG, ngj = NG, ngk = NG  ! ghost layers per direction (1 in a degenerate direction)
  real(dp), parameter :: PI   = 3.141592653589793238462643_dp
  real(dp), parameter :: MU0  = 4.0e-7_dp*PI
  real(dp), parameter :: KB   = 1.380649e-23_dp
  real(dp), parameter :: AMU  = 1.66053906660e-27_dp
  real(dp), parameter :: QE   = 1.602176634e-19_dp

  !------------------ $Case_title ------------------------------------
  character(len=128) :: casename = 'Case1'
  integer  :: run_mode = 0                 ! 0 standard, 1 restart from the latest <case>_NNNNNNNN.rsto

  !------------------ $Control ---------------------------------------
  real(dp) :: Lch = 1.0_dp, Bch = 1.0_dp, Rhoch = 1.0_dp
  real(dp) :: Ion_mass_amu = 2.014_dp     ! deuterium
  real(dp) :: dt_dim = 0.0_dp             ! Del_T [s]; 0 = automatic from CFL
  real(dp) :: End_T  = 1.0_dp             ! [s]
  real(dp) :: cfl_max = 0.4_dp, vn_max = 0.4_dp
  real(dp) :: gam = 5.0_dp/3.0_dp
  integer  :: nsteps_max = 100000000
  logical  :: nondim_inputs = .false.     ! Units = 0 (SI) / 1 (code units: Char_* = 1 and all inputs non-dimensional)

  !------------------ $Grid ------------------------------------------
  character(len=256) :: gridfile = 'grid.xyz'
  integer  :: grid_format = 1             ! 1 PLOT3D ASCII (probCFD), 2 PLOT3D Fortran unformatted float64 (KSTAR file), 3 generated box
  integer  :: grid_multiblock = 0         ! 1: leading nblocks record (must be 1)
  integer  :: ncell_box(3) = [64, 64, 1]  ! generated box: cells
  real(dp) :: box_min(3) = 0.0_dp, box_max(3) = 1.0_dp
  integer  :: periodic(3) = 0             ! 1: periodic in i/j/k
  integer  :: periodic_dup(3) = 0         ! 1: the periodic closing node is repeated in the file
  real(dp) :: periodic_rot(3) = 0.0_dp    ! rotational periodicity about the z axis [deg] (toroidal wedge); 0 = translation
  integer  :: bc_type(2,3) = 2            ! (min/max, dir): 0 periodic, 1 conducting wall, 2 zero gradient, 3 fixed (initial values held), 4 axis hole (conducting core)
  integer  :: nproc_dir(3) = 0            ! MPI decomposition (0 = automatic)

  !------------------ $Initial_Condition -----------------------------
  integer  :: ic_type = 6                 ! 1 Brio-Wu, 2 Orszag-Tang, 3 CP Alfven wave, 4 cylindrical tearing, 5 tokamak (Solov'ev), 6 uniform
  real(dp) :: rho_init = 1.0_dp, p_init = 1.0_dp, v_init(3) = 0.0_dp, B_init(3) = 0.0_dp
  real(dp) :: ic_param(10) = 0.0_dp        ! case-dependent parameters (see equilibrium.F90)
  integer  :: ic_dir = 1                   ! 1-D tests: direction of variation
  real(dp) :: perturb_amp = 0.0_dp         ! velocity perturbation amplitude (code units)
  integer  :: perturb_m = 2, perturb_n = 1 ! poloidal / toroidal (axial) mode numbers of the perturbation
  integer  :: perturb_seed = 12345         ! random perturbation: seed (perturb_m = 0)

  !------------------ $Tokamak ---------------------------------------
  real(dp) :: tok_R0 = 1.8_dp, tok_a = 0.5_dp, tok_kappa = 1.8_dp, tok_delta = 0.6_dp
  real(dp) :: tok_B0 = 3.5_dp, tok_Asol = -0.155_dp, tok_psi0 = 8.318766676137974_dp   ! KSTAR grid of the delivery
  real(dp) :: tok_p_floor = 0.01_dp        ! p_floor / p_axis
  real(dp) :: tok_n_axis = 5.0e19_dp       ! ion density on axis [m-3]
  real(dp) :: tok_n_exp = 0.5_dp           ! n = n_axis (p/p_axis)^n_exp (floored)
  real(dp) :: tok_n_floor = 0.05_dp        ! n_floor / n_axis
  real(dp) :: tok_phi_extent = 360.0_dp    ! toroidal extent of the grid [deg]

  !------------------ $Numerical_Scheme ------------------------------
  integer  :: Cnvt_scheme = 4              ! 1 1st-order upwind, 2 2nd-order upwind, 3 QUICK, 4 2nd-order ENO (probCFD), 5 3rd-order ENO
  integer  :: Riemann_solver = 3           ! 1 Rusanov (local Lax-Friedrichs), 2 HLL, 3 HLLD
  integer  :: Time_scheme = 3              ! 1 Euler, 2 SSP RK2 (Heun), 3 SSP RK3
  integer  :: Recon_vars = 1               ! 1 primitive variables, 2 conserved variables
  integer  :: CT_scheme = 1                ! 1 Gardiner-Stone upwind CT, 2 Balsara-Spicer averaging, 0 cell-centred B only (no CT; div B grows)
  real(dp) :: rho_floor = 1.0e-6_dp, p_floor = 1.0e-8_dp
  integer  :: SW_Eq_Subtract = 0           ! 1: subtract the discrete residual of the initial state (well-balanced equilibrium)
  integer  :: SW_Energy = 1                ! 1 total energy (conservative), 2 internal-energy (pressure) equation, 3 dual (total where the thermal fraction > Dual_eps)
  real(dp) :: dual_eps = 0.01_dp           ! SW_Energy = 3: thermal fraction below which the internal-energy equation sets p

  !------------------ $Physics ---------------------------------------
  real(dp) :: eta_dim = 0.0_dp             ! resistivity [Ohm m] (or code units)
  integer  :: Spitzer = 0                  ! 1: eta = eta_ref (T/T_ref)^-3/2
  real(dp) :: T_ref_eV = 1000.0_dp         ! reference temperature of the Spitzer law [eV]
  real(dp) :: eta_max_factor = 100.0_dp    ! cap of the Spitzer resistivity: eta <= eta_ref * factor
  real(dp) :: nu_dim = 0.0_dp              ! kinematic viscosity [m2/s]
  real(dp) :: chi_par_dim = 0.0_dp, chi_perp_dim = 0.0_dp   ! heat diffusivities [m2/s] (q = -rho chi grad T)
  real(dp) :: chi_ratio_max = 1.0e6_dp     ! cap of chi_par / chi_perp
  real(dp) :: D_perp_dim = 0.0_dp          ! particle diffusivity [m2/s]
  real(dp) :: P_heat_dim = 0.0_dp          ! heating power density on axis [W/m3]
  real(dp) :: heat_width = 0.5_dp          ! Gaussian heating width / a (tokamak: distance from the axis in the poloidal plane)
  real(dp) :: heat_start = 0.0_dp          ! [s]
  integer  :: SW_radiation = 0             ! 1: bremsstrahlung P = C_br Z^2 n_e n_i sqrt(T_keV)
  real(dp) :: Z_eff = 1.0_dp
  integer  :: SW_Hall = 0                  ! reserved (not implemented in this version)
  real(dp) :: grav(3) = 0.0_dp             ! body force [m/s2] (tests)

  !------------------ $Turbulence ------------------------------------
  integer  :: SW_Turb_MHD = 0              ! 0 none, 1 MHD Smagorinsky (nu_sgs, eta_sgs), 2 dynamic (reserved), 3 anomalous transport
  real(dp) :: Cs_sgs = 0.1_dp, Ceta_sgs = 0.1_dp, Pr_sgs = 0.7_dp
  integer  :: Anom_model = 1               ! 1 constant chi_anom / D_anom / nu_anom, 2 Bohm: chi = Bohm_coef T/(16 e B)
  real(dp) :: chi_anom_dim = 0.0_dp, D_anom_dim = 0.0_dp, nu_anom_dim = 0.0_dp
  real(dp) :: Bohm_coef = 1.0_dp
  real(dp) :: anom_max_dim = 1.0e3_dp      ! cap of the Bohm diffusivity [m2/s]

  !------------------ $Export_Interval --------------------------------
  integer  :: out_interval = 100, log_interval = 10, monitor_interval = 1, restart_interval = 0
  integer  :: out_slice = 0                ! 1: also write the k = mid toroidal plane (2-D ASCII)

  !------------------ $Stop ------------------------------------------
  real(dp) :: stop_cfl_max = 1.0e4_dp, stop_vmax = 1.0e3_dp
  logical  :: debug_balance = .false.      ! Debug_balance = 1: print the energy balance of every RHS evaluation

  !------------------ derived scales ----------------------------------
  real(dp) :: vA0 = 1.0_dp, t0 = 1.0_dp, p0 = 1.0_dp, n0 = 1.0_dp, T0_K = 1.0_dp, mion = 1.0_dp
  real(dp) :: eta_c = 0.0_dp, nu_c = 0.0_dp, chi_par_c = 0.0_dp, chi_perp_c = 0.0_dp, D_perp_c = 0.0_dp
  real(dp) :: P_heat_c = 0.0_dp, chi_anom_c = 0.0_dp, D_anom_c = 0.0_dp, nu_anom_c = 0.0_dp, anom_max_c = 0.0_dp
  real(dp) :: T_ref_c = 1.0_dp, grav_c(3) = 0.0_dp, heat_start_c = 0.0_dp
  logical  :: diffusion_on = .false.

  !------------------ grid (local sub-domain) --------------------------
  integer  :: ngx, ngy, ngz                ! global cells
  integer  :: nx, ny, nz                   ! local cells
  integer  :: ioff, joff, koff             ! global index of local cell 1 is ioff+1
  ! node coordinates xn(3, 0-NG:nx+ngi, ...) : node (i,j,k) is the (minus,minus,minus) corner of cell (i,j,k)
  real(dp), allocatable :: xn(:,:,:,:)
  real(dp), allocatable :: vol(:,:,:)      ! cell volume
  real(dp), allocatable :: xc(:,:,:,:)     ! cell centre
  real(dp), allocatable :: sf(:,:,:,:,:)   ! sf(3, c, i,j,k): area vector of the minus face of cell (i,j,k) in direction c
  real(dp), allocatable :: xf(:,:,:,:,:)   ! face centres
  real(dp), allocatable :: le(:,:,:,:,:)   ! le(3, c, i,j,k): edge vector of edge family c (along xi_c) at the (minus,minus) corner
  real(dp), allocatable :: xe(:,:,:,:,:)   ! edge midpoints
  real(dp), allocatable :: hc(:,:,:,:)     ! hc(c,i,j,k): characteristic cell length in direction c
  real(dp), allocatable :: dvec(:,:,:,:,:) ! dvec(3,c,i,j,k): vector from cell (i-1) centre to cell i centre across the minus face
  real(dp) :: vol_total = 0.0_dp

  !------------------ fields ------------------------------------------
  real(dp), allocatable :: u(:,:,:,:)      ! conserved (NVAR, cells with ghosts)
  real(dp), allocatable :: w(:,:,:,:)      ! primitive  rho, vx, vy, vz, p, Bx, By, Bz
  real(dp), allocatable :: bf(:,:,:,:)     ! bf(c, i,j,k): magnetic flux through the minus face in direction c
  real(dp), allocatable :: u0(:,:,:,:), bf0(:,:,:,:)   ! RK stage storage (level n)
  real(dp), allocatable :: ueq(:,:,:,:), bfeq(:,:,:,:) ! initial / equilibrium state (fixed BCs, perturbation energy)
  real(dp), allocatable :: rhs(:,:,:,:)    ! dU/dt (NVAR)
  real(dp), allocatable :: emf(:,:,:,:)    ! emf(c, i,j,k): E . L_e on edge family c
  real(dp), allocatable :: ef(:,:,:,:,:)   ! ef(3, c, i,j,k): Cartesian E at face centres (Riemann)
  real(dp), allocatable :: gradw(:,:,:,:,:)! gradw(3, 1:5+3, cells): cell gradients of rho, v, p, B (Green-Gauss)
  real(dp), allocatable :: jc(:,:,:,:)     ! cell current density (3)
  real(dp), allocatable :: eta_cell(:,:,:), nu_cell(:,:,:), chi_par_cell(:,:,:), chi_perp_cell(:,:,:), D_cell(:,:,:)
  real(dp), allocatable :: nu_sgs(:,:,:), eta_sgs(:,:,:)
  real(dp), allocatable :: src_heat(:,:,:) ! heating power density profile (code units, without the time switch)
  real(dp), allocatable :: divb(:,:,:)
  real(dp), allocatable :: eint(:,:,:,:), eint0(:,:,:,:), rhs_e(:,:,:)   ! internal energy p/(gam-1) (SW_Energy >= 2), (1, cells)
  real(dp), allocatable :: einteq(:,:,:,:)

  !------------------ run state ---------------------------------------
  real(dp) :: t = 0.0_dp, dt = 0.0_dp
  integer  :: istep = 0
  real(dp) :: cfl_now = 0.0_dp
  real(dp) :: e_kin = 0.0_dp, e_mag = 0.0_dp, e_th = 0.0_dp, e_tot = 0.0_dp, mass_tot = 0.0_dp
  real(dp) :: vmax = 0.0_dp, divb_max = 0.0_dp, e_pert = 0.0_dp, e_kin0 = 0.0_dp, rho_min = 0.0_dp, p_min = 0.0_dp
  real(dp) :: p_axis_c = 0.0_dp, R_axis_c = 0.0_dp   ! tokamak: axis pressure and radius in code units

contains

  pure integer function c2(c)   ! cyclic next direction
    integer, intent(in) :: c
    c2 = mod(c, 3) + 1
  end function c2

  pure integer function c3(c)   ! cyclic second next direction
    integer, intent(in) :: c
    c3 = mod(c + 1, 3) + 1
  end function c3

  ! primitive from conserved (one cell)
  pure subroutine cons2prim(uc, wc)
    real(dp), intent(in)  :: uc(NVAR)
    real(dp), intent(out) :: wc(NVAR)
    real(dp) :: r, ek, eb
    r = max(uc(1), rho_floor)
    wc(1) = r
    wc(2:4) = uc(2:4)/r
    wc(6:8) = uc(6:8)
    ek = 0.5_dp*r*sum(wc(2:4)**2)
    eb = 0.5_dp*sum(uc(6:8)**2)
    wc(5) = max((gam - 1.0_dp)*(uc(5) - ek - eb), p_floor)
  end subroutine cons2prim

  pure subroutine prim2cons(wc, uc)
    real(dp), intent(in)  :: wc(NVAR)
    real(dp), intent(out) :: uc(NVAR)
    uc(1) = wc(1)
    uc(2:4) = wc(1)*wc(2:4)
    uc(6:8) = wc(6:8)
    uc(5) = wc(5)/(gam - 1.0_dp) + 0.5_dp*wc(1)*sum(wc(2:4)**2) + 0.5_dp*sum(wc(6:8)**2)
  end subroutine prim2cons

end module global_vars
