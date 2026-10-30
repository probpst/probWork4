!=======================================================================
! equilibrium.F90 -- initial conditions / equilibria
!
!   IC_type 1  Brio-Wu shock tube (gamma = 2)            [code units]
!           2  Orszag-Tang vortex (periodic unit square)  [code units]
!           3  circularly polarised Alfven wave (Toth / Gardiner-Stone)
!           4  periodic cylinder, cylindrical tearing-mode equilibrium
!           5  tokamak: Cerfon-Freidberg analytic Solov'ev equilibrium
!              (the equilibrium of the KSTAR flux-aligned grid)   [SI]
!           6  uniform state (rho_init, p_init, v_init, B_init)
!
!   The magnetic field is initialised on the FACES (constrained
!   transport): wherever a vector potential is known the face flux is
!   the circulation of A around the face (3-point Gauss quadrature per
!   edge), which is divergence-free to round-off; otherwise B(face
!   centre) . S_f.  The cell-centred B is then averaged from the faces.
!=======================================================================
module equilibrium
  use global_vars
  use mpi_vars
  use geometry, only: cross
  implicit none
  private
  public :: initialize_fields, set_units, cell_B_from_faces, solovev_setup, solovev_psiN, solovev_axis_R, solovev_coeff

  ! Solov'ev data (normalised coordinates x = R/R0, y = Z/R0)
  real(dp) :: sv_c(7), sv_eps, sv_psi_axis, sv_x_axis
  real(dp) :: tear_r_s = 0.0_dp
  integer, parameter :: NT = 4001
  real(dp) :: tab_r(NT), tab_bth(NT), tab_bz(NT), tab_ath(NT), tab_az(NT), tab_q(NT)
  real(dp) :: tear_p = 0.0_dp

contains

  !--------------------------------------------------------------------
  subroutine set_units
    real(dp) :: fac
    if (nondim_inputs) then
       Lch = 1.0_dp;  Bch = 1.0_dp;  Rhoch = 1.0_dp
       vA0 = 1.0_dp;  t0 = 1.0_dp;  p0 = 1.0_dp;  n0 = 1.0_dp
       eta_c = eta_dim;  nu_c = nu_dim;  chi_par_c = chi_par_dim;  chi_perp_c = chi_perp_dim;  D_perp_c = D_perp_dim
       P_heat_c = P_heat_dim;  chi_anom_c = chi_anom_dim;  D_anom_c = D_anom_dim;  nu_anom_c = nu_anom_dim
       anom_max_c = anom_max_dim;  grav_c = grav;  T_ref_c = T_ref_eV;  heat_start_c = heat_start
       mion = Ion_mass_amu*AMU;  T0_K = 1.0_dp
    else
       mion = Ion_mass_amu*AMU
       vA0 = Bch/sqrt(MU0*Rhoch)
       t0  = Lch/vA0
       p0  = Bch**2/MU0
       n0  = Rhoch/mion
       T0_K = mion*p0/(2.0_dp*KB*Rhoch)          ! T [K] = T0_K * (p*/rho*)
       eta_c = eta_dim/(MU0*Lch*vA0)
       fac = 1.0_dp/(Lch*vA0)
       nu_c = nu_dim*fac;  chi_par_c = chi_par_dim*fac;  chi_perp_c = chi_perp_dim*fac;  D_perp_c = D_perp_dim*fac
       chi_anom_c = chi_anom_dim*fac;  D_anom_c = D_anom_dim*fac;  nu_anom_c = nu_anom_dim*fac;  anom_max_c = anom_max_dim*fac
       P_heat_c = P_heat_dim*Lch/(p0*vA0)
       grav_c = grav*Lch/vA0**2
       T_ref_c = T_ref_eV*QE/KB/T0_K
       heat_start_c = heat_start/t0
    end if
    diffusion_on = (eta_c > 0.0_dp .or. nu_c > 0.0_dp .or. chi_par_c > 0.0_dp .or. chi_perp_c > 0.0_dp .or. &
                    D_perp_c > 0.0_dp .or. SW_Turb_MHD /= 0 .or. Spitzer == 1)
  end subroutine set_units

  !--------------------------------------------------------------------
  subroutine initialize_fields
    integer :: i, j, k
    real(dp) :: wc(NVAR), x(3)
    logical :: useA

    allocate(u(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), w(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(bf(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(u0(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), bf0(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(ueq(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), bfeq(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(rhs(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(emf(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), ef(3, 3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(gradw(3, 8, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), jc(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(eta_cell(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), nu_cell(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(chi_par_cell(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), chi_perp_cell(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(D_cell(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), nu_sgs(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), eta_sgs(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(src_heat(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), divb(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(eint(1, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), eint0(1, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(einteq(1, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), rhs_e(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    eint = 0.0_dp;  eint0 = 0.0_dp;  rhs_e = 0.0_dp
    u = 0.0_dp;  w = 0.0_dp;  bf = 0.0_dp;  rhs = 0.0_dp;  emf = 0.0_dp;  ef = 0.0_dp;  jc = 0.0_dp;  gradw = 0.0_dp
    eta_cell = eta_c;  nu_cell = nu_c;  chi_par_cell = chi_par_c;  chi_perp_cell = chi_perp_c;  D_cell = D_perp_c
    nu_sgs = 0.0_dp;  eta_sgs = 0.0_dp;  src_heat = 0.0_dp;  divb = 0.0_dp

    if (ic_type == 4) call tearing_setup
    if (ic_type == 5) call solovev_setup

    ! cell-centred primitive state (all cells incl. ghosts; periodic ghosts are overwritten by the exchange)
    do k = 1-ngk, nz+ngk
       do j = 1-ngj, ny+ngj
          do i = 1-ngi, nx+ngi
             x = xc(:, i, j, k)
             call ic_state(x, wc)
             w(:, i, j, k) = wc
          end do
       end do
    end do

    ! face fluxes
    useA = (ic_type == 2 .or. ic_type == 3 .or. ic_type == 4 .or. ic_type == 5)
    if (useA) then
       call faces_from_potential
    else
       call faces_from_B
    end if
    if (ic_type == 5) call add_toroidal_faces
    call exchange(bf, 3)
    call cell_B_from_faces(bf, w)
    call heating_profile
    ! the unperturbed state: fixed boundary values and the equilibrium residual (SW_Eq_Subtract)
    do k = 1-ngk, nz+ngk
       do j = 1-ngj, ny+ngj
          do i = 1-ngi, nx+ngi
             call prim2cons(w(:, i, j, k), ueq(:, i, j, k))
          end do
       end do
    end do
    call exchange(ueq, NVAR)
    einteq(1, :, :, :) = w(5, :, :, :)/(gam - 1.0_dp)
    call exchange(einteq, 1)
    bfeq = bf
    ! the initial state: equilibrium + perturbation
    call add_perturbation
    do k = 1-ngk, nz+ngk
       do j = 1-ngj, ny+ngj
          do i = 1-ngi, nx+ngi
             call prim2cons(w(:, i, j, k), u(:, i, j, k))
          end do
       end do
    end do
    call exchange(u, NVAR)
    eint(1, :, :, :) = w(5, :, :, :)/(gam - 1.0_dp)
    call exchange(eint, 1)
    call report_boundary_flux
  end subroutine initialize_fields

  ! largest normal field on the non-periodic boundaries (relative to |B|): zero on a flux-surface boundary
  subroutine report_boundary_flux
    integer :: c, side, i, j, k
    real(dp) :: bmax(2,3), bl, gm
    bmax = 0.0_dp
    do c = 1, 3
       if (periodic(c) == 1) cycle
       do side = 1, 2
          if (side == 1 .and. coords(c) /= 0) cycle
          if (side == 2 .and. coords(c) /= dims(c) - 1) cycle
          bl = 0.0_dp
          select case (c)
          case (1); i = merge(1, nx+1, side == 1)
             do k = 1, nz; do j = 1, ny
                bl = max(bl, abs(bf(1,i,j,k))/max(norm2(sf(:,1,i,j,k))*norm2(w(6:8,min(i,nx),j,k)), 1e-30_dp))
             end do; end do
          case (2); j = merge(1, ny+1, side == 1)
             do k = 1, nz; do i = 1, nx
                bl = max(bl, abs(bf(2,i,j,k))/max(norm2(sf(:,2,i,j,k))*norm2(w(6:8,i,min(j,ny),k)), 1e-30_dp))
             end do; end do
          case (3); k = merge(1, nz+1, side == 1)
             do j = 1, ny; do i = 1, nx
                bl = max(bl, abs(bf(3,i,j,k))/max(norm2(sf(:,3,i,j,k))*norm2(w(6:8,i,j,min(k,nz))), 1e-30_dp))
             end do; end do
          end select
          bmax(side, c) = bl
       end do
    end do
    do c = 1, 3
       do side = 1, 2
          gm = gmax(bmax(side, c))
          if (myrank == 0 .and. periodic(c) /= 1) write(*,'(A,I0,A,I0,A,ES10.3)') '  boundary B_n/|B| (dir ', c, ', side ', side, '): ', gm
       end do
    end do
  end subroutine report_boundary_flux

  !--------------------------------------------------------------------
  ! primitive state at a point (code units)
  !--------------------------------------------------------------------
  subroutine ic_state(x, wc)
    real(dp), intent(in)  :: x(3)
    real(dp), intent(out) :: wc(NVAR)
    real(dp) :: s, ca, sa, x1, amp, r, th, bth, bz, pr, psiN, R_, Z_, phi, bRR, bZZ, bPh, pp, nn
    wc = 0.0_dp
    select case (ic_type)
    case (1)   ! Brio-Wu: direction ic_dir, jump at ic_param(1) (default 0.5)
       s = x(ic_dir)
       if (s < merge(ic_param(1), 0.5_dp, ic_param(1) /= 0.0_dp)) then
          wc(1) = 1.0_dp;  wc(5) = 1.0_dp;  wc(5+c2(ic_dir)) = 1.0_dp
       else
          wc(1) = 0.125_dp;  wc(5) = 0.1_dp;  wc(5+c2(ic_dir)) = -1.0_dp
       end if
       wc(5+ic_dir) = 0.75_dp
    case (2)   ! Orszag-Tang on [0,1]^2
       wc(1) = 25.0_dp/(36.0_dp*PI);  wc(5) = 5.0_dp/(12.0_dp*PI)
       wc(2) = -sin(2.0_dp*PI*x(2));  wc(3) = sin(2.0_dp*PI*x(1))
       wc(6) = -sin(2.0_dp*PI*x(2))/sqrt(4.0_dp*PI);  wc(7) = sin(4.0_dp*PI*x(1))/sqrt(4.0_dp*PI)
    case (3)   ! CP Alfven wave: k = 2 pi (cos a, sin a), tan a = ic_param(1) (default 2), amplitude ic_param(2) (0.1)
       call alfven_frame(ca, sa, amp)
       x1 = x(1)*ca + x(2)*sa
       wc(1) = 1.0_dp;  wc(5) = 0.1_dp
       ! rotated-frame components: B1 = 1, B2 = amp sin(2 pi x1), B3 = amp cos(2 pi x1); v = (0, B2, B3) (wave travelling towards -x1)
       wc(6) = ca - amp*sin(2.0_dp*PI*x1)*sa;   wc(7) = sa + amp*sin(2.0_dp*PI*x1)*ca;   wc(8) = amp*cos(2.0_dp*PI*x1)
       wc(2) = -amp*sin(2.0_dp*PI*x1)*sa;       wc(3) = amp*sin(2.0_dp*PI*x1)*ca;        wc(4) = amp*cos(2.0_dp*PI*x1)
    case (4)   ! cylinder
       r = sqrt(x(1)**2 + x(2)**2);  th = atan2(x(2), x(1))
       call tearing_profile(r, bth, bz, pr)
       wc(1) = 1.0_dp;  wc(5) = pr
       wc(6) = -bth*sin(th);  wc(7) = bth*cos(th);  wc(8) = bz
    case (5)   ! tokamak (SI inside, converted)
       R_ = sqrt(x(1)**2 + x(2)**2)*Lch;  Z_ = x(3)*Lch;  phi = atan2(x(2), x(1))
       call solovev_fields(R_, Z_, bRR, bZZ, bPh, pp, nn, psiN)
       wc(1) = nn*mion/Rhoch
       wc(5) = pp/p0
       wc(6) = (bRR*cos(phi) - bPh*sin(phi))/Bch
       wc(7) = (bRR*sin(phi) + bPh*cos(phi))/Bch
       wc(8) = bZZ/Bch
    case default
       wc(1) = rho_init;  wc(5) = p_init;  wc(2:4) = v_init;  wc(6:8) = B_init
    end select
  end subroutine ic_state

  subroutine alfven_frame(ca, sa, amp)
    real(dp), intent(out) :: ca, sa, amp
    real(dp) :: ta
    ta = merge(ic_param(1), 2.0_dp, ic_param(1) /= 0.0_dp)
    ca = 1.0_dp/sqrt(1.0_dp + ta**2);  sa = ta*ca
    amp = merge(ic_param(2), 0.1_dp, ic_param(2) /= 0.0_dp)
  end subroutine alfven_frame

  !--------------------------------------------------------------------
  ! vector potential at a point (code units) for the IC types that have one
  !--------------------------------------------------------------------
  function ic_potential(x) result(a)
    real(dp), intent(in) :: x(3)
    real(dp) :: a(3)
    real(dp) :: ca, sa, amp, x1, x2, a2, a3, r, th, ath, az, R_, Z_, psi, phi
    a = 0.0_dp
    select case (ic_type)
    case (2)
       a(3) = cos(4.0_dp*PI*x(1))/(4.0_dp*PI*sqrt(4.0_dp*PI)) + cos(2.0_dp*PI*x(2))/(2.0_dp*PI*sqrt(4.0_dp*PI))
    case (3)
       call alfven_frame(ca, sa, amp)
       x1 = x(1)*ca + x(2)*sa;  x2 = -x(1)*sa + x(2)*ca
       a2 = amp*sin(2.0_dp*PI*x1)/(2.0_dp*PI)
       a3 = x2 + amp*cos(2.0_dp*PI*x1)/(2.0_dp*PI)
       a(1) = -a2*sa;  a(2) = a2*ca;  a(3) = a3
    case (4)
       r = sqrt(x(1)**2 + x(2)**2);  th = atan2(x(2), x(1))
       call tearing_potential(r, ath, az)
       a(1) = -ath*sin(th);  a(2) = ath*cos(th);  a(3) = az
    case (5)   ! poloidal field from A_phi = psi/R (SI); code units: A* = A / (B0 L0)
       R_ = sqrt(x(1)**2 + x(2)**2)*Lch;  Z_ = x(3)*Lch;  phi = atan2(x(2), x(1))
       psi = tok_psi0*solovev_Psi(R_/tok_R0, Z_/tok_R0)
       a(1) = -psi/R_*sin(phi)/(Bch*Lch);  a(2) = psi/R_*cos(phi)/(Bch*Lch);  a(3) = 0.0_dp
    end select
  end function ic_potential

  ! face flux = circulation of A around the face (edges in the right-hand order of S_f), 3-point Gauss per edge
  subroutine faces_from_potential
    integer :: i, j, k, c, e
    real(dp) :: p(3, 4), circ, xa(3), xb(3), gl(3), gw(3)
    gl = [-sqrt(0.6_dp), 0.0_dp, sqrt(0.6_dp)];  gw = [5.0_dp, 8.0_dp, 5.0_dp]/18.0_dp
    do k = 1-ngk, nz+ngk
       do j = 1-ngj, ny+ngj
          do i = 1-ngi, nx+ngi
             do c = 1, 3
                select case (c)
                case (1)
                   if (j+1 > ny+ngj+1 .or. k+1 > nz+ngk+1) cycle
                   p(:,1) = xn(:,i,j,k); p(:,2) = xn(:,i,j+1,k); p(:,3) = xn(:,i,j+1,k+1); p(:,4) = xn(:,i,j,k+1)
                case (2)
                   if (i+1 > nx+ngi+1 .or. k+1 > nz+ngk+1) cycle
                   p(:,1) = xn(:,i,j,k); p(:,2) = xn(:,i,j,k+1); p(:,3) = xn(:,i+1,j,k+1); p(:,4) = xn(:,i+1,j,k)
                case (3)
                   if (i+1 > nx+ngi+1 .or. j+1 > ny+ngj+1) cycle
                   p(:,1) = xn(:,i,j,k); p(:,2) = xn(:,i+1,j,k); p(:,3) = xn(:,i+1,j+1,k); p(:,4) = xn(:,i,j+1,k)
                end select
                circ = 0.0_dp
                do e = 1, 4
                   xa = p(:, e);  xb = p(:, mod(e, 4) + 1)
                   circ = circ + edge_integral(xa, xb, gl, gw)
                end do
                bf(c, i, j, k) = circ
             end do
          end do
       end do
    end do
  end subroutine faces_from_potential

  real(dp) function edge_integral(xa, xb, gl, gw)
    real(dp), intent(in) :: xa(3), xb(3), gl(3), gw(3)
    integer :: q
    real(dp) :: xq(3), a(3)
    edge_integral = 0.0_dp
    do q = 1, 3
       xq = 0.5_dp*(xa + xb) + 0.5_dp*gl(q)*(xb - xa)
       a = ic_potential(xq)
       edge_integral = edge_integral + gw(q)*dot_product(a, xb - xa)
    end do
  end function edge_integral

  subroutine faces_from_B
    integer :: i, j, k, c
    real(dp) :: wc(NVAR)
    do k = 1-ngk, nz+ngk
       do j = 1-ngj, ny+ngj
          do i = 1-ngi, nx+ngi
             do c = 1, 3
                call ic_state(xf(:, c, i, j, k), wc)
                bf(c, i, j, k) = dot_product(wc(6:8), sf(:, c, i, j, k))
             end do
          end do
       end do
    end do
  end subroutine faces_from_B

  ! tokamak: the toroidal field F(psi)/R through every face (B_phi . S_f, 2x2 Gauss points on the face);
  ! exactly divergence-free for an axisymmetric grid (it only threads the k faces)
  subroutine add_toroidal_faces
    integer :: i, j, k, c, q1, q2
    real(dp) :: p(3,4), xq(3), R_, Z_, phi, bRR, bZZ, bPh, pp, nn, psiN, bvec(3), s(3), flux, g2(2), xi, et
    g2 = [-1.0_dp/sqrt(3.0_dp), 1.0_dp/sqrt(3.0_dp)]
    do k = 1-ngk, nz+ngk
       do j = 1-ngj, ny+ngj
          do i = 1-ngi, nx+ngi
             do c = 1, 3
                select case (c)
                case (1)
                   if (j+1 > ny+ngj+1 .or. k+1 > nz+ngk+1) cycle
                   p(:,1) = xn(:,i,j,k); p(:,2) = xn(:,i,j+1,k); p(:,3) = xn(:,i,j+1,k+1); p(:,4) = xn(:,i,j,k+1)
                case (2)
                   if (i+1 > nx+ngi+1 .or. k+1 > nz+ngk+1) cycle
                   p(:,1) = xn(:,i,j,k); p(:,2) = xn(:,i,j,k+1); p(:,3) = xn(:,i+1,j,k+1); p(:,4) = xn(:,i+1,j,k)
                case (3)
                   if (i+1 > nx+ngi+1 .or. j+1 > ny+ngj+1) cycle
                   p(:,1) = xn(:,i,j,k); p(:,2) = xn(:,i+1,j,k); p(:,3) = xn(:,i+1,j+1,k); p(:,4) = xn(:,i,j+1,k)
                end select
                flux = 0.0_dp
                do q1 = 1, 2
                   do q2 = 1, 2
                      xi = g2(q1);  et = g2(q2)
                      ! bilinear map of the quadrilateral
                      xq = 0.25_dp*((1-xi)*(1-et)*p(:,1) + (1+xi)*(1-et)*p(:,2) + (1+xi)*(1+et)*p(:,3) + (1-xi)*(1+et)*p(:,4))
                      s  = 0.25_dp*cross( 0.25_dp*(-(1-et)*p(:,1) + (1-et)*p(:,2) + (1+et)*p(:,3) - (1+et)*p(:,4)), &
                                          0.25_dp*(-(1-xi)*p(:,1) - (1+xi)*p(:,2) + (1+xi)*p(:,3) + (1-xi)*p(:,4)) )*4.0_dp
                      R_ = sqrt(xq(1)**2 + xq(2)**2)*Lch;  Z_ = xq(3)*Lch;  phi = atan2(xq(2), xq(1))
                      call solovev_fields(R_, Z_, bRR, bZZ, bPh, pp, nn, psiN)
                      bvec = [-bPh*sin(phi), bPh*cos(phi), 0.0_dp]/Bch
                      flux = flux + dot_product(bvec, s)
                   end do
                end do
                bf(c, i, j, k) = bf(c, i, j, k) + flux
             end do
          end do
       end do
    end do
  end subroutine add_toroidal_faces

  !--------------------------------------------------------------------
  ! cell-centred Cartesian B from the six face fluxes:
  !   B = (1/V) sum_c 0.5 (Phi_c- + Phi_c+) a_c,  a_c = x_f(c+) - x_f(c-)
  ! (exact for a parallelepiped; the dual-basis identity)
  !--------------------------------------------------------------------
  subroutine cell_B_from_faces(b, wp)
    real(dp), intent(in)    :: b(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(inout) :: wp(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: i, j, k
    real(dp) :: bb(3)
    do k = 1-ngk, nz+ngk-1
       do j = 1-ngj, ny+ngj-1
          do i = 1-ngi, nx+ngi-1
             bb = 0.5_dp*(b(1,i,j,k) + b(1,i+1,j,k))*(xf(:,1,i+1,j,k) - xf(:,1,i,j,k)) + &
                  0.5_dp*(b(2,i,j,k) + b(2,i,j+1,k))*(xf(:,2,i,j+1,k) - xf(:,2,i,j,k)) + &
                  0.5_dp*(b(3,i,j,k) + b(3,i,j,k+1))*(xf(:,3,i,j,k+1) - xf(:,3,i,j,k))
             wp(6:8, i, j, k) = bb/vol(i, j, k)
          end do
       end do
    end do
  end subroutine cell_B_from_faces

  !--------------------------------------------------------------------
  subroutine add_perturbation
    integer :: i, j, k, m, n
    real(dp) :: x(3), r, th, ph, f, vr, rs, wdt, rnd(3), R_, Z_, rho_pol, psiN, ang
    integer, allocatable :: seed(:)
    integer :: ns
    if (perturb_amp == 0.0_dp) return
    m = perturb_m;  n = perturb_n
    call random_seed(size=ns);  allocate(seed(ns));  seed = perturb_seed + 97*myrank;  call random_seed(put=seed)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             x = xc(:, i, j, k)
             if (m == 0) then
                call random_number(rnd)
                w(2:4, i, j, k) = w(2:4, i, j, k) + perturb_amp*(2.0_dp*rnd - 1.0_dp)
                cycle
             end if
             select case (ic_type)
             case (4)
                r = sqrt(x(1)**2 + x(2)**2);  th = atan2(x(2), x(1))
                rs = tear_r_s;  wdt = merge(ic_param(5), 0.2_dp, ic_param(5) /= 0.0_dp)
                f = exp(-((r - rs)/wdt)**2)
                vr = perturb_amp*f*cos(real(m,dp)*th - real(n,dp)*2.0_dp*PI*x(3)/cyl_length())
                w(2, i, j, k) = w(2, i, j, k) + vr*cos(th);  w(3, i, j, k) = w(3, i, j, k) + vr*sin(th)
             case (5)
                R_ = sqrt(x(1)**2 + x(2)**2)*Lch;  Z_ = x(3)*Lch;  ph = atan2(x(2), x(1))
                psiN = solovev_psiN(R_/tok_R0, Z_/tok_R0)
                rho_pol = sqrt(max(psiN, 0.0_dp))
                th = atan2(Z_, R_ - sv_x_axis*tok_R0)
                f = exp(-((rho_pol - 0.5_dp)/0.2_dp)**2)*max(1.0_dp - psiN, 0.0_dp)
                vr = perturb_amp*f*cos(real(m,dp)*th - real(n,dp)*ph)
                ! radial in the poloidal plane (unit vector from the axis)
                ang = th
                w(2, i, j, k) = w(2, i, j, k) + vr*cos(ang)*cos(ph)
                w(3, i, j, k) = w(3, i, j, k) + vr*cos(ang)*sin(ph)
                w(4, i, j, k) = w(4, i, j, k) + vr*sin(ang)
             case default
                w(2, i, j, k) = w(2, i, j, k) + perturb_amp*sin(2.0_dp*PI*real(m,dp)*x(1))*cos(2.0_dp*PI*real(n,dp)*x(2))
             end select
          end do
       end do
    end do
  end subroutine add_perturbation

  real(dp) function cyl_length()
    cyl_length = merge(ic_param(4), 2.0_dp*PI*3.0_dp, ic_param(4) /= 0.0_dp)
  end function cyl_length

  ! heating profile (code units): tokamak - Gaussian in the poloidal distance from the axis; others - uniform
  subroutine heating_profile
    integer :: i, j, k
    real(dp) :: x(3), R_, Z_, d, wdt
    if (P_heat_c == 0.0_dp) return
    do k = 1-ngk, nz+ngk
       do j = 1-ngj, ny+ngj
          do i = 1-ngi, nx+ngi
             x = xc(:, i, j, k)
             if (ic_type == 5) then
                R_ = sqrt(x(1)**2 + x(2)**2)*Lch;  Z_ = x(3)*Lch
                d = sqrt((R_ - sv_x_axis*tok_R0)**2 + Z_**2)
                wdt = heat_width*tok_a
                src_heat(i, j, k) = P_heat_c*exp(-(d/wdt)**2)
             else
                src_heat(i, j, k) = P_heat_c
             end if
          end do
       end do
    end do
  end subroutine heating_profile

  !====================================================================
  ! cylindrical tearing equilibrium (code units):
  !   B_z0 = 1 on axis, B_theta = r / (R0 q(r)),  q = q0 (1 + (r/r0)^(2 lam))^(1/lam)
  !   p = beta/2 uniform, B_z(r) from radial force balance
  !   ic_param: (1) q0 [1.2], (2) r0 [0.6], (3) lam [2], (4) L_z = 2 pi R0 [2 pi * 3], (5) perturbation width
  !====================================================================

  subroutine tearing_setup
    real(dp) :: q0, r0, lam, R0_, dr, r, rmax, bth2_int, s
    integer :: n
    q0 = merge(ic_param(1), 1.2_dp, ic_param(1) /= 0.0_dp)
    r0 = merge(ic_param(2), 0.6_dp, ic_param(2) /= 0.0_dp)
    lam = merge(ic_param(3), 2.0_dp, ic_param(3) /= 0.0_dp)
    R0_ = cyl_length()/(2.0_dp*PI)
    tear_p = merge(p_init, 0.1_dp, p_init > 0.0_dp)
    rmax = 1.5_dp
    dr = rmax/real(NT-1, dp)
    bth2_int = 0.0_dp
    do n = 1, NT
       r = real(n-1, dp)*dr
       tab_r(n) = r
       tab_q(n) = q0*(1.0_dp + (r/r0)**(2.0_dp*lam))**(1.0_dp/lam)
       tab_bth(n) = r/(R0_*tab_q(n))
       if (n > 1) bth2_int = bth2_int + 0.5_dp*(tab_bth(n)**2/max(r, 1e-12_dp) + tab_bth(n-1)**2/max(tab_r(n-1), 1e-12_dp))*dr
       s = 1.0_dp - tab_bth(n)**2 - 2.0_dp*bth2_int
       tab_bz(n) = sqrt(max(s, 1e-6_dp))
    end do
    ! potentials: A_z = -int B_theta dr,  A_theta = (1/r) int r B_z dr
    tab_az(1) = 0.0_dp;  tab_ath(1) = 0.0_dp;  s = 0.0_dp
    do n = 2, NT
       tab_az(n) = tab_az(n-1) - 0.5_dp*(tab_bth(n) + tab_bth(n-1))*dr
       s = s + 0.5_dp*(tab_r(n)*tab_bz(n) + tab_r(n-1)*tab_bz(n-1))*dr
       tab_ath(n) = s/tab_r(n)
    end do
    ! rational surface of the perturbation: q = m/n
    tear_r_s = 0.0_dp
    do n = 2, NT
       if ((tab_q(n-1) - real(perturb_m,dp)/real(max(perturb_n,1),dp))*(tab_q(n) - real(perturb_m,dp)/real(max(perturb_n,1),dp)) <= 0.0_dp) then
          tear_r_s = tab_r(n);  exit
       end if
    end do
    if (myrank == 0) then
       write(*,'(A,F6.3,A,F6.3,A,F5.2,A,F8.4)') '  tearing equilibrium: q0 ', q0, '  r0 ', r0, '  lambda ', lam, '  L_z ', cyl_length()
       write(*,'(A,F8.4,A,I0,A,I0,A)') '                       rational surface r_s = ', tear_r_s, '  (m/n = ', perturb_m, '/', perturb_n, ')'
       write(*,'(A,F8.4,A,F8.4)')  '                       B_z(axis) 1.0000  B_z(1) ', tab_bz(nint(1.0_dp/dr)+1), '  q(1) ', tab_q(nint(1.0_dp/dr)+1)
    end if
  end subroutine tearing_setup

  subroutine tearing_profile(r, bth, bz, pr)
    real(dp), intent(in)  :: r
    real(dp), intent(out) :: bth, bz, pr
    real(dp) :: f
    integer :: n
    f = r/(tab_r(2) - tab_r(1));  n = min(max(int(f) + 1, 1), NT - 1);  f = f - real(n-1, dp)
    bth = (1.0_dp - f)*tab_bth(n) + f*tab_bth(n+1)
    bz  = (1.0_dp - f)*tab_bz(n)  + f*tab_bz(n+1)
    pr  = tear_p
  end subroutine tearing_profile

  subroutine tearing_potential(r, ath, az)
    real(dp), intent(in)  :: r
    real(dp), intent(out) :: ath, az
    real(dp) :: f
    integer :: n
    f = r/(tab_r(2) - tab_r(1));  n = min(max(int(f) + 1, 1), NT - 1);  f = f - real(n-1, dp)
    ath = (1.0_dp - f)*tab_ath(n) + f*tab_ath(n+1)
    az  = (1.0_dp - f)*tab_az(n)  + f*tab_az(n+1)
  end subroutine tearing_potential

  !====================================================================
  ! Cerfon & Freidberg (2010) Solov'ev equilibrium, up-down symmetric:
  !   Delta* Psi = (1-A) x^2 + A,  Psi = Psi_p + sum c_i Psi_i
  !   7 coefficients from 3 boundary points, 3 curvature conditions and
  !   one extremum condition (the same construction as the Python
  !   generator of the KSTAR grid).  Dimensional: psi = psi0 Psi.
  !====================================================================
  subroutine solovev_setup
    real(dp) :: M(7,7), b(7), x, dx, f, fp
    integer :: it
    sv_eps = tok_a/tok_R0
    call solovev_coeff(sv_eps, tok_kappa, tok_delta, tok_Asol, sv_c)
    ! magnetic axis: extremum of Psi on y = 0 (Newton on dPsi/dx)
    x = 1.0_dp + 0.3_dp*sv_eps
    do it = 1, 100
       f  = Psi_x(x, 0.0_dp)
       dx = 1.0e-6_dp
       fp = (Psi_x(x + dx, 0.0_dp) - Psi_x(x - dx, 0.0_dp))/(2.0_dp*dx)
       x = x - f/fp
       if (abs(f/fp) < 1.0e-14_dp) exit
    end do
    sv_x_axis = x
    sv_psi_axis = solovev_Psi(x, 0.0_dp)
    p_axis_c = -(1.0_dp - tok_Asol)*tok_psi0**2*sv_psi_axis/(MU0*tok_R0**4)/p0
    R_axis_c = sv_x_axis*tok_R0/Lch
    if (myrank == 0) then
       write(*,'(A,ES12.4,A,ES12.4)') "  Solov'ev: Psi_axis ", sv_psi_axis, '  R_axis [m] ', sv_x_axis*tok_R0
       write(*,'(A,ES12.4,A,ES12.4)') '            p_axis [Pa] ', p_axis_c*p0, '  Shafranov shift [m] ', (sv_x_axis - 1.0_dp)*tok_R0
       write(*,'(A,7ES12.4)') '            coefficients ', sv_c
    end if
  end subroutine solovev_setup

  subroutine solovev_coeff(eps, kap, dlt, A, c)
    real(dp), intent(in)  :: eps, kap, dlt, A
    real(dp), intent(out) :: c(7)
    real(dp) :: M(7,7), b(7), alpha, N1, N2, N3, xo, xi, xt, yt
    real(dp) :: val(0:7, 5)     ! rows: particular + 7 basis; cols: f, fx, fxx, fy, fyy at a point
    integer :: n
    alpha = asin(dlt)
    N1 = -(1.0_dp + alpha)**2/(eps*kap**2)
    N2 =  (1.0_dp - alpha)**2/(eps*kap**2)
    N3 = -kap/(eps*cos(alpha)**2)
    xo = 1.0_dp + eps;  xi = 1.0_dp - eps;  xt = 1.0_dp - dlt*eps;  yt = kap*eps
    ! constraint rows: for each function n: [f(P_out), f(P_in), f(P_top), fx(P_top), fyy(P_out)+N1 fx(P_out), fyy(P_in)+N2 fx(P_in), fxx(P_top)+N3 fy(P_top)]
    do n = 0, 7
       call basis_all(n, xo, 0.0_dp, A, val(n, :))
       M(1, max(n,1)) = val(n, 1)
       M(5, max(n,1)) = val(n, 5) + N1*val(n, 2)
       if (n == 0) then; b(1) = val(0,1); b(5) = val(0,5) + N1*val(0,2); end if
       call basis_all(n, xi, 0.0_dp, A, val(n, :))
       M(2, max(n,1)) = val(n, 1)
       M(6, max(n,1)) = val(n, 5) + N2*val(n, 2)
       if (n == 0) then; b(2) = val(0,1); b(6) = val(0,5) + N2*val(0,2); end if
       call basis_all(n, xt, yt, A, val(n, :))
       M(3, max(n,1)) = val(n, 1)
       M(4, max(n,1)) = val(n, 2)
       M(7, max(n,1)) = val(n, 3) + N3*val(n, 4)
       if (n == 0) then; b(3) = val(0,1); b(4) = val(0,2); b(7) = val(0,3) + N3*val(0,4); end if
    end do
    b = -b
    call solve7(M, b, c)
  end subroutine solovev_coeff

  ! f, fx, fxx, fy, fyy of the particular solution (n = 0) and the basis functions (n = 1..7)
  subroutine basis_all(n, x, y, A, v)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: x, y, A
    real(dp), intent(out) :: v(5)
    real(dp) :: L
    L = log(x)
    select case (n)
    case (0)
       v(1) = x**4/8 + A*(x**2*L/2 - x**4/8)
       v(2) = x**3/2 + A*(x*L + x/2 - x**3/2)
       v(3) = 1.5_dp*x**2 + A*(L + 1.5_dp - 1.5_dp*x**2)
       v(4) = 0.0_dp;  v(5) = 0.0_dp
    case (1)
       v = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    case (2)
       v = [x**2, 2*x, 2.0_dp, 0.0_dp, 0.0_dp]
    case (3)
       v = [y**2 - x**2*L, -2*x*L - x, -2*L - 3, 2*y, 2.0_dp]
    case (4)
       v = [x**4 - 4*x**2*y**2, 4*x**3 - 8*x*y**2, 12*x**2 - 8*y**2, -8*x**2*y, -8*x**2]
    case (5)
       v(1) = 2*y**4 - 9*y**2*x**2 + 3*x**4*L - 12*x**2*y**2*L
       v(2) = -30*x*y**2 + 12*x**3*L + 3*x**3 - 24*x*y**2*L
       v(3) = -54*y**2 + 36*x**2*L + 21*x**2 - 24*y**2*L
       v(4) = 8*y**3 - 18*y*x**2 - 24*x**2*y*L
       v(5) = 24*y**2 - 18*x**2 - 24*x**2*L
    case (6)
       v(1) = x**6 - 12*x**4*y**2 + 8*x**2*y**4
       v(2) = 6*x**5 - 48*x**3*y**2 + 16*x*y**4
       v(3) = 30*x**4 - 144*x**2*y**2 + 16*y**4
       v(4) = -24*x**4*y + 32*x**2*y**3
       v(5) = -24*x**4 + 96*x**2*y**2
    case (7)
       v(1) = 8*y**6 - 140*y**4*x**2 + 75*y**2*x**4 - 15*x**6*L + 180*x**4*y**2*L - 120*x**2*y**4*L
       v(2) = -400*x*y**4 + 480*x**3*y**2 - 90*x**5*L - 15*x**5 + 720*x**3*y**2*L - 240*x*y**4*L
       v(3) = -640*y**4 + 2160*x**2*y**2 - 450*x**4*L - 165*x**4 + 2160*x**2*y**2*L - 240*y**4*L
       v(4) = 48*y**5 - 560*y**3*x**2 + 150*y*x**4 + 360*x**4*y*L - 480*x**2*y**3*L
       v(5) = 240*y**4 - 1680*y**2*x**2 + 150*x**4 + 360*x**4*L - 1440*x**2*y**2*L
    end select
  end subroutine basis_all

  subroutine solve7(M, b, c)
    real(dp), intent(inout) :: M(7,7), b(7)
    real(dp), intent(out)   :: c(7)
    integer :: i, j, k, p
    real(dp) :: f, tmp(7), tb
    do k = 1, 7
       p = k
       do i = k+1, 7
          if (abs(M(i,k)) > abs(M(p,k))) p = i
       end do
       if (p /= k) then
          tmp = M(k,:);  M(k,:) = M(p,:);  M(p,:) = tmp
          tb = b(k);  b(k) = b(p);  b(p) = tb
       end if
       do i = k+1, 7
          f = M(i,k)/M(k,k)
          M(i,:) = M(i,:) - f*M(k,:)
          b(i) = b(i) - f*b(k)
       end do
    end do
    do i = 7, 1, -1
       c(i) = b(i)
       do j = i+1, 7
          c(i) = c(i) - M(i,j)*c(j)
       end do
       c(i) = c(i)/M(i,i)
    end do
  end subroutine solve7

  real(dp) function solovev_Psi(x, y)
    real(dp), intent(in) :: x, y
    real(dp) :: v(5)
    integer :: n
    call basis_all(0, x, y, tok_Asol, v);  solovev_Psi = v(1)
    do n = 1, 7
       call basis_all(n, x, y, tok_Asol, v);  solovev_Psi = solovev_Psi + sv_c(n)*v(1)
    end do
  end function solovev_Psi

  real(dp) function Psi_x(x, y)
    real(dp), intent(in) :: x, y
    real(dp) :: v(5)
    integer :: n
    call basis_all(0, x, y, tok_Asol, v);  Psi_x = v(2)
    do n = 1, 7
       call basis_all(n, x, y, tok_Asol, v);  Psi_x = Psi_x + sv_c(n)*v(2)
    end do
  end function Psi_x

  real(dp) function Psi_y(x, y)
    real(dp), intent(in) :: x, y
    real(dp) :: v(5)
    integer :: n
    call basis_all(0, x, y, tok_Asol, v);  Psi_y = v(4)
    do n = 1, 7
       call basis_all(n, x, y, tok_Asol, v);  Psi_y = Psi_y + sv_c(n)*v(4)
    end do
  end function Psi_y

  real(dp) function solovev_psiN(x, y)
    real(dp), intent(in) :: x, y
    solovev_psiN = (solovev_Psi(x, y) - sv_psi_axis)/(0.0_dp - sv_psi_axis)
  end function solovev_psiN

  real(dp) function solovev_axis_R()
    solovev_axis_R = sv_x_axis*tok_R0
  end function solovev_axis_R

  ! SI fields at (R, Z): B_R, B_Z, B_phi [T], p [Pa], n [m-3], psiN
  subroutine solovev_fields(R_, Z_, bRR, bZZ, bPh, pp, nn, psiN)
    real(dp), intent(in)  :: R_, Z_
    real(dp), intent(out) :: bRR, bZZ, bPh, pp, nn, psiN
    real(dp) :: x, y, Ps, Px, Py, F2, pax, prel
    x = R_/tok_R0;  y = Z_/tok_R0
    Ps = solovev_Psi(x, y);  Px = Psi_x(x, y);  Py = Psi_y(x, y)
    psiN = (Ps - sv_psi_axis)/(0.0_dp - sv_psi_axis)
    ! outside the boundary (psiN > 1): hold the boundary values of p and n, keep the field
    bRR = -tok_psi0*Py/(tok_R0**2*x)
    bZZ =  tok_psi0*Px/(tok_R0**2*x)
    F2 = (tok_R0*tok_B0)**2 - 2.0_dp*tok_Asol*tok_psi0**2*min(Ps, 0.0_dp)/tok_R0**2
    bPh = sqrt(max(F2, 0.0_dp))/R_
    pax = -(1.0_dp - tok_Asol)*tok_psi0**2*sv_psi_axis/(MU0*tok_R0**4)
    prel = max(-(1.0_dp - tok_Asol)*tok_psi0**2*min(Ps, 0.0_dp)/(MU0*tok_R0**4)/pax, 0.0_dp)   ! p/p_axis in 0..1
    pp = pax*(prel + tok_p_floor)     ! uniform offset: the force balance J x B = grad p is kept exactly
    nn = tok_n_axis*max(prel**tok_n_exp, tok_n_floor)
  end subroutine solovev_fields

end module equilibrium
