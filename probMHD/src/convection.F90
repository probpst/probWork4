!=======================================================================
! convection.F90 -- convective (ideal MHD) face fluxes
!
!   Godunov finite-volume form in the transformed space: for every face
!   the left / right states are reconstructed along the index direction
!   (unit spacing in xi, as in probCFD) and an approximate Riemann
!   solver (Rusanov / HLL / HLLD) gives the flux in the face-normal
!   frame (n, t1, t2); the momentum flux is rotated back to Cartesian.
!
!   Reconstruction (Cnvt_scheme, probCFD convection.F90 function upw):
!     1  1st-order upwind      phi_f = phi_D
!     2  2nd-order upwind      phi_f = 3/2 phi_D - 1/2 phi_U
!     3  QUICK (3rd order)     phi_f = 3/8 phi_A + 3/4 phi_D - 1/8 phi_U
!     4  2nd-order ENO         phi_f = phi_D + 1/2 m(phi_A - phi_D, phi_D - phi_U)
!        (D = donor, A = acceptor, U = upstream of D; m = minmod: the
!         smoother of the central and the upwind slope, 0 at an
!         extremum -> essentially non-oscillatory)
!     5  3rd-order ENO (Harten-Engquist-Osher-Chakravarthy): the
!        smoothest of the three 3-cell stencils by hierarchical
!        divided differences, cell-average to face-value weights
!   The left state is the donor-side value of upw with vel > 0, the
!   right state the one with vel < 0, so the scheme reduces to the
!   probCFD ENO flux for a scalar advected quantity.
!
!   The face-normal magnetic field is the constrained-transport face
!   value (identical on both sides); the Riemann flux of the tangential
!   field components gives the face electric field for induction.F90.
!=======================================================================
module convection
  use global_vars
  use mpi_vars
  implicit none
  private
  public :: convective_fluxes, upw, minmod, max_signal

  real(dp), allocatable, public :: flx(:,:,:,:,:)   ! flx(NVAR, c, i,j,k): flux through the minus face of cell (i,j,k), direction c (times area)
  real(dp), allocatable, public :: flx_e(:,:,:,:)    ! upwind internal-energy flux e_int v . S (SW_Energy >= 2)
  real(dp), allocatable, public :: vflux(:,:,:,:)    ! volume flux v . S of the face (mass flux / upwind density)
  real(dp) :: dbg_fmax = 0.0_dp, dbg_ql(NVAR), dbg_qr(NVAR), dbg_ff(NVAR)
  integer :: dbg_loc(3) = 0

contains

  !--------------------------------------------------------------------
  ! probCFD_부산 convection.F90: upwind interpolation of the advected value
  !--------------------------------------------------------------------
  pure real(dp) function upw(vel, fmm, fm, fp, fpp)
    real(dp), intent(in) :: vel, fmm, fm, fp, fpp
    select case (Cnvt_scheme)
    case (2)                                   ! 2nd-order upwind
       if (vel >= 0.0_dp) then
          upw = 1.5_dp*fm - 0.5_dp*fmm
       else
          upw = 1.5_dp*fp - 0.5_dp*fpp
       end if
    case (3)                                   ! QUICK
       if (vel >= 0.0_dp) then
          upw = 0.375_dp*fp + 0.75_dp*fm - 0.125_dp*fmm
       else
          upw = 0.375_dp*fm + 0.75_dp*fp - 0.125_dp*fpp
       end if
    case (4)                                   ! 2nd-order ENO (minmod)
       if (vel >= 0.0_dp) then
          upw = fm + 0.5_dp*minmod(fp - fm, fm - fmm)
       else
          upw = fp - 0.5_dp*minmod(fp - fm, fpp - fp)
       end if
    case default                               ! 1st-order upwind
       if (vel >= 0.0_dp) then
          upw = fm
       else
          upw = fp
       end if
    end select
  end function upw

  ! m(x,y) = x if |x| <= |y|, y if |x| > |y|, 0 if x y < 0
  pure real(dp) function minmod(x, y)
    real(dp), intent(in) :: x, y
    if (x*y <= 0.0_dp) then
       minmod = 0.0_dp
    else if (abs(x) <= abs(y)) then
       minmod = x
    else
       minmod = y
    end if
  end function minmod

  !--------------------------------------------------------------------
  ! 3rd-order ENO: left face value at the interface between f(-1) and f(0)
  !   f(-3:2) = cell values; stencils {-3,-2,-1}, {-2,-1,0}, {-1,0,1}
  !--------------------------------------------------------------------
  pure real(dp) function eno3_left(f)
    real(dp), intent(in) :: f(-3:2)
    integer :: s
    ! start from cell -1; add the neighbour with the smaller first difference
    if (abs(f(-1) - f(-2)) <= abs(f(0) - f(-1))) then
       ! stencil {-2,-1}: second differences of {-3,-2,-1} and {-2,-1,0}
       if (abs(f(-1) - 2.0_dp*f(-2) + f(-3)) <= abs(f(0) - 2.0_dp*f(-1) + f(-2))) then
          s = 0
       else
          s = 1
       end if
    else
       ! stencil {-1,0}: {-2,-1,0} and {-1,0,1}
       if (abs(f(0) - 2.0_dp*f(-1) + f(-2)) <= abs(f(1) - 2.0_dp*f(0) + f(-1))) then
          s = 1
       else
          s = 2
       end if
    end if
    select case (s)
    case (0); eno3_left = (2.0_dp*f(-3) - 7.0_dp*f(-2) + 11.0_dp*f(-1))/6.0_dp
    case (1); eno3_left = (-f(-2) + 5.0_dp*f(-1) + 2.0_dp*f(0))/6.0_dp
    case default; eno3_left = (2.0_dp*f(-1) + 5.0_dp*f(0) - f(1))/6.0_dp
    end select
  end function eno3_left

  ! left and right face values from the 6-cell stencil f(-3:2) (face between -1 and 0)
  pure subroutine recon(f, fl, fr)
    real(dp), intent(in)  :: f(-3:2)
    real(dp), intent(out) :: fl, fr
    real(dp) :: g(-3:2)
    integer :: m
    if (Cnvt_scheme == 5) then
       fl = eno3_left(f)
       do m = -3, 2
          g(m) = f(-1 - m)          ! mirror: the right state is the left state of the reflected sequence
       end do
       fr = eno3_left(g)
    else
       fl = upw( 1.0_dp, f(-2), f(-1), f(0), f(1))
       fr = upw(-1.0_dp, f(-2), f(-1), f(0), f(1))
    end if
  end subroutine recon

  !--------------------------------------------------------------------
  ! Riemann fluxes on the minus faces of cells 0..n+1 (all directions)
  !--------------------------------------------------------------------
  subroutine convective_fluxes(wp, bb, smax)
    real(dp), intent(in)  :: wp(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)  :: bb(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(out) :: smax(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)   ! max signal speed per face (for dt)
    integer :: i, j, k, c, m, q
    real(dp) :: f(-3:2, NVAR), wl(NVAR), s(3), area, n(3), t1(3), t2(3), bn, fl(NVAR), ff(NVAR), sm
    real(dp) :: ql(NVAR), qr(NVAR), e1, e2
    logical :: degenerate
    if (.not. allocated(flx)) then
       allocate(flx(NVAR, 3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
       allocate(flx_e(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), vflux(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
       flx_e = 0.0_dp;  vflux = 0.0_dp
    end if
    smax = 0.0_dp
    do c = 1, 3
       do k = 0, nz+1
          do j = 0, ny+1
             do i = 0, nx+1
                degenerate = (c == 1 .and. ngx == 1) .or. (c == 2 .and. ngy == 1) .or. (c == 3 .and. ngz == 1)
                if (degenerate) then
                   ! a single cell plane in this direction (2-D / axisymmetric wedge): the two sides of the
                   ! face are the cell and its periodic image (one ghost layer) - no reconstruction needed,
                   ! but the flux through the face is kept (the hoop force of a wedge lives in it)
                   if ((c == 1 .and. i == 0) .or. (c == 2 .and. j == 0) .or. (c == 3 .and. k == 0)) then
                      ! minus face of the single ghost layer: not needed (no second ghost cell behind it)
                      flx(:, c, i, j, k) = 0.0_dp;  flx_e(c, i, j, k) = 0.0_dp;  vflux(c, i, j, k) = 0.0_dp
                      wl = wp(:, i, j, k)
                      ef(:, c, i, j, k) = -cross3(wl(2:4), wl(6:8))
                      cycle
                   end if
                   select case (c)
                   case (1); f(-1, :) = wp(:, i-1, j, k);  f(0, :) = wp(:, i, j, k)
                   case (2); f(-1, :) = wp(:, i, j-1, k);  f(0, :) = wp(:, i, j, k)
                   case (3); f(-1, :) = wp(:, i, j, k-1);  f(0, :) = wp(:, i, j, k)
                   end select
                else
                   ! stencil along c
                   select case (c)
                   case (1); do m = -3, 2; f(m, :) = wp(:, i+m, j, k); end do
                   case (2); do m = -3, 2; f(m, :) = wp(:, i, j+m, k); end do
                   case (3); do m = -3, 2; f(m, :) = wp(:, i, j, k+m); end do
                   end select
                end if
                s = sf(:, c, i, j, k);  area = norm2(s)
                if (area <= 0.0_dp) then
                   flx(:, c, i, j, k) = 0.0_dp;  ef(:, c, i, j, k) = 0.0_dp;  flx_e(c, i, j, k) = 0.0_dp;  vflux(c, i, j, k) = 0.0_dp
                   cycle
                end if
                n = s/area
                ! tangential basis from the face's own edge along the next grid direction: intrinsic to
                ! the grid, so it rotates with a toroidal wedge and the limiter sees the same components
                call tangents(n, le(:, c2(c), i, j, k), t1, t2)
                bn = bb(c, i, j, k)/area
                ! the stencil vectors in the face frame (n, t1, t2): the limiter then acts on the
                ! normal / tangential components, which are invariant under the rotation of a
                ! toroidal wedge and keep a mirrored wall state exactly symmetric
                do m = merge(-1, -3, degenerate), merge(0, 2, degenerate)
                   wl = f(m, :)
                   call to_frame(wl, n, t1, t2, wl(6)*n(1) + wl(7)*n(2) + wl(8)*n(3), wl)
                   f(m, :) = wl
                end do
                if (Recon_vars == 2) then
                   do m = -3, 2
                      wl = f(m, :);  call prim2cons(wl, wl);  f(m, :) = wl
                   end do
                end if
                if (degenerate) then
                   ql = f(-1, :);  qr = f(0, :)
                else
                   do q = 1, NVAR
                      call recon(f(:, q), ql(q), qr(q))
                   end do
                end if
                if (Recon_vars == 2) then
                   call cons2prim(ql, ql);  call cons2prim(qr, qr)
                end if
                ql(1) = max(ql(1), rho_floor);  qr(1) = max(qr(1), rho_floor)
                ql(5) = max(ql(5), p_floor);    qr(5) = max(qr(5), p_floor)
                ql(6) = bn;  qr(6) = bn               ! the constrained-transport face value
                select case (Riemann_solver)
                case (1); call rusanov(ql, qr, ff, sm)
                case (2); call hll(ql, qr, ff, sm)
                case default; call hlld(ql, qr, ff, sm)
                end select
                if (debug_balance .and. c == 1 .and. i == nx+1 .and. coords(1) == dims(1)-1) then
                   if (abs(ff(5))*area > dbg_fmax) then
                      dbg_fmax = abs(ff(5))*area;  dbg_ql = ql;  dbg_qr = qr;  dbg_ff = ff;  dbg_loc = [i, j, k]
                   end if
                end if
                ! back to Cartesian, times the area
                fl(1) = ff(1)
                fl(2:4) = ff(2)*n + ff(3)*t1 + ff(4)*t2
                fl(5) = ff(5)
                fl(6:8) = ff(7)*t1 + ff(8)*t2            ! (F_Bt1, F_Bt2) rotated: cell-centred B flux (unused with CT)
                flx(:, c, i, j, k) = fl*area
                ! internal-energy transport: the upwind side by the sign of the mass flux
                if (ff(1) >= 0.0_dp) then
                   vflux(c, i, j, k) = ff(1)*area/ql(1);  flx_e(c, i, j, k) = vflux(c, i, j, k)*ql(5)/(gam - 1.0_dp)
                else
                   vflux(c, i, j, k) = ff(1)*area/qr(1);  flx_e(c, i, j, k) = vflux(c, i, j, k)*qr(5)/(gam - 1.0_dp)
                end if
                ! face electric field E = -v x B: E_t1 = F_Bt2, E_t2 = -F_Bt1 (n, t1, t2 right-handed)
                e1 = ff(8);  e2 = -ff(7)
                ef(:, c, i, j, k) = e1*t1 + e2*t2
                smax(c, i, j, k) = sm
             end do
          end do
       end do
    end do
    if (debug_balance .and. coords(1) == dims(1)-1) then
       write(*,'(A,3I5,A,ES11.3)') '  [wall face] largest |F_E|: face ', dbg_loc, '  F_E A = ', dbg_fmax
       write(*,'(A,8ES12.4)') '  [wall face] qL = ', dbg_ql
       write(*,'(A,8ES12.4)') '  [wall face] qR = ', dbg_qr
       write(*,'(A,8ES12.4)') '  [wall face] F  = ', dbg_ff
       dbg_fmax = 0.0_dp
    end if
  end subroutine convective_fluxes

  pure function cross3(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2);  c(2) = a(3)*b(1) - a(1)*b(3);  c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross3

  pure subroutine tangents(n, a, t1, t2)
    real(dp), intent(in)  :: n(3), a(3)
    real(dp), intent(out) :: t1(3), t2(3)
    real(dp) :: e(3)
    integer :: imin
    t1 = a - dot_product(a, n)*n
    if (norm2(t1) < 1.0e-12_dp*max(norm2(a), 1.0e-300_dp)) then      ! degenerate edge: fall back to a Cartesian axis
       imin = minloc(abs(n), 1)
       e = 0.0_dp;  e(imin) = 1.0_dp
       t1 = e - dot_product(e, n)*n
    end if
    t1 = t1/norm2(t1)
    t2(1) = n(2)*t1(3) - n(3)*t1(2)
    t2(2) = n(3)*t1(1) - n(1)*t1(3)
    t2(3) = n(1)*t1(2) - n(2)*t1(1)
  end subroutine tangents

  pure subroutine to_frame(wc, n, t1, t2, bn, q)
    real(dp), intent(in)  :: wc(NVAR), n(3), t1(3), t2(3), bn
    real(dp), intent(out) :: q(NVAR)
    real(dp) :: v(3), b(3)
    v = wc(2:4);  b = wc(6:8)
    q(1) = wc(1)
    q(2) = dot_product(v, n);  q(3) = dot_product(v, t1);  q(4) = dot_product(v, t2)
    q(5) = wc(5)
    q(6) = bn
    q(7) = dot_product(b, t1);  q(8) = dot_product(b, t2)
  end subroutine to_frame

  ! fast magnetosonic speed along the normal
  pure real(dp) function cfast(q)
    real(dp), intent(in) :: q(NVAR)
    real(dp) :: a2, b2, bn2, d
    a2 = gam*q(5)/q(1)
    b2 = (q(6)**2 + q(7)**2 + q(8)**2)/q(1)
    bn2 = q(6)**2/q(1)
    d = sqrt(max((a2 + b2)**2 - 4.0_dp*a2*bn2, 0.0_dp))
    cfast = sqrt(0.5_dp*(a2 + b2 + d))
  end function cfast

  ! physical flux in the face frame and the conserved state
  pure subroutine phys_flux(q, U, F)
    real(dp), intent(in)  :: q(NVAR)
    real(dp), intent(out) :: U(NVAR), F(NVAR)
    real(dp) :: pt, vb, e
    pt = q(5) + 0.5_dp*(q(6)**2 + q(7)**2 + q(8)**2)
    vb = q(2)*q(6) + q(3)*q(7) + q(4)*q(8)
    e  = q(5)/(gam - 1.0_dp) + 0.5_dp*q(1)*(q(2)**2 + q(3)**2 + q(4)**2) + 0.5_dp*(q(6)**2 + q(7)**2 + q(8)**2)
    U(1) = q(1);  U(2) = q(1)*q(2);  U(3) = q(1)*q(3);  U(4) = q(1)*q(4);  U(5) = e;  U(6:8) = q(6:8)
    F(1) = q(1)*q(2)
    F(2) = q(1)*q(2)**2 + pt - q(6)**2
    F(3) = q(1)*q(2)*q(3) - q(6)*q(7)
    F(4) = q(1)*q(2)*q(4) - q(6)*q(8)
    F(5) = (e + pt)*q(2) - q(6)*vb
    F(6) = 0.0_dp
    F(7) = q(7)*q(2) - q(6)*q(3)
    F(8) = q(8)*q(2) - q(6)*q(4)
  end subroutine phys_flux

  pure subroutine rusanov(ql, qr, F, sm)
    real(dp), intent(in)  :: ql(NVAR), qr(NVAR)
    real(dp), intent(out) :: F(NVAR), sm
    real(dp) :: UL(NVAR), UR(NVAR), FL(NVAR), FR(NVAR)
    call phys_flux(ql, UL, FL);  call phys_flux(qr, UR, FR)
    sm = max(abs(ql(2)) + cfast(ql), abs(qr(2)) + cfast(qr))
    F = 0.5_dp*(FL + FR) - 0.5_dp*sm*(UR - UL)
  end subroutine rusanov

  pure subroutine hll(ql, qr, F, sm)
    real(dp), intent(in)  :: ql(NVAR), qr(NVAR)
    real(dp), intent(out) :: F(NVAR), sm
    real(dp) :: UL(NVAR), UR(NVAR), FL(NVAR), FR(NVAR), SL, SR, cl, cr
    call phys_flux(ql, UL, FL);  call phys_flux(qr, UR, FR)
    cl = cfast(ql);  cr = cfast(qr)
    SL = min(ql(2), qr(2)) - max(cl, cr)
    SR = max(ql(2), qr(2)) + max(cl, cr)
    sm = max(abs(SL), abs(SR))
    if (SL >= 0.0_dp) then
       F = FL
    else if (SR <= 0.0_dp) then
       F = FR
    else
       F = (SR*FL - SL*FR + SL*SR*(UR - UL))/(SR - SL)
    end if
  end subroutine hll

  !--------------------------------------------------------------------
  ! HLLD (Miyoshi & Kusano, J. Comput. Phys. 208, 315 (2005))
  !--------------------------------------------------------------------
  pure subroutine hlld(ql, qr, F, sm)
    real(dp), intent(in)  :: ql(NVAR), qr(NVAR)
    real(dp), intent(out) :: F(NVAR), sm
    real(dp) :: UL(NVAR), UR(NVAR), FL(NVAR), FR(NVAR), SL, SR, SM_, cl, cr, bx, sgnb
    real(dp) :: ptL, ptR, pts, rL, rR, vnL, vnR
    real(dp) :: rsL, rsR, vsL, vsR, wsL, wsR, bysL, bysR, bzsL, bzsR, esL, esR, den, fac
    real(dp) :: SsL, SsR, sqL, sqR, vss, wss, byss, bzss, essL, essR
    real(dp) :: UsL(NVAR), UsR(NVAR), UssL(NVAR), UssR(NVAR)
    real(dp), parameter :: eps = 1.0e-12_dp

    call phys_flux(ql, UL, FL);  call phys_flux(qr, UR, FR)
    cl = cfast(ql);  cr = cfast(qr)
    vnL = ql(2);  vnR = qr(2);  rL = ql(1);  rR = qr(1)
    SL = min(vnL, vnR) - max(cl, cr)
    SR = max(vnL, vnR) + max(cl, cr)
    sm = max(abs(SL), abs(SR))
    if (SL >= 0.0_dp) then
       F = FL;  return
    end if
    if (SR <= 0.0_dp) then
       F = FR;  return
    end if
    bx = ql(6)
    ptL = ql(5) + 0.5_dp*(ql(6)**2 + ql(7)**2 + ql(8)**2)
    ptR = qr(5) + 0.5_dp*(qr(6)**2 + qr(7)**2 + qr(8)**2)
    den = (SR - vnR)*rR - (SL - vnL)*rL
    SM_ = ((SR - vnR)*rR*vnR - (SL - vnL)*rL*vnL - ptR + ptL)/den
    pts = ((SR - vnR)*rR*ptL - (SL - vnL)*rL*ptR + rL*rR*(SR - vnR)*(SL - vnL)*(vnR - vnL))/den

    ! left star state
    rsL = rL*(SL - vnL)/(SL - SM_)
    den = rL*(SL - vnL)*(SL - SM_) - bx**2
    if (abs(den) < eps*max(rL*(SL - vnL)**2, bx**2, eps)) then
       vsL = ql(3);  wsL = ql(4);  bysL = ql(7);  bzsL = ql(8)
    else
       fac = (SM_ - vnL)/den
       vsL = ql(3) - bx*ql(7)*fac
       wsL = ql(4) - bx*ql(8)*fac
       fac = (rL*(SL - vnL)**2 - bx**2)/den
       bysL = ql(7)*fac;  bzsL = ql(8)*fac
    end if
    esL = ((SL - vnL)*UL(5) - ptL*vnL + pts*SM_ + bx*((ql(2)*ql(6) + ql(3)*ql(7) + ql(4)*ql(8)) - (SM_*bx + vsL*bysL + wsL*bzsL)))/(SL - SM_)
    UsL = [rsL, rsL*SM_, rsL*vsL, rsL*wsL, esL, bx, bysL, bzsL]
    ! right star state
    rsR = rR*(SR - vnR)/(SR - SM_)
    den = rR*(SR - vnR)*(SR - SM_) - bx**2
    if (abs(den) < eps*max(rR*(SR - vnR)**2, bx**2, eps)) then
       vsR = qr(3);  wsR = qr(4);  bysR = qr(7);  bzsR = qr(8)
    else
       fac = (SM_ - vnR)/den
       vsR = qr(3) - bx*qr(7)*fac
       wsR = qr(4) - bx*qr(8)*fac
       fac = (rR*(SR - vnR)**2 - bx**2)/den
       bysR = qr(7)*fac;  bzsR = qr(8)*fac
    end if
    esR = ((SR - vnR)*UR(5) - ptR*vnR + pts*SM_ + bx*((qr(2)*qr(6) + qr(3)*qr(7) + qr(4)*qr(8)) - (SM_*bx + vsR*bysR + wsR*bzsR)))/(SR - SM_)
    UsR = [rsR, rsR*SM_, rsR*vsR, rsR*wsR, esR, bx, bysR, bzsR]

    sqL = sqrt(rsL);  sqR = sqrt(rsR)
    SsL = SM_ - abs(bx)/sqL
    SsR = SM_ + abs(bx)/sqR

    if (SsL >= 0.0_dp) then
       F = FL + SL*(UsL - UL);  return
    end if
    if (SsR <= 0.0_dp) then
       F = FR + SR*(UsR - UR);  return
    end if
    ! double-star states
    sgnb = sign(1.0_dp, bx)
    vss  = (sqL*vsL + sqR*vsR + (bysR - bysL)*sgnb)/(sqL + sqR)
    wss  = (sqL*wsL + sqR*wsR + (bzsR - bzsL)*sgnb)/(sqL + sqR)
    byss = (sqL*bysR + sqR*bysL + sqL*sqR*(vsR - vsL)*sgnb)/(sqL + sqR)
    bzss = (sqL*bzsR + sqR*bzsL + sqL*sqR*(wsR - wsL)*sgnb)/(sqL + sqR)
    essL = esL - sqL*((SM_*bx + vsL*bysL + wsL*bzsL) - (SM_*bx + vss*byss + wss*bzss))*sgnb
    essR = esR + sqR*((SM_*bx + vsR*bysR + wsR*bzsR) - (SM_*bx + vss*byss + wss*bzss))*sgnb
    UssL = [rsL, rsL*SM_, rsL*vss, rsL*wss, essL, bx, byss, bzss]
    UssR = [rsR, rsR*SM_, rsR*vss, rsR*wss, essR, bx, byss, bzss]
    if (SM_ >= 0.0_dp) then
       F = FL + SsL*UssL - (SsL - SL)*UsL - SL*UL
    else
       F = FR + SsR*UssR - (SsR - SR)*UsR - SR*UR
    end if
  end subroutine hlld

  ! maximum signal speed over the own cells (for the time step)
  real(dp) function max_signal(wp)
    real(dp), intent(in) :: wp(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: i, j, k
    real(dp) :: q(NVAR), cf, a2, b2
    max_signal = 0.0_dp
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             q = wp(:, i, j, k)
             a2 = gam*q(5)/q(1);  b2 = sum(q(6:8)**2)/q(1)
             cf = sqrt(a2 + b2)
             max_signal = max(max_signal, norm2(q(2:4)) + cf)
          end do
       end do
    end do
  end function max_signal

end module convection
