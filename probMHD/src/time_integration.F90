!=======================================================================
! time_integration.F90 -- right-hand side assembly and the explicit
!   SSP Runge-Kutta time stepping (Time_scheme 1 Euler, 2 SSP-RK2 Heun
!   as in probCFD, 3 SSP-RK3 Shu-Osher).
!
!   Every stage: primitive variables and cell-centred B (from the faces),
!   gradients / current / transport coefficients, Riemann fluxes, non-
!   ideal fluxes, sources, edge EMFs; then the conserved cells and the
!   face fluxes advance together.
!=======================================================================
module time_integration
  use global_vars
  use mpi_vars
  use boundary
  use convection
  use induction
  use diffusion
  use equilibrium, only: cell_B_from_faces
  implicit none
  private
  public :: compute_rhs, rk_step, compute_dt, refresh_state

  real(dp), allocatable :: smax(:,:,:,:), u1(:,:,:,:), b1(:,:,:,:), u2(:,:,:,:), b2(:,:,:,:)
  real(dp), allocatable :: rhs_eq(:,:,:,:), emf_eq(:,:,:,:), rhs_e_eq(:,:,:)
  real(dp), allocatable :: e1int(:,:,:,:), e2int(:,:,:,:)
  logical :: eq_residual_stored = .false.

contains

  ! ghost layers of (uu, bb): exchange + physical BCs
  subroutine refresh_state(uu, bb, ee)
    real(dp), intent(inout) :: uu(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(inout) :: bb(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(inout), optional :: ee(1, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    call exchange(bb, 3)
    call apply_bc_faces(bb)
    call exchange(uu, NVAR)
    call apply_bc(uu, bb)
    if (present(ee) .and. SW_Energy >= 2) then
       call exchange(ee, 1)
       call apply_bc_scalar(ee, einteq)
    end if
  end subroutine refresh_state

  ! SW_Energy >= 2: the pressure from the internal energy (or, dual, from the total energy where the
  ! thermal fraction is large); the total energy is then re-synchronised with it
  subroutine sync_energy(uu, ee, wp)
    real(dp), intent(inout) :: uu(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(inout) :: ee(1, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(inout) :: wp(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: i, j, k
    real(dp) :: ek, eb, eth
    do k = 1-ngk, nz+ngk
       do j = 1-ngj, ny+ngj
          do i = 1-ngi, nx+ngi
             ek = 0.5_dp*uu(1,i,j,k)*sum(wp(2:4,i,j,k)**2);  eb = 0.5_dp*sum(wp(6:8,i,j,k)**2)
             eth = uu(5,i,j,k) - ek - eb
             if (SW_Energy == 3 .and. eth > dual_eps*uu(5,i,j,k)) then
                ee(1,i,j,k) = eth
             end if
             ee(1,i,j,k) = max(ee(1,i,j,k), p_floor/(gam - 1.0_dp))
             wp(5,i,j,k) = (gam - 1.0_dp)*ee(1,i,j,k)
             uu(5,i,j,k) = ee(1,i,j,k) + ek + eb
          end do
       end do
    end do
  end subroutine sync_energy

  !--------------------------------------------------------------------
  subroutine compute_rhs(uu, bb, ee, tnow)
    real(dp), intent(inout) :: uu(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)    :: bb(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(inout) :: ee(1, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)    :: tnow
    integer :: i, j, k
    logical :: resist
    real(dp) :: g1, g2
    if (.not. allocated(smax)) allocate(smax(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    ! primitive state
    do k = 1-ngk, nz+ngk
       do j = 1-ngj, ny+ngj
          do i = 1-ngi, nx+ngi
             call cons2prim(uu(:,i,j,k), w(:,i,j,k))
          end do
       end do
    end do
    if (CT_scheme >= 1) then
       call cell_B_from_faces(bb, w)
       uu(6:8, :, :, :) = w(6:8, :, :, :)
       ! the ghost cells carry the (rotated) face-averaged field of their images, and the mirrored
       ! interior field on the physical boundaries (a ghost cell averaged from its own extrapolated
       ! faces would differ from the mirror by the grid curvature and radiate waves from the wall)
       call exchange(uu, NVAR)
       call apply_bc(uu, bb)
       ! pressure with the face-averaged B
       do k = 1-ngk, nz+ngk
          do j = 1-ngj, ny+ngj
             do i = 1-ngi, nx+ngi
                call cons2prim(uu(:,i,j,k), w(:,i,j,k))
             end do
          end do
       end do
    end if
    if (SW_Energy >= 2) call sync_energy(uu, ee, w)
    resist = .false.
    if (diffusion_on) then
       call cell_gradients(w)
       call cell_current(w)
       call transport_coefficients(w)
       resist = (eta_c > 0.0_dp .or. SW_Turb_MHD == 1 .or. Spitzer == 1)
    end if
    call convective_fluxes(w, bb, smax)
    if (diffusion_on) call diffusive_fluxes(w, flx)
    ! cell balance
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             rhs(:,i,j,k) = -( flx(:,1,i+1,j,k) - flx(:,1,i,j,k) + flx(:,2,i,j+1,k) - flx(:,2,i,j,k) &
                             + flx(:,3,i,j,k+1) - flx(:,3,i,j,k) )/vol(i,j,k)
          end do
       end do
    end do
    if (SW_Energy >= 2) then
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                rhs_e(i,j,k) = -( flx_e(1,i+1,j,k) - flx_e(1,i,j,k) + flx_e(2,i,j+1,k) - flx_e(2,i,j,k) &
                                + flx_e(3,i,j,k+1) - flx_e(3,i,j,k) )/vol(i,j,k) &
                               - w(5,i,j,k)*( vflux(1,i+1,j,k) - vflux(1,i,j,k) + vflux(2,i,j+1,k) - vflux(2,i,j,k) &
                                            + vflux(3,i,j,k+1) - vflux(3,i,j,k) )/vol(i,j,k)
                if (diffusion_on) rhs_e(i,j,k) = rhs_e(i,j,k) - ( flx_q(1,i+1,j,k) - flx_q(1,i,j,k) + flx_q(2,i,j+1,k) - flx_q(2,i,j,k) &
                                                                + flx_q(3,i,j,k+1) - flx_q(3,i,j,k) )/vol(i,j,k)
             end do
          end do
       end do
       if (diffusion_on) call internal_energy_terms(w, rhs_e)
    end if
    if (debug_balance) call energy_balance_report
    call add_sources(w, rhs, rhs_e, tnow)
    if (CT_scheme >= 1) call edge_emf(resist)
    ! well-balancing: the discrete residual of the initial equilibrium is subtracted, so the
    ! equilibrium is an exact fixed point of the scheme (SW_Eq_Subtract = 1)
    if (SW_Eq_Subtract == 1) then
       if (.not. eq_residual_stored) then
          allocate(rhs_eq(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), emf_eq(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
          allocate(rhs_e_eq(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
          rhs_eq = rhs;  emf_eq = emf;  rhs_e_eq = rhs_e;  eq_residual_stored = .true.
          g1 = gmax(maxval(abs(rhs(2:4,1:nx,1:ny,1:nz))));  g2 = gmax(maxval(abs(emf(:,1:nx+1,1:ny+1,1:nz+1))))   ! collective: every rank
          if (myrank == 0) write(*,'(A,ES10.3,A,ES10.3)') '  equilibrium subtraction: max |rhs_eq| (momentum) ', g1, '  max |emf_eq| ', g2
       end if
       rhs = rhs - rhs_eq
       emf = emf - emf_eq
       rhs_e = rhs_e - rhs_e_eq
    end if
  end subroutine compute_rhs

  ! debugging aid: net energy rate of the own cells against the boundary / seam fluxes
  subroutine energy_balance_report
    integer :: i, j, k
    real(dp) :: tot, fi1, fi2, fj1, fj2, fk1, fk2
    tot = 0.0_dp;  fi1 = 0.0_dp;  fi2 = 0.0_dp;  fj1 = 0.0_dp;  fj2 = 0.0_dp;  fk1 = 0.0_dp;  fk2 = 0.0_dp
    do k = 1, nz; do j = 1, ny; do i = 1, nx
       tot = tot + rhs(5,i,j,k)*vol(i,j,k)
    end do; end do; end do
    do k = 1, nz; do j = 1, ny
       fi1 = fi1 + flx(5,1,1,j,k);  fi2 = fi2 + flx(5,1,nx+1,j,k)
    end do; end do
    do k = 1, nz; do i = 1, nx
       fj1 = fj1 + flx(5,2,i,1,k);  fj2 = fj2 + flx(5,2,i,ny+1,k)
    end do; end do
    do j = 1, ny; do i = 1, nx
       fk1 = fk1 + flx(5,3,i,j,1);  fk2 = fk2 + flx(5,3,i,j,nz+1)
    end do; end do
    tot = gsum(tot);  fi1 = gsum(fi1);  fi2 = gsum(fi2);  fj1 = gsum(fj1);  fj2 = gsum(fj2);  fk1 = gsum(fk1);  fk2 = gsum(fk2)
    if (myrank == 0) write(*,'(A,ES11.3,A,2ES11.3,A,2ES11.3,A,2ES11.3)') '  [balance] sum rhs_E V = ', tot, &
         '  E-flux i-faces (1, n+1): ', fi1, fi2, '  j: ', fj1, fj2, '  k: ', fk1, fk2
  end subroutine energy_balance_report

  !--------------------------------------------------------------------
  ! time step from the signal speeds of the last RHS evaluation
  !--------------------------------------------------------------------
  subroutine compute_dt
    integer :: i, j, k, c
    real(dp) :: sc, dtc, dtd, dmax, h2, dloc
    if (dt_dim > 0.0_dp) then
       dt = dt_dim/t0
       call cfl_of_dt
       return
    end if
    dtc = huge(1.0_dp)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             sc = 0.0_dp
             if (ngx > 1) sc = sc + max(smax(1,i,j,k), smax(1,i+1,j,k))/hc(1,i,j,k)
             if (ngy > 1) sc = sc + max(smax(2,i,j,k), smax(2,i,j+1,k))/hc(2,i,j,k)
             if (ngz > 1) sc = sc + max(smax(3,i,j,k), smax(3,i,j,k+1))/hc(3,i,j,k)
             if (sc > 0.0_dp) dtc = min(dtc, cfl_max/sc)
          end do
       end do
    end do
    dtd = huge(1.0_dp)
    if (diffusion_on) then
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                dloc = max(nu_cell(i,j,k), eta_cell(i,j,k) + eta_sgs(i,j,k), chi_par_cell(i,j,k), chi_perp_cell(i,j,k), D_cell(i,j,k))
                if (dloc <= 0.0_dp) cycle
                h2 = 0.0_dp
                do c = 1, 3
                   if ((c == 1 .and. ngx == 1) .or. (c == 2 .and. ngy == 1) .or. (c == 3 .and. ngz == 1)) cycle
                   h2 = h2 + 1.0_dp/hc(c,i,j,k)**2
                end do
                dtd = min(dtd, vn_max/(2.0_dp*dloc*h2))
             end do
          end do
       end do
    end if
    dt = gmin(min(dtc, dtd))
    cfl_now = cfl_max*dt/gmin(dtc)
  end subroutine compute_dt

  subroutine cfl_of_dt
    integer :: i, j, k
    real(dp) :: sc, cmax
    cmax = 0.0_dp
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             sc = 0.0_dp
             if (ngx > 1) sc = sc + max(smax(1,i,j,k), smax(1,i+1,j,k))/hc(1,i,j,k)
             if (ngy > 1) sc = sc + max(smax(2,i,j,k), smax(2,i,j+1,k))/hc(2,i,j,k)
             if (ngz > 1) sc = sc + max(smax(3,i,j,k), smax(3,i,j,k+1))/hc(3,i,j,k)
             cmax = max(cmax, sc*dt)
          end do
       end do
    end do
    cfl_now = gmax(cmax)
  end subroutine cfl_of_dt

  !--------------------------------------------------------------------
  ! one time step (dt known); on entry rhs/emf hold L(U^n) from compute_rhs
  !--------------------------------------------------------------------
  subroutine rk_step
    integer :: i, j, k
    if (.not. allocated(u1)) then
       allocate(u1(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), b1(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
       allocate(u2(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), b2(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
       allocate(e1int(1, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), e2int(1, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
       u1 = 0.0_dp;  b1 = 0.0_dp;  u2 = 0.0_dp;  b2 = 0.0_dp;  e1int = 0.0_dp;  e2int = 0.0_dp
    end if
    u0 = u;  bf0 = bf;  eint0 = eint
    select case (Time_scheme)
    case (1)
       call stage(u0, bf0, eint0, u, bf, eint, dt)
       call refresh_state(u, bf, eint)
    case (2)
       call stage(u0, bf0, eint0, u1, b1, e1int, dt)
       call refresh_state(u1, b1, e1int)
       call compute_rhs(u1, b1, e1int, t + dt)
       call stage(u1, b1, e1int, u2, b2, e2int, dt)
       u(:,1:nx,1:ny,1:nz) = 0.5_dp*(u0(:,1:nx,1:ny,1:nz) + u2(:,1:nx,1:ny,1:nz))
       bf(:,1:nx+1,1:ny+1,1:nz+1) = 0.5_dp*(bf0(:,1:nx+1,1:ny+1,1:nz+1) + b2(:,1:nx+1,1:ny+1,1:nz+1))
       eint(:,1:nx,1:ny,1:nz) = 0.5_dp*(eint0(:,1:nx,1:ny,1:nz) + e2int(:,1:nx,1:ny,1:nz))
       call refresh_state(u, bf, eint)
    case default
       call stage(u0, bf0, eint0, u1, b1, e1int, dt)
       call refresh_state(u1, b1, e1int)
       call compute_rhs(u1, b1, e1int, t + dt)
       call stage(u1, b1, e1int, u2, b2, e2int, dt)
       u2(:,1:nx,1:ny,1:nz) = 0.75_dp*u0(:,1:nx,1:ny,1:nz) + 0.25_dp*u2(:,1:nx,1:ny,1:nz)
       b2(:,1:nx+1,1:ny+1,1:nz+1) = 0.75_dp*bf0(:,1:nx+1,1:ny+1,1:nz+1) + 0.25_dp*b2(:,1:nx+1,1:ny+1,1:nz+1)
       e2int(:,1:nx,1:ny,1:nz) = 0.75_dp*eint0(:,1:nx,1:ny,1:nz) + 0.25_dp*e2int(:,1:nx,1:ny,1:nz)
       call refresh_state(u2, b2, e2int)
       call compute_rhs(u2, b2, e2int, t + 0.5_dp*dt)
       call stage(u2, b2, e2int, u, bf, eint, dt)
       u(:,1:nx,1:ny,1:nz) = (u0(:,1:nx,1:ny,1:nz) + 2.0_dp*u(:,1:nx,1:ny,1:nz))/3.0_dp
       bf(:,1:nx+1,1:ny+1,1:nz+1) = (bf0(:,1:nx+1,1:ny+1,1:nz+1) + 2.0_dp*bf(:,1:nx+1,1:ny+1,1:nz+1))/3.0_dp
       eint(:,1:nx,1:ny,1:nz) = (eint0(:,1:nx,1:ny,1:nz) + 2.0_dp*eint(:,1:nx,1:ny,1:nz))/3.0_dp
       call refresh_state(u, bf, eint)
    end select
    call floor_state(u)
  end subroutine rk_step

  ! Euler stage: (un, bn, en) + dtt * (rhs, -curl E, rhs_e) -> (uo, bo, eo) on the own cells / faces
  subroutine stage(un, bn, en, uo, bo, eo, dtt)
    real(dp), intent(in)    :: un(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), bn(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)    :: en(1, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(inout) :: uo(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk), bo(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(inout) :: eo(1, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)    :: dtt
    uo(:,1:nx,1:ny,1:nz) = un(:,1:nx,1:ny,1:nz) + dtt*rhs(:,1:nx,1:ny,1:nz)
    if (SW_Energy >= 2) eo(1,1:nx,1:ny,1:nz) = en(1,1:nx,1:ny,1:nz) + dtt*rhs_e(1:nx,1:ny,1:nz)
    if (CT_scheme >= 1) then
       call update_faces(bo, bn, dtt)
    else
       bo = bn
    end if
  end subroutine stage

  ! density / pressure floors on the own cells
  subroutine floor_state(uu)
    real(dp), intent(inout) :: uu(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: i, j, k
    real(dp) :: wc(NVAR)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             if (uu(1,i,j,k) < rho_floor) uu(1,i,j,k) = rho_floor
             call cons2prim(uu(:,i,j,k), wc)
             if (wc(5) <= p_floor*1.0000001_dp) then
                wc(5) = p_floor
                call prim2cons(wc, uu(:,i,j,k))
             end if
          end do
       end do
    end do
  end subroutine floor_state

end module time_integration
