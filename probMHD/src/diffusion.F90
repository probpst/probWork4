!=======================================================================
! diffusion.F90 -- non-ideal terms, turbulence models and sources
!
!   Viscous stress     tau = rho nu_eff (grad v + grad v^T - 2/3 div v I)
!   Heat flux          q = -rho [chi_par b b + chi_perp (I - b b)] . grad T,  T = p/rho
!   Particle diffusion F_rho = -D_perp (I - b b) . grad rho   (anomalous transport)
!   Resistive terms    induction.F90 (eta J on the edges); the energy
!                      flux (eta J x B) . S_f here keeps the total
!                      energy conservative (Ohmic heating implied)
!   Face gradients: average of the Green-Gauss cell gradients with the
!   over-relaxed correction along the centre-to-centre vector.
!
!   TURBULENCE SWITCH  SW_Turb_MHD
!     0  none           (the ENO dissipation acts as an implicit LES)
!     1  MHD Smagorinsky (Mueller & Carati 2002 form)
!          nu_sgs  = (Cs   Delta)^2 |S|,     |S| = sqrt(2 S_ij S_ij)
!          eta_sgs = (Ceta Delta)^2 |J|,     Delta = V^(1/3)
!          chi_sgs = nu_sgs / Pr_sgs (isotropic)
!     2  dynamic Smagorinsky: reserved (stops with a message)
!     3  anomalous perpendicular transport (the transport-model view of
!        the micro-turbulence a fluid code cannot resolve):
!          chi_perp += chi_anom, D_perp += D_anom, nu += nu_anom
!          Anom_model 1: constants;  2: Bohm, chi = Bohm_coef T/(16 e B)
!   Sources: heating src_heat(x) (t >= heat_start), bremsstrahlung
!   radiation P_br = 5.35e-37 Z_eff n_e n_i sqrt(T_keV) [W/m3], gravity.
!=======================================================================
module diffusion
  use global_vars
  use mpi_vars
  use geometry, only: cross
  implicit none
  private
  public :: cell_gradients, transport_coefficients, diffusive_fluxes, add_sources, max_diffusivity, internal_energy_terms
  real(dp), allocatable, public :: flx_q(:,:,:,:)   ! heat-conduction (+ diffusing-mass enthalpy) flux through the faces

contains

  ! Green-Gauss gradients of rho, v, p, B (gradw(3, 1:8, cells 2-NG..n+NG-1))
  subroutine cell_gradients(wp)
    real(dp), intent(in) :: wp(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: i, j, k, q
    real(dp) :: g(3, NVAR), fv(NVAR)
    do k = 2-ngk, nz+ngk-1
       do j = 2-ngj, ny+ngj-1
          do i = 2-ngi, nx+ngi-1
             g = 0.0_dp
             fv = 0.5_dp*(wp(:,i,j,k) + wp(:,i+1,j,k)); do q = 1, NVAR; g(:,q) = g(:,q) + fv(q)*sf(:,1,i+1,j,k); end do
             fv = 0.5_dp*(wp(:,i,j,k) + wp(:,i-1,j,k)); do q = 1, NVAR; g(:,q) = g(:,q) - fv(q)*sf(:,1,i,j,k);   end do
             fv = 0.5_dp*(wp(:,i,j,k) + wp(:,i,j+1,k)); do q = 1, NVAR; g(:,q) = g(:,q) + fv(q)*sf(:,2,i,j+1,k); end do
             fv = 0.5_dp*(wp(:,i,j,k) + wp(:,i,j-1,k)); do q = 1, NVAR; g(:,q) = g(:,q) - fv(q)*sf(:,2,i,j,k);   end do
             fv = 0.5_dp*(wp(:,i,j,k) + wp(:,i,j,k+1)); do q = 1, NVAR; g(:,q) = g(:,q) + fv(q)*sf(:,3,i,j,k+1); end do
             fv = 0.5_dp*(wp(:,i,j,k) + wp(:,i,j,k-1)); do q = 1, NVAR; g(:,q) = g(:,q) - fv(q)*sf(:,3,i,j,k);   end do
             gradw(:, :, i, j, k) = g/vol(i, j, k)
          end do
       end do
    end do
  end subroutine cell_gradients

  !--------------------------------------------------------------------
  ! transport coefficients per cell (resistivity, viscosity, heat and
  ! particle diffusivities), including the turbulence model
  !--------------------------------------------------------------------
  subroutine transport_coefficients(wp)
    real(dp), intent(in) :: wp(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: i, j, k
    real(dp) :: T, delta, smag, g(3,3), sij, bmag, chiB, Tdim, Bdim
    do k = 2-ngk, nz+ngk-1
       do j = 2-ngj, ny+ngj-1
          do i = 2-ngi, nx+ngi-1
             T = wp(5,i,j,k)/wp(1,i,j,k)
             ! resistivity
             if (Spitzer == 1) then
                eta_cell(i,j,k) = min(eta_c*(max(T, 1.0e-30_dp)/T_ref_c)**(-1.5_dp), eta_c*eta_max_factor)
             else
                eta_cell(i,j,k) = eta_c
             end if
             nu_cell(i,j,k) = nu_c
             chi_par_cell(i,j,k) = chi_par_c
             chi_perp_cell(i,j,k) = chi_perp_c
             D_cell(i,j,k) = D_perp_c
             nu_sgs(i,j,k) = 0.0_dp;  eta_sgs(i,j,k) = 0.0_dp
             select case (SW_Turb_MHD)
             case (1)
                delta = vol(i,j,k)**(1.0_dp/3.0_dp)
                g = gradw(:, 2:4, i, j, k)          ! g(a, b) = d v_b / d x_a
                sij = 0.0_dp
                sij = sij + 2.0_dp*(g(1,1)**2 + g(2,2)**2 + g(3,3)**2) &
                          + (g(1,2) + g(2,1))**2 + (g(1,3) + g(3,1))**2 + (g(2,3) + g(3,2))**2
                smag = sqrt(sij)                      ! |S| = sqrt(2 S_ij S_ij)
                nu_sgs(i,j,k)  = (Cs_sgs*delta)**2*smag
                eta_sgs(i,j,k) = (Ceta_sgs*delta)**2*norm2(jc(:,i,j,k))
                nu_cell(i,j,k) = nu_cell(i,j,k) + nu_sgs(i,j,k)
                chi_par_cell(i,j,k)  = chi_par_cell(i,j,k)  + nu_sgs(i,j,k)/Pr_sgs
                chi_perp_cell(i,j,k) = chi_perp_cell(i,j,k) + nu_sgs(i,j,k)/Pr_sgs
             case (3)
                if (Anom_model == 2) then
                   ! Bohm: chi = T_e / (16 e B)  [m2/s] with T_e = T/2 ... convention: T_e = T_i = T_code*T0_K
                   Tdim = T*T0_K*KB/QE                      ! eV
                   Bdim = max(norm2(wp(6:8,i,j,k))*Bch, 1.0e-6_dp)
                   chiB = Bohm_coef*Tdim/(16.0_dp*Bdim)     ! m2/s
                   chiB = min(chiB/(Lch*vA0), anom_max_c)   ! code units
                   if (nondim_inputs) chiB = min(Bohm_coef*T/(16.0_dp*max(norm2(wp(6:8,i,j,k)), 1.0e-6_dp)), anom_max_c)
                   chi_perp_cell(i,j,k) = chi_perp_cell(i,j,k) + chiB
                   D_cell(i,j,k) = D_cell(i,j,k) + chiB*merge(D_anom_c/max(chi_anom_c, 1.0e-30_dp), 0.2_dp, chi_anom_c > 0.0_dp)
                else
                   chi_perp_cell(i,j,k) = chi_perp_cell(i,j,k) + chi_anom_c
                   D_cell(i,j,k) = D_cell(i,j,k) + D_anom_c
                   nu_cell(i,j,k) = nu_cell(i,j,k) + nu_anom_c
                end if
             end select
             chi_par_cell(i,j,k) = min(chi_par_cell(i,j,k), chi_ratio_max*max(chi_perp_cell(i,j,k), 1.0e-30_dp))
          end do
       end do
    end do
  end subroutine transport_coefficients

  ! largest diffusivity (for the time step)
  real(dp) function max_diffusivity()
    max_diffusivity = max(maxval(nu_cell(1:nx,1:ny,1:nz)), maxval(eta_cell(1:nx,1:ny,1:nz) + eta_sgs(1:nx,1:ny,1:nz)), &
                          maxval(chi_perp_cell(1:nx,1:ny,1:nz)), maxval(D_cell(1:nx,1:ny,1:nz)), &
                          maxval(chi_par_cell(1:nx,1:ny,1:nz)))
  end function max_diffusivity

  !--------------------------------------------------------------------
  ! diffusive fluxes added to flx (minus faces of cells 1..n+1 in each direction)
  !--------------------------------------------------------------------
  subroutine diffusive_fluxes(wp, fx)
    real(dp), intent(in)    :: wp(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(inout) :: fx(NVAR, 3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: i, j, k, c, im, jm, km, a, b
    real(dp) :: s(3), gf(3, NVAR), d(3), dn, wf(NVAR), nuf, chpar, chper, Df, etaf, T1, T2, gT(3), bf3(3), bhat(3), bmag
    real(dp) :: gv(3,3), divv, tau(3,3), tv(3), q(3), rhof, jf(3), ejb(3), gr(3), fr(3), bb(3), btg
    if (.not. allocated(flx_q)) then
       allocate(flx_q(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    end if
    flx_q = 0.0_dp
    do c = 1, 3
       if (c == 1 .and. ngx == 1) cycle
       if (c == 2 .and. ngy == 1) cycle
       if (c == 3 .and. ngz == 1) cycle
       do k = 1, nz + merge(1, 0, c == 3)
          do j = 1, ny + merge(1, 0, c == 2)
             do i = 1, nx + merge(1, 0, c == 1)
                im = i;  jm = j;  km = k
                if (c == 1) im = i - 1
                if (c == 2) jm = j - 1
                if (c == 3) km = k - 1
                s = sf(:, c, i, j, k)
                d = dvec(:, c, i, j, k)                       ! from cell (im,jm,km) to cell (i,j,k)
                wf = 0.5_dp*(wp(:,i,j,k) + wp(:,im,jm,km))
                ! face gradients: average + correction along d
                gf = 0.5_dp*(gradw(:,:,i,j,k) + gradw(:,:,im,jm,km))
                dn = dot_product(d, d)
                do a = 1, NVAR
                   gf(:,a) = gf(:,a) + ((wp(a,i,j,k) - wp(a,im,jm,km)) - dot_product(gf(:,a), d))*d/dn
                end do
                nuf   = 0.5_dp*(nu_cell(i,j,k) + nu_cell(im,jm,km))
                chpar = 0.5_dp*(chi_par_cell(i,j,k) + chi_par_cell(im,jm,km))
                chper = 0.5_dp*(chi_perp_cell(i,j,k) + chi_perp_cell(im,jm,km))
                Df    = 0.5_dp*(D_cell(i,j,k) + D_cell(im,jm,km))
                etaf  = 0.5_dp*(eta_cell(i,j,k) + eta_sgs(i,j,k) + eta_cell(im,jm,km) + eta_sgs(im,jm,km))
                rhof  = wf(1)
                bb = wf(6:8);  bmag = norm2(bb)
                if (bmag > 1.0e-12_dp) then
                   bhat = bb/bmag
                else
                   bhat = 0.0_dp
                end if
                ! ---- viscous stress ----
                if (nuf > 0.0_dp) then
                   gv = gf(:, 2:4)                            ! gv(a,b) = d v_b / d x_a
                   divv = gv(1,1) + gv(2,2) + gv(3,3)
                   do a = 1, 3
                      do b = 1, 3
                         tau(a,b) = rhof*nuf*(gv(a,b) + gv(b,a))
                      end do
                      tau(a,a) = tau(a,a) - rhof*nuf*2.0_dp/3.0_dp*divv
                   end do
                   tv = matmul(tau, s)                        ! tau . S  (momentum flux out of the minus cell: -tau.S)
                   fx(2:4, c, i, j, k) = fx(2:4, c, i, j, k) - tv
                   fx(5, c, i, j, k)   = fx(5, c, i, j, k)   - dot_product(tv, wf(2:4))
                end if
                ! ---- heat conduction: T = p / rho ----
                if (chpar > 0.0_dp .or. chper > 0.0_dp) then
                   ! grad T = (grad p - T grad rho) / rho
                   T1 = wf(5)/wf(1)
                   gT = (gf(:,5) - T1*gf(:,1))/wf(1)
                   ! exact face-normal temperature difference for the compact part
                   gT = gT + (((wp(5,i,j,k)/wp(1,i,j,k)) - (wp(5,im,jm,km)/wp(1,im,jm,km))) - dot_product(gT, d))*d/dn
                   btg = dot_product(bhat, gT)
                   q = -rhof*(chpar*btg*bhat + chper*(gT - btg*bhat))
                   fx(5, c, i, j, k) = fx(5, c, i, j, k) + dot_product(q, s)
                   flx_q(c, i, j, k) = flx_q(c, i, j, k) + dot_product(q, s)
                end if
                ! ---- particle diffusion (perpendicular) ----
                if (Df > 0.0_dp) then
                   gr = gf(:,1)
                   btg = dot_product(bhat, gr)
                   fr = -Df*(gr - btg*bhat)
                   fx(1, c, i, j, k) = fx(1, c, i, j, k) + dot_product(fr, s)
                   ! the diffusing mass carries its momentum and enthalpy (consistency of the conserved variables)
                   fx(2:4, c, i, j, k) = fx(2:4, c, i, j, k) + dot_product(fr, s)*wf(2:4)
                   fx(5, c, i, j, k)   = fx(5, c, i, j, k)   + dot_product(fr, s)*(0.5_dp*sum(wf(2:4)**2) + gam/(gam-1.0_dp)*wf(5)/wf(1))
                   flx_q(c, i, j, k)   = flx_q(c, i, j, k)   + dot_product(fr, s)*wf(5)/wf(1)/(gam-1.0_dp)
                end if
                ! ---- resistive Poynting flux: (eta J x B) . S ----
                if (etaf > 0.0_dp) then
                   jf = 0.5_dp*(jc(:,i,j,k) + jc(:,im,jm,km))
                   ejb = etaf*cross(jf, bb)
                   fx(5, c, i, j, k) = fx(5, c, i, j, k) + dot_product(ejb, s)
                end if
             end do
          end do
       end do
    end do
  end subroutine diffusive_fluxes

  ! cell heating of the internal-energy equation: Ohmic eta J^2 and viscous dissipation
  subroutine internal_energy_terms(wp, rr)
    real(dp), intent(in)    :: wp(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(inout) :: rr(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: i, j, k, a, b
    real(dp) :: g(3,3), divv, phi, nueff
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             rr(i,j,k) = rr(i,j,k) + (eta_cell(i,j,k) + eta_sgs(i,j,k))*sum(jc(:,i,j,k)**2)
             nueff = nu_cell(i,j,k)
             if (nueff > 0.0_dp) then
                g = gradw(:, 2:4, i, j, k)
                divv = g(1,1) + g(2,2) + g(3,3)
                phi = 0.0_dp
                do a = 1, 3
                   do b = 1, 3
                      phi = phi + 0.5_dp*(g(a,b) + g(b,a))**2
                   end do
                end do
                phi = phi - 2.0_dp/3.0_dp*divv**2
                rr(i,j,k) = rr(i,j,k) + wp(1,i,j,k)*nueff*max(phi, 0.0_dp)
             end if
          end do
       end do
    end do
  end subroutine internal_energy_terms

  !--------------------------------------------------------------------
  ! volumetric sources: heating, radiation, gravity
  !--------------------------------------------------------------------
  subroutine add_sources(wp, rr, rre, tnow)
    real(dp), intent(in)    :: wp(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(inout) :: rr(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(inout) :: rre(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)    :: tnow
    integer :: i, j, k
    real(dp) :: n_m3, T_keV, pbr, T
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             if (P_heat_c /= 0.0_dp .and. tnow >= heat_start_c) then
                rr(5,i,j,k) = rr(5,i,j,k) + src_heat(i,j,k);  rre(i,j,k) = rre(i,j,k) + src_heat(i,j,k)
             end if
             if (SW_radiation == 1) then
                T = wp(5,i,j,k)/wp(1,i,j,k)
                if (nondim_inputs) then
                   pbr = Z_eff*wp(1,i,j,k)**2*sqrt(max(T, 0.0_dp))*1.0e-3_dp   ! code-unit proxy
                else
                   n_m3 = wp(1,i,j,k)*Rhoch/mion
                   T_keV = T*T0_K*KB/QE/1000.0_dp
                   pbr = 5.35e-37_dp*Z_eff*n_m3**2*sqrt(max(T_keV, 0.0_dp))    ! W/m3
                   pbr = pbr*Lch/(p0*vA0)
                end if
                pbr = min(pbr, 0.5_dp*wp(5,i,j,k)/(gam-1.0_dp)/max(dt, 1.0e-30_dp))
                rr(5,i,j,k) = rr(5,i,j,k) - pbr;  rre(i,j,k) = rre(i,j,k) - pbr
             end if
             if (any(grav_c /= 0.0_dp)) then
                rr(2:4,i,j,k) = rr(2:4,i,j,k) + wp(1,i,j,k)*grav_c
                rr(5,i,j,k)   = rr(5,i,j,k)   + wp(1,i,j,k)*dot_product(grav_c, wp(2:4,i,j,k))
             end if
          end do
       end do
    end do
  end subroutine add_sources

end module diffusion
