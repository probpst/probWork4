!=======================================================================
! induction.F90 -- constrained transport of the face magnetic fluxes
!
!   d Phi_f / dt = - circulation of E around the face,
!   E on the edges (E . L_e) from
!     (a) the face electric fields of the Riemann solver, averaged over
!         the four faces that share the edge with the upwind corrections
!         of Gardiner & Stone (2005) (CT_scheme = 1; = 2: the plain
!         Balsara & Spicer (1999) average);
!     (b) the resistive part eta J . L_e with the cell current (Green-
!         Gauss curl of the cell-centred B) averaged to the edge, plus
!         the sub-grid resistivity of the turbulence model.
!   The sum of the six face fluxes of every cell is then conserved to
!   round-off: div B stays at its initial (zero) value.
!   Edges on conducting-wall / fixed boundaries carry E = 0.
!=======================================================================
module induction
  use global_vars
  use mpi_vars
  use geometry, only: cross
  use convection, only: flx
  implicit none
  private
  public :: edge_emf, update_faces, cell_current, divergence_B

contains

  ! cell-centred current J = curl B = (1/V) sum_f S_f x B_f (B_f = face average of the cells)
  subroutine cell_current(wp)
    real(dp), intent(in) :: wp(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: i, j, k
    real(dp) :: s(3), bfc(3), cur(3)
    do k = 2-ngk, nz+ngk-1
       do j = 2-ngj, ny+ngj-1
          do i = 2-ngi, nx+ngi-1
             cur = 0.0_dp
             ! minus faces (outward normal = -S), plus faces (+S)
             bfc = 0.5_dp*(wp(6:8,i,j,k) + wp(6:8,i-1,j,k));  cur = cur - cross(sf(:,1,i,j,k), bfc)
             bfc = 0.5_dp*(wp(6:8,i,j,k) + wp(6:8,i+1,j,k));  cur = cur + cross(sf(:,1,i+1,j,k), bfc)
             bfc = 0.5_dp*(wp(6:8,i,j,k) + wp(6:8,i,j-1,k));  cur = cur - cross(sf(:,2,i,j,k), bfc)
             bfc = 0.5_dp*(wp(6:8,i,j,k) + wp(6:8,i,j+1,k));  cur = cur + cross(sf(:,2,i,j+1,k), bfc)
             bfc = 0.5_dp*(wp(6:8,i,j,k) + wp(6:8,i,j,k-1));  cur = cur - cross(sf(:,3,i,j,k), bfc)
             bfc = 0.5_dp*(wp(6:8,i,j,k) + wp(6:8,i,j,k+1));  cur = cur + cross(sf(:,3,i,j,k+1), bfc)
             jc(:, i, j, k) = cur/vol(i, j, k)
          end do
       end do
    end do
  end subroutine cell_current

  !--------------------------------------------------------------------
  ! edge EMFs: emf(c,i,j,k) = E . L_e for the edge along xi_c at the
  ! (minus, minus) corner of cell (i,j,k) in the other two directions.
  !   edge 1 at (i,j,k): shared by faces 2 at (i,j,k), (i,j,k-1) and faces 3 at (i,j,k), (i,j-1,k)
  !   edge 2 at (i,j,k): faces 3 at (i,j,k), (i-1,j,k) and faces 1 at (i,j,k), (i,j,k-1)
  !   edge 3 at (i,j,k): faces 1 at (i,j,k), (i,j-1,k) and faces 2 at (i,j,k), (i-1,j,k)
  !--------------------------------------------------------------------
  subroutine edge_emf(resistive)
    logical, intent(in) :: resistive
    integer :: i, j, k, iu, ju, ku
    real(dp) :: e(3), jl(3), etae, d1a(3), d1b(3), d2a(3), d2b(3), d3a(3), d3b(3)
    real(dp), allocatable, save :: er(:,:,:,:)      ! cell reference field -v x B
    if (.not. allocated(er)) allocate(er(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    do k = 0, nz+1
       do j = 0, ny+1
          do i = 0, nx+1
             er(:,i,j,k) = -cross(w(2:4,i,j,k), w(6:8,i,j,k))
          end do
       end do
    end do
    do k = 1, nz+1
       do j = 1, ny+1
          do i = 1, nx+1
             ! ---------------- edge 1 (along xi_1 at j-1/2, k-1/2) ----------------
             e = 0.25_dp*(ef(:,2,i,j,k) + ef(:,2,i,j,k-1) + ef(:,3,i,j,k) + ef(:,3,i,j-1,k))
             if (CT_scheme == 1) then
                ! Gardiner-Stone upwind corrections (index space, unit spacing)
                ju = upwind_cell(flx(1,2,i,j,k), j-1, j);    d3a = 2.0_dp*(er(:,i,ju,k) - ef(:,3,i,ju,k))
                if (flx(1,2,i,j,k) == 0.0_dp) d3a = (er(:,i,j-1,k) - ef(:,3,i,j-1,k)) + (er(:,i,j,k) - ef(:,3,i,j,k))
                ju = upwind_cell(flx(1,2,i,j,k-1), j-1, j);  d3b = 2.0_dp*(ef(:,3,i,ju,k) - er(:,i,ju,k-1))
                if (flx(1,2,i,j,k-1) == 0.0_dp) d3b = (ef(:,3,i,j-1,k) - er(:,i,j-1,k-1)) + (ef(:,3,i,j,k) - er(:,i,j,k-1))
                ku = upwind_cell(flx(1,3,i,j,k), k-1, k);    d2a = 2.0_dp*(er(:,i,j,ku) - ef(:,2,i,j,ku))
                if (flx(1,3,i,j,k) == 0.0_dp) d2a = (er(:,i,j,k-1) - ef(:,2,i,j,k-1)) + (er(:,i,j,k) - ef(:,2,i,j,k))
                ku = upwind_cell(flx(1,3,i,j-1,k), k-1, k);  d2b = 2.0_dp*(ef(:,2,i,j,ku) - er(:,i,j-1,ku))
                if (flx(1,3,i,j-1,k) == 0.0_dp) d2b = (ef(:,2,i,j,k-1) - er(:,i,j-1,k-1)) + (ef(:,2,i,j,k) - er(:,i,j-1,k))
                e = e - 0.125_dp*(d3a - d3b) - 0.125_dp*(d2a - d2b)
             end if
             emf(1, i, j, k) = dot_product(e, le(:,1,i,j,k))
             ! ---------------- edge 2 (along xi_2 at k-1/2, i-1/2) ----------------
             e = 0.25_dp*(ef(:,3,i,j,k) + ef(:,3,i-1,j,k) + ef(:,1,i,j,k) + ef(:,1,i,j,k-1))
             if (CT_scheme == 1) then
                ku = upwind_cell(flx(1,3,i,j,k), k-1, k);    d1a = 2.0_dp*(er(:,i,j,ku) - ef(:,1,i,j,ku))
                if (flx(1,3,i,j,k) == 0.0_dp) d1a = (er(:,i,j,k-1) - ef(:,1,i,j,k-1)) + (er(:,i,j,k) - ef(:,1,i,j,k))
                ku = upwind_cell(flx(1,3,i-1,j,k), k-1, k);  d1b = 2.0_dp*(ef(:,1,i,j,ku) - er(:,i-1,j,ku))
                if (flx(1,3,i-1,j,k) == 0.0_dp) d1b = (ef(:,1,i,j,k-1) - er(:,i-1,j,k-1)) + (ef(:,1,i,j,k) - er(:,i-1,j,k))
                iu = upwind_cell(flx(1,1,i,j,k), i-1, i);    d3a = 2.0_dp*(er(:,iu,j,k) - ef(:,3,iu,j,k))
                if (flx(1,1,i,j,k) == 0.0_dp) d3a = (er(:,i-1,j,k) - ef(:,3,i-1,j,k)) + (er(:,i,j,k) - ef(:,3,i,j,k))
                iu = upwind_cell(flx(1,1,i,j,k-1), i-1, i);  d3b = 2.0_dp*(ef(:,3,iu,j,k) - er(:,iu,j,k-1))
                if (flx(1,1,i,j,k-1) == 0.0_dp) d3b = (ef(:,3,i-1,j,k) - er(:,i-1,j,k-1)) + (ef(:,3,i,j,k) - er(:,i,j,k-1))
                e = e - 0.125_dp*(d1a - d1b) - 0.125_dp*(d3a - d3b)
             end if
             emf(2, i, j, k) = dot_product(e, le(:,2,i,j,k))
             ! ---------------- edge 3 (along xi_3 at i-1/2, j-1/2) ----------------
             e = 0.25_dp*(ef(:,1,i,j,k) + ef(:,1,i,j-1,k) + ef(:,2,i,j,k) + ef(:,2,i-1,j,k))
             if (CT_scheme == 1) then
                iu = upwind_cell(flx(1,1,i,j,k), i-1, i);    d2a = 2.0_dp*(er(:,iu,j,k) - ef(:,2,iu,j,k))
                if (flx(1,1,i,j,k) == 0.0_dp) d2a = (er(:,i-1,j,k) - ef(:,2,i-1,j,k)) + (er(:,i,j,k) - ef(:,2,i,j,k))
                iu = upwind_cell(flx(1,1,i,j-1,k), i-1, i);  d2b = 2.0_dp*(ef(:,2,iu,j,k) - er(:,iu,j-1,k))
                if (flx(1,1,i,j-1,k) == 0.0_dp) d2b = (ef(:,2,i-1,j,k) - er(:,i-1,j-1,k)) + (ef(:,2,i,j,k) - er(:,i,j-1,k))
                ju = upwind_cell(flx(1,2,i,j,k), j-1, j);    d1a = 2.0_dp*(er(:,i,ju,k) - ef(:,1,i,ju,k))
                if (flx(1,2,i,j,k) == 0.0_dp) d1a = (er(:,i,j-1,k) - ef(:,1,i,j-1,k)) + (er(:,i,j,k) - ef(:,1,i,j,k))
                ju = upwind_cell(flx(1,2,i-1,j,k), j-1, j);  d1b = 2.0_dp*(ef(:,1,i,ju,k) - er(:,i-1,ju,k))
                if (flx(1,2,i-1,j,k) == 0.0_dp) d1b = (ef(:,1,i,j-1,k) - er(:,i-1,j-1,k)) + (ef(:,1,i,j,k) - er(:,i-1,j,k))
                e = e - 0.125_dp*(d2a - d2b) - 0.125_dp*(d1a - d1b)
             end if
             emf(3, i, j, k) = dot_product(e, le(:,3,i,j,k))
             if (resistive) then
                jl = 0.25_dp*(jc(:,i,j,k) + jc(:,i,j-1,k) + jc(:,i,j,k-1) + jc(:,i,j-1,k-1))
                etae = 0.25_dp*(eta_tot(i,j,k) + eta_tot(i,j-1,k) + eta_tot(i,j,k-1) + eta_tot(i,j-1,k-1))
                emf(1, i, j, k) = emf(1, i, j, k) + etae*dot_product(jl, le(:,1,i,j,k))
                jl = 0.25_dp*(jc(:,i,j,k) + jc(:,i-1,j,k) + jc(:,i,j,k-1) + jc(:,i-1,j,k-1))
                etae = 0.25_dp*(eta_tot(i,j,k) + eta_tot(i-1,j,k) + eta_tot(i,j,k-1) + eta_tot(i-1,j,k-1))
                emf(2, i, j, k) = emf(2, i, j, k) + etae*dot_product(jl, le(:,2,i,j,k))
                jl = 0.25_dp*(jc(:,i,j,k) + jc(:,i-1,j,k) + jc(:,i,j-1,k) + jc(:,i-1,j-1,k))
                etae = 0.25_dp*(eta_tot(i,j,k) + eta_tot(i-1,j,k) + eta_tot(i,j-1,k) + eta_tot(i-1,j-1,k))
                emf(3, i, j, k) = emf(3, i, j, k) + etae*dot_product(jl, le(:,3,i,j,k))
             end if
          end do
       end do
    end do
    call wall_edges
  end subroutine edge_emf

  pure integer function upwind_cell(mflux, iminus, iplus)
    real(dp), intent(in) :: mflux
    integer,  intent(in) :: iminus, iplus
    if (mflux > 0.0_dp) then
       upwind_cell = iminus
    else
       upwind_cell = iplus
    end if
  end function upwind_cell

  pure real(dp) function eta_tot(i, j, k)
    integer, intent(in) :: i, j, k
    eta_tot = eta_cell(i, j, k) + eta_sgs(i, j, k)
  end function eta_tot

  ! E = 0 on the edges lying in a conducting-wall / fixed / axis boundary plane
  subroutine wall_edges
    integer :: c, side, typ
    do c = 1, 3
       if (periodic(c) == 1) cycle
       do side = 1, 2
          if (side == 1 .and. coords(c) /= 0) cycle
          if (side == 2 .and. coords(c) /= dims(c) - 1) cycle
          typ = bc_type(side, c)
          if (typ == 2) cycle
          select case (c)
          case (1)
             if (side == 1) then; emf(2:3, 1, :, :) = 0.0_dp; else; emf(2:3, nx+1, :, :) = 0.0_dp; end if
          case (2)
             if (side == 1) then; emf(1, :, 1, :) = 0.0_dp; emf(3, :, 1, :) = 0.0_dp
             else; emf(1, :, ny+1, :) = 0.0_dp; emf(3, :, ny+1, :) = 0.0_dp; end if
          case (3)
             if (side == 1) then; emf(1:2, :, :, 1) = 0.0_dp; else; emf(1:2, :, :, nz+1) = 0.0_dp; end if
          end select
       end do
    end do
  end subroutine wall_edges

  !--------------------------------------------------------------------
  ! face flux rate: dPhi/dt = - circulation of E (right-hand rule of S_f)
  !   face 1 (i,j,k): + E2(i,j,k) + E3(i,j+1,k) - E2(i,j,k+1) - E3(i,j,k)
  !   face 2 (i,j,k): + E3(i,j,k) + E1(i,j,k+1) - E3(i+1,j,k) - E1(i,j,k)
  !   face 3 (i,j,k): + E1(i,j,k) + E2(i+1,j,k) - E1(i,j+1,k) - E2(i,j,k)
  !--------------------------------------------------------------------
  subroutine update_faces(bnew, bold, dtt)
    real(dp), intent(inout) :: bnew(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)    :: bold(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)    :: dtt
    integer :: i, j, k
    do k = 1, nz+1
       do j = 1, ny+1
          do i = 1, nx+1
             if (j <= ny .and. k <= nz) &
                bnew(1,i,j,k) = bold(1,i,j,k) - dtt*(emf(2,i,j,k) + emf(3,i,j+1,k) - emf(2,i,j,k+1) - emf(3,i,j,k))
             if (i <= nx .and. k <= nz) &
                bnew(2,i,j,k) = bold(2,i,j,k) - dtt*(emf(3,i,j,k) + emf(1,i,j,k+1) - emf(3,i+1,j,k) - emf(1,i,j,k))
             if (i <= nx .and. j <= ny) &
                bnew(3,i,j,k) = bold(3,i,j,k) - dtt*(emf(1,i,j,k) + emf(2,i+1,j,k) - emf(1,i,j+1,k) - emf(2,i,j,k))
          end do
       end do
    end do
    if (debug_balance .and. periodic(3) == 1 .and. dims(3) == 1) then
       write(*,'(A,ES11.3,A,ES11.3)') '  [seam k] max |dPhi3(nz+1) - dPhi3(1)| = ', &
            maxval(abs((bnew(3,1:nx,1:ny,nz+1) - bold(3,1:nx,1:ny,nz+1)) - (bnew(3,1:nx,1:ny,1) - bold(3,1:nx,1:ny,1)))), &
            '   max |dPhi3| = ', maxval(abs(bnew(3,1:nx,1:ny,1) - bold(3,1:nx,1:ny,1)))
       write(*,'(A,ES11.3,A,ES11.3)') '  [seam k] max |E1(nz+1)-E1(1)| = ', maxval(abs(emf(1,1:nx,1:ny+1,nz+1) - emf(1,1:nx,1:ny+1,1))), &
            '  max |E2(nz+1)-E2(1)| = ', maxval(abs(emf(2,1:nx+1,1:ny,nz+1) - emf(2,1:nx+1,1:ny,1)))
    end if
  end subroutine update_faces

  ! normalised divergence: |sum_f Phi_f| / (B_ref V / h) per cell, maximum over the domain,
  ! with B_ref the larger of the local |B| and 1 % of the domain maximum (a cell where the
  ! field vanishes must not read as a divergence error)
  subroutine divergence_B(bb, wp, dmax)
    real(dp), intent(in)  :: bb(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)  :: wp(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(out) :: dmax
    integer :: i, j, k
    real(dp) :: d, bmag, h, bfloor
    bfloor = 0.0_dp
    do k = 1, nz; do j = 1, ny; do i = 1, nx
       bfloor = max(bfloor, norm2(wp(6:8,i,j,k)))
    end do; end do; end do
    bfloor = 0.01_dp*max(gmax(bfloor), 1.0e-30_dp)
    dmax = 0.0_dp
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             d = bb(1,i+1,j,k) - bb(1,i,j,k) + bb(2,i,j+1,k) - bb(2,i,j,k) + bb(3,i,j,k+1) - bb(3,i,j,k)
             bmag = max(norm2(wp(6:8,i,j,k)), bfloor)
             h = vol(i,j,k)**(1.0_dp/3.0_dp)
             divb(i,j,k) = d/vol(i,j,k)*h/bmag
             dmax = max(dmax, abs(divb(i,j,k)))
          end do
       end do
    end do
  end subroutine divergence_B

end module induction
