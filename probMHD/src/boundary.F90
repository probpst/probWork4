!=======================================================================
! boundary.F90 -- ghost cells of the physical (non-periodic) boundaries
!
!   BC types (bc_type(side, dir)):
!     0 periodic       handled by the MPI Cartesian exchange
!     1 conducting wall (slip): rho, p, B_t copied; v_n and B_n mirrored;
!                      the wall-edge EMF is zero (induction.F90), so the
!                      flux through the wall faces stays at its initial
!                      value (zero for a flux-surface boundary)
!     2 zero gradient  all variables and face fluxes copied outward
!     3 fixed          ghost cells hold the initial state (ueq, bfeq);
!                      the boundary-edge EMF is zero
!     4 axis hole      conducting core: as 1 (the inner boundary of a
!                      flux-aligned grid is a hole, not the axis)
!   Periodic directions are exchanged by mpi_vars (also across ranks).
!=======================================================================
module boundary
  use global_vars
  use mpi_vars
  implicit none
  private
  public :: apply_bc, apply_bc_faces, apply_bc_scalar

contains

  ! ghost cells of the conserved array (and the primitive w for the ghosts)
  subroutine apply_bc(uu, bb)
    real(dp), intent(inout) :: uu(NVAR, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)    :: bb(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: c, side, g, i, j, k, ib, ig, typ, ng3(3)
    ng3 = [ngi, ngj, ngk]
    do c = 1, 3
       if (periodic(c) == 1) cycle
       do side = 1, 2
          if (side == 1 .and. coords(c) /= 0) cycle
          if (side == 2 .and. coords(c) /= dims(c) - 1) cycle
          typ = bc_type(side, c)
          do g = 1, ng3(c)
             select case (c)
             case (1)
                do k = 1-ngk, nz+ngk
                   do j = 1-ngj, ny+ngj
                      call ghost_pair(side, g, nx, ib, ig)
                      call set_ghost(typ, uu(:, ib, j, k), uu(:, ig, j, k), ueq(:, ig, j, k), face_normal(side, 1, ib, j, k))
                   end do
                end do
             case (2)
                do k = 1-ngk, nz+ngk
                   do i = 1-ngi, nx+ngi
                      call ghost_pair(side, g, ny, ib, ig)
                      call set_ghost(typ, uu(:, i, ib, k), uu(:, i, ig, k), ueq(:, i, ig, k), face_normal(side, 2, i, ib, k))
                   end do
                end do
             case (3)
                do j = 1-ngj, ny+ngj
                   do i = 1-ngi, nx+ngi
                      call ghost_pair(side, g, nz, ib, ig)
                      call set_ghost(typ, uu(:, i, j, ib), uu(:, i, j, ig), ueq(:, i, j, ig), face_normal(side, 3, i, j, ib))
                   end do
                end do
             end select
          end do
       end do
    end do
  end subroutine apply_bc

  ! ghost layers of a cell scalar (1, cells): copied from the mirror cell, or held (fixed BC)
  subroutine apply_bc_scalar(a, ae)
    real(dp), intent(inout) :: a(1, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)    :: ae(1, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: c, side, g, ib, ig, typ, ng3(3)
    ng3 = [ngi, ngj, ngk]
    do c = 1, 3
       if (periodic(c) == 1) cycle
       do side = 1, 2
          if (side == 1 .and. coords(c) /= 0) cycle
          if (side == 2 .and. coords(c) /= dims(c) - 1) cycle
          typ = bc_type(side, c)
          do g = 1, ng3(c)
             select case (c)
             case (1); call ghost_pair(side, g, nx, ib, ig)
                if (typ == 3) then; a(1, ig, :, :) = ae(1, ig, :, :); else; a(1, ig, :, :) = a(1, ib, :, :); end if
             case (2); call ghost_pair(side, g, ny, ib, ig)
                if (typ == 3) then; a(1, :, ig, :) = ae(1, :, ig, :); else; a(1, :, ig, :) = a(1, :, ib, :); end if
             case (3); call ghost_pair(side, g, nz, ib, ig)
                if (typ == 3) then; a(1, :, :, ig) = ae(1, :, :, ig); else; a(1, :, :, ig) = a(1, :, :, ib); end if
             end select
          end do
       end do
    end do
  end subroutine apply_bc_scalar

  ! interior mirror cell ib and ghost ig for layer g (mirror about the boundary face)
  pure subroutine ghost_pair(side, g, n, ib, ig)
    integer, intent(in)  :: side, g, n
    integer, intent(out) :: ib, ig
    if (side == 1) then
       ig = 1 - g;  ib = g
    else
       ig = n + g;  ib = n + 1 - g
    end if
    ib = max(min(ib, n), 1)
  end subroutine ghost_pair

  function face_normal(side, c, i, j, k) result(nrm)
    integer, intent(in) :: side, c, i, j, k
    real(dp) :: nrm(3), s(3)
    integer :: ii, jj, kk
    ii = i;  jj = j;  kk = k
    if (side == 2) then        ! the plus face of the last cell = minus face of cell n+1
       if (c == 1) ii = nx + 1
       if (c == 2) jj = ny + 1
       if (c == 3) kk = nz + 1
    else
       if (c == 1) ii = 1
       if (c == 2) jj = 1
       if (c == 3) kk = 1
    end if
    s = sf(:, c, ii, jj, kk)
    nrm = s/max(norm2(s), tiny(1.0_dp))
  end function face_normal

  subroutine set_ghost(typ, ub, ug, ue, nrm)
    integer,  intent(in)    :: typ
    real(dp), intent(in)    :: ub(NVAR), ue(NVAR), nrm(3)
    real(dp), intent(inout) :: ug(NVAR)
    real(dp) :: wc(NVAR), vn, bn
    select case (typ)
    case (3)
       ug = ue
    case (1, 4)
       call cons2prim(ub, wc)
       vn = dot_product(wc(2:4), nrm);  wc(2:4) = wc(2:4) - 2.0_dp*vn*nrm
       bn = dot_product(wc(6:8), nrm);  wc(6:8) = wc(6:8) - 2.0_dp*bn*nrm
       call prim2cons(wc, ug)
    case default
       ug = ub
    end select
  end subroutine set_ghost

  ! ghost face fluxes
  subroutine apply_bc_faces(bb)
    real(dp), intent(inout) :: bb(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: c, side, g, typ, ng3(3)
    ng3 = [ngi, ngj, ngk]
    do c = 1, 3
       if (periodic(c) == 1) cycle
       do side = 1, 2
          if (side == 1 .and. coords(c) /= 0) cycle
          if (side == 2 .and. coords(c) /= dims(c) - 1) cycle
          typ = bc_type(side, c)
          do g = 1, ng3(c)
             select case (typ)
             case (3)
                call copy_face_layer(bb, bfeq, c, side, g)
             case (1, 4)
                call mirror_face_layer(bb, c, side, g)
             case default
                call extrap_face_layer(bb, c, side, g)
             end select
          end do
       end do
    end do
  end subroutine apply_bc_faces

  ! fixed: ghost faces = initial
  subroutine copy_face_layer(bb, be, c, side, g)
    real(dp), intent(inout) :: bb(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)    :: be(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer, intent(in) :: c, side, g
    integer :: ig
    select case (c)
    case (1); ig = merge(1 - g, nx + g, side == 1); bb(:, ig, :, :) = be(:, ig, :, :)
              if (side == 2 .and. g == 1) bb(1, nx+1, :, :) = be(1, nx+1, :, :)
    case (2); ig = merge(1 - g, ny + g, side == 1); bb(:, :, ig, :) = be(:, :, ig, :)
              if (side == 2 .and. g == 1) bb(2, :, ny+1, :) = be(2, :, ny+1, :)
    case (3); ig = merge(1 - g, nz + g, side == 1); bb(:, :, :, ig) = be(:, :, :, ig)
              if (side == 2 .and. g == 1) bb(3, :, :, nz+1) = be(3, :, :, nz+1)
    end select
  end subroutine copy_face_layer

  ! zero gradient: ghost faces copy the nearest interior faces
  subroutine extrap_face_layer(bb, c, side, g)
    real(dp), intent(inout) :: bb(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer, intent(in) :: c, side, g
    integer :: ig, ib
    select case (c)
    case (1)
       if (side == 1) then
          ig = 1 - g;  ib = 1
          bb(:, ig, :, :) = bb(:, ib, :, :)
       else
          ig = nx + g;  ib = nx
          bb(2:3, ig, :, :) = bb(2:3, ib, :, :)
          if (ig + 1 <= nx + ngi) bb(1, ig+1, :, :) = bb(1, nx+1, :, :)
       end if
    case (2)
       if (side == 1) then
          ig = 1 - g;  ib = 1
          bb(:, :, ig, :) = bb(:, :, ib, :)
       else
          ig = ny + g;  ib = ny
          bb(1, :, ig, :) = bb(1, :, ib, :);  bb(3, :, ig, :) = bb(3, :, ib, :)
          if (ig + 1 <= ny + ngj) bb(2, :, ig+1, :) = bb(2, :, ny+1, :)
       end if
    case (3)
       if (side == 1) then
          ig = 1 - g;  ib = 1
          bb(:, :, :, ig) = bb(:, :, :, ib)
       else
          ig = nz + g;  ib = nz
          bb(1:2, :, :, ig) = bb(1:2, :, :, ib)
          if (ig + 1 <= nz + ngk) bb(3, :, :, ig+1) = bb(3, :, :, nz+1)
       end if
    end select
  end subroutine extrap_face_layer

  ! conducting wall: normal flux mirrored (the wall face keeps its value), tangential faces copied
  subroutine mirror_face_layer(bb, c, side, g)
    real(dp), intent(inout) :: bb(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer, intent(in) :: c, side, g
    integer :: ig, ib
    select case (c)
    case (1)
       if (side == 1) then
          ig = 1 - g;  ib = g
          bb(2:3, ig, :, :) = bb(2:3, ib, :, :)
          ! normal faces: face ig (minus face of ghost ig) mirrors face ib+1 (plus face of interior ib) about face 1
          bb(1, ig, :, :) = 2.0_dp*bb(1, 1, :, :) - bb(1, ib+1, :, :)
       else
          ig = nx + g;  ib = nx + 1 - g
          bb(2:3, ig, :, :) = bb(2:3, ib, :, :)
          if (ig + 1 <= nx + ngi) bb(1, ig+1, :, :) = 2.0_dp*bb(1, nx+1, :, :) - bb(1, ib, :, :)
       end if
    case (2)
       if (side == 1) then
          ig = 1 - g;  ib = g
          bb(1, :, ig, :) = bb(1, :, ib, :);  bb(3, :, ig, :) = bb(3, :, ib, :)
          bb(2, :, ig, :) = 2.0_dp*bb(2, :, 1, :) - bb(2, :, ib+1, :)
       else
          ig = ny + g;  ib = ny + 1 - g
          bb(1, :, ig, :) = bb(1, :, ib, :);  bb(3, :, ig, :) = bb(3, :, ib, :)
          if (ig + 1 <= ny + ngj) bb(2, :, ig+1, :) = 2.0_dp*bb(2, :, ny+1, :) - bb(2, :, ib, :)
       end if
    case (3)
       if (side == 1) then
          ig = 1 - g;  ib = g
          bb(1:2, :, :, ig) = bb(1:2, :, :, ib)
          bb(3, :, :, ig) = 2.0_dp*bb(3, :, :, 1) - bb(3, :, :, ib+1)
       else
          ig = nz + g;  ib = nz + 1 - g
          bb(1:2, :, :, ig) = bb(1:2, :, :, ib)
          if (ig + 1 <= nz + ngk) bb(3, :, :, ig+1) = 2.0_dp*bb(3, :, :, nz+1) - bb(3, :, :, ib)
       end if
    end select
  end subroutine mirror_face_layer

end module boundary
