!=======================================================================
! geometry.F90 -- finite-volume metrics of the curvilinear grid
!
!   Node (i,j,k) is the (minus,minus,minus) corner of cell (i,j,k).
!   Face c of cell (i,j,k): the minus face in direction c, area vector
!     sf(:,c,i,j,k) oriented towards +xi_c (right-handed grid required:
!     the cell volume computed from the face vectors must be positive).
!   Edge c at (i,j,k): the edge along xi_c at the (minus,minus) corner
!     of the other two directions, le = x(node end) - x(node start).
!   Volumes by the divergence theorem, V = (1/3) sum_f x_f . S_f, with
!   the face area vector of a (possibly non-planar) quadrilateral as
!   half the cross product of its diagonals.
!=======================================================================
module geometry
  use global_vars
  use mpi_vars
  implicit none
  private
  public :: build_metrics, cross

contains

  pure function cross(a, b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
  end function cross

  subroutine build_metrics
    integer :: i, j, k, c
    real(dp) :: p1(3), p2(3), p3(3), p4(3), s(3), xcf(3), v, vmin, vl
    real(dp) :: dmax, sfm(3), sfp(3), xm(3), xp(3)

    allocate(vol(1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(xc(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(sf(3, 3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(xf(3, 3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(le(3, 3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(xe(3, 3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(hc(3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))
    allocate(dvec(3, 3, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk))

    ! faces: minus face of cell (i,j,k) in direction c
    do k = 1-ngk, nz+ngk
       do j = 1-ngj, ny+ngj
          do i = 1-ngi, nx+ngi
             ! face 1 (constant i): corners (i,j,k),(i,j+1,k),(i,j+1,k+1),(i,j,k+1): normal ~ xi_2 x xi_3 = +xi_1
             p1 = xn(:, i, j,   k);   p2 = xn(:, i, j+1, k)
             p3 = xn(:, i, j+1, k+1); p4 = xn(:, i, j,   k+1)
             sf(:, 1, i, j, k) = 0.5_dp*cross(p3 - p1, p4 - p2)
             xf(:, 1, i, j, k) = 0.25_dp*(p1 + p2 + p3 + p4)
             ! face 2 (constant j): corners (i,j,k),(i,j,k+1),(i+1,j,k+1),(i+1,j,k): xi_3 x xi_1 = +xi_2
             p1 = xn(:, i,   j, k);   p2 = xn(:, i,   j, k+1)
             p3 = xn(:, i+1, j, k+1); p4 = xn(:, i+1, j, k)
             sf(:, 2, i, j, k) = 0.5_dp*cross(p3 - p1, p4 - p2)
             xf(:, 2, i, j, k) = 0.25_dp*(p1 + p2 + p3 + p4)
             ! face 3 (constant k): corners (i,j,k),(i+1,j,k),(i+1,j+1,k),(i,j+1,k): xi_1 x xi_2 = +xi_3
             p1 = xn(:, i,   j,   k); p2 = xn(:, i+1, j,   k)
             p3 = xn(:, i+1, j+1, k); p4 = xn(:, i,   j+1, k)
             sf(:, 3, i, j, k) = 0.5_dp*cross(p3 - p1, p4 - p2)
             xf(:, 3, i, j, k) = 0.25_dp*(p1 + p2 + p3 + p4)
             ! edges
             le(:, 1, i, j, k) = xn(:, i+1, j, k) - xn(:, i, j, k);  xe(:, 1, i, j, k) = 0.5_dp*(xn(:, i+1, j, k) + xn(:, i, j, k))
             le(:, 2, i, j, k) = xn(:, i, j+1, k) - xn(:, i, j, k);  xe(:, 2, i, j, k) = 0.5_dp*(xn(:, i, j+1, k) + xn(:, i, j, k))
             le(:, 3, i, j, k) = xn(:, i, j, k+1) - xn(:, i, j, k);  xe(:, 3, i, j, k) = 0.5_dp*(xn(:, i, j, k+1) + xn(:, i, j, k))
          end do
       end do
    end do

    ! volumes and centres (cells 1-NG .. n+NG-1: the plus faces of the last ghost layer are not stored)
    vol = 0.0_dp;  xc = 0.0_dp
    vmin = huge(1.0_dp)
    do k = 1-ngk, nz+ngk-1
       do j = 1-ngj, ny+ngj-1
          do i = 1-ngi, nx+ngi-1
             v = 0.0_dp;  xcf = 0.0_dp
             do c = 1, 3
                select case (c)
                case (1); sfm = sf(:,1,i,j,k); xm = xf(:,1,i,j,k); sfp = sf(:,1,i+1,j,k); xp = xf(:,1,i+1,j,k)
                case (2); sfm = sf(:,2,i,j,k); xm = xf(:,2,i,j,k); sfp = sf(:,2,i,j+1,k); xp = xf(:,2,i,j+1,k)
                case (3); sfm = sf(:,3,i,j,k); xm = xf(:,3,i,j,k); sfp = sf(:,3,i,j,k+1); xp = xf(:,3,i,j,k+1)
                end select
                v = v + (dot_product(xp, sfp) - dot_product(xm, sfm))/3.0_dp
                xcf = xcf + 0.5_dp*(xm + xp)/3.0_dp
             end do
             vol(i, j, k) = v
             ! centre: average of the 8 nodes (adequate for smooth grids)
             xc(:, i, j, k) = 0.125_dp*(xn(:,i,j,k) + xn(:,i+1,j,k) + xn(:,i,j+1,k) + xn(:,i+1,j+1,k) + &
                                        xn(:,i,j,k+1) + xn(:,i+1,j,k+1) + xn(:,i,j+1,k+1) + xn(:,i+1,j+1,k+1))
             if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny .and. k >= 1 .and. k <= nz) vmin = min(vmin, v)
          end do
       end do
    end do
    ! last ghost layer: copy neighbours (never used for fluxes)
    vol(nx+ngi, :, :) = vol(nx+ngi-1, :, :);  vol(:, ny+ngj, :) = vol(:, ny+ngj-1, :);  vol(:, :, nz+ngk) = vol(:, :, nz+ngk-1)
    xc(:, nx+ngi, :, :) = xc(:, nx+ngi-1, :, :);  xc(:, :, ny+ngj, :) = xc(:, :, ny+ngj-1, :);  xc(:, :, :, nz+ngk) = xc(:, :, :, nz+ngk-1)

    vl = gmin(vmin)
    if (vl <= 0.0_dp) then
       if (myrank == 0) write(*,'(A,ES12.4)') ' ERROR: non-positive cell volume (left-handed or folded grid), min V = ', vl
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    vol_total = gsum(sum(vol(1:nx, 1:ny, 1:nz)))

    ! characteristic lengths and centre-to-centre vectors
    do k = 2-ngk, nz+ngk-1
       do j = 2-ngj, ny+ngj-1
          do i = 2-ngi, nx+ngi-1
             do c = 1, 3
                select case (c)
                case (1); dmax = 0.5_dp*(norm2(sf(:,1,i,j,k)) + norm2(sf(:,1,i+1,j,k))); dvec(:,1,i,j,k) = xc(:,i,j,k) - xc(:,i-1,j,k)
                case (2); dmax = 0.5_dp*(norm2(sf(:,2,i,j,k)) + norm2(sf(:,2,i,j+1,k))); dvec(:,2,i,j,k) = xc(:,i,j,k) - xc(:,i,j-1,k)
                case (3); dmax = 0.5_dp*(norm2(sf(:,3,i,j,k)) + norm2(sf(:,3,i,j,k+1))); dvec(:,3,i,j,k) = xc(:,i,j,k) - xc(:,i,j,k-1)
                end select
                hc(c, i, j, k) = vol(i, j, k)/max(dmax, tiny(1.0_dp))
             end do
          end do
       end do
    end do
    hc(:, 1-ngi, :, :) = hc(:, 2-ngi, :, :);  hc(:, :, 1-ngj, :) = hc(:, :, 2-ngj, :);  hc(:, :, :, 1-ngk) = hc(:, :, :, 2-ngk)
    hc(:, nx+ngi, :, :) = hc(:, nx+ngi-1, :, :);  hc(:, :, ny+ngj, :) = hc(:, :, ny+ngj-1, :);  hc(:, :, :, nz+ngk) = hc(:, :, :, nz+ngk-1)
    dvec(:, :, 1-ngi, :, :) = dvec(:, :, 2-ngi, :, :);  dvec(:, :, :, 1-ngj, :) = dvec(:, :, :, 2-ngj, :);  dvec(:, :, :, :, 1-ngk) = dvec(:, :, :, :, 2-ngk)
    dvec(:, :, nx+ngi, :, :) = dvec(:, :, nx+ngi-1, :, :);  dvec(:, :, :, ny+ngj, :) = dvec(:, :, :, ny+ngj-1, :);  dvec(:, :, :, :, nz+ngk) = dvec(:, :, :, :, nz+ngk-1)

    if (myrank == 0) then
       write(*,'(A,ES12.4,A,ES12.4)') '  grid metrics     : min cell volume* ', vl, '  total volume* ', vol_total
    end if
  end subroutine build_metrics

end module geometry
