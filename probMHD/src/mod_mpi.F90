!=======================================================================
! mod_mpi.F90 -- MPI set-up, 3-D Cartesian decomposition of the (single)
!   grid block, ghost-layer exchange and global reductions / gathers.
!
!   Decomposition: MPI_Cart_create with the periods of the grid, so the
!   ghost layers across periodic boundaries are filled by the same
!   exchange as the interior sub-domain interfaces.  Physical (non-
!   periodic) boundaries have MPI_PROC_NULL neighbours; their ghost
!   cells are set by boundary.F90.
!=======================================================================
module mpi_vars
  use mpi
  use global_vars
  implicit none
  integer :: myrank = 0, nprocs = 1, ierr
  integer :: comm_cart
  integer :: dims(3), coords(3)
  logical :: periods(3)
  integer :: nbr_lo(3), nbr_hi(3)
  real(dp) :: wall_t0

contains

  subroutine mpi_setup
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
    wall_t0 = MPI_Wtime()
  end subroutine mpi_setup

  ! build the Cartesian topology once the global cell counts are known
  subroutine mpi_decompose
    integer :: c
    logical :: reorder
    dims = nproc_dir
    ! directions with a single cell plane cannot be split
    do c = 1, 3
       if (dims(c) == 0) then
          if ((c == 1 .and. ngx == 1) .or. (c == 2 .and. ngy == 1) .or. (c == 3 .and. ngz == 1)) dims(c) = 1
       end if
    end do
    call MPI_Dims_create(nprocs, 3, dims, ierr)
    if (product(dims) /= nprocs) then
       if (myrank == 0) write(*,'(A,3I4,A,I4)') ' ERROR: MPI decomposition ', dims, ' does not match the process count ', nprocs
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
    periods = (periodic == 1)
    reorder = .false.
    call MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, comm_cart, ierr)
    call MPI_Cart_coords(comm_cart, myrank, 3, coords, ierr)
    do c = 1, 3
       call MPI_Cart_shift(comm_cart, c - 1, 1, nbr_lo(c), nbr_hi(c), ierr)
    end do
    call split_range(ngx, dims(1), coords(1), ioff, nx)
    call split_range(ngy, dims(2), coords(2), joff, ny)
    call split_range(ngz, dims(3), coords(3), koff, nz)
    ngi = merge(1, NG, ngx == 1);  ngj = merge(1, NG, ngy == 1);  ngk = merge(1, NG, ngz == 1)
    if ((nx < ngi .and. dims(1) > 1) .or. (ny < ngj .and. dims(2) > 1) .or. (nz < ngk .and. dims(3) > 1)) then
       if (myrank == 0) write(*,'(A)') ' ERROR: a sub-domain is thinner than the ghost layer; use fewer processes in that direction'
       call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
  end subroutine mpi_decompose

  subroutine split_range(n, np, ip, off, nl)
    integer, intent(in)  :: n, np, ip
    integer, intent(out) :: off, nl
    integer :: base, rem
    base = n/np;  rem = mod(n, np)
    nl  = base + merge(1, 0, ip < rem)
    off = ip*base + min(ip, rem)
  end subroutine split_range

  !--------------------------------------------------------------------
  ! ghost exchange of a (nv, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk) array
  !--------------------------------------------------------------------
  subroutine exchange(a, nv)
    integer,  intent(in)    :: nv
    real(dp), intent(inout) :: a(nv, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: c, n(3), lo(3), hi(3), cnt, ng(3)
    real(dp), allocatable, save :: sbuf(:,:,:,:), rbuf(:,:,:,:)
    integer :: st(MPI_STATUS_SIZE)
    real(dp) :: ang
    n = [nx, ny, nz];  ng = [ngi, ngj, ngk]
    do c = 1, 3
       lo = 1 - ng;  hi = n + ng       ! full extent in the other directions (corners filled progressively)
       ang = periodic_rot(c)*PI/180.0_dp
       if (n(c) < ng(c) .and. periods(c) .and. dims(c) == 1) then
          ! fewer cells than ghost layers in a periodic direction held by one rank: wrap by copying
          call fill_self(a, nv, c)
          if (ang /= 0.0_dp .and. nv == NVAR) then
             call rotate_ghosts(a, nv, c, 1, -ang)
             call rotate_ghosts(a, nv, c, 2,  ang)
          end if
          cycle
       end if
       cnt = nv*ng(c)*product(hi - lo + 1)/(hi(c) - lo(c) + 1)
       select case (c)
       case (1)
          allocate(sbuf(nv, ngi, lo(2):hi(2), lo(3):hi(3)), rbuf(nv, ngi, lo(2):hi(2), lo(3):hi(3)))
          ! send to hi neighbour the last NG cells, receive into the low ghosts
          sbuf = a(:, nx-ngi+1:nx, :, :)
          call MPI_Sendrecv(sbuf, cnt, MPI_DOUBLE_PRECISION, nbr_hi(1), 1, rbuf, cnt, MPI_DOUBLE_PRECISION, nbr_lo(1), 1, comm_cart, st, ierr)
          if (nbr_lo(1) /= MPI_PROC_NULL) a(:, 1-ngi:0, :, :) = rbuf
          sbuf = a(:, 1:ngi, :, :)
          call MPI_Sendrecv(sbuf, cnt, MPI_DOUBLE_PRECISION, nbr_lo(1), 2, rbuf, cnt, MPI_DOUBLE_PRECISION, nbr_hi(1), 2, comm_cart, st, ierr)
          if (nbr_hi(1) /= MPI_PROC_NULL) a(:, nx+1:nx+ngi, :, :) = rbuf
       case (2)
          allocate(sbuf(nv, lo(1):hi(1), ngj, lo(3):hi(3)), rbuf(nv, lo(1):hi(1), ngj, lo(3):hi(3)))
          sbuf = a(:, :, ny-ngj+1:ny, :)
          call MPI_Sendrecv(sbuf, cnt, MPI_DOUBLE_PRECISION, nbr_hi(2), 3, rbuf, cnt, MPI_DOUBLE_PRECISION, nbr_lo(2), 3, comm_cart, st, ierr)
          if (nbr_lo(2) /= MPI_PROC_NULL) a(:, :, 1-ngj:0, :) = rbuf
          sbuf = a(:, :, 1:ngj, :)
          call MPI_Sendrecv(sbuf, cnt, MPI_DOUBLE_PRECISION, nbr_lo(2), 4, rbuf, cnt, MPI_DOUBLE_PRECISION, nbr_hi(2), 4, comm_cart, st, ierr)
          if (nbr_hi(2) /= MPI_PROC_NULL) a(:, :, ny+1:ny+ngj, :) = rbuf
       case (3)
          allocate(sbuf(nv, lo(1):hi(1), lo(2):hi(2), ngk), rbuf(nv, lo(1):hi(1), lo(2):hi(2), ngk))
          sbuf = a(:, :, :, nz-ngk+1:nz)
          call MPI_Sendrecv(sbuf, cnt, MPI_DOUBLE_PRECISION, nbr_hi(3), 5, rbuf, cnt, MPI_DOUBLE_PRECISION, nbr_lo(3), 5, comm_cart, st, ierr)
          if (nbr_lo(3) /= MPI_PROC_NULL) a(:, :, :, 1-ngk:0) = rbuf
          sbuf = a(:, :, :, 1:ngk)
          call MPI_Sendrecv(sbuf, cnt, MPI_DOUBLE_PRECISION, nbr_lo(3), 6, rbuf, cnt, MPI_DOUBLE_PRECISION, nbr_hi(3), 6, comm_cart, st, ierr)
          if (nbr_hi(3) /= MPI_PROC_NULL) a(:, :, :, nz+1:nz+ngk) = rbuf
       end select
       deallocate(sbuf, rbuf)      ! (shape depends on nv and the direction)
       ! rotational periodicity: vectors of the ghost images are rotated
       if (ang /= 0.0_dp .and. nv == NVAR .and. periods(c)) then
          if (coords(c) == 0)         call rotate_ghosts(a, nv, c, 1, -ang)
          if (coords(c) == dims(c)-1) call rotate_ghosts(a, nv, c, 2,  ang)
       end if
    end do
  end subroutine exchange

  ! rotate the vector components (2:4 velocity, 6:8 field) of the ghost layer on side 1 (low) / 2 (high) about z
  subroutine rotate_ghosts(a, nv, c, side, ang)
    integer,  intent(in)    :: nv, c, side
    real(dp), intent(inout) :: a(nv, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(in)    :: ang
    integer :: i, j, k, g, n(3), idx, ng(3)
    real(dp) :: ca, sa
    ca = cos(ang);  sa = sin(ang)
    n = [nx, ny, nz];  ng = [ngi, ngj, ngk]
    do g = 1, ng(c)
       idx = merge(1 - g, n(c) + g, side == 1)
       select case (c)
       case (1)
          do k = 1-ngk, nz+ngk; do j = 1-ngj, ny+ngj
             call rot2(a(2:4, idx, j, k));  call rot2(a(6:8, idx, j, k))
          end do; end do
       case (2)
          do k = 1-ngk, nz+ngk; do i = 1-ngi, nx+ngi
             call rot2(a(2:4, i, idx, k));  call rot2(a(6:8, i, idx, k))
          end do; end do
       case (3)
          do j = 1-ngj, ny+ngj; do i = 1-ngi, nx+ngi
             call rot2(a(2:4, i, j, idx));  call rot2(a(6:8, i, j, idx))
          end do; end do
       end select
    end do
  contains
    subroutine rot2(v)
      real(dp), intent(inout) :: v(3)
      real(dp) :: x, y
      x = v(1);  y = v(2)
      v(1) = ca*x - sa*y;  v(2) = sa*x + ca*y
    end subroutine rot2
  end subroutine rotate_ghosts

  ! periodic direction with fewer local cells than NG on one rank: wrap by copying (n(c) >= 1)
  subroutine fill_self(a, nv, c)
    integer,  intent(in)    :: nv, c
    real(dp), intent(inout) :: a(nv, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: g, n(3), src, ng(3)
    n = [nx, ny, nz];  ng = [ngi, ngj, ngk]
    do g = 1, ng(c)
       src = modulo(n(c) - g, n(c)) + 1          ! image of ghost 1-g
       select case (c)
       case (1); a(:, 1-g, :, :) = a(:, src, :, :)
       case (2); a(:, :, 1-g, :) = a(:, :, src, :)
       case (3); a(:, :, :, 1-g) = a(:, :, :, src)
       end select
       src = modulo(g - 1, n(c)) + 1             ! image of ghost n+g
       select case (c)
       case (1); a(:, nx+g, :, :) = a(:, src, :, :)
       case (2); a(:, :, ny+g, :) = a(:, :, src, :)
       case (3); a(:, :, :, nz+g) = a(:, :, :, src)
       end select
    end do
  end subroutine fill_self

  !--------------------------------------------------------------------
  real(dp) function gmax(x)
    real(dp), intent(in) :: x
    call MPI_Allreduce(x, gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  end function gmax

  real(dp) function gmin(x)
    real(dp), intent(in) :: x
    call MPI_Allreduce(x, gmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
  end function gmin

  real(dp) function gsum(x)
    real(dp), intent(in) :: x
    call MPI_Allreduce(x, gsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  end function gsum

  logical function gany(l)
    logical, intent(in) :: l
    call MPI_Allreduce(l, gany, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
  end function gany

  !--------------------------------------------------------------------
  ! gather a (nv, interior cells) array to rank 0 as (nv, ngx, ngy, ngz)
  !--------------------------------------------------------------------
  subroutine gather_cells(a, nv, g)
    integer,  intent(in)  :: nv
    real(dp), intent(in)  :: a(nv, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    real(dp), intent(out) :: g(:,:,:,:)      ! allocated on rank 0 as (nv, ngx, ngy, ngz)
    integer :: r, info(6), rinfo(6), st(MPI_STATUS_SIZE), cnt
    real(dp), allocatable :: buf(:,:,:,:)
    info = [ioff, joff, koff, nx, ny, nz]
    if (myrank == 0) then
       g(:, ioff+1:ioff+nx, joff+1:joff+ny, koff+1:koff+nz) = a(:, 1:nx, 1:ny, 1:nz)
       do r = 1, nprocs - 1
          call MPI_Recv(rinfo, 6, MPI_INTEGER, r, 10, MPI_COMM_WORLD, st, ierr)
          allocate(buf(nv, rinfo(4), rinfo(5), rinfo(6)))
          cnt = nv*rinfo(4)*rinfo(5)*rinfo(6)
          call MPI_Recv(buf, cnt, MPI_DOUBLE_PRECISION, r, 11, MPI_COMM_WORLD, st, ierr)
          g(:, rinfo(1)+1:rinfo(1)+rinfo(4), rinfo(2)+1:rinfo(2)+rinfo(5), rinfo(3)+1:rinfo(3)+rinfo(6)) = buf
          deallocate(buf)
       end do
    else
       call MPI_Send(info, 6, MPI_INTEGER, 0, 10, MPI_COMM_WORLD, ierr)
       allocate(buf(nv, nx, ny, nz))
       buf = a(:, 1:nx, 1:ny, 1:nz)
       cnt = nv*nx*ny*nz
       call MPI_Send(buf, cnt, MPI_DOUBLE_PRECISION, 0, 11, MPI_COMM_WORLD, ierr)
       deallocate(buf)
    end if
  end subroutine gather_cells

  ! scatter the inverse: rank 0 holds g(nv, ngx, ngy, ngz); every rank receives its interior
  subroutine scatter_cells(g, nv, a)
    integer,  intent(in)    :: nv
    real(dp), intent(in)    :: g(:,:,:,:)
    real(dp), intent(inout) :: a(nv, 1-ngi:nx+ngi, 1-ngj:ny+ngj, 1-ngk:nz+ngk)
    integer :: r, info(6), rinfo(6), st(MPI_STATUS_SIZE), cnt
    real(dp), allocatable :: buf(:,:,:,:)
    info = [ioff, joff, koff, nx, ny, nz]
    if (myrank == 0) then
       a(:, 1:nx, 1:ny, 1:nz) = g(:, ioff+1:ioff+nx, joff+1:joff+ny, koff+1:koff+nz)
       do r = 1, nprocs - 1
          call MPI_Recv(rinfo, 6, MPI_INTEGER, r, 12, MPI_COMM_WORLD, st, ierr)
          allocate(buf(nv, rinfo(4), rinfo(5), rinfo(6)))
          buf = g(:, rinfo(1)+1:rinfo(1)+rinfo(4), rinfo(2)+1:rinfo(2)+rinfo(5), rinfo(3)+1:rinfo(3)+rinfo(6))
          cnt = nv*rinfo(4)*rinfo(5)*rinfo(6)
          call MPI_Send(buf, cnt, MPI_DOUBLE_PRECISION, r, 13, MPI_COMM_WORLD, ierr)
          deallocate(buf)
       end do
    else
       call MPI_Send(info, 6, MPI_INTEGER, 0, 12, MPI_COMM_WORLD, ierr)
       allocate(buf(nv, nx, ny, nz))
       cnt = nv*nx*ny*nz
       call MPI_Recv(buf, cnt, MPI_DOUBLE_PRECISION, 0, 13, MPI_COMM_WORLD, st, ierr)
       a(:, 1:nx, 1:ny, 1:nz) = buf
       deallocate(buf)
    end if
  end subroutine scatter_cells

  subroutine mpi_abort_msg(msg)
    character(len=*), intent(in) :: msg
    write(*,'(A)') ' ERROR: '//trim(msg)
    flush(6)
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end subroutine mpi_abort_msg

end module mpi_vars
