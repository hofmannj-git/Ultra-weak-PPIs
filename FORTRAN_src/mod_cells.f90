module mod_cells

  use mod_particle, only : particle_t
  use mod_domain, only : domain_t
  use mod_group, only : group_t
  implicit none
  private

  type, public :: cells_t

     ! handles to domain, particle, and group

     type(domain_t), pointer :: domain => null()
     type(particle_t), pointer :: particle => null()
     type(group_t), pointer :: group => null()

     ! state flags

     logical :: cells_exist
     logical :: current ! necessary?

     ! inclusion and exclusion

     integer :: includegroup

     ! cell parameters
     ! need minimum size specification for real-space stencil
     ! need maximum size specification for reciprocal-space grid
     ! stay_under says to keep grid spacing below a certain value
     ! (sphere with diameter 'value' can enclose triclinic cell)
     ! and is .true. when using reciprocal space (.false. for real space)

     real :: target_size
     logical :: stay_under

     ! (dynamic) number of cells in each dimension

     integer :: ncells
     integer, dimension(3) :: ncell
     integer :: ncellx,ncelly,ncellz

     ! pointers to the first particles in each cell (cellhead)
     ! and the pointer to the next particle in the cell (next)
     ! and the number of (distinct) cells in stencil around current

     integer, dimension(:), allocatable :: cellhead
     integer, dimension(:), allocatable :: next
     integer :: nstencil

     ! sizes and inverse sizes of the regularly distributed bins

     real :: binsizex,binsizey,binsizez
     real :: bininvx,bininvy,bininvz

     ! rectangular or lamda extreme coordinates
     ! have to set what these point to

     real, dimension(:), pointer :: boxlo,boxhi

   contains

     procedure, public :: init
     procedure, public :: reset_grid
     procedure :: coord2cell
     procedure, public :: assign_particles
     procedure, public :: stencil_neighbor
     procedure, public :: cell2triplet
     procedure, public :: triplet2cell

  end type cells_t

contains

  ! ----------------------------------------------------------------------
  ! allocate and set parameters
  ! ----------------------------------------------------------------------

  subroutine init(this, group, target_size, stay_under, includegroup)

    class(cells_t), intent(inout)     :: this
    type(group_t), intent(in), target :: group
    real, intent(in)                  :: target_size
    logical, intent(in)               :: stay_under
    integer, intent(in), optional     :: includegroup

    ! check inputs

    if (target_size <= 0.) then
       print *, 'Minimum size in cells must be > 0'
       stop
    end if

    ! set handles and defaults

    this%group => group
    this%domain => group%domain
    this%particle => group%particle
    this%target_size = target_size
    this%stay_under = stay_under
    this%ncells = 0
    this%nstencil = 13

    ! optional arguments

    if (present(includegroup)) then
       if (includegroup < 1) then
          print *, 'Invalid group ID for cell building'
          stop
       end if
    else
       this%includegroup = 0
    end if

    ! calculate the cells grid, allocating the cellhead

    this%current = .false.
    call this%reset_grid()

    ! allocate the 'next' array and default to zero (no 'next' particles)

    allocate(this%next(this%particle%nparticles))
    this%next(:) = 0

    ! choose whether triclinic or orthogonal here

    if (this%domain%triclinic) then
       this%boxlo => this%domain%boxlo_lamda
       this%boxhi => this%domain%boxhi_lamda
    else
       this%boxlo => this%domain%boxlo
       this%boxhi => this%domain%boxhi
    end if

    ! cells now exist

    this%cells_exist = .true.

  end subroutine init

  ! ----------------------------------------------------------------------
  ! reset the grid, assuming same cutoff sizes
  ! this may result in a different number of cells per dimension
  ! perform allocation if the current cell dimensions are too small
  ! if nx, ny, nz specified, just use those values
  ! ----------------------------------------------------------------------

  subroutine reset_grid(this, nx, ny, nz)

    class(cells_t), intent(inout) :: this
    integer, intent(in), optional :: nx, ny, nz
    integer                       :: cell_product
    real, dimension(3)            :: binsize,bininv

    if (.not. (present(nx) .and. present(ny) .and. present(nz))) then
       call determine_ncell(this%domain, this%target_size, &
            this%stay_under, this%ncell)
    else
       if ((nx > 2) .and. (ny > 2) .and. (nz > 2)) then
          this%ncell = (/ nx, ny, nz /)
       else
          write(*, '(a)') 'error, need to specify nx,ny,nz > 2 in reset_grid'
          stop
       end if
    end if

    ! reallocate if what we have is too small

    cell_product = product(this%ncell(:))
    if (this%ncells < cell_product) then
       if (allocated(this%cellhead)) then
          deallocate(this%cellhead)
       end if
       allocate(this%cellhead(cell_product))
    end if

    ! set the variables
    ! if triclinic, the bin size is in lamda coords

    if (this%domain%triclinic) then
       binsize(:) = this%domain%prd_lamda(:) / this%ncell(:)
    else
       binsize(:) = this%domain%prd(:) / this%ncell(:)
    end if

    bininv(:) = 1.0 / binsize(:)

    this%ncellx = this%ncell(1)
    this%ncelly = this%ncell(2)
    this%ncellz = this%ncell(3)
    this%binsizex = binsize(1)
    this%binsizey = binsize(2)
    this%binsizez = binsize(3)
    this%bininvx = bininv(1)
    this%bininvy = bininv(2)
    this%bininvz = bininv(3)
    this%ncells = product(this%ncell(:))

    this%current = .true.

  end subroutine reset_grid

  ! ----------------------------------------------------------------------
  ! convert particle coords into cell #
  ! coords assumed to match (regular, triclinic) bininv variables
  ! ----------------------------------------------------------------------

  integer function coord2cell(this, x)

    class(cells_t), intent(in)     :: this
    real, dimension(:), intent(in) :: x
    integer                        :: ix,iy,iz
    real                           :: bininvx,bininvy,bininvz
    integer                        :: ncellx,ncelly,ncellz
    real, dimension(3)             :: boxlo,boxhi

    boxlo(:) = this%boxlo(:)
    boxhi(:) = this%boxhi(:)
    bininvx = this%bininvx
    bininvy = this%bininvy
    bininvz = this%bininvz
    ncellx = this%ncellx
    ncelly = this%ncelly
    ncellz = this%ncellz

    ! note: indices start from 1, unlike in LAMMPS (0)
    ! TODO: is min necessary here?

    ix = int((x(1)-boxlo(1))*bininvx) + 1
    ix = min(ix,ncellx)
    iy = int((x(2)-boxlo(2))*bininvy) + 1
    iy = min(iy, ncelly)
    iz = int((x(3)-boxlo(3))*bininvz) + 1
    iz = min(iz,ncellz)

    coord2cell = (iz-1) * ncelly * ncellx &
               + (iy-1) * ncellx &
               + ix

  end function coord2cell

  ! ----------------------------------------------------------------------
  ! assign particles to cells
  ! works with orthogonal and triclinic boxes
  ! ----------------------------------------------------------------------

  subroutine assign_particles(this)

    class(cells_t), intent(inout)  :: this
    real, dimension(:,:), pointer  :: x => null()
    integer, dimension(:), pointer :: mask => null()
    integer                        :: i,icell,nparticles,bitmask

    ! verify that particles are wrapped inside periodic boundaries

    if (.not. this%domain%particles_inside()) then
       call this%domain%remap()
    end if

    ! if we don't have more than 3 cells per dimension, don't use cells
    ! and return

    do i = 1, 3
       if (this%ncell(i) < 3) then
          this%cells_exist = .false.
          return
       end if
    end do

    x => this%particle%x
    mask => this%group%mask
    nparticles = this%particle%nparticles

    ! cellhead gets a zero (less than minimum particle index of 1)
    ! if no particle is inside

    this%cellhead(:) = 0
    this%next(:) = 0

    ! convert to lamda coords for assigning cells if triclinic
    ! TODO too many remaps

    if (this%domain%triclinic) then
       call this%domain%x2lamda()
    end if

    ! bin in reverse order, so that binhead->next->next increases in particle id

    if (this%includegroup > 0) then
       bitmask = this%group%bitmask(this%includegroup)

       do i = nparticles, 1, -1
          if (iand(mask(i), bitmask) > 0) then
             icell = this%coord2cell(x(:,i))
             this%next(i) = this%cellhead(icell)
             this%cellhead(icell) = i
          end if
       end do
    else ! loop over all particles, 
       do i = nparticles, 1, -1
          icell = this%coord2cell(x(:,i))
          this%next(i) = this%cellhead(icell)
          this%cellhead(icell) = i
       end do
    end if

    ! convert back if triclinic

    if (this%domain%triclinic) then
       call this%domain%lamda2x()
    end if

  end subroutine assign_particles

  ! ----------------------------------------------------------------------
  ! Determine the number of cells used for a calculation.
  !
  ! In real space, we want to select cell dimensions such that
  ! we don't have to look beyond the first shell for particles.
  ! That means that cell planes have to be bigger than some spherical
  ! diameter apart in all directions.  This means that we have to find
  ! the minimum triclinic cell s.t. it can fit a circle with diameter
  ! equal to the target size.
  !
  ! In reciprocal space, we want to find the smallest triclinic cell
  ! that can fit inside a sphere with diameter equal to at most
  ! (1/2)*(2*pi/q_{max}).  I didn't carry this out exactly.  What I did
  ! was ensure that each triclinic cell vector has less than a given
  ! length (see plot_ssf.f90) associated with the Nyquist criterion
  ! and used fudge factors that gave extra sampling, enough that
  ! there were no visible aliasing artifacts visible during shear flow.
  ! (No artifacts produced by the limited box distortion that LAMMPS allows.)
  ! This may require more work to obtain a more optimal number of grid
  ! points.
  ! ----------------------------------------------------------------------

  subroutine determine_ncell(domain, target_size, stay_under, ncell)

    type(domain_t), intent(in)         :: domain
    real, intent(in)                   :: target_size
    logical, intent(in)                :: stay_under
    integer, dimension(3), intent(out) :: ncell
    real, dimension(3)                 :: v
    real                               :: xprd,yprd,zprd,xy,xz,yz

    xprd = domain%xprd
    yprd = domain%yprd
    zprd = domain%zprd
    xy = domain%xy
    xz = domain%xz
    yz = domain%yz

    if (stay_under) then

       ! Get grid spacings along a,b,c with lengths smaller than target size.
       ! Not sure if this is the most appropriate way of doing things!
       ! It might be better to find the triclinic cell with maximum volume
       ! that fits inside a sphere with diameter = target_size.
       ! Might be good to consult the imaging literature on a multi-dimensional
       ! Nyquist criterion with a distorted grid.

       v(1) = xprd
       v(2) = sqrt(xy*xy + yprd*yprd)
       v(3) = sqrt(xz*xz + yz*yz + zprd*zprd)
       ncell(:) = ceiling(v(:) / target_size)
    else

       ! Use cross products between triclinic box vectors to get spacing
       ! between the planes spanned by these vectors (e.g. the length
       ! between bounding planes spanned by b and c vectors).

       v(1) = (xy*yz - xz*yprd)
       v(1) = v(1)*v(1)
       v(1) = v(1) + (yprd*yprd + xy*xy)*zprd*zprd
       v(1) = sqrt(v(1))
       v(2) = xprd*sqrt(zprd*zprd + yz*yz)
       v(3) = xprd*yprd
       v(:) = xprd*yprd*zprd / target_size / v(:)
       ncell(:) = floor(v(:))
    end if

  end subroutine determine_ncell

  ! ----------------------------------------------------------------------
  ! Return the ith stencil neighbor for cell b
  ! Neighbors are all in the +x,+y,+z directions, to avoid double counting
  ! The sequence is
  ! 1 (+x)
  ! 2 (+y)
  ! 3 (+z)
  ! 4 (+x)(+y)
  ! 5 (-x)(+y)
  ! 6 (+x)(+z)
  ! 7 (-x)(+z)
  ! 8 (+y)(+z)
  ! 9 (-y)(+z)
  ! 10 (+x)(+y)(+z)
  ! 11 (-x)(+y)(+z)
  ! 12 (+x)(-y)(+z)
  ! 13 (+x)(+y)(-z)
  ! Returns 0 if this would cross a non-periodic boundary
  ! ----------------------------------------------------------------------

  integer function stencil_neighbor(this, b, i)

    class(cells_t), intent(in) :: this
    integer, intent(in)        :: b,i
    integer                    :: ncellx,ncelly,ncellz
    logical                    :: xperiodic,yperiodic,zperiodic
    integer, dimension(3)      :: triplet,offset

    ncellx = this%ncellx
    ncelly = this%ncelly
    ncellz = this%ncellz
    xperiodic = this%domain%xperiodic
    yperiodic = this%domain%yperiodic
    zperiodic = this%domain%zperiodic

    triplet = this%cell2triplet(b)

    select case (i)
    case (1)
       offset = (/ +1,  0,  0 /)
    case(2)
       offset = (/  0, +1,  0 /)
    case(3)
       offset = (/  0,  0, +1 /)
    case(4)
       offset = (/ +1, +1,  0 /)
    case(5)
       offset = (/ -1, +1,  0 /)
    case(6)
       offset = (/ +1,  0, +1 /)
    case(7)
       offset = (/ -1,  0, +1 /)
    case(8)
       offset = (/  0, +1, +1 /)
    case(9)
       offset = (/  0, -1, +1 /)
    case(10)
       offset = (/ +1, +1, +1 /)
    case(11)
       offset = (/ -1, +1, +1 /)
    case(12)
       offset = (/ +1, -1, +1 /)
    case(13)
       offset = (/ +1, +1, -1 /)
    end select

    triplet(:) = modulo((triplet(:) + offset(:) - 1), &
         (/ ncellx, ncelly, ncellz /)) + 1

    stencil_neighbor = this%triplet2cell(triplet)

  end function stencil_neighbor

  function cell2triplet (this, n) result (trip)

    class(cells_t), intent(in) :: this
    integer, intent(in) :: n
    integer, dimension(3) :: trip
    integer :: ncellx,ncelly

    ncellx = this%ncellx
    ncelly = this%ncelly

    trip(3) = (n-1)/ncellx/ncelly + 1
    trip(2) = modulo(n-1,ncellx*ncelly) / ncellx + 1
    trip(1) = modulo(modulo(n-1,ncellx*ncelly),ncellx) + 1    

  end function cell2triplet

  function triplet2cell (this, trip) result (n)

    class(cells_t), intent(in) :: this
    integer, dimension(3), intent(in) :: trip
    integer :: n,ncellx,ncelly

    ncellx = this%ncellx
    ncelly = this%ncelly
    n = (trip(3)-1)*ncellx*ncelly &
                     + (trip(2)-1)*ncellx &
                     + trip(1)
  end function triplet2cell

end module mod_cells
