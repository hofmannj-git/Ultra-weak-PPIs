! groups are static, unlike in LAMMPS
! mask for particles has been moved here

module mod_group

  use mod_particle, only : particle_t
  use mod_domain, only : domain_t
  implicit none
  private

  integer, parameter, public :: MAX_GROUP = 32

  type, public :: excludes_t
     character(:), dimension(:), allocatable :: style
     integer, dimension(:), allocatable :: m
     integer, dimension(:), allocatable :: n
     integer :: nex
  end type excludes_t

  type, public :: group_t
     integer :: ngroup
     integer, dimension(:), allocatable :: bitmask
     type(particle_t), pointer :: particle => null()
     type(domain_t), pointer :: domain => null()
     integer, dimension(:), allocatable :: mask
   contains
     procedure, public :: count => count_group
     procedure, public :: init
     procedure, public :: build_ex_bits
     procedure, public :: build_ex_type
     procedure, public :: create
     procedure, public :: update
  end type group_t

contains

  ! ----------------------------------------------------------------------
  ! allocate and set flags
  ! ----------------------------------------------------------------------

  subroutine init(this, domain)

    class(group_t), intent(out)         :: this
    class(domain_t), intent(in), target :: domain
    integer :: nparticles,i

    this%domain => domain
    this%particle => domain%particle
    nparticles = this%particle%nparticles

    allocate(this%bitmask(MAX_GROUP))
    do i = 1, MAX_GROUP
       this%bitmask(i) = ishft(1,i-1)
    end do

    allocate(this%mask(nparticles))

    ! need to have all particles in group 1 (all)

    this%ngroup = 1
    this%mask(:) = 1

  end subroutine init

  ! ----------------------------------------------------------------------
  ! build type-based exclusion list
  ! ----------------------------------------------------------------------

  subroutine build_ex_type(this, pairs, ex_type)

    class(group_t), intent(in)           :: this
    integer, dimension(:,:), intent(in)  :: pairs
    logical, dimension(:,:), intent(out) :: ex_type
    integer                              :: nex_type,i,j

    ! check sizes

    nex_type = size(pairs,2)
    if (size(pairs,1) /= 2) then
       print *, 'exclusion needs pairs'
       stop
    end if
    if (nex_type < 1) then
       print *, 'exclusion needs at least one pair'
       stop
    end if
    if ((size(ex_type,1) /= size(ex_type,2)) .or. &
         (size(ex_type,1) /= nex_type)) then
       print *, 'bad dimensions for type exclusion'
       stop
    end if

    ! check type existence, and assign ex_type if legit

    ex_type = .false.
    do j = 1, nex_type
       do i = 1, 2
          if ((pairs(i,j) < 1) .or. (pairs(i,j) > this%particle%ntypes)) then
             print *, 'bad type for exclusion'
             stop
          end if
       end do
       ex_type(pairs(1,j),pairs(2,j)) = .true.
    end do

  end subroutine build_ex_type

  ! ----------------------------------------------------------------------
  ! build group-based exclusions
  ! one pair of groupbits per exclusion pair
  ! ----------------------------------------------------------------------

  subroutine build_ex_bits(this, pairs, ex1_bit, ex2_bit)

    class(group_t), intent(in)          :: this
    integer, dimension(:,:), intent(in) :: pairs
    integer, dimension(:), intent(out)  :: ex1_bit,ex2_bit
    integer                             :: nex_group,i,j

    ! check dimensions

    nex_group = size(pairs,2)
    if (size(pairs,1) /= 2) then
       print *, 'exclusion needs pairs'
       stop
    end if
    if (nex_group < 1) then
       print *, 'exclusion needs at least one pair'
       stop
    end if
    if ((size(ex1_bit) /= size(ex2_bit)) .or. (size(ex1_bit)) /= nex_group) then
       print *, 'bad dimensions for group exclusion'
       stop
    end if

    ! check group existence, and find bitmask for each pair

    do j = 1, nex_group
       do i = 1, 2
          if ((pairs(i,j) < 1) .or. (pairs(i,j) > this%ngroup)) then
             print *, 'could not find group'
             stop
          end if
       end do
       ex1_bit(j) = this%bitmask(pairs(1,j))
       ex2_bit(j) = this%bitmask(pairs(2,j))
    end do

  end subroutine build_ex_bits

  ! ----------------------------------------------------------------------
  ! Add flagged particles to a group
  ! ----------------------------------------------------------------------

  integer function create(this, flag)

    class(group_t), intent(inout) :: this
    logical, dimension(:), intent(in) :: flag
    integer :: i,bit,igroup,nparticles

    nparticles = this%particle%nparticles

    ! exit with error if there are too many groups

    if (this%ngroup == MAX_GROUP) then
       print *, 'too many groups'
       stop
    end if

    ! incrementing gets current group id

    this%ngroup = this%ngroup + 1
    igroup = this%ngroup

    ! add particles to group whose flags are set

    bit = this%bitmask(igroup)
    do i = 1, nparticles
       if (flag(i)) this%mask(i) = ior(this%mask(i),bit)
    end do

    create = igroup

  end function create

  ! ----------------------------------------------------------------------
  ! Update members of an existing group
  ! ----------------------------------------------------------------------

  subroutine update(this, igroup, flag)

    class(group_t), intent(inout) :: this
    integer, intent(in) :: igroup
    logical, dimension(:), intent(in) :: flag
    integer :: nparticles,bit,i

    nparticles = this%particle%nparticles
    bit = this%bitmask(igroup)

    ! return error if dimension of flag is not equal to nparticles

    if (size(flag) .ne. nparticles) then
       write(*, '(a)') 'error: flag for group has improper size'
       stop
    end if

    ! verify that the group in question already exists

    if (igroup > this%ngroup) then
       write(*, '(a)') 'error: igroup not found'
       stop
    end if
    if (igroup < 1) then
       write(*, '(a)') 'error: igroup must be greater than 1'
       stop
    end if

    ! clear all bits (set to 0) for current group in mask

    do i = 1, nparticles
       this%mask(i) = ibclr(this%mask(i),igroup-1)
    end do

    ! add the particles as in the create function
    ! by setting the bits recently cleared to 1 if the particle
    ! is in the group

    do i = 1, nparticles
       if (flag(i)) this%mask(i) = ior(this%mask(i),bit)
    end do

  end subroutine update

  ! ----------------------------------------------------------------------
  ! Return number of particles in a group
  ! ----------------------------------------------------------------------

  integer function count_group(this, igroup)

    class(group_t), intent(in) :: this
    integer, intent(in) :: igroup
    integer :: groupbit,nparticles,i

    groupbit = this%bitmask(igroup)
    nparticles = this%particle%nparticles

    count_group = 0
    do i = 1, nparticles
       if (iand(this%mask(i),groupbit) /= 0) count_group = count_group + 1
    end do

  end function count_group

end module mod_group
