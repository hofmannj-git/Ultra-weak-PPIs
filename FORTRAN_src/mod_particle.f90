module mod_particle

  implicit none
  private

  integer, public, parameter :: IMGMASK = 1023
  integer, public, parameter :: IMGMAX = 512
  integer, public, parameter :: IMGBITS = 10
  integer, public, parameter :: IMG2BITS = 20

  type, public :: particle_t

     integer :: nparticles
     integer :: ntypes

     integer, dimension(:), allocatable :: type
     integer, dimension(:), allocatable :: molec
     integer, dimension(:), allocatable :: image
     real, dimension(:,:), allocatable  :: x, v, quat, wbody, euler
     real, dimension(:), allocatable    :: radius

     ! optional arrays

     real, dimension(:), allocatable :: typeradius
     real, dimension(:,:), allocatable  :: omega
     real, dimension(:),   allocatable  :: pe
     real, dimension(:,:), allocatable  :: stress

     logical :: omega_flag
     logical :: pe_flag
     logical :: stress_flag
     logical :: typeradius_flag

   contains

     procedure :: init => particle_init
     procedure :: free => particle_free
     procedure :: copy_one
     procedure :: copy
     procedure :: zero_image

  end type particle_t

contains

  ! ----------------------------------------------------------------------
  ! allocate and set flags
  ! ----------------------------------------------------------------------

  subroutine particle_init(this, nparticles, &
       omega_flag, pe_flag, stress_flag, typeradius_flag)

    class(particle_t), intent(inout)  :: this
    integer, intent(in)               :: nparticles
    logical, intent(in), optional     :: omega_flag, pe_flag, stress_flag
    logical, intent(in), optional     :: typeradius_flag

    this%nparticles = nparticles
    this%omega_flag = .false.
    this%pe_flag = .false.
    this%stress_flag = .false.
    this%typeradius_flag = .false.
    if (present(omega_flag)) then
       this%omega_flag = omega_flag
    end if
    if (present(pe_flag)) then
       this%pe_flag = pe_flag
    end if
    if (present(stress_flag)) then
       this%stress_flag = stress_flag
    end if
    if (present(typeradius_flag)) then
       this%typeradius_flag = typeradius_flag
    end if

    allocate(this%type(nparticles))
    this%type(:) = 1
    allocate(this%molec(nparticles))
    this%molec(:) = 0
    this%ntypes = 1

    ! For image, not initializing to zero but the zero image

    allocate(this%image(nparticles))
    call this%zero_image()

    allocate(this%x(3,nparticles))
    this%x(:,:) = 0.0
    allocate(this%v(3,nparticles))
    this%v(:,:) = 0.0
    allocate(this%wbody(3,nparticles))
    this%wbody(:,:) = 0.0
    allocate(this%quat(4,nparticles))
    this%quat(:,:) = 0.0
    allocate(this%euler(3,nparticles))
    this%euler(:,:) = 0.0
    allocate(this%radius(nparticles))
    this%radius(:) = 0.0
    if (this%omega_flag) then
       allocate(this%omega(3,nparticles))
       this%omega(:,:) = 0.0
    end if
    if (this%pe_flag) then
       allocate(this%pe(nparticles))
       this%pe(:) = 0.0
    end if
    if (this%stress_flag) then
       allocate(this%stress(6,nparticles))
       this%stress(:,:) = 0.0
    end if

  end subroutine particle_init

  ! ----------------------------------------------------------------------
  ! free the particle type (necessary if have allocation?)
  ! ----------------------------------------------------------------------

  elemental subroutine particle_free(this)

    class(particle_t), intent(inout) :: this

    deallocate(this%image)
    deallocate(this%molec)
    deallocate(this%x)
    deallocate(this%v)
    deallocate(this%wbody)
    deallocate(this%quat)
    deallocate(this%euler)
    deallocate(this%radius)
    if (allocated(this%omega)) then
       deallocate(this%omega)
    end if
    if (allocated(this%pe)) then
       deallocate(this%pe)
    end if
    if (allocated(this%stress)) then
       deallocate(this%stress)
    end if

  end subroutine particle_free

  ! ----------------------------------------------------------------------
  ! copy one particle's info to another particle index
  ! ----------------------------------------------------------------------

  subroutine copy_one(this, i, j)

    class(particle_t), intent(inout) :: this
    integer, intent(in)              :: i,j

    this%type(j) = this%type(i)
    this%molec(j) = this%molec(i)
    this%image(j) = this%image(i)
    this%x(:,j) = this%x(:,i)
    this%v(:,j) = this%v(:,i)
    this%wbody(:,j) = this%wbody(:,i)
    this%quat(:,j) = this%quat(:,i)
    this%euler(:,j) = this%euler(:,i)
    this%radius(j) = this%radius(i)
    if (this%omega_flag) then
       this%omega(:,j) = this%omega(:,i)
    end if
    if (this%pe_flag) then
       this%pe(j) = this%pe(i)
    end if
    if (this%stress_flag) then
       this%stress(:,j) = this%stress(:,i)
    end if

  end subroutine copy_one

  ! ----------------------------------------------------------------------
  ! use this particle type to create an identical particle type
  ! if the particle type has allocated variables and the same
  ! ntypes and nparticles, then don't reallocate
  ! ----------------------------------------------------------------------

  subroutine copy(this, that)

    class(particle_t), intent(in) :: this
    type(particle_t), intent(inout) :: that
    integer :: nparticles,ntypes
    logical :: ntypes_different,nparticles_different

    ntypes_different = .false.
    nparticles_different = .false.

    ! copy scalar variables from this to that

    nparticles = this%nparticles
    if (that%nparticles /= nparticles) then
       that%nparticles = nparticles
       nparticles_different = .true.
    end if
    ntypes = this%ntypes
    if (that%ntypes /= ntypes) then
       that%ntypes = ntypes
       ntypes_different = .true.
    end if
    that%omega_flag = this%omega_flag
    that%pe_flag = this%pe_flag
    that%stress_flag = this%stress_flag
    that%typeradius_flag = this%typeradius_flag

    ! allocate that's arrays only if different numbers of types/particles

    if (nparticles_different) then
       if (allocated(that%type)) deallocate(that%type)
       allocate(that%type(nparticles))
       if (allocated(that%molec)) deallocate(that%molec)
       allocate(that%molec(nparticles))
       if (allocated(that%image)) deallocate(that%image)
       allocate(that%image(nparticles))
       if (allocated(that%x)) deallocate(that%x)
       allocate(that%x(3,nparticles))
       if (allocated(that%v)) deallocate(that%v)
       allocate(that%v(3,nparticles))
       if (allocated(that%quat)) deallocate(that%quat)
       allocate(that%quat(4,nparticles))
       if (allocated(that%wbody)) deallocate(that%wbody)
       allocate(that%wbody(3,nparticles))
       if (allocated(that%euler)) deallocate(that%euler)
       allocate(that%euler(3,nparticles))
       if (allocated(that%radius)) deallocate(that%radius)
       allocate(that%radius(nparticles))
       if (that%omega_flag) then
          if (allocated(that%omega)) deallocate(that%omega)
          allocate(that%omega(3,nparticles))
       end if
       if (that%pe_flag) then
          if (allocated(that%pe)) deallocate(that%pe)
          allocate(that%pe(nparticles))
       end if
       if (that%stress_flag) then
          if (allocated(that%stress)) deallocate(that%stress)
          allocate(that%stress(6,nparticles))
       end if
    end if
    if (ntypes_different) then
       if (that%typeradius_flag) then
          if (allocated(that%typeradius)) deallocate(that%typeradius)
          allocate(that%typeradius(ntypes))
       end if
    end if

    ! copy data

    that%type(:) = this%type(:)
    that%molec(:) = this%molec(:)
    that%image(:) = this%image(:)
    that%x(:,:) = this%x(:,:)
    that%v(:,:) = this%v(:,:)
    that%quat(:,:) = this%quat(:,:)
    that%wbody(:,:) = this%wbody(:,:)
    that%euler(:,:) = this%euler(:,:)
    that%radius(:) = this%radius(:)

    if (that%typeradius_flag) that%typeradius(:) = this%typeradius(:)
    if (that%omega_flag) that%omega(:,:) = this%omega(:,:)
    if (that%pe_flag) that%pe(:) = this%pe(:)
    if (that%stress_flag) that%stress(:,:) = this%stress(:,:)

  end subroutine copy

  ! ----------------------------------------------------------------------
  ! zero the images
  ! useful after unmapping, as another remap with increment the images
  ! ----------------------------------------------------------------------

  subroutine zero_image(this)

    class(particle_t), intent(inout) :: this

    this%image(:) = ior( &
                      ior(ishft(IMGMAX,IMG2BITS), ishft(IMGMAX,IMGBITS)), &
                      IMGMAX)

  end subroutine zero_image

end module mod_particle
