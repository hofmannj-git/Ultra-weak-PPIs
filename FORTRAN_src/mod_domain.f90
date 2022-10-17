module mod_domain

  use mod_particle, only : particle_t, IMGMASK, IMGMAX, IMGBITS, IMG2BITS
  implicit none
  private

  ! overloading .and. for logical and integer types
  ! in order to use it like LAMMPS's domain.cpp

  interface operator (.and.)
     module procedure log_and_int
     module procedure int_and_log
  end interface operator (.and.)

  public :: get_image_flip, image_flip, propagate_in_flow

  type, public :: domain_t
 
     logical :: box_exist                      ! true = domain created
     integer :: dimension                      ! 2 = 2d, 3 = 3d
                                               ! 2d not supported
     logical :: triclinic                      ! true for triclinic

     logical :: xperiodic,yperiodic,zperiodic  ! true = periodic
     logical, dimension(3) :: periodicity      ! xyz periodicity as array

     real :: xprd,yprd,zprd                    ! global box dimensions
     real :: xprd_half,yprd_half,zprd_half     ! half dimensions

     real, dimension(3) :: prd                 ! array form of dimensions
     real, dimension(3) :: prd_half
     real, dimension(3) :: prd_lamda
     real, dimension(3) :: boxlo,boxhi         ! orthogonal box global bounds
     real, dimension(3) :: boxlo_lamda
     real, dimension(3) :: boxhi_lamda
     real, dimension(3) :: boxlo_bound
     real, dimension(3) :: boxhi_bound

     real               :: xy,xz,yz            ! 3 tilt factors
     real, dimension(6) :: h,h_inv             ! shape matrix in Voigt notation
     real, dimension(6) :: h_rate              ! rate of box size/shape change

     type(particle_t), pointer :: particle

   contains

     procedure, public :: init
     procedure, public :: set_global_box
     procedure, public :: remap
     procedure, public :: unmap
     procedure, public :: x2lamda
     procedure, public :: lamda2x
     procedure, public :: minimum_image
     procedure, public :: particles_inside
     procedure, public :: get_flow_unmapped
     procedure, public :: volume

     procedure :: remap_one
     procedure :: unmap_one
     procedure :: lamda2x_one
     procedure :: x2lamda_one
     procedure :: swap_particle_pointer

  end type domain_t

contains

  ! ----------------------------------------------------------------------
  ! Set defaults for the domain type
  ! Domain does not "exist" unless set_global_box is called
  ! ----------------------------------------------------------------------

  subroutine init(this,particle)

    class(domain_t), intent(out)          :: this
    class(particle_t), intent(in), target :: particle

    this%particle => particle

    this%dimension = 3

    this%box_exist = .false.
    this%xperiodic = .true.
    this%yperiodic = .true.
    this%zperiodic = .true.
    this%periodicity = (/ .true., .true., .true. /)
    this%triclinic = .false.

    this%boxlo(1) = -0.5
    this%boxlo(2) = -0.5
    this%boxlo(3) = -0.5
    this%boxhi(1) = 0.5
    this%boxhi(2) = 0.5
    this%boxhi(3) = 0.5
    this%xy = 0.0
    this%xz = 0.0
    this%yz = 0.0

    this%h(4) = 0.0
    this%h(5) = 0.0
    this%h(6) = 0.0
    this%h_inv(4) = 0.0
    this%h_inv(5) = 0.0
    this%h_inv(6) = 0.0
    this%h_rate(:) = 0.0

    this%prd_lamda = (/ 1.0, 1.0, 1.0 /)
    this%boxlo_lamda = (/ 0.0, 0.0, 0.0 /)
    this%boxhi_lamda = (/ 1.0, 1.0, 1.0 /)

  end subroutine init

  ! ----------------------------------------------------------------------
  ! Set or reset the global box (domain)
  ! ----------------------------------------------------------------------

  subroutine set_global_box(this,boxlo,boxhi,xy,xz,yz,h_rate)

    class(domain_t), intent(inout)           :: this
    real, dimension(3), intent(in), optional :: boxlo,boxhi
    real, intent(in), optional               :: xy,xz,yz
    real, dimension(6), intent(in), optional :: h_rate
    integer :: i

    ! check if inputs are legit

    if (present(boxlo) .or. present(boxhi)) then
       do i = 1, 3
          if (boxlo(i) > boxhi(i)) then
             print *, 'boxhi must be greater than boxlo'
             stop
          end if
       end do
    end if

    ! change to triclinic if not already and xy,xz,yz specified

    if ((.not. this%triclinic) .and. &
        (present(xy) .or. present(xz) .or. present(yz)) &
       ) then
       this%triclinic = .true.
    end if

    ! set box variables

    if (present(boxlo)) this%boxlo = boxlo
    if (present(boxhi)) this%boxhi = boxhi
    if (present(xy)) this%xy = xy
    if (present(xz)) this%xz = xz
    if (present(yz)) this%yz = yz
    if (present(h_rate)) this%h_rate = h_rate

    ! calculate all secondary box variables
    ! slightly inefficient, but I didn't want to introduce a ton of switches

    this%xprd = this%boxhi(1) - this%boxlo(1)
    this%yprd = this%boxhi(2) - this%boxlo(2)
    this%zprd = this%boxhi(3) - this%boxlo(3)
    this%prd = (/ this%xprd, this%yprd, this%zprd /)

    this%h(1) = this%xprd
    this%h(2) = this%yprd
    this%h(3) = this%zprd
    this%h_inv(1) = 1.0/this%h(1)
    this%h_inv(2) = 1.0/this%h(2)
    this%h_inv(3) = 1.0/this%h(3)

    this%xprd_half = 0.5*this%xprd
    this%yprd_half = 0.5*this%yprd
    this%zprd_half = 0.5*this%zprd
    this%prd_half = (/ this%xprd_half, this%yprd_half, this%zprd_half /)

    if (this%triclinic) then
       this%h(4) = this%yz
       this%h(5) = this%xz
       this%h(6) = this%xy
       this%h_inv(4) = -this%h(4) / (this%h(2)*this%h(3))
       this%h_inv(5) = (this%h(4)*this%h(6) - this%h(2)*this%h(5)) &
                       / (this%h(1)*this%h(2)*this%h(3))
       this%h_inv(6) = -this%h(6) / (this%h(1)*this%h(2))

       this%boxlo_bound(1) = min(this%boxlo(1), &
                                 this%boxlo(1)+this%xy)
       this%boxlo_bound(1) = min(this%boxlo_bound(1), &
                                 this%boxlo_bound(1)+this%xz)
       this%boxlo_bound(2) = min(this%boxlo(2), &
                                 this%boxlo(2)+this%yz)
       this%boxlo_bound(3) = this%boxlo(3)

       this%boxhi_bound(1) = max(this%boxhi(1), &
                                 this%boxhi(1)+this%xy)
       this%boxhi_bound(1) = max(this%boxhi_bound(1), &
                                 this%boxhi_bound(1)+this%xz)
       this%boxhi_bound(2) = max(this%boxhi(2), &
                                 this%boxhi(2)+this%yz)
       this%boxhi_bound(3) = this%boxhi(3)
    end if

    ! if "not existing", the box now "exists"

    if (.not. this%box_exist) this%box_exist = .true.

  end subroutine set_global_box

  ! ----------------------------------------------------------------------
  ! remap the point into the periodic box no matter how far away
  ! adjust 3 image flags encoded in image accordingly
  ! resulting coord must satisfy lo <= coord < hi
  ! MAX is important since coord - prd < lo can happen when coord = hi
  ! for triclinic, point is converted to lamda coords (0-1) before doing remap
  ! image = 10 bits for each dimension
  ! increment/decrement in wrap-around fashion
  ! ----------------------------------------------------------------------

  subroutine remap_one(this, x, image)

    class(domain_t), intent(in), target       :: this
    real, dimension(3), intent(inout), target :: x
    integer, intent(inout)                    :: image
    real, dimension(3)                        :: lo,hi,period
    real, dimension(:), pointer               :: coord
    real, dimension(3), target                :: lamda
    integer                                   :: idim,otherdims ! imageint?

    if (.not. this%triclinic) then
       lo = this%boxlo
       hi = this%boxhi
       period = this%prd
       coord => x
    else
       lo = this%boxlo_lamda
       hi = this%boxhi_lamda
       period = this%prd_lamda
       call this%x2lamda_one(x,lamda)
       coord => lamda
    end if

    yes_xperiodic: if (this%xperiodic) then
       do while (coord(1) < lo(1))
          coord(1) = coord(1) + period(1)
          idim = iand(image,IMGMASK)
          otherdims = ieor(image,idim)
          idim = idim - 1
          idim = iand(idim,IMGMASK)
          image = ior(otherdims,idim)
       end do
       do while (coord(1) >= hi(1))
          coord(1) = coord(1) - period(1)
          idim = iand(image,IMGMASK)
          otherdims = ieor(image,idim)
          idim = idim + 1
          idim = iand(idim,IMGMASK)
          image = ior(otherdims,idim)
       end do
       coord(1) = max(coord(1),lo(1))
    end if yes_xperiodic

    yes_yperiodic: if (this%yperiodic) then
       do while (coord(2) < lo(2))
          coord(2) = coord(2) + period(2)
          idim = iand(ishft(image,-IMGBITS),IMGMASK)
          otherdims = ieor(image,ishft(idim,IMGBITS))
          idim = idim - 1
          idim = iand(idim,IMGMASK)
          image = ior(otherdims,ishft(idim,IMGBITS))
       end do
       do while (coord(2) >= hi(2))
          coord(2) = coord(2) - period(2)
          idim = iand(ishft(image,-IMGBITS),IMGMASK)
          otherdims = ieor(image,ishft(idim,IMGBITS))
          idim = idim + 1
          idim = iand(idim,IMGMASK)
          image = ior(otherdims,ishft(idim,IMGBITS))
       end do
       coord(2) = max(coord(2),lo(2))
    end if yes_yperiodic

    yes_zperiodic: if (this%zperiodic) then
       do while (coord(3) < lo(3))
          coord(3) = coord(3) + period(3)
          idim = ishft(image,-IMG2BITS)
          otherdims = ieor(image,ishft(idim,IMG2BITS))
          idim = idim - 1
          idim = iand(idim,IMGMASK)
          image = ior(otherdims,ishft(idim,IMG2BITS))
       end do
       do while (coord(3) >= hi(3))
          coord(3) = coord(3) - period(3)
          idim = ishft(image,-IMG2BITS)
          otherdims = ieor(image,ishft(idim,IMG2BITS))
          idim = idim + 1
          idim = iand(idim,IMGMASK)
          image = ior(otherdims,ishft(idim,IMG2BITS))
       end do
       coord(3) = max(coord(3),lo(3))
    end if yes_zperiodic

    if (this%triclinic) call this%lamda2x_one(coord,x)

  end subroutine remap_one

  ! ----------------------------------------------------------------------
  ! remap all particles back into simulation box, setting images
  ! ----------------------------------------------------------------------

  subroutine remap(this)

    class(domain_t), intent(inout) :: this
    integer                        :: i

    do i = 1, this%particle%nparticles
       call this%remap_one(this%particle%x(:,i), this%particle%image(i))
    end do

  end subroutine remap

  ! ----------------------------------------------------------------------
  ! unmap the point via image flags
  ! x overwritten with result, don't reset image flag
  ! for triclinic, use h[] to add in tilt factors in other dims as needed
  ! ----------------------------------------------------------------------

  subroutine unmap_one(this, x, image)

    class(domain_t), intent(in), target :: this
    real, dimension(3), intent(inout)   :: x
    integer, intent(in)                 :: image
    integer                             :: xbox,ybox,zbox

    xbox = iand(image,IMGMASK) - IMGMAX
    ybox = iand(ishft(image,-IMGBITS),IMGMASK) - IMGMAX
    zbox = ishft(image,-IMG2BITS) - IMGMAX

    if (.not. this%triclinic) then
       x(1) = x(1) + xbox * this%xprd
       x(2) = x(2) + ybox * this%yprd
       x(3) = x(3) + zbox * this%zprd
    else
       x(1) = x(1) + this%h(1)*xbox + this%h(6)*ybox + this%h(5)*zbox
       x(2) = x(2) + this%h(2)*ybox + this%h(4)*zbox
       x(3) = x(3) + this%h(3)*zbox
    end if

  end subroutine unmap_one

  ! ----------------------------------------------------------------------
  ! Unmap (unwrap) all the particles
  ! IMPORTANT: unlike in LAMMPS, this zeros all particle images
  ! ----------------------------------------------------------------------

  subroutine unmap(this)

    class(domain_t), intent(inout) :: this
    integer :: i

    do i = 1, this%particle%nparticles
       call this%unmap_one(this%particle%x(:,i), this%particle%image(i))
    end do

    call this%particle%zero_image()

  end subroutine unmap

  ! ----------------------------------------------------------------------
  ! convert box coords to triclinic 0-1 lamda coords for all N particles
  ! lamda = H^-1 (x - x0)
  ! ----------------------------------------------------------------------

  subroutine x2lamda(this)

    class(domain_t), intent(in), target :: this
    real, dimension(3)                  :: delta
    integer                             :: i
    real, dimension(:,:), pointer       :: x => null()
    real, dimension(:), pointer         :: h_inv => null()
    real, dimension(:), pointer         :: boxlo => null()

    x => this%particle%x
    h_inv => this%h_inv
    boxlo => this%boxlo

    do i = 1, this%particle%nparticles
       delta(:) = x(:,i) - boxlo(:)
       x(1,i) = h_inv(1)*delta(1) + h_inv(6)*delta(2) + h_inv(5)*delta(3)
       x(2,i) = h_inv(2)*delta(2) + h_inv(4)*delta(3)
       x(3,i) = h_inv(3)*delta(3)
    end do

  end subroutine x2lamda

  ! ----------------------------------------------------------------------
  ! convert box coords to triclinic 0-1 lamda coords for one particle
  ! use my_boxlo & my_h_inv stored by caller for previous state of box
  ! lamda = H^-1 (x - x0)
  ! x and lamda can point to same 3-vector
  ! ----------------------------------------------------------------------

  pure subroutine x2lamda_one(this, x, lamda)

    class(domain_t), intent(in), target :: this
    real, dimension(3), intent(in) :: x
    real, dimension(3), intent(out) :: lamda
    real, dimension(3) :: delta

    delta(1) = x(1) - this%boxlo(1)
    delta(2) = x(2) - this%boxlo(2)
    delta(3) = x(3) - this%boxlo(3)

    lamda(1) = &
         this%h_inv(1)*delta(1) + this%h_inv(6)*delta(2) + this%h_inv(5)*delta(3)
    lamda(2) = &
         this%h_inv(2)*delta(2) + this%h_inv(4)*delta(3)
    lamda(3) = &
         this%h_inv(3)*delta(3)

  end subroutine x2lamda_one

  ! ----------------------------------------------------------------------
  ! convert triclinic 0-1 lamda coords to box coords for all particles
  ! x = H lamda + x0
  ! ----------------------------------------------------------------------

  subroutine lamda2x(this)

    class(domain_t), intent(in), target :: this
    integer                             :: i
    real, dimension(:,:), pointer       :: x => null()
    real, dimension(:), pointer         :: boxlo => null()
    real, dimension(:), pointer         :: h => null()

    x => this%particle%x
    h => this%h
    boxlo => this%boxlo

    do i = 1, this%particle%nparticles
       x(1,i) = h(1)*x(1,i) + h(6)*x(2,i) + h(5)*x(3,i) + boxlo(1)
       x(2,i) = h(2)*x(2,i) + h(4)*x(3,i) + boxlo(2)
       x(3,i) = h(3)*x(3,i) + boxlo(3)
    end do

  end subroutine lamda2x

  ! ----------------------------------------------------------------------
  ! convert triclinic 0-1 lamda coords to box coords for a single particle
  ! ----------------------------------------------------------------------

  pure subroutine lamda2x_one(this, lamda, x)

    class(domain_t), intent(in), target :: this
    real, dimension(3), intent(in)      :: lamda
    real, dimension(3), intent(out)     :: x
    real, dimension(6) :: h
    real, dimension(3) :: boxlo

    h = this%h
    boxlo = this%boxlo

    x(1) = h(1)*lamda(1) + h(6)*lamda(2) + h(5)*lamda(3) + boxlo(1)
    x(2) = h(2)*lamda(2) + h(4)*lamda(3) + boxlo(2)
    x(3) = h(3)*lamda(3) + boxlo(3)

  end subroutine lamda2x_one

  ! ----------------------------------------------------------------------
  ! minimum image convention
  ! use 1/2 of box size as test
  ! for triclinic, also add/subtract tilt factors in other dims as needed
  ! ----------------------------------------------------------------------

  pure subroutine minimum_image(this, dx, dy, dz)

    class(domain_t), intent(in) :: this
    real, intent(inout)         :: dx,dy,dz

    not_triclinic: if (.not. this%triclinic) then
       if (this%xperiodic) then
          if (abs(dx) > this%xprd_half) then
             if (dx < 0.0) then
                dx = dx + this%xprd
             else
                dx = dx - this%xprd
             end if
          end if
       end if
       if (this%yperiodic) then
          if (abs(dy) > this%yprd_half) then
             if (dy < 0.0) then
                dy = dy + this%yprd
             else
                dy = dy - this%yprd
             end if
          end if
       end if
       if (this%zperiodic) then
          if (abs(dz) > this%zprd_half) then
             if (dz < 0.0) then
                dz = dz + this%zprd
             else
                dz = dz - this%zprd
             end if
          end if
       end if
    else ! triclinic
       if (this%zperiodic) then
          if (abs(dz) > this%zprd_half) then
             if (dz < 0.0) then
                dz = dz + this%zprd
                dy = dy + this%yz
                dx = dx + this%xz
             else
                dz = dz - this%zprd
                dy = dy - this%yz
                dx = dx - this%xz
             end if
          end if
       end if
       if (this%yperiodic) then
          if (abs(dy) > this%yprd_half) then
             if (dy < 0.0) then
                dy = dy + this%yprd
                dx = dx + this%xy
             else
                dy = dy - this%yprd
                dx = dx - this%xy
             end if
          end if
       end if
       if (this%xperiodic) then
          if (abs(dx) > this%xprd_half) then
             if (dx < 0.0) then
                dx = dx + this%xprd
             else
                dx = dx - this%xprd
             end if
          end if
       end if
    end if not_triclinic

  end subroutine minimum_image

  ! ----------------------------------------------------------------------
  ! test if all particles are inside the box
  ! need to run when not in lamda coords!
  ! ----------------------------------------------------------------------

  elemental logical function particles_inside(this)

    class(domain_t), intent(in) :: this
    integer                     :: i,j
    real, dimension(3)          :: coord

    if (this%triclinic) then
       do i = 1, this%particle%nparticles
          call this%x2lamda_one(this%particle%x(:,i),coord)
          do j = 1, 3
             if ((coord(j) < this%boxlo_lamda(j)) .or. &
                 (coord(j) >= this%boxhi_lamda(j))) then
                particles_inside = .false.
                return                
             end if
          end do
       end do
    else ! not triclinic
       do i = 1, this%particle%nparticles
          do j = 1, 3
             if ((this%particle%x(j,i) <  this%boxlo(j)) .or. &
                 (this%particle%x(j,i) >= this%boxhi(j))) then
                particles_inside = .false.
                return
             end if
          end do
       end do
    end if

    particles_inside = .true.

  end function particles_inside

  ! ----------------------------------------------------------------------
  ! Get the likeliest image flip between snapshots i and j, given their
  ! box bounds.  Should return the correct value if the snapshots are close
  ! together, but not if the box is changing rapidly between snapshots
  ! ----------------------------------------------------------------------

  subroutine get_image_flip(bounds_a, bounds_b, flip)

    real, dimension(9), intent(in) :: bounds_a,bounds_b
    integer, dimension(3), intent(out) :: flip
    real :: rxy_a,rxz_a,ryz_a,rxy_b,rxz_b,ryz_b
    real, dimension(6) :: h_a,h_b
    integer :: flipxy,flipxz,flipyz

    flipxy = 0
    flipxz = 0
    flipyz = 0

    ! get the tilt ratios

    call bounds2h(bounds_a,h_a)
    call bounds2h(bounds_b,h_b)
    rxy_a = h_a(6) / h_a(1)
    rxy_b = h_b(6) / h_b(1)
    rxz_a = h_a(5) / h_a(1)
    rxz_b = h_b(5) / h_b(1)
    ryz_a = h_a(4) / h_a(2)
    ryz_b = h_b(4) / h_b(2)

    if (abs(rxy_b - rxy_a) > 0.5) then
       if (rxy_b > rxy_a) then
          flipxy = 1
       else
          flipxy = -1
       end if
    end if

    if (abs(ryz_b - ryz_a) > 0.5) then
       if (ryz_b > ryz_a) then
          flipyz = 1
       else
          flipyz = -1
       end if
    end if

    if (abs(rxz_b - rxz_a - rxy_a*flipyz) > 0.5) then
       if (rxz_b > rxz_a + rxy_a*flipyz) then
          flipxz = 1
       else
          flipxz = -1
       end if
    end if

    flip = (/ flipxy, flipxz, flipyz /)

  end subroutine get_image_flip

  ! propagate particles from domain_init to domain_fini through box_flips
  ! coords start unwrapped and are output as unwrapped

  subroutine propagate_in_flow(domain_init, domain_fini, box_flips, particle_flow)

    type(domain_t), intent(inout) :: domain_init,domain_fini
    integer, dimension(:), intent(in) :: box_flips
    type(particle_t), pointer :: particle_init,particle_fini
    type(particle_t), intent(inout) :: particle_flow
    real, dimension(3) :: boxlo_init,boxlo_fini,boxhi_fini
    real :: xy,xz,yz,xy_fini,xz_fini,yz_fini,xprd_fini,yprd_fini
    integer :: nparticles,flipxy,flipxz,flipyz

    nparticles = domain_init%particle%nparticles
    boxlo_fini(:) = domain_fini%boxlo(:)
    boxlo_init(:) = domain_init%boxlo(:)
    boxhi_fini(:) = domain_fini%boxhi(:)
    xy_fini = domain_fini%xy
    xz_fini = domain_fini%xz
    yz_fini = domain_fini%yz
    flipxy = box_flips(1)
    flipxz = box_flips(2)
    flipyz = box_flips(3)
    xprd_fini = domain_fini%xprd
    yprd_fini = domain_fini%yprd

    ! first, simply copy the particle coords to the flow particles

    call domain_init%particle%copy(particle_flow)

    ! juggle the particle pointer and use x2lamda with
    ! domain_init on particle_flow
    ! make sure to return domain_init to original state

    particle_init => domain_init%particle
    call domain_init%swap_particle_pointer(particle_flow)
    call domain_init%x2lamda()
    call domain_init%swap_particle_pointer(particle_init)

    ! add (negative) flips to adjust domain_fini
    ! (will switch it back)
    ! A vec is (xhi-xlo,0,0)
    ! B vec is (xy,yhi-ylo,0)
    ! C vec is (xz,yz,xhi-zlo)
    ! normal flip process is below (flip signs reversed below)
    ! A' = A
    ! B' = B + mA
    ! C' = C + pB + nA
    ! xy changes due to flipxy with xprd
    ! xz changes due to flipxz with xprd and flipyz with xy
    ! yz changes due to flipyz with yprd

    yz = yz_fini - flipyz*yprd_fini
    xz = xz_fini - flipxz*xprd_fini - flipyz*xy_fini
    xy = xy_fini - flipxy*xprd_fini

    ! now, juggle domain_fini's particle pointer
    ! and lamda2x using domain_fini's altered geometry
    ! then, set domain_fini back to normal

    particle_fini => domain_fini%particle
    call domain_fini%swap_particle_pointer(particle_flow)
    call domain_fini%set_global_box(boxlo_fini, boxhi_fini, &
         xy=xy, xz=xz, yz=yz)
    call domain_fini%lamda2x()
    call domain_fini%set_global_box(boxlo_fini, boxhi_fini, &
         xy=xy_fini, xz=xz_fini, yz=yz_fini)
    call domain_fini%swap_particle_pointer(particle_fini)

  end subroutine propagate_in_flow

  ! ----------------------------------------------------------------------
  ! Get flow-advanced positions in current box
  ! Requires a domain type (and associated particle type) from prev timestep
  ! And the change in img flags just due to box flips between the previous
  ! snapshot and this one
  ! ----------------------------------------------------------------------

  subroutine get_flow_unmapped(this, total_flips, domain_init, particle_flow)

    class(domain_t), target, intent(inout) :: this
    integer, dimension(3), intent(in) :: total_flips
    type(domain_t), target, intent(inout) :: domain_init
    type(particle_t), target, intent(inout) :: particle_flow
    integer :: nparticles,m,n,p
    type(particle_t), pointer :: particle_init => null()
    type(particle_t), pointer :: particle_fini => null()
    real, dimension(3) :: boxlo_init,boxhi_init
    real :: xy_init,xz_init,yz_init,xprd,yprd,xy,xz,yz
    logical :: xperiodic,yperiodic

    nparticles = this%particle%nparticles
    if (nparticles /= particle_flow%nparticles) then
       print *, 'unequal particle numbers encountered in get_flow_unmapped'
       stop
    end if

    ! not handling nonperiodic for now

    xperiodic = domain_init%xperiodic
    yperiodic = domain_init%yperiodic
    if ((.not. xperiodic) .or. (.not. yperiodic)) then
       print *, 'get_flow_unmapped needs periodic for now'
       stop
    end if

    ! copy previous UNMAPPED coords to particle_flow's coords

    particle_flow%x(:,:) = domain_init%particle%x(:,:)

    ! temporarily associate domain_init with particle_flow
    ! does this point to pointer or data?

    particle_init => domain_init%particle
    call domain_init%swap_particle_pointer(particle_flow)

    ! use remap to get the correct INITIAL images in particle_flow

    call domain_init%remap()

    ! now, particle_flow is 'in' domain_init
    ! as in fix_deform
    ! first set an updated box (here could be very deformed)
    ! first use image_flip to prepare images for flipping
    ! then remap the particles back into the box (uses x2lamda)

    m = total_flips(1)
    n = total_flips(2)
    p = total_flips(3)
    call image_flip(particle_flow, m, n, p)

    ! store current values first

    boxlo_init(:) = domain_init%boxlo(:)
    boxhi_init(:) = domain_init%boxhi(:)
    xy_init = domain_init%xy
    xz_init = domain_init%xz
    yz_init = domain_init%yz

    ! shift box just due to flips between initial and final states
    ! suspect???

    xprd = domain_init%xprd
    yprd = domain_init%yprd
    yz = domain_init%yz
    xz = domain_init%xz
    xy = domain_init%xy
    yz = yz + p*yprd
    xz = xz + p*yz
    xz = xz + n*xprd
    xy = xy + m*xprd

    call domain_init%set_global_box(domain_init%boxlo, domain_init%boxhi, &
         xy=xy, xz=xz, yz=yz)

    ! remap to advance particle images

    call domain_init%remap()

    ! need to shift into new box
    ! change to lamda, and then place into new box

    call domain_init%x2lamda()

    ! return domain_init to initial state

    domain_init%particle => particle_init
    domain_init%boxlo(:) = boxlo_init(:)
    domain_init%boxhi(:) = boxhi_init(:)
    domain_init%xy = xy_init
    domain_init%xz = xz_init
    domain_init%yz = yz_init
    call domain_init%set_global_box(boxlo_init, boxhi_init, &
         xy=xy_init, xz=xz_init, yz=yz_init)

    ! particle images are as originally in previous timestep
    ! apply delta image to flip all images

    ! we now have lamda coords with correct images
    ! temporarily associate current domain with these lamda coords
    ! then unmap them
    ! this zeros the images by default

    particle_fini => this%particle
    call this%swap_particle_pointer(particle_flow)

    call this%lamda2x()
    call this%unmap()

    ! reset the pointer for this domain

    this%particle => particle_fini

  end subroutine get_flow_unmapped

  ! ----------------------------------------------------------------------
  ! Update images due to box flips
  ! ----------------------------------------------------------------------

  subroutine image_flip(particle, m, n, p)

    type(particle_t), target, intent(inout) :: particle
    integer, intent(in) :: m,n,p
    integer, dimension(:), pointer :: image => null()
    integer :: nparticles,i,xbox,ybox,zbox

    image => particle%image
    nparticles = particle%nparticles

    do i = 1, nparticles
       xbox = iand(image(i), IMGMASK) - IMGMAX
       ybox = iand(ishft(image(i), -IMGBITS), IMGMASK) - IMGMAX
       zbox = ishft(image(i), -IMG2BITS) - IMGMAX

       ybox = ybox - p*zbox
       xbox = xbox - m*ybox - n*zbox

       image(i) = iand(xbox + IMGMAX, IMGMASK)
       image(i) = ior(image(i), &
            ishft(iand(ybox + IMGMAX, IMGMASK), IMGBITS))
       image(i) = ior(image(i), &
            ishft(iand(zbox + IMGMAX, IMGMASK), IMG2BITS))
    end do

  end subroutine image_flip

  ! ----------------------------------------------------------------------
  ! Convert box bounds to h array (used between mod data and here)
  ! ----------------------------------------------------------------------

  subroutine bounds2h(bounds, h)

    real, dimension(9), intent(in) :: bounds
    real, dimension(6), intent(out) :: h
    real :: xlo_bound,xhi_bound,ylo_bound,yhi_bound,zlo_bound,zhi_bound
    real :: xy,xz,yz,xlo,xhi,ylo,yhi,zlo,zhi

    xlo_bound = bounds(1)
    xhi_bound = bounds(2)
    xy = bounds(3)
    ylo_bound = bounds(4)
    yhi_bound = bounds(5)
    xz = bounds(6)
    zlo_bound = bounds(7)
    zhi_bound = bounds(8)
    yz = bounds(9)

    xlo = min(0.0,xy)
    xlo = min(xlo,xz)
    xlo = min(xlo,xy+xz)
    xlo = xlo_bound - xlo
    xhi = max(0.0,xy)
    xhi = max(xlo,xz)
    xhi = max(xlo,xy+xz)
    xhi = xhi_bound - xhi
    ylo = min(0.0,yz)
    ylo = ylo_bound - ylo
    yhi = max(0.0,yz)
    yhi = yhi_bound - yhi
    zlo = zlo_bound
    zhi = zhi_bound

    h(1) = xhi - xlo
    h(2) = yhi - ylo
    h(3) = zhi - zlo
    h(4) = yz
    h(5) = xz
    h(6) = xy

  end subroutine bounds2h

!!$  ! ----------------------------------------------------------------------
!!$  ! Get inverse of h matrix
!!$  ! ----------------------------------------------------------------------
!!$
!!$  subroutine h2h_inv(h, h_inv)
!!$
!!$    real, dimension(6), intent(in) :: h
!!$    real, dimension(6), intent(out) :: h_inv
!!$
!!$    h_inv(1) = 1.0/h(1)
!!$    h_inv(2) = 1.0/h(2)
!!$    h_inv(3) = 1.0/h(3)
!!$    h_inv(4) = -h(4) / (h(2)*h(3))
!!$    h_inv(5) = (h(4)*h(6) - h(2)*h(5)) &
!!$               / (h(1)*h(2)*h(3))
!!$    h_inv(6) = -h(6) / (h(1)*h(2))
!!$
!!$  end subroutine h2h_inv

  ! ----------------------------------------------------------------------
  ! Switch particle pointer
  ! ----------------------------------------------------------------------

  subroutine swap_particle_pointer(this, particle)

    class(domain_t), intent(inout) :: this
    type(particle_t), target, intent(in) :: particle

    if (.not. associated(this%particle)) then
       print *, 'error, particle not associated'
       stop
    end if

    this%particle => particle

  end subroutine swap_particle_pointer

  ! ----------------------------------------------------------------------
  ! Get box volume
  ! ----------------------------------------------------------------------

  function volume(this)

    class(domain_t), intent(in) :: this
    real :: volume

    volume = this%xprd * this%yprd * this%zprd

  end function volume

  ! ----------------------------------------------------------------------
  ! AND operators between logical and integer (like in C++)
  ! ----------------------------------------------------------------------

  elemental logical function int_and_log(i,l)

    integer, intent(in) :: i
    logical, intent(in) :: l
    int_and_log = .false.
    if (l .and. (i /= 0)) then
       int_and_log = .true.
    end if

  end function int_and_log

  ! ----------------------------------------------------------------------

  elemental logical function log_and_int(l,i)

    logical, intent(in) :: l
    integer, intent(in) :: i
    log_and_int = .false.
    if (l .and. (i /= 0)) then
       log_and_int = .true.
    end if

  end function log_and_int

end module mod_domain
