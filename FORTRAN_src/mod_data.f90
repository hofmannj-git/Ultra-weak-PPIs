!
! manage hdf5 and text input (dynamic reading, etc.)
!

! TODO find a way to read in faster: no extra buffer?

module mod_data

  use iso_fortran_env
  use hdf5
  use mod_particle, only : particle_t
  use mod_domain, only : domain_t, get_image_flip
  implicit none
  private

  type :: dumpfile_t
     character(len=235) :: path
     integer(int64)     :: tstart, tend, tinc
  end type dumpfile_t

  type, public :: source_data_t

     logical :: regular ! true if equally spaced data without gaps
     integer :: nsnaps,ndumps
     integer, dimension(:), allocatable :: dump_nsnaps
     integer(int64), dimension(:), allocatable :: timestep ! timestep by index
     character(len=:), allocatable :: directory ! directory
     type(dumpfile_t), dimension(:), allocatable :: dumps
     real, dimension(:,:), allocatable :: box_bounds ! time-dependent bounds
     integer, dimension(:,:), allocatable :: image_flips

   contains

     procedure, public :: init
     procedure, public :: set_box
     procedure, public :: read_snap
     procedure, public :: step_exists
     procedure, public :: step_index
     procedure, public :: free => free_source_data

     procedure :: get_timesteps
     procedure :: get_dump_snap

  end type source_data_t

contains

  ! ----------------------------------------------------------------------
  ! read the *.txt files in the given directory and verify source data
  ! also set up a particle_t
  ! can set additional particle_t at a later time with init_particle
  ! ----------------------------------------------------------------------

  subroutine init(this, directory, particle, domain)

    class(source_data_t), intent(inout)     :: this
    character(len=*), intent(in)            :: directory
    type(particle_t), intent(inout), target :: particle
    type(domain_t), intent(inout), target   :: domain
    integer                                 :: i,nparticles,hdferr,count
    integer(int64)                          :: tstep
    character(len=2000000)                  :: buf
    real, dimension(9)                      :: bounds_i,bounds_j

    this%directory = trim(directory)

    ! open hdf5

    call h5open_f(hdferr)

    ! read in the paths.txt data (locations of hdf5 dump files and steps)

    open(7, file=trim(this%directory)//'paths.txt')
    read(7, *) ! ITEM: path tstart tend tinc
    i = 0
    do
       read(7, *, end=10)
       i = i + 1
    end do
10  close(7)
    this%ndumps = i
    allocate(this%dumps(this%ndumps))
    allocate(this%dump_nsnaps(this%ndumps))

    open(7, file=directory//'paths.txt')
    read(7, *) ! ITEM: path tstart tend tinc
    do i = 1, this%ndumps
       read(7, *) buf, this%dumps(i)%tstart, this%dumps(i)%tend, &
            this%dumps(i)%tinc
       write(this%dumps(i)%path, *) &
            this%directory//trim(adjustl(buf))
    end do
    close(7)

    ! get master list of timesteps from info read from paths.txt
    ! verify if these are regular
    ! this produces the nsnaps variable needed for allocation

    call this%get_timesteps()

    ! read in the particle.txt data (# particles, types, and size by type)
    ! init a particle_t (using defaults) to store the data

    open(7, file=this%directory//'particle.txt')
    read(7, *) ! ITEM: NUMBER OF PARTICLES
    read(7, *) nparticles
    call particle%init(nparticles, .false., .true., .true., .true.)
    read(7, *) ! ITEM: NUMBER OF TYPES
    read(7, *) particle%ntypes

    allocate(particle%typeradius(particle%ntypes))
    read(7, *) ! ITEM: SIZES IN TYPE ORDER
    do i = 1, particle%ntypes
       read(7, *) particle%typeradius(i)
    end do
    particle%typeradius(:) = particle%typeradius(:) / 2.0

    if(allocated(particle%type)) then
       deallocate(particle%type)
    end if
    allocate(particle%type(particle%nparticles))

    read(7, *) ! ITEM: LIST OF TYPES IN LAMMPS ATOM ORDER
    count = 0
    do i = 1, particle%nparticles
       read(7, '(I10)', end=20) particle%type(i)
       if (particle%type(i) == 1) then
          count = count + 1
       else if (particle%type(i) == 2) then
          particle%molec(i) = count
       end if
    end do
20  close(7)

    ! verify that box.txt exists and has data for all timesteps
    ! read all of this in
    ! if encountering a repeat timestep, tell the user

    allocate(this%box_bounds(9,this%nsnaps))
    open(7, file=this%directory//'box.txt')
    read(7, *) ! ITEM: TIMESTEP XLO_BOUND XHI_BOUND XY ..
    i = 1
    do
       read(7, *, end=30) tstep, this%box_bounds(:,i)

       ! first check is the timestep just isn't represented in paths.txt

       if (.not. this%step_exists(tstep)) then
          print *, 'inconsistent timestep in box.txt: ', tstep
          stop
       end if

       ! if the timestep read is beyond the timestep that we'd need,
       ! then produce an error
       ! repeat timesteps are okay, but most recent taken?

       if (tstep > this%timestep(i)) then
          print *, 'Error: no step ', this%timestep(i), ' in box.txt'
          print *, 'or the timesteps are out of order'
          stop
       end if

       i = i + 1
    end do
30  close(7)

    ! get all image flips between snapshots
    ! first flip is 0,0,0

    allocate(this%image_flips(3,this%nsnaps))
    this%image_flips(:,:) = 0
    do i = 1, this%nsnaps-1
       bounds_i(:) = this%box_bounds(:,i)
       bounds_j(:) = this%box_bounds(:,i+1)
       call get_image_flip(bounds_i, bounds_j, this%image_flips(:,i+1))
    end do

    ! set up the box and init to first snapshot data

    call domain%init(particle)
    call this%set_box(domain, 1)

  end subroutine init

  ! ----------------------------------------------------------------------
  ! free resources
  ! ----------------------------------------------------------------------

  subroutine free_source_data(this)

    class(source_data_t), intent(inout) :: this
    integer                             :: hdferr

    ! not really necessary for allocatable

    deallocate(this%dumps)
    deallocate(this%dump_nsnaps)
    deallocate(this%timestep)
    deallocate(this%box_bounds)
    deallocate(this%directory)

    call h5close_f(hdferr)

  end subroutine free_source_data

  ! ----------------------------------------------------------------------
  ! count nsnaps & build list of timesteps from what was listed in paths.txt
  ! the tstarts must be monotonically increasing,
  !   and modulo(tstart-tend,tinc) must be zero for each dump
  ! for regularity, all tinc must be the same (equal spacing within dumps),
  !   all modulo(tstart-torigin,tinc) must be zero (reasonable spacing between),
  !   all tstart-last_tend <= tinc (no gaps)
  ! ----------------------------------------------------------------------

  subroutine get_timesteps(this)

    class(source_data_t), intent(inout) :: this
    integer                             :: i,j,k,ndumps
    integer(int64)                      :: tstart,tend,tinc
    integer(int64)                      :: last_tstart,last_tend,origin
    integer(int64)                      :: next_tstart,last_in_this
    logical                             :: regular

    regular = .true.
    ndumps = size(this%dumps)
    origin = this%dumps(1)%tstart

    ! first, check if all of the dump files have legitimate tstart, tend, tinc
    ! each tstart must be > last tstart

    do i = 1, ndumps
       tstart = this%dumps(i)%tstart
       tend = this%dumps(i)%tend
       tinc = this%dumps(i)%tinc
       if (i > 1) then
          if (tstart <= last_tstart) then
             print *, 'paths.txt tstarts not monotonically increasing ', tstart
             stop
          end if
       end if
       if (tend < tstart) then
          print *, 'tend is less than tstart in paths.txt ', tstart, tend
          stop
       end if
       if (tinc <= 0) then
          print *, 'tinc must be > 0: ', tinc
          stop
       end if
       if (modulo(tend-tstart,tinc) .ne. 0) then
          print *, 'bad tstart-tend not multiple of tinc ', tstart, tend, tinc
          stop
       end if
       last_tstart = tstart
       last_tend = tend
    end do

    ! now, check if all the tinc are the same number
    ! if not, then cannot be regular

    tinc = this%dumps(1)%tinc
    do i = 2, ndumps
       if (this%dumps(i)%tinc .ne. tinc) then
          regular = .false.
          exit
       end if
    end do

    ! check if tstarts are multiples of tinc away from the origin

    if (regular) then
       do i = 2, ndumps
          tstart = this%dumps(i)%tstart
          if (modulo(tstart-origin,tinc) .ne. 0) then
             regular = .false.
             exit
          end if
       end do
    end if

    ! finally, if none of the tests above failed, check to see
    ! if there are any gaps between the snapshots

    if (regular) then
       last_tend = this%dumps(1)%tend
       do i = 2, ndumps
          tstart = this%dumps(i)%tstart
          if (tstart-last_tend > tinc) then
             print *, 'gaps in paths.txt: ', last_tend, tstart
             stop
          end if
          last_tend = this%dumps(i)%tend
       end do
    end if

    ! at this point, it is known whether the dumps are regular or not
    ! count the number of snaps

    this%nsnaps = 0
    do i = 1, ndumps-1
       tstart = this%dumps(i)%tstart
       tend = this%dumps(i)%tend
       tinc = this%dumps(i)%tinc
       next_tstart = this%dumps(i+1)%tstart

       ! get the maximum timestep in this dump
       ! that doesn't go beyond first of next dump

       last_in_this = tend
       do while (last_in_this > next_tstart)
          last_in_this = last_in_this - tinc
       end do

       ! divide to find the number of snaps
       ! if we don't already have the last_in_this covered in next dump,
       ! add one

       this%dump_nsnaps(i) = int((last_in_this-tstart)/tinc)
       if (last_in_this /= next_tstart) then
          this%dump_nsnaps(i) = this%dump_nsnaps(i) + 1
       end if
    end do
    this%dump_nsnaps(ndumps) = &
         int((this%dumps(ndumps)%tend-this%dumps(ndumps)%tstart) &
         /this%dumps(ndumps)%tinc) &
         + 1
    this%nsnaps = sum(this%dump_nsnaps)
    allocate(this%timestep(this%nsnaps))

    ! fill the timestep array

    k = 1
    do i = 1, ndumps
       tstart = this%dumps(i)%tstart
       tinc = this%dumps(i)%tinc
       do j = 1, this%dump_nsnaps(i)
          this%timestep(k) = tstart + (j-1)*tinc
          k = k + 1
       end do
    end do

    this%regular = regular

  end subroutine get_timesteps

  ! ----------------------------------------------------------------------
  ! Check if a step exists (no assumption on regularity)
  ! ----------------------------------------------------------------------

  elemental logical function step_exists(this, step)

    class(source_data_t), intent(in) :: this
    integer(int64), intent(in)       :: step
    integer                          :: i

    step_exists = .false.
    do i = 1, this%nsnaps
       if (this%timestep(i) .eq. step) then
          step_exists = .true.
          return
       end if
    end do

  end function step_exists

  ! ----------------------------------------------------------------------
  ! Get step index
  ! Stop program if does not exist
  ! ----------------------------------------------------------------------

  integer function step_index(this, step)

    class(source_data_t), intent(in) :: this
    integer(int64), intent(in)       :: step
    integer                          :: i

    do i = 1, this%nsnaps
       if (this%timestep(i) .eq. step) then
          step_index = i
          return
       end if
    end do

    print *, 'bad ttarg'
    stop

  end function step_index

  ! ----------------------------------------------------------------------
  ! Pass a snapshot's box dimensions to a domain type, to make it current
  ! IMPORTANT: LAMMPS dump file uses xlo_bound,xhi_bound variables,
  !            not xlo,xhi ... have to convert here
  ! ----------------------------------------------------------------------

  subroutine set_box(this, domain, i)

    class(source_data_t), intent(in), target :: this
    type(domain_t), intent(inout)            :: domain
    integer, intent(in)                      :: i
    real, dimension(3)                       :: boxlo,boxhi
    real                                     :: xy,xz,yz
    real                                     :: xlo,xhi,ylo,yhi,zlo,zhi
    real                                     :: xlo_bound,xhi_bound
    real                                     :: ylo_bound,yhi_bound
    real                                     :: zlo_bound,zhi_bound
    real, dimension(:), pointer              :: box => null()

    box => this%box_bounds(:,i)

    xlo_bound = box(1)
    xhi_bound = box(2)
    xy = box(3)
    ylo_bound = box(4)
    yhi_bound = box(5)
    xz = box(6)
    zlo_bound = box(7)
    zhi_bound = box(8)
    yz = box(9)

    xlo = min(0.0,xy)
    xlo = min(xlo,xz)
    xlo = min(xlo,xy+xz)
    xlo = xlo_bound - xlo
    xhi = max(0.0,xy)
    xhi = max(xhi,xz)
    xhi = max(xhi,xy+xz)
    xhi = xhi_bound - xhi
    ylo = min(0.0,yz)
    ylo = ylo_bound - ylo
    yhi = max(0.0,yz)
    yhi = yhi_bound - yhi
    zlo = zlo_bound
    zhi = zhi_bound

    boxlo = (/ xlo, ylo, zlo /)
    boxhi = (/ xhi, yhi, zhi /)

    call domain%set_global_box(boxlo,boxhi,xy,xz,yz)

  end subroutine set_box

  ! ----------------------------------------------------------------------
  ! Get dump and snap ids for given snapshot
  ! ----------------------------------------------------------------------

  pure subroutine get_dump_snap(this, n, dump_id, snap_id)

    class(source_data_t), intent(in) :: this
    integer, intent(in)              :: n
    integer, intent(out)             :: dump_id, snap_id
    integer                          :: i,j

    j = n
    do i = 1, this%ndumps
       j = j - this%dump_nsnaps(i)
       if (j <= 0) then
          dump_id = i
          snap_id = this%dump_nsnaps(i) + j
          return
       end if
    end do

  end subroutine get_dump_snap

  ! ----------------------------------------------------------------------
  ! Read a snapshot from the source data into the particle type
  ! ----------------------------------------------------------------------

  subroutine read_snap(this, particle, domain, n)

    class(source_data_t), intent(in) :: this
    type(particle_t), intent(inout) :: particle
    type(domain_t), intent(inout) :: domain
    integer, intent(in) :: n

    real, dimension(:,:,:), allocatable :: buf_pos, buf_vel, buf_pe, buf_sts
    integer :: i, dump, nparticles, snap
    integer(HSIZE_T), dimension(3) :: start
    integer :: hdferr
    integer(HID_T) :: d_pos, d_vel, d_pe, d_sts
    integer(HID_T) :: s_pos, s_vel, s_pe, s_sts
    integer(HID_T) :: s_m_pos, s_m_vel, s_m_pe, s_m_sts
    integer(HID_T) :: file_id
    integer(HSIZE_T), dimension(3) :: ct_pos, ct_vel, ct_pe, ct_sts

    nparticles = particle%nparticles

    call this%get_dump_snap(n, dump, snap)

    ! can I get rid of the trivial third dimension?

    allocate(buf_pos(3, nparticles, 1))
    allocate(buf_vel(3, nparticles, 1))
    allocate(buf_pe (1, nparticles, 1))
    allocate(buf_sts(6, nparticles, 1))

    ! prepare the file for reading

    call h5fopen_f(trim(adjustl(this%dumps(dump)%path)), H5F_ACC_RDONLY_F, &
         file_id, hdferr)

    ! get dataset handles

    call h5dopen_f(file_id, 'pos'         , d_pos, hdferr)
    call h5dopen_f(file_id, 'vel'         , d_vel, hdferr)
    call h5dopen_f(file_id, 'pe'          , d_pe , hdferr)
    call h5dopen_f(file_id, 'virialstress', d_sts, hdferr)

    ! get file dataspace handles

    call h5dget_space_f(d_pos, s_pos, hdferr)
    call h5dget_space_f(d_vel, s_vel, hdferr)
    call h5dget_space_f(d_pe , s_pe , hdferr)
    call h5dget_space_f(d_sts, s_sts, hdferr)

    ! get starting position for selection

    start = (/ 0, 0, snap-1 /)

    ! set count for hyperslab selection

    ct_pos = (/ 3, nparticles, 1 /)
    ct_vel = (/ 3, nparticles, 1 /)
    ct_pe  = (/ 1, nparticles, 1 /)
    ct_sts = (/ 6, nparticles, 1 /)

    ! get dataspace handle for memory

    call h5screate_simple_f(3, ct_pos, s_m_pos, hdferr)
    call h5screate_simple_f(3, ct_vel, s_m_vel, hdferr)
    call h5screate_simple_f(3, ct_pe , s_m_pe , hdferr)
    call h5screate_simple_f(3, ct_sts, s_m_sts, hdferr)

    ! select by hyperslabs

    call h5sselect_hyperslab_f(s_pos, H5S_SELECT_SET_F, &
         start, ct_pos, hdferr)
    call h5sselect_hyperslab_f(s_vel, H5S_SELECT_SET_F, &
         start, ct_vel, hdferr)
    call h5sselect_hyperslab_f(s_pe , H5S_SELECT_SET_F, &
         start, ct_pe , hdferr)
    call h5sselect_hyperslab_f(s_sts, H5S_SELECT_SET_F, &
         start, ct_sts, hdferr)

    ! read file data to buffers (most expensive part of subroutine by far)

    call h5dread_f(d_pos, H5T_NATIVE_REAL, buf_pos, ct_pos, hdferr, &
         mem_space_id=s_m_pos, file_space_id=s_pos)
    call h5dread_f(d_vel, H5T_NATIVE_REAL, buf_vel, ct_vel, hdferr, &
         mem_space_id=s_m_vel, file_space_id=s_vel)
    call h5dread_f(d_pe , H5T_NATIVE_REAL, buf_pe , ct_pe , hdferr, &
         mem_space_id=s_m_pe , file_space_id=s_pe )
    call h5dread_f(d_sts, H5T_NATIVE_REAL, buf_sts, ct_sts, hdferr, &
         mem_space_id=s_m_sts, file_space_id=s_sts)

    ! read buffers to particle structure

    do i = 1, particle%nparticles
       particle%x(:,i) =  buf_pos(:, i, 1)
       particle%v(:,i) = buf_vel(:, i, 1)
       particle%wbody(:,i) = buf_vel(:, i, 1)
       particle%quat(:,i) = buf_sts(1:4, i, 1)
       particle%euler(1,i) = atan2(2*(buf_sts(1,i,1)*buf_sts(2,i,1) + &
                                   buf_sts(3,i,1)*buf_sts(4,i,1)), &
                                   1 - 2*(buf_sts(2,i,1)*buf_sts(2,i,1) + &
                                   buf_sts(3,i,1)*buf_sts(3,i,1)))
       particle%euler(2,i) = asin(2*(buf_sts(1,i,1)*buf_sts(3,i,1) - &
                                  buf_sts(4,i,1)*buf_sts(2,i,1)))
       particle%euler(3,i) = atan2(2*(buf_sts(1,i,1)*buf_sts(4,i,1) - &
                                   buf_sts(2,i,1)*buf_sts(3,i,1)), &
                                   1 - 2*(buf_sts(3,i,1)*buf_sts(3,i,1) + &
                                   buf_sts(4,i,1)*buf_sts(4,i,1)))
       if (particle%pe_flag) then
          particle%pe(i) = buf_pe (1, i, 1)
       end if
       if (particle%stress_flag) then
          particle%stress(:,i) = buf_sts(:, i, 1)
       end if
    end do

    ! free resources

    call h5sclose_f(s_pos, hdferr)
    call h5sclose_f(s_vel, hdferr)
    call h5sclose_f(s_pe , hdferr)
    call h5sclose_f(s_sts, hdferr)
    call h5sclose_f(s_m_pos, hdferr)
    call h5sclose_f(s_m_vel, hdferr)
    call h5sclose_f(s_m_pe , hdferr)
    call h5sclose_f(s_m_sts, hdferr)
    call h5dclose_f(d_pos, hdferr)
    call h5dclose_f(d_vel, hdferr)
    call h5dclose_f(d_pe , hdferr)
    call h5dclose_f(d_sts, hdferr)
    call h5fclose_f(file_id, hdferr)

    ! reset the box using stored box dimensions

    call this%set_box(domain, n)

  end subroutine read_snap

end module mod_data
