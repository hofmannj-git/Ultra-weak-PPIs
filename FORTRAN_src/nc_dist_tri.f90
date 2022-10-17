program nc_dist_tri

  use, intrinsic :: iso_fortran_env
  use mod_particle, only : particle_t
  use mod_domain, only : domain_t
  use mod_group, only : group_t
  use mod_cells, only : cells_t
  use mod_stat_func, only : eval_nc, eval_nc_parttype
  use mod_data, only : source_data_t
  implicit none

  integer :: i,j,k,nparticles,nargs,snap,nsnaps,type_a,type_b,group_a,group_b,na,nb,npartcalc
  type(particle_t) :: particle
  type(domain_t) :: domain
  type(group_t) :: group
  type(cells_t) :: cells
  type(source_data_t) :: source_data
  character(len=300) :: buf,directory,mode
  integer, dimension(:), allocatable :: nc
  real :: rsep,minsize,meanval
  integer(int64) :: tstart,tdur,tinc,tstep
  integer, parameter :: ncmax = 20
  real, dimension(ncmax+1) :: ncdist
  logical, dimension(:), allocatable :: group_ids
  logical :: mean_only

  ! read inputs

  nargs = command_argument_count()
  if (nargs .lt. 6) then
     write(*, '(A)') 'nc_dist_tri [mode] [directory] [tstart] [tdur] [tinc] [rsep] {type_a} {type_b}'
     stop
  end if
  call get_command_argument(1, buf); mode = trim(adjustl(buf))
  call get_command_argument(2, buf); directory = trim(adjustl(buf))
  call get_command_argument(3, buf); read(buf, *) tstart
  call get_command_argument(4, buf); read(buf, *) tdur
  call get_command_argument(5, buf); read(buf, *) tinc
  call get_command_argument(6, buf); read(buf, *) rsep
  type_a = 0
  type_b = 0
  if (nargs > 6) then
    call get_command_argument(7, buf); read(buf, *) type_a !particle i, at center of calculation
    call get_command_argument(8, buf); read(buf, *) type_b !particle j, those in contact with type i
     if ((type_a == 0) .neqv. (type_b == 0)) then
        write(*, '(a)') 'only both types a and b can be zero, not one'
        stop
     end if
  end if

  if (trim(adjustl(mode)) == 'mean') then
     mean_only = .true.
  else if (trim(adjustl(mode)) == 'full') then
     mean_only = .false.
  else if (trim(adjustl(mode)) == 'band') then
     mean_only = .false.
  else
     print *, 'mode must be mean or full or band (full and band)'
     stop
  end if
  call source_data%init(trim(adjustl(directory)), particle, domain)
  if (.not. source_data%step_exists(tstart)) then
     print *, 'bad tstart value'
     stop
  end if
  if (.not. source_data%regular) then
     print *, 'hdf5 timesteps appear irregular'
     print *, 'not supporting nonzero tdur for irregular'
     stop
  end if
  if (modulo(tdur,tinc) /= 0) then
     print *, 'tdur is not a multiple of tinc'
     stop
  end if
  nsnaps = int(tdur/tinc) + 1

  ! init group and cells
  nparticles = particle%nparticles
  call group%init(domain)
  if (type_a /= 0) then
     if (type_a > particle%ntypes) then
        write(*, '(a)') 'type a must be <= ntypes'
        stop
     end if
     if (type_b > particle%ntypes) then
        write(*, '(a)') 'type b must be <= ntypes'
        stop
     end if
     if (type_a < 1) then
        write(*, '(a)') 'type a must be >= 1'
        stop
     end if
     if (type_b < 1) then
        write(*, '(a)') 'type b must be >= 1'
        stop
     end if
     allocate(group_ids(nparticles))
     group_ids(:) = .false.
     do i = 1, nparticles
        if (particle%type(i) == type_a) group_ids(i) = .true.
     end do
     group_a = group%create(group_ids(:))
     if (type_b /= type_a) then
        group_ids(:) = .false.
        do i = 1, nparticles
           if (particle%type(i) == type_b) group_ids(i) = .true.
        end do
        group_b = group%create(group_ids(:))
     else
        group_b = group_a
     end if
  else
     group_a = 1
     group_b = 1
  end if
  npartcalc = nparticles
  na = 0
  nb = 0
  if (type_a /= 0) then
     na = group%count(group_a)
     nb = group%count(group_b)
     !npartcalc = na
  end if
  allocate(nc(npartcalc))

  minsize = 2*maxval(particle%typeradius(:)) + rsep
  call cells%init(group, minsize, .false.)

  ! loop over snaps

  do i = 1, nsnaps

     tstep = tstart + (i-1)*tinc
     snap = source_data%step_index(tstep)

     ! get contact numbers

     nc(:) = 0
     call source_data%read_snap(particle, domain, snap)
     call cells%reset_grid()
     call cells%assign_particles()
     if (type_a /= 0) then
        call eval_nc_parttype(cells, group, group_a, group_b, type_a, type_b, &
            (particle%typeradius(type_a)+particle%typeradius(type_b))*(1+rsep), nc)
     else
        call eval_nc(cells, .true., rsep*2*(sum(particle%typeradius(:))/ &
                    size(particle%typeradius(:))), nc)
     end if
     !print *, nc
     ! get the distribution
     ncdist(:) = 0.0
     do j = 1, npartcalc
        k = nc(j)
        if (type_a /= 0) then
           if ((k <= ncmax) .and. (particle%type(j) == type_a)) then
              ncdist(k+1) = ncdist(k+1) + 1
              !print *, 'calc'
           end if
        else
           if (k <= ncmax) then
              ncdist(k+1) = ncdist(k+1) + 1
           end if
        end if
     end do
     if (type_a /= 0) then
        ncdist = ncdist / na
     else
        ncdist = ncdist / npartcalc
     end if

     ! print output

     write(*, '(i10)', advance='no') tstep
     if (.not. mean_only) then
        do j = 0, ncmax
           write(*, '(f8.5)', advance='no') ncdist(j+1)
        end do
        write (*, '(a)') ''
     else
        meanval = 0.0
        do j = 0, ncmax
           meanval = meanval + j*ncdist(j+1)
        end do
        write (*, '(f8.5)') meanval
     end if

  end do

end program nc_dist_tri
