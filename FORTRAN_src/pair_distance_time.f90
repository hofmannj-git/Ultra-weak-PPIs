program pair_distance_time

  use, intrinsic :: iso_fortran_env
  use mod_particle, only : particle_t
  use mod_domain, only : domain_t
  use mod_group, only : group_t
  use mod_cells, only : cells_t
  use mod_stat_func, only : eval_pairdisttime
  use mod_data, only : source_data_t
  use mod_statistics, only : running_stat_t
  implicit none

  character(len=300) :: buf,directory,subtr_str,subtr_g0_str
  integer(int64) :: tstart,tdur,tinc,tstep,g0
  real :: rmax,rmin,minsize,rval,rgrid,rsep
  real(real64) :: gmin,gmax
  integer :: nbin,i,j,k,m,n,nsnaps,nargs,snap,nparticles,type_a,type_b,group_a,group_b,num_to_calc
  character(len=2) :: plane
  character(len=3) :: mode
  character(len=300) :: fname
  type(source_data_t) :: source_data
  type(particle_t) :: particle
  type(group_t) :: group
  type(domain_t) :: domain
  type(cells_t) :: cells
  real, dimension(:,:), allocatable :: pair_sep_time, pair_sep_id
  real, dimension(:), allocatable :: bin_bindunbtimes
  logical, dimension(:), allocatable :: group_ids
  logical :: subtr,subtr_g0

  ! read inputs

  nargs = command_argument_count()
  if (nargs < 6) then
     write(*, '(a)') 'pair_distance_time [directory] [tstart] &
                     &[tdur] [tinc] [num_to_calc] [rsep] {type_a} {type_b}'
     stop
  end if
  call get_command_argument(1, buf); directory = trim(adjustl(buf))
  call get_command_argument(2, buf); read(buf, *) tstart
  call get_command_argument(3, buf); read(buf, *) tdur
  call get_command_argument(4, buf); read(buf, *) tinc
  call get_command_argument(5, buf); read(buf, *) num_to_calc
  call get_command_argument(6, buf); read(buf, *) rsep
  print *, "read first set of input arguments"
  type_a = 0
  type_b = 0
  if (nargs > 6) then
     call get_command_argument(7, buf); read(buf, *) type_a
     call get_command_argument(8, buf); read(buf, *) type_b
     if ((type_a == 0) .neqv. (type_b == 0)) then
        write(*, '(a)') 'only both types a and b can be zero, not one'
        stop
     end if
  end if
  print *, "read in all inputs"

  call source_data%init(trim(adjustl(directory)), particle, domain)
  if (.not. source_data%step_exists(tstart)) then
     write(*, '(a)') 'bad tstart value'
     stop
  end if
  print *, "read in timestep for measurement"
  if (.not. source_data%regular) then
     write(*, '(a)') 'hdf5 timesteps appear irregular'
     write(*, '(a)') 'not supporting nonzero tdur for irregular'
     stop
  end if
  nsnaps = int(tdur/tinc) + 1
  nparticles = particle%nparticles

  ! init group

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

  allocate(pair_sep_time(nsnaps,num_to_calc))
  pair_sep_time(:,:)=0
  allocate(pair_sep_id(2,num_to_calc))
  pair_sep_id(:,:)=0
  ! TODO is this an appropriate minsize?
  minsize = 2*maxval(particle%typeradius(:)) + rsep
  call cells%init(group, minsize, .false.)

  ! loop over snaps, accumulating to average

  print *, group_a, group_b
  print *, 'looping over snaps'
  do i = 1, nsnaps

     tstep = tstart + (i-1)*tinc
     snap = source_data%step_index(tstep)

     call source_data%read_snap(particle, domain, snap)
     call cells%reset_grid()
     call cells%assign_particles()

     call eval_pairdisttime(cells, group, group_a, group_b, rsep, pair_sep_time, pair_sep_id, num_to_calc, i)

  end do

  ! output data to files
  open( unit = 103, file = 'pair_distance_time.out')
  do k = 1, nsnaps
     write( 103, * ) pair_sep_time(k,:)
  end do
  close( 103 )

end program pair_distance_time
