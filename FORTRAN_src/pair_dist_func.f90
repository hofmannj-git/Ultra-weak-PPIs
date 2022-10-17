program pair_dist_func

  use, intrinsic :: iso_fortran_env
  use mod_particle, only : particle_t
  use mod_domain, only : domain_t
  use mod_group, only : group_t
  use mod_cells, only : cells_t
  use mod_stat_func, only : eval_rdf, eval_pair_dist_func
  use mod_data, only : source_data_t
  use mod_statistics, only : running_stat_t
  implicit none

  character(len=300) :: buf,directory,subtr_str,subtr_g0_str
  integer(int64) :: tstart,tdur,tinc,tstep,g0
  real :: rmax,rmin,minsize,rval,rgrid
  real(real64) :: gmin,gmax
  integer :: nbin,i,j,k,m,n,nsnaps,nargs,snap,nparticles,type_a,type_b,group_a,group_b
  character(len=2) :: plane
  character(len=3) :: mode
  character(len=300) :: fname
  type(source_data_t) :: source_data
  type(particle_t) :: particle
  type(group_t) :: group
  type(domain_t) :: domain
  type(cells_t) :: cells
  type(running_stat_t), dimension(:), allocatable :: rs_array_rad,rs0_array_rad
  type(running_stat_t), dimension(:,:,:), allocatable :: rs_array_xyz,rs0_array_xyz
  real, dimension(:), allocatable :: g_rad,g0_rad
  real, dimension(:,:,:), allocatable :: g_xyz,g0_xyz
  real(real64), dimension(:,:), allocatable :: g_xyz_flat
  logical, dimension(:), allocatable :: group_ids
  logical :: subtr,subtr_g0

  ! read inputs

  nargs = command_argument_count()
  if (nargs < 9) then
     write(*, '(a)') 'pair_dist_func [mode] [directory] [tstart] &
                     &[tdur] [tinc] [rmax] [rmin] [nbin] [subtr] {gmin} {gmax} &
                     &{plane} {subtr_g0} {g0} {type_a} {type_b}'
     stop
  end if
  call get_command_argument(1, buf); mode = trim(adjustl(buf))
  call get_command_argument(2, buf); directory = trim(adjustl(buf))
  call get_command_argument(3, buf); read(buf, *) tstart
  call get_command_argument(4, buf); read(buf, *) tdur
  call get_command_argument(5, buf); read(buf, *) tinc
  call get_command_argument(6, buf); read(buf, *) rmax
  call get_command_argument(7, buf); read(buf, *) rmin
  call get_command_argument(8, buf); read(buf, *) nbin
  call get_command_argument(9, buf); subtr_str = trim(adjustl(buf))
  if (subtr_str == 'yes') then
     subtr = .true.
  else if (subtr_str == 'no') then
     subtr = .false.
  else
     write(*, '(a)') 'subtr must be yes or no'
     stop
  end if
  if ((mode /= 'xyz') .and. (mode /= 'rad')) then
     write(*, '(a)') 'mode must be xyz or rad'
     stop
  end if   
  if (rmax <= 0.0) then
     write(*, '(a)') 'rmax must be > 0.0'
     stop
  end if
  ! TODO does 1 actually work?
  if (nbin < 1) then
     write(*, '(a)') 'nbin must be >= 1'
     stop
  end if
  if ((mode == 'xyz') .and. (mod(nbin,2) == 0)) then
     write(*, '(a)') 'must have odd nbin for mode xyz'
     stop
  end if
  print *, "read first set of input arguments"
  if (nargs > 9) then
     call get_command_argument(10, buf); read(buf, *) gmin
     call get_command_argument(11, buf); read(buf, *) gmax
     call get_command_argument(12, buf); read(buf, *) plane
     ! Expected input for subtr_g0 is timestep of 
     !  snap present in same folder
     call get_command_argument(13, buf); subtr_g0_str = trim(adjustl(buf))
     print *, "read in g0 information"
     if (gmax <= 0.0) then
        write(*, '(a)') 'requiring gmax be > 0.0'
        stop
     end if
     if ((plane /= 'xy') .and. (plane /= 'xz') .and. (plane /= 'yz')) then
        write(*, '(a)') 'plane must be xy, xz, or yz'
        stop
     end if
     if (subtr_g0_str == 'yes') then
        subtr_g0 = .true.
        call get_command_argument(14, buf); read(buf, *) g0
     else if (subtr_g0_str == 'no') then
        subtr_g0 = .false.
     else
        write(*, '(a)') 'subtr_g0 must be yes or no'
        stop
     end if
  else
     if (mode == 'xyz') then
        write(*, '(a)') 'must supply gmax and plane arguments for mode xyz'
        stop
     end if
  end if
  type_a = 0
  type_b = 0
  print *, "sucessfully read in all inputs up to type inputs"
  print *, subtr_g0
  if (nargs > 13) then
     if ((subtr_g0)) then
        print *, "read subtr_g0 variable as true"
        call get_command_argument(15, buf); read(buf, *) type_a
        call get_command_argument(16, buf); read(buf, *) type_b
     elseif ((.not. subtr_g0)) then
        call get_command_argument(14, buf); read(buf, *) type_a
        call get_command_argument(15, buf); read(buf, *) type_b
     else
        write(*, '(a)') 'must input correct number of arguments at this time'
     end if
     if ((type_a == 0) .neqv. (type_b == 0)) then
        write(*, '(a)') 'only both types a and b can be zero, not one'
        stop
     end if
  end if
  print *, "read in all inputs"

  ! get bin spacing
  ! TODO does this get things right with center square?
  
  if (mode == 'xyz') then
     rgrid = 2 * rmax / nbin
  else
     rgrid = (rmax - rmin) / nbin
  end if

  call source_data%init(trim(adjustl(directory)), particle, domain)
  if (.not. source_data%step_exists(tstart)) then
     write(*, '(a)') 'bad tstart value'
     stop
  end if
  print *, "read in timestep for measurement, now trying to read in g0"
  if (subtr_g0) then ! check to see if subtr_g0 is used
     print *, "trying to see if g0 data exists"
     if (.not. source_data%step_exists(g0)) then
        write(*, '(a)') 'g0 value'
        stop
     end if
  end if
  if (.not. source_data%regular) then
     write(*, '(a)') 'hdf5 timesteps appear irregular'
     write(*, '(a)') 'not supporting nonzero tdur for irregular'
     stop
  end if
  nsnaps = int(tdur/tinc) + 1
  nparticles = particle%nparticles

  ! init group
  ! using all-all for group pair in pair distribution function

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

  ! TODO is this an appropriate minsize?

  minsize = sqrt(3.0) * rmax
  call cells%init(group, minsize, .false.)

  ! allocate arrays for statistical functions

  if (mode == 'xyz') then
     allocate(g_xyz(nbin,nbin,nbin))
     allocate(rs_array_xyz(nbin,nbin,nbin))
     allocate(g_xyz_flat(nbin,nbin))
     if (subtr_g0) then
        allocate(g0_xyz(nbin,nbin,nbin))
        allocate(rs0_array_xyz(nbin,nbin,nbin))
     end if
  else
     allocate(g_rad(nbin))
     allocate(rs_array_rad(nbin))
     if (subtr_g0) then
        allocate(g0_rad(nbin))
        allocate(rs0_array_rad(nbin))
     end if
  end if

  ! loop over snaps, accumulating to average

  print *, group_a, group_b
  print *, 'looping over snaps'
  do i = 1, nsnaps

     tstep = tstart + (i-1)*tinc
     snap = source_data%step_index(tstep)

     call source_data%read_snap(particle, domain, snap)
     call cells%reset_grid()
     call cells%assign_particles()

     if (mode == 'xyz') then
        print *, 'eval pdf snap ', i
        call eval_pair_dist_func(cells, group, group_a, group_b, rmax, nbin, g_xyz)
        do j = 1, nbin
           do k = 1, nbin
              do m = 1, nbin
                 call rs_array_xyz(j,k,m)%push(real(g_xyz(j,k,m), real64))
              end do
           end do
        end do
     else
        call eval_rdf(cells, group, group_a, group_b, rmax, rmin, nbin, g_rad)
        do j = 1, nbin
           call rs_array_rad(j)%push(real(g_rad(j), real64))
        end do
     end if

  end do

  ! perform calculations on timestep selected for g0 if called
 
  if (subtr_g0) then
     print *, 'calculating g0'
     
     snap = source_data%step_index(g0)

     call source_data%read_snap(particle, domain, snap)
     call cells%reset_grid()
     call cells%assign_particles()

     if (mode == 'xyz') then
        print *, 'eval pdf snap ', i
        call eval_pair_dist_func(cells, group, group_a, group_b, rmax, nbin, g0_xyz)
        do j = 1, nbin
           do k = 1, nbin
              do m = 1, nbin
                 call rs0_array_xyz(j,k,m)%push(real(g0_xyz(j,k,m), real64))
              end do
           end do
        end do
     else
        call eval_rdf(cells, group, group_a, group_b, rmax, rmin, nbin, g0_rad)
        do j = 1, nbin
           call rs0_array_rad(j)%push(real(g_rad(j), real64))
        end do
     end if
  end if

  ! output data
  ! make sure to use averaged results

  ! output data if subtr_g0 is present properly
 
  if (subtr_g0) then
     if (mode == 'xyz') then
        do k = 1, nbin
           do j = 1, nbin
              do i = 1, nbin
                 g_xyz(i,j,k) = real(rs_array_xyz(i,j,k)%mean())-real(rs0_array_xyz(i,j,k)%mean())
              end do
           end do
        end do
        n = size(g_xyz, 1)
        m = (n + 1) / 2
        if (plane == 'xy') then
           do i = 1, n
              do j = 1, n
                 g_xyz_flat(i,j) = real(g_xyz(i,j,m))
              end do
           end do
        else if (plane == 'xz') then
           do i = 1, n
              do j = 1, n
                 g_xyz_flat(i,j) = real(g_xyz(i,m,j))
              end do
           end do
        else
           do i = 1, n
              do j = 1, n
                 g_xyz_flat(i,j) = real(g_xyz(m,i,j))
              end do
           end do
        end if
        write(buf,'(i10)') tstep
        write(fname,'(5a)') 'pdf-', plane, '-', trim(adjustl(buf)), '.ppm'

     else
        write(*, '(a)') 'r, g(r)'
        do i = 1, nbin
           g_rad(i) = real(rs_array_rad(i)%mean())-real(rs0_array_rad(i)%mean())
           rval = rmin + (real(i)-0.5) * rgrid
           write(*, '(2g14.5)') rval, g_rad(i)
        end do
     end if
  else 
     if (mode == 'xyz') then
        do k = 1, nbin
           do j = 1, nbin
              do i = 1, nbin
                 g_xyz(i,j,k) = real(rs_array_xyz(i,j,k)%mean())
              end do
           end do
        end do
        n = size(g_xyz, 1)
        m = (n + 1) / 2
        if (plane == 'xy') then
           do i = 1, n
              do j = 1, n
                 g_xyz_flat(i,j) = real(g_xyz(i,j,m))
              end do
           end do
        else if (plane == 'xz') then
           do i = 1, n
              do j = 1, n
                 g_xyz_flat(i,j) = real(g_xyz(i,m,j))
              end do
           end do
        else
           do i = 1, n
              do j = 1, n
                 g_xyz_flat(i,j) = real(g_xyz(m,i,j))
              end do
           end do
        end if
        write(buf,'(i10)') tstep
        write(fname,'(5a)') 'pdf-', plane, '-', trim(adjustl(buf)), '.ppm'

     else
        write(*, '(a)') 'r, g(r)'
        do i = 1, nbin
           g_rad(i) = real(rs_array_rad(i)%mean())
           rval = rmin + (real(i)-0.5) * rgrid
           write(*, '(2g14.5)') rval, g_rad(i)
        end do
     end if
  end if
  

end program pair_dist_func
