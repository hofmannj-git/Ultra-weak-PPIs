program cluster_dist_tri

  use, intrinsic :: iso_fortran_env
  use mod_particle, only : particle_t
  use mod_domain, only : domain_t
  use mod_group, only : group_t
  use mod_cells, only : cells_t
  use mod_stat_func, only : eval_cluster_id
  use mod_data, only : source_data_t
  implicit none

  integer :: i,j,k,nparticles,nargs,snap,nsnaps
  type(particle_t) :: particle
  type(domain_t) :: domain
  type(group_t) :: group
  type(cells_t) :: cells
  type(source_data_t) :: source_data
  character(len=300) :: buf,directory,mode
  integer, dimension(:), allocatable :: cluster_id,cluster_count
  real :: rsep,minsize
  integer(int64) :: tstart,tdur,tinc,tstep
  real, dimension(:), allocatable :: cluster_dist
  logical :: mean_only

  ! read inputs

  nargs = command_argument_count()
  if (nargs .lt. 6) then
     write(*, '(A)') 'cluster_dist_tri [mode] [directory] [tstart] [tdur] [tinc] [rsep]'
     stop
  end if
  call get_command_argument(1, buf); mode = trim(adjustl(buf))
  call get_command_argument(2, buf); directory = trim(adjustl(buf))
  call get_command_argument(3, buf); read(buf, *) tstart
  call get_command_argument(4, buf); read(buf, *) tdur
  call get_command_argument(5, buf); read(buf, *) tinc
  call get_command_argument(6, buf); read(buf, *) rsep

  if (trim(adjustl(mode)) == 'mean') then
     mean_only = .true.
  else if (trim(adjustl(mode)) == 'full') then
     mean_only = .false.
  else
     print *, 'mode must be mean or full'
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
  nparticles = particle%nparticles
  allocate(cluster_id(nparticles))
  allocate(cluster_count(nparticles))
  allocate(cluster_dist(nparticles))

  ! init group and cells

  call group%init(domain)
  minsize = 2*maxval(particle%typeradius(:)) + rsep
  call cells%init(group, minsize, .false.)

  ! loop over snaps

  do i = 1, nsnaps

     tstep = tstart + (i-1)*tinc
     snap = source_data%step_index(tstep)

     ! get cluster ids

     cluster_id(:) = 0
     call source_data%read_snap(particle, domain, snap)
     call cells%reset_grid()
     call cells%assign_particles()
     call eval_cluster_id(cells, .true., rsep, cluster_id)

     ! get the number of particles in each cluster

     cluster_count(:) = 0
     do j = 1, nparticles
        k = cluster_id(j)
        cluster_count(k) = cluster_count(k) + 1
     end do

     ! get the number of clusters with size k

     cluster_dist(:) = 0.0
     do j = 1, nparticles
        k = cluster_count(j)
        cluster_dist(k) = cluster_dist(k) + 1
     end do

     ! change these into the fraction of particles in each cluster size

     do j = 1, nparticles
        cluster_dist(j) = j * cluster_dist(j)
     end do
     cluster_dist(:) = cluster_dist(:) / nparticles

     ! print output

     if (trim(adjustl(mode)) == 'full') then

        ! if 'full', will print columns of cluster size and fraction
        ! of particles in those cluster sizes
        ! there will be a header separating each time step

        write(*, '(a)') 'ITEM: step'
        write(*, '(i14)') tstep
        write(*, '(a)') 'ITEM: cluster_nparticles fraction_of_particles'
        do j = 1, nparticles
           if (cluster_dist(j) > 0.0) then
              write(*, '(i14,e14.7)') j, cluster_dist(j)
           end if
        end do

     else

        ! if mean, then write not supported for no

        print *, 'mean not supported yet'
        stop

     end if

  end do

end program cluster_dist_tri
