program cluster_dist_tri

  use, intrinsic :: iso_fortran_env
  use mod_particle, only : particle_t
  use mod_domain, only : domain_t
  use mod_group, only : group_t
  use mod_cells, only : cells_t
  use mod_stat_func, only : eval_cluster_id_time, eval_cluster_breakup, eval_clustersize_decay,eval_final_cluster
  use mod_data, only : source_data_t
  implicit none

  integer :: i,j,k,nparticles,nargs,snap,nsnaps
  type(particle_t) :: particle
  type(domain_t) :: domain
  type(group_t) :: group
  type(cells_t) :: cells
  type(source_data_t) :: source_data
  character(len=300) :: buf,directory,mode
  integer, dimension(:), allocatable :: cluster_count, cluster_duration, cluster_decaytime
  integer, dimension(:,:), allocatable :: cluster_id, cluster_id_init
  real, dimension(:), allocatable :: cluster_dist, clustersize_decay
  real :: rsep,minsize,timestep,BTs
  integer(int64) :: tstart,tdur,tinc,tstep
  logical :: mean_only

  ! read inputs

  nargs = command_argument_count()
  if (nargs .lt. 8) then
     write(*, '(A)') 'cluster_dist_tri_tine [mode] [directory] [tstart] [tdur] [tinc] [rsep] [timestep] [BTs]'
     stop
  end if
  call get_command_argument(1, buf); mode = trim(adjustl(buf))
  call get_command_argument(2, buf); directory = trim(adjustl(buf))
  call get_command_argument(3, buf); read(buf, *) tstart
  call get_command_argument(4, buf); read(buf, *) tdur
  call get_command_argument(5, buf); read(buf, *) tinc
  call get_command_argument(6, buf); read(buf, *) rsep
  call get_command_argument(7, buf); read(buf, *) timestep
  call get_command_argument(8, buf); read(buf, *) BTs

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
  allocate(cluster_id(nparticles,nparticles))
  allocate(cluster_id_init(nparticles,nparticles))
  allocate(cluster_count(nparticles))
  allocate(cluster_duration(nparticles))
  allocate(cluster_decaytime(nparticles))
  allocate(clustersize_decay(nparticles))

  ! init group and cells

  call group%init(domain)
  minsize = 2*maxval(particle%typeradius(:)) + rsep
  call cells%init(group, minsize, .false.)
  cluster_id_init(:,:) = 0
  cluster_duration(:) = 0
  cluster_decaytime(:) = 0
  clustersize_decay(:) = 0
  cluster_id(:,:) = 0

  ! loop over snaps

  do i = 1, nsnaps

     tstep = tstart + (i-1)*tinc
     snap = source_data%step_index(tstep)

     cluster_id(:,:) = 0
     call source_data%read_snap(particle, domain, snap)
     call cells%reset_grid()
     call cells%assign_particles()
    
    ! store cluster list as previous if not
     if (i==1) then
        ! EVALUATE PARTICLES IN EACH UNIQUE CLUSTER, store unique cluster id for each
        call eval_cluster_id_time(cells, .true., rsep, cluster_id)
        ! get the number of particles in each cluster
        cluster_count(:) = 0
        do j = 1, nparticles
           cluster_count(j) = count(cluster_id(j,:)/=0)
        end do
        !print*,cluster_count
        do j = 1, nparticles
           if (cluster_count(j) .lt. 2) then           ! checks for cluster sizes less than 2
              cluster_id(j,:) = 0
              !where(cluster_id .eq. j) cluster_id = 0  ! and replaces monomers with ID 0, therefore not clusters
           end if
        end do
        !print*,cluster_id(2,1)
        cluster_id_init = cluster_id
     else
    ! IF I/=1, EVALUATE IF ANY PARTICLES IN CLUSTER ARE STILL PRESENT (IF SO, ADD +1TS OF DURATION). IF NOT, OUTPUT CURRENT CLUSTER LIFETIME TO HIST WITH DURATIONS FOR CLUSTER OF INITIAL SIZE
        call eval_cluster_breakup(cells, .true., rsep, cluster_id_init, cluster_duration, cluster_decaytime)
        !print*,cluster_duration
     end if
     if (i == nsnaps) then
        call eval_final_cluster(cells, cluster_duration, cluster_decaytime)
     end if
  end do
  !print*,cluster_decaytime
! CALCULATE MEAN CLUSTER DECAY TIME FOR EACH CLUSTER SIZE
  call eval_clustersize_decay(nparticles, cluster_count, cluster_decaytime, clustersize_decay)
! Print
  do j = 1, nparticles
!      if (clustersize_decay(j) > 0.0) then
!         write(*, '(i14,e14.7)') j, clustersize_decay(j)
!      end if
     if (cluster_decaytime(j) /= 0.0) then
        write(*, '(i14,i14)') cluster_count(j), cluster_decaytime(j)
     end if
  end do

end program cluster_dist_tri
