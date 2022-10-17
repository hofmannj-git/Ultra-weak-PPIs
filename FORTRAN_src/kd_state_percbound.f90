! IMPORTANT: assumes unwrapped coord input!
! TODO: add flag to data read that says unwrapped or not
! TODO: add switch to unwrap coords if wrapped and has image

program kd_state_percbound

  use, intrinsic :: iso_fortran_env
  use mod_particle, only : particle_t
  use mod_domain, only : domain_t
  use mod_group, only : group_t
  use mod_cells, only : cells_t
  use mod_stat_func, only : eval_binding, add_numbound, add_tempbound, &
                            add_state, elim_tempbound, elim_tempunbound, &
                            check_tempbound, check_tempunbound, &
                            add_numbound_time, final_check_bound, &
                            final_check_search
  use mod_data, only : source_data_t
  implicit none

  character(len=300) :: buf,directory
  integer(int64) :: tstart,tdur,tinc,tstep
  real :: rmax,rmin,minsize,rval,rgrid
  integer :: nbin,i,j,k,m,n,nsnaps,nargs,snap,nparticles,group_a,group_b, &
             type_a,type_b,num_bound,num_a,num_b,c,groupbit_a,groupbit_b, &
             num_bindevents, num_unbindevents, num_times, num_bound_time, &
             num_times_cut, count_num, num_bindevents_alt
  type(source_data_t) :: source_data
  type(particle_t) :: particle
  type(group_t) :: group
  type(domain_t) :: domain
  type(cells_t) :: cells
  integer, dimension(:), allocatable :: state_space, bindtimes_hist_alt, &
                                        unbindtimes_hist, state_space_time
  real, dimension(:), allocatable :: bin_bindunbtimes
  integer, dimension(:,:), allocatable :: temp_unbound, temp_bound_alt
  logical, dimension(:), allocatable :: group_ids
  real :: cut_dist, BTs, R_s_um, D_A, D_B, rad_a_um, rad_b_um, kon_smol, &
          kon_NBT, koff, Kd_NBT, Kd_smol, Kd_conc, perc_bound_avg, timestep, &
          avg_bind_time, avg_unbind_time, mult_c, cut_time, &
          perc_bound_time_avg, Kd_conc_time, avg_bind_time_alt
  real(16), parameter :: PI = 4 * atan (1.0_16)
!  integer, dimension(:), pointer :: mask => null()

  ! read inputs
  nargs = command_argument_count()
  if (nargs .lt. 13) then
     write(*, '(A)') 'kd_state_percbound [directory] [tstart] [tdur] [tinc] [cut_dist]'
     write(*, '(A)') '[type_a] [type_b] [BTs] [R_s_um] [D_A] [D_B] [ts] [cut_time]'
     stop
  end if
  call get_command_argument(1, buf); directory = trim(adjustl(buf))
  call get_command_argument(2, buf); read(buf, *) tstart
  call get_command_argument(3, buf); read(buf, *) tdur
  call get_command_argument(4, buf); read(buf, *) tinc
  call get_command_argument(5, buf); read(buf, *) cut_dist
  call get_command_argument(6, buf); read(buf, *) type_a
  call get_command_argument(7, buf); read(buf, *) type_b
  call get_command_argument(8, buf); read(buf, *) BTs
  call get_command_argument(9, buf); read(buf, *) R_s_um
  call get_command_argument(10, buf); read(buf, *) D_A
  call get_command_argument(11, buf); read(buf, *) D_B
  call get_command_argument(12, buf); read(buf, *) timestep
  call get_command_argument(13, buf); read(buf, *) cut_time

! FROM DATA: box, Nl, Nm, Ns, Rl, Rm

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
  if (modulo(tdur,tinc) /= 0) then
     print *, 'tdur is not a multiple of tinc'
     stop
  end if
  nsnaps = int(tdur/tinc)
  nparticles = particle%nparticles

 ! init group for group pair in Kd calculation

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

  minsize = 2*maxval(particle%typeradius(:)) + cut_dist
  call cells%init(group, minsize, .false.)
  ! count number of each particle group, a & b
  groupbit_a = group%bitmask(group_a)
  groupbit_b = group%bitmask(group_b)
  num_a = group%count(group_a)
  num_b = group%count(group_b)

  ! allocating necessary arrays for bin, bin_unb, num_bound
  allocate(state_space(nparticles))
  state_space(:) = 0
  allocate(state_space_time(nparticles))
  state_space_time(:) = 0
  allocate(temp_unbound(nparticles,2))
  temp_unbound(:,:) = 0
  allocate(temp_bound_alt(nparticles,2))
  temp_bound_alt(:,:) = 0
  allocate(bin_bindunbtimes(tdur/tinc))
  bin_bindunbtimes(:) = 0
  allocate(bindtimes_hist_alt(tdur/tinc))
  bindtimes_hist_alt(:) = 0
  allocate(unbindtimes_hist(tdur/tinc))
  unbindtimes_hist(:) = 0
  num_bound = 0
  num_bound_time = 0

  ! loop over snaps, accumulating to average
  print *, group_a, group_b
  print *, 'looping over snaps'
  do i = 1, nsnaps
     tstep = tstart + (i-1)*tinc
     snap = source_data%step_index(tstep)
     call source_data%read_snap(particle, domain, snap)
     call cells%reset_grid()
     call cells%assign_particles()
     call eval_binding(cells, group, group_a, group_b, cut_dist, state_space) ! calculate state_space @ ts=i

     call add_numbound(cells, group, group_b, state_space, num_bound) ! add number of bound B at ts=i to total bound
!     print *, num_bound
     !print *, num_bound
     if (i == 1) then
        ! add all bound_pairs(A,B) to temp_bound w/ dur=1ts
        call add_tempbound(cells, group, group_b, &
                           temp_bound_alt, state_space, tinc)
        ! add all state_space(B=0) to temp_unbound w/ dur=1ts
        call add_state(cells, group, group_b, state_space, temp_unbound, tinc)
     else
        ! eliminate newly unbound pairs, if (A,B) in temp_bound not in bound_pairs for ts=i
        call elim_tempbound(tinc, state_space_time, cut_time, &
                            cells, group, group_b, state_space, temp_bound_alt, &
                            bindtimes_hist_alt)
        ! for (A,B) in bound_pairs for ts=i, check if in temp_bound & either add (A,B) w/ dur=1ts or add 1ts to existing dur
        call check_tempbound(tinc, cells, group, &
                             group_b, state_space, temp_bound_alt)
        if (i >= int(cut_time)) then
           call add_numbound_time(cells, group, group_b, state_space_time, num_bound_time) ! add number of bound B at ts=i to total bound
        end if
        ! eliminate newly bound B, if (B) in temp_unbound but state_space(B/=0) for ts=i
        call elim_tempunbound(cells, group, group_b, state_space, temp_unbound, tinc, unbindtimes_hist)
        ! for state_space(B=0) for ts=i, check if B in temp_unbound & either add B w/ dur=1ts or add 1ts to existing dur
        call check_tempunbound(cells, group, group_b, state_space, temp_unbound, tinc)
     end if
     if (i == nsnaps) then
        ! add current length of time in temp_bound to part_bound_times
        call final_check_bound(tinc, &
                                state_space_time, cut_time, &
                                cells, group, group_b, state_space, &
                                temp_bound_alt, bindtimes_hist_alt)
        ! add current length of time in temp_unbound to part_search_times
        call final_check_search(cells, group, group_b, state_space, temp_unbound, tinc, unbindtimes_hist)
     end if

  end do

  num_times = tdur/tinc
  do k = 1, num_times
     bin_bindunbtimes(k) = (k-1/2)*timestep*tinc*BTs*4
  end do

  ! use (num_bound) to calc. Kd from concentration
  perc_bound_avg = real(num_bound) / real(num_b) / real(num_times+1)
  write(*, '(a)') 'The average percent of B bound is :'
  write(*, '(e12.3)') perc_bound_avg
  Kd_conc = 1/(6.022*(10**8))/((2*R_s_um*cells%domain%xprd)**3) * &
            (1 - perc_bound_avg)**2 / perc_bound_avg * num_b
  write(*, '(a)') 'Kd calculated with concentration is [M]:'
  write(*, '(e12.3)') Kd_conc

  ! output data to files
  open( unit = 103, file = 'bindunbtimebin_state.out' )
  open( unit = 1002, file = 'unbindtimehist_state.out' )
  open( unit = 1004, file = 'bindtimehist_state.out' )
  do k = 1, num_times
     write( 103, * ) bin_bindunbtimes(k)
     write( 1002, * ) unbindtimes_hist(k)
     write( 1004, * ) bindtimes_hist_alt(k)
  end do
  close( 1004 )
  close( 1002 )
  close( 103 )

end program kd_state_percbound
