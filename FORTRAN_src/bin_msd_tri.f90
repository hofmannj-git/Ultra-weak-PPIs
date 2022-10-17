! IMPORTANT: assumes unwrapped coord input!
! TODO: add flag to data read that says unwrapped or not
! TODO: add switch to unwrap coords if wrapped and has image

program bin_msd_tri

  use, intrinsic :: iso_fortran_env
  use mod_particle, only : particle_t
  use mod_domain, only : domain_t, propagate_in_flow
  use mod_group, only : group_t
  use mod_cells, only : cells_t
  use mod_stat_func, only : eval_nc, get_displ_var
  use mod_data, only : source_data_t
  implicit none

  integer :: i,j,nparticles,nargs,snap,nsnaps
  type(particle_t) :: particle_init,particle_flow,particle
  type(domain_t) :: domain_init,domain
  type(group_t) :: group
  type(cells_t) :: cells
  type(source_data_t) :: source_data
  character(len=300) :: buf,directory,mode_str,subtr
  character(len=2) :: mode
  integer, dimension(:), allocatable :: nc
  real :: rsep,minsize
  integer(int64) :: tstart,tdur,tinc,tstep
  integer, parameter :: ncmax = 20
  integer, parameter :: part_numtypes = 3
  real, dimension(6) :: var
  real, dimension(ncmax+1) :: var_by_nc
  integer :: modeint
  integer, dimension(ncmax+1) :: nc_igroup
  integer, dimension(part_numtypes+1) :: parttype_igroup
  logical, dimension(:), allocatable :: group_membership
  integer, dimension(3) :: delta_img
  logical :: nobin, subtr_mean, parttype

  ! read inputs

  nargs = command_argument_count()
  if (nargs .lt. 8) then
  write(*, '(A)') 'bin_msd_tri_restart [mode] [directory] [tstart] [tdur] [tinc] [rsep] [subtr_mean] {nobin} {parttype}'
     stop
  end if
  call get_command_argument(1, buf); mode_str = trim(adjustl(buf))
  call get_command_argument(2, buf); directory = trim(adjustl(buf))
  call get_command_argument(3, buf); read(buf, *) tstart
  call get_command_argument(4, buf); read(buf, *) tdur
  call get_command_argument(5, buf); read(buf, *) tinc
  call get_command_argument(6, buf); read(buf, *) rsep
  call get_command_argument(7, buf); subtr = trim(adjustl(buf))
  if (trim(adjustl(subtr)) == 'yes') then
     subtr_mean = .true.
  else if (trim(adjustl(subtr)) == 'no') then
     subtr_mean = .false.
  else
     print *, 'argument for [subtr_mean] can only be yes or no'
  end if
  nobin = .false.
  if (nargs .gt. 7) then
     call get_command_argument(8, buf)
     if (trim(adjustl(buf)) == 'nobin') then
        nobin = .true.
     else
        print *, 'eighth argument can only be nobin'
     end if
  end if
  parttype = .false.
  if (nargs .gt. 8) then
     call get_command_argument(9, buf)
     if (trim(adjustl(buf)) == 'parttype') then
        parttype = .true.
     else
        print *, 'ninth argument can only be parttype'
     end if
  end if

  ! can output just msd (trace) or components of the displacement variance tensor
  ! for the trace option, output should be 6Dt for non-interacting particles
  ! note: using Voigt notation, like in H in LAMMPS

  if (trim(adjustl(mode_str)) == 'tr') then
     mode = 'tr'
     modeint = 0
  else if (trim(adjustl(mode_str)) == 'xx') then
     mode = 'xx'
     modeint = 1
  else if (trim(adjustl(mode_str)) == 'yy') then
     mode = 'yy'
     modeint = 2
  else if (trim(adjustl(mode_str)) == 'zz') then
     mode = 'zz'
     modeint = 3
  else if (trim(adjustl(mode_str)) == 'yz') then
     mode = 'yz'
     modeint = 4
  else if (trim(adjustl(mode_str)) == 'xz') then
     mode = 'xz'
     modeint = 5
  else if (trim(adjustl(mode_str)) == 'xy') then
     mode = 'xy'
     modeint = 6
  else
     print *, 'mode must be tr (trace), xx, yy, zz, xy, xz, or yz'
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
  if (.not. nobin) then
     allocate(nc(nparticles))
  end if
  allocate(group_membership(nparticles))

  ! need two domain types and two particle types
  ! one for original coords and box
  ! one for time-advanced coords and box
  ! init by copying data read from data files

  call particle%copy(particle_init)
  call particle%copy(particle_flow)
  call domain_init%init(particle_init)
  snap = source_data%step_index(tstart)
  call source_data%set_box(domain_init, snap)

  ! init group, dynam and cells

  call group%init(domain_init)
  if (.not. nobin) then
     minsize = 2*maxval(particle%typeradius(:)) + rsep
     call cells%init(group, minsize,.false.)
  end if

  ! determine contact numbers at tstart

  call source_data%read_snap(particle_init, domain_init, snap)
  if (.not. nobin) then
     call cells%reset_grid()
     call cells%assign_particles()
     call eval_nc(cells, .true., rsep*2.0*(sum(particle%typeradius(:))/size(particle%typeradius(:))), nc)
  end if

  ! make groups one-at-a-time
  ! collect igroups as groups are created

  if (.not. nobin) then
     do i = 0, ncmax
        group_membership(:) = .false.
        do j = 1, nparticles
           if (nc(j) == i) then
              group_membership(j) = .true.
           end if
        end do
        nc_igroup(i+1) = group%create(group_membership)
     end do
  end if
  if (parttype) then
     do i = 0, part_numtypes
        group_membership(:) = .false.
        do j = 1, nparticles
           if (particle%type(j) == (i+1)) then
              group_membership(j) = .true.
           end if
        end do
        parttype_igroup(i+1) = group%create(group_membership)
     end do
  end if

  ! need to unmap particles, since assign_particles remaps them

  call domain_init%unmap()

  ! print zero-line of output

  write(*, '(i14)', advance='no') 0
  if (.not. nobin) then
     do j = 0, ncmax
        write(*, '(g14.5)', advance='no') 0.0
     end do
     write (*, '(a)') ''
  else if (parttype) then
     do j = 0, part_numtypes
        write(*, '(g14.5)', advance='no') 0.0
     end do
     write (*, '(a)') ''
  else
     write(*, '(g14.5)') 0.0
  end if

  ! loop over snaps

  delta_img(:) = 0

  do i = 2, nsnaps

     tstep = tstart + (i-1)*tinc
     snap = source_data%step_index(tstep)

     call source_data%read_snap(particle, domain, snap)

     ! update flow-advanced mean displacement particle type

     delta_img(:) = delta_img(:) + source_data%image_flips(:,i)
     call propagate_in_flow(domain_init, domain, delta_img, particle_flow)

     ! determine displacement variance
     ! fill output var_by_nc array by contact number

     if (.not. nobin) then
        do j = 0, ncmax
           call get_displ_var(particle_flow, particle, group, var, subtr_mean, igroup=nc_igroup(j+1))
           select case (modeint)
           case(0)
              var_by_nc(j+1) = var(1) + var(2) + var(3)
           case(1:6)
              var_by_nc(j+1) = var(modeint)
           end select
        end do
     else if (parttype) then
        do j = 0, part_numtypes
           call get_displ_var(particle_flow, particle, group, var,  subtr_mean, igroup=parttype_igroup(j+1))
           select case (modeint)
           case(0)
              var_by_nc(j+1) = var(1) + var(2) + var(3)
           case(1:6)
              var_by_nc(j+1) = var(modeint)
           end select
        end do
     else
        call get_displ_var(particle_flow, particle, group, var, subtr_mean)
        select case (modeint)
        case(0)
           var_by_nc(1) = var(1) + var(2) + var(3)
        case(1:6)
           var_by_nc(1) = var(modeint)
        end select
     end if

     ! print output

     write(*, '(i14)', advance='no') tstep-tstart
     if (.not. nobin) then
        do j = 0, ncmax
           write(*, '(g14.5)', advance='no') var_by_nc(j+1)
        end do
        write (*, '(a)') ''
     else if (parttype) then
        do j = 0, part_numtypes
           write(*, '(g14.5)', advance='no') var_by_nc(j+1)
        end do
        write (*, '(a)') ''
     else
        write(*, '(g14.5)') var_by_nc(1)
     end if

  end do

end program bin_msd_tri
