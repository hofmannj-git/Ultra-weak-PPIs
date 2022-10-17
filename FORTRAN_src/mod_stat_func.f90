! TODO add group selection and exclusion

module mod_stat_func

  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  use mod_constants, only : pi, twopi
  use mod_particle, only : particle_t
  use mod_group, only : group_t
  use mod_cells, only : cells_t
  implicit none
  private

  public :: eval_rdf, eval_nc, eval_nc_parttype, get_displ_var
  public :: get_mean_pos, get_displ_var_sq, get_displ_var_cr
  public :: eval_binding, add_numbound, add_tempbound, add_state
  public :: add_numbound_time, check_tempunbound
  public :: elim_tempbound, elim_tempunbound, check_tempbound
  public :: final_check_bound, final_check_search, symmetrize_3d
  public :: eval_pairdisttime, eval_bindsearch_traj, eval_pairdisttime_search
  public :: eval_cluster_id_time,eval_cluster_breakup, eval_cluster_id
  public :: eval_bindsearch_traj_strictloose, eval_cluster_id_parttype
  public :: eval_clustersize_decay,eval_final_cluster, eval_pair_dist_func

contains

  ! ----------------------------------------------------------------------
  ! Calculate RDF in periodic/nonperiodic and rectangular/triclinic boxes
  ! TODO add group specification
  ! ----------------------------------------------------------------------

  subroutine eval_rdf(cells, group, agroup, bgroup, rcut, rmin, nbin, rdf)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    real, intent(in) :: rcut,rmin
    integer, intent(in) :: agroup,bgroup,nbin
    real, dimension(:), intent(out) :: rdf
    integer(int64), dimension(nbin) :: hist
    real :: rsq,rcutsq,rminsq,normfac,rgrid,inv_rgrid,delx,dely,delz,xtmp,ytmp,ztmp
    real :: normdum_a, normdum_b, normdum_c
    integer :: i,j,c,d,n,groupbit_a,groupbit_b,othergroup,neighcell,na,nb
    integer :: nparticles
    integer, dimension(:), pointer :: mask => null()
    logical :: groups_identical

    ! If group ids are not equal, cannot have a particle in both groups.
    ! TODO insert this conditional with a new group function

    na = group%count(agroup)
    nb = group%count(bgroup)
    if (na < 1) then
       write(*, '(a)') 'group A in pair distribution function has bad count'
       stop
    end if
    if (nb < 1) then
       write(*, '(a)') 'group B in pair distribution function has bad count'
       stop
    end if

    groups_identical = .false.
    if (agroup == bgroup) groups_identical = .true.

    groupbit_a = group%bitmask(agroup)
    groupbit_b = group%bitmask(bgroup)
    mask => group%mask
    rcutsq = rcut * rcut
    rminsq = rmin * rmin
    rgrid = (rcut - rmin) / nbin
    inv_rgrid = 1.0 / rgrid
    hist(:) = 0

    ! check if particle coordinates need wrapping, and wrap if they do
    ! can do this here, because no reason to be in lamda coords

    if (.not. cells%domain%particles_inside()) then
       call cells%domain%remap()
    end if

    ! if cells are functional, then use them to calculate RDF
    ! if not, loop over particle pairs

    if_cells: if (cells%cells_exist) then

       cell_loop: do c = 1, cells%ncells

          ! below, all in the ith cell
          ! loop over all particles i,j in this cell c
          ! cellhead(c) is the first i particle
          ! the first (distinct) j is the next of i

          i = cells%cellhead(c)
          self_loop: do while (i > 0)

             if (iand(mask(i),groupbit_a) /= 0) then
                othergroup = groupbit_b
             else if (iand(mask(i),groupbit_b) /= 0) then
                othergroup = groupbit_a
             else
                i = cells%next(i)
                cycle
             end if

             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)
             j = cells%next(i)
             do while (j > 0)
                if (iand(mask(j),othergroup) == 0) then
                   j = cells%next(j)
                   cycle
                end if
                delx = cells%particle%x(1,j) - xtmp
                dely = cells%particle%x(2,j) - ytmp
                delz = cells%particle%x(3,j) - ztmp
                call cells%domain%minimum_image(delx, dely, delz)
                rsq = delx*delx + dely*dely + delz*delz
                if ((rsq < rcutsq) .and. (rsq > rminsq)) then
                   n = floor((sqrt(rsq) - rmin) * inv_rgrid) + 1
                   n = min(n,nbin)
                   hist(n) = hist(n) + 1
                end if
                j = cells%next(j)
             end do
             i = cells%next(i)
          end do self_loop

          ! loop over distinct cells
          ! should I make this part of previous loop?

          i = cells%cellhead(c)
          distinct_loop: do while (i > 0)

             if (iand(mask(i),groupbit_a) /= 0) then
                othergroup = groupbit_b
             else if (iand(mask(i),groupbit_b) /= 0) then
                othergroup = groupbit_a
             else
                i = cells%next(i)
                cycle
             end if

             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)
             do d = 1, cells%nstencil
                neighcell = cells%stencil_neighbor(c,d)
                if (neighcell == 0) cycle
                j = cells%cellhead(neighcell)
                do while (j > 0)
                   if (iand(mask(j),othergroup) == 0) then
                      j = cells%next(j)
                      cycle
                   end if
                   delx = cells%particle%x(1,j) - xtmp
                   dely = cells%particle%x(2,j) - ytmp
                   delz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delx, dely, delz)
                   rsq = delx*delx + dely*dely + delz*delz
                   if ((rsq < rcutsq) .and. (rsq > rminsq)) then
                      n = floor((sqrt(rsq) - rmin) * inv_rgrid) + 1
                      n = min(n,nbin)
                      hist(n) = hist(n) + 1
                   end if
                   j = cells%next(j)
                end do
             end do
             i = cells%next(i)
          end do distinct_loop
       end do cell_loop

    else ! not using cells, so loop over all pairs

       nparticles = cells%particle%nparticles
       do i = 1, nparticles-1

          if (iand(mask(i),groupbit_a) /= 0) then
             othergroup = groupbit_b
          else if (iand(mask(i),groupbit_b) /= 0) then
             othergroup = groupbit_a
          else
             cycle
          end if

          xtmp = cells%particle%x(1,i)
          ytmp = cells%particle%x(2,i)
          ztmp = cells%particle%x(3,i)
          do j = i+1, nparticles
             if (iand(mask(j),othergroup) == 0) cycle
             delx = cells%particle%x(1,j) - xtmp
             dely = cells%particle%x(2,j) - ytmp
             delz = cells%particle%x(3,j) - ztmp
             call cells%domain%minimum_image(delx, dely, delz)
             rsq = delx*delx + dely*dely + delz*delz
             if ((rsq < rcutsq) .and. (rsq > rminsq)) then
                n = floor((sqrt(rsq) - rmin) * inv_rgrid) + 1
                n = min(n,nbin)
                hist(n) = hist(n) + 1
             end if
          end do
       end do

    end if if_cells

    ! Normalize the histogram.

    normfac = 3 * product(cells%domain%prd(:)) &
         / (4 * pi * na * nb)
    normdum_a = 3*rmin*rmin*rgrid
    normdum_b = 3*rmin*rgrid*rgrid
    normdum_c = rgrid*rgrid*rgrid
    if (groups_identical) normfac = normfac * 2
    do i = 1, nbin
       rdf(i) = hist(i) * normfac &
                /(normdum_a + normdum_b*(2*i - 1) + normdum_c*(3*i*i - 3*i + 1))
    end do
  end subroutine eval_rdf

  ! ----------------------------------------------------------------------
  ! Calculate pair distribution function
  ! nbin is the # of bins in one direction (used for all)
  ! rmax is the maximum distance along one axis
  ! in reality, the maximum distance is sqrt(3)*rmax in corner of cube
  ! agroup is the group of the center particle
  ! bgroup is the group of surrounding particles
  ! ----------------------------------------------------------------------

  subroutine eval_pair_dist_func(cells, group, agroup, bgroup, rmax, nbin, g)

    class(cells_t), intent(inout) :: cells
    class(group_t), intent(in), target :: group
    real, intent(in) :: rmax
    integer, intent(in) :: agroup,bgroup,nbin
    real, dimension(:,:,:), intent(out) :: g
    integer(int64), dimension(nbin,nbin,nbin) :: hist
    real :: normfac,delx,dely,delz,xtmp,ytmp,ztmp,rgrid,inv_rgrid
    integer :: i,j,c,d,xbin,ybin,zbin,neighcell,na,nb,groupbit_a,groupbit_b
    integer :: nparticles,othergroup,vec_sign
    logical :: groups_identical
    integer, dimension(:), pointer :: mask => null()

    if (mod(nbin,2) == 0) then
       write(*, '(a)') 'nbin must be odd in eval_pair_dist_func'
       stop
    end if

    na = group%count(agroup)
    nb = group%count(bgroup)
    if (na < 1) then
       write(*, '(a)') 'group A in pair distribution function has bad count'
       stop
    end if
    if (nb < 1) then
       write(*, '(a)') 'group B in pair distribution function has bad count'
       stop
    end if

    groups_identical = .false.
    if (agroup == bgroup) groups_identical = .true.
    mask => group%mask
    hist(:,:,:) = 0_int64
    rgrid = 2 * rmax / nbin
    inv_rgrid = 1.0 / rgrid

    ! check if particle coordinates need wrapping, and wrap if they do
    ! can do this here, because no reason to be in lamda coords

    if (.not. cells%domain%particles_inside()) then
       call cells%domain%remap()
    end if

    ! if cells are functional, then use them to calculate RDF
    ! if not, loop over particle pairs

    groupbit_a = group%bitmask(agroup)
    groupbit_b = group%bitmask(bgroup)

    if_cells: if (cells%cells_exist) then

       cell_loop: do c = 1, cells%ncells

          ! below, all in the ith cell
          ! loop over all particles i,j in this cell c
          ! cellhead(c) is the first i particle
          ! the first (distinct) j is the next of i

          i = cells%cellhead(c)
          self_loop: do while (i > 0)

             if (iand(mask(i),groupbit_a) /= 0) then
                vec_sign = 1
                othergroup = groupbit_b
             else if (iand(mask(i),groupbit_b) /= 0) then
                vec_sign = -1
                othergroup = groupbit_a
             else
                i = cells%next(i)
                cycle
             end if

             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)

             j = cells%next(i)
             do while (j > 0)

                if (iand(mask(j),othergroup) == 0) then
                   j = cells%next(j)
                   cycle
                end if

                delx = cells%particle%x(1,j) - xtmp
                dely = cells%particle%x(2,j) - ytmp
                delz = cells%particle%x(3,j) - ztmp
                delx = vec_sign * delx
                dely = vec_sign * dely
                delz = vec_sign * delz

                call cells%domain%minimum_image(delx, dely, delz)

                xbin = ceiling( (delx+rmax) * inv_rgrid )
                ybin = ceiling( (dely+rmax) * inv_rgrid )
                zbin = ceiling( (delz+rmax) * inv_rgrid )

                if ((xbin < 1) .or. (xbin > nbin)) then
                   j = cells%next(j)
                   cycle
                end if
                if ((ybin < 1) .or. (ybin > nbin)) then
                   j = cells%next(j)
                   cycle
                end if
                if ((zbin < 1) .or. (zbin > nbin)) then
                   j = cells%next(j)
                   cycle
                end if
                
                hist(xbin,ybin,zbin) = hist(xbin,ybin,zbin) + 1
                j = cells%next(j)
             end do
             i = cells%next(i)
          end do self_loop

          ! loop over distinct cells
          ! should I make this part of previous loop?

          i = cells%cellhead(c)
          distinct_loop: do while (i > 0)

             if (iand(mask(i),groupbit_a) /= 0) then
                vec_sign = 1
                othergroup = groupbit_b
             else if (iand(mask(i),groupbit_b) /= 0) then
                vec_sign = -1
                othergroup = groupbit_a
             else
                i = cells%next(i)
                cycle
             end if

             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)

             do d = 1, cells%nstencil
                neighcell = cells%stencil_neighbor(c,d)
                if (neighcell == 0) cycle
                j = cells%cellhead(neighcell)
                do while (j > 0)

                   if (iand(mask(j),othergroup) == 0) then
                      j = cells%next(j)
                      cycle
                   end if

                   delx = cells%particle%x(1,j) - xtmp
                   dely = cells%particle%x(2,j) - ytmp
                   delz = cells%particle%x(3,j) - ztmp
                   delx = vec_sign * delx
                   dely = vec_sign * dely
                   delz = vec_sign * delz
                   call cells%domain%minimum_image(delx, dely, delz)

                   xbin = ceiling( (delx+rmax) * inv_rgrid )
                   ybin = ceiling( (dely+rmax) * inv_rgrid )
                   zbin = ceiling( (delz+rmax) * inv_rgrid )

                   if ((xbin < 1) .or. (xbin > nbin)) then
                      j = cells%next(j)
                      cycle
                   end if
                   if ((ybin < 1) .or. (ybin > nbin)) then
                      j = cells%next(j)
                      cycle
                   end if
                   if ((zbin < 1) .or. (zbin > nbin)) then
                      j = cells%next(j)
                      cycle
                   end if
                
                   hist(xbin,ybin,zbin) = hist(xbin,ybin,zbin) + 1

                   j = cells%next(j)
                end do
             end do
             i = cells%next(i)
          end do distinct_loop
       end do cell_loop

    else ! not using cells, so loop over all pairs

       nparticles = cells%particle%nparticles
       do i = 1, nparticles-1

          if (iand(mask(i),groupbit_a) /= 0) then
             vec_sign = 1
             othergroup = groupbit_b
          else if (iand(mask(i),groupbit_b) /= 0) then
             vec_sign = -1
             othergroup = groupbit_a
          else
             cycle
          end if

          xtmp = cells%particle%x(1,i)
          ytmp = cells%particle%x(2,i)
          ztmp = cells%particle%x(3,i)
          do j = i+1, nparticles
             if (iand(mask(j),othergroup) == 0) cycle
             delx = cells%particle%x(1,j) - xtmp
             dely = cells%particle%x(2,j) - ytmp
             delz = cells%particle%x(3,j) - ztmp
             delx = vec_sign * delx
             dely = vec_sign * dely
             delz = vec_sign * delz
             call cells%domain%minimum_image(delx, dely, delz)

             xbin = ceiling( (delx+rmax) * inv_rgrid )
             ybin = ceiling( (dely+rmax) * inv_rgrid )
             zbin = ceiling( (delz+rmax) * inv_rgrid )

             if (xbin < 1) cycle
             if (ybin < 1) cycle
             if (zbin < 1) cycle
             if (xbin > nbin) cycle
             if (ybin > nbin) cycle
             if (zbin > nbin) cycle
                
             hist(xbin,ybin,zbin) = hist(xbin,ybin,zbin) + 1

          end do
       end do

    end if if_cells

    ! Normalize the histogram.

    normfac = product(cells%domain%prd(:)) &
         / (rgrid*rgrid*rgrid * na * nb)
    if (groups_identical) normfac = 2 * normfac
    g(:,:,:) = normfac * hist(:,:,:)
    if (groups_identical) call symmetrize_3d(g(:,:,:))

  end subroutine eval_pair_dist_func


  ! ----------------------------------------------------------------------
  ! Calculate contact numbers, with/without cells in
  ! rectangular/triclinic box in periodic/nonperiodic boundaries
  ! Can apply a cutoff based on a_i+a_j+delta (TRUE) or a uniform
  ! cutoff (FALSE) in the yes_delta argument.  The distance r that follows
  ! is either the delta or the uniform cutoff.
  ! Get the number of bins from the dimensions of the ncdist array
  ! Assuming that r>0 if yes_delta is false
  ! And that r+a_i+a_j > 0 for all i,j if yes_delta is true
  ! ----------------------------------------------------------------------

  subroutine eval_nc(cells, yes_delta, r, nc)

    type(cells_t), intent(inout) :: cells
    logical, intent(in) :: yes_delta
    real, intent(in) :: r
    integer, dimension(:), intent(inout) :: nc
    integer :: c,d,i,j,neighcell,itype,jtype,nparticles,ntypes
    logical :: typeradius_flag
    real :: rcutsq,rsq,delx,dely,delz,xtmp,ytmp,ztmp
    real, dimension(:,:), allocatable :: rsqmat

    nparticles = cells%particle%nparticles
    typeradius_flag = cells%particle%typeradius_flag
    ntypes = cells%particle%ntypes

    ! default rcutsq for yes_delta is false
    rcutsq = r*r

    if (yes_delta .and. typeradius_flag) then
       allocate(rsqmat(ntypes,ntypes))
       do i = 1, ntypes
          do j = 1, ntypes
             rsqmat(i,j) = cells%particle%typeradius(i) &
                         + cells%particle%typeradius(j) &
                         + r
          end do
       end do
       rsqmat(:,:) = rsqmat(:,:)*rsqmat(:,:)
    end if

    nc(:) = 0

    ! keeping logical evaluation out of inner loop
    ! first case: uniform cutoff

    option_switch: if (.not. yes_delta) then

       if_cells_uniform: if (cells%cells_exist) then

          cell_loop_uniform: do c = 1, cells%ncells

             i = cells%cellhead(c)
             do while (i > 0)
                xtmp = cells%particle%x(1,i)
                ytmp = cells%particle%x(2,i)
                ztmp = cells%particle%x(3,i)
                j = cells%next(i)
                do while (j > 0)
                   delx = cells%particle%x(1,j) - xtmp
                   dely = cells%particle%x(2,j) - ytmp
                   delz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delx, dely, delz)
                   rsq = delx*delx + dely*dely + delz*delz
                   if (rsq < rcutsq) then
                      nc(i) = nc(i) + 1
                      nc(j) = nc(j) + 1
                   end if
                   j = cells%next(j)
                end do
                i = cells%next(i)
             end do

             ! loop over distinct cells
             ! should I make this part of previous loop?

             i = cells%cellhead(c)
             do while (i > 0)
                xtmp = cells%particle%x(1,i)
                ytmp = cells%particle%x(2,i)
                ztmp = cells%particle%x(3,i)
                do d = 1, cells%nstencil
                   neighcell = cells%stencil_neighbor(c,d)
                   if (neighcell == 0) cycle
                   j = cells%cellhead(neighcell)
                   do while (j > 0)
                      delx = cells%particle%x(1,j) - xtmp
                      dely = cells%particle%x(2,j) - ytmp
                      delz = cells%particle%x(3,j) - ztmp
                      call cells%domain%minimum_image(delx, dely, delz)
                      rsq = delx*delx + dely*dely + delz*delz
                      if (rsq < rcutsq) then
                         nc(i) = nc(i) + 1
                         nc(j) = nc(j) + 1
                      end if
                      j = cells%next(j)
                   end do
                end do
                i = cells%next(i)
             end do
                   
          end do cell_loop_uniform

       else ! not using cells, so loop over all pairs

          do i = 1, nparticles-1
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)
             do j = i+1, nparticles
                delx = cells%particle%x(1,j) - xtmp
                dely = cells%particle%x(2,j) - ytmp
                delz = cells%particle%x(3,j) - ztmp
                call cells%domain%minimum_image(delx, dely, delz)
                rsq = delx*delx + dely*dely + delz*delz
                if (rsq < rcutsq) then
                   nc(i) = nc(i) + 1
                   nc(j) = nc(j) + 1
                end if
             end do
          end do

       end if if_cells_uniform

    else if (typeradius_flag) then ! can use matrix cutoffs

       if_cells_typeradius: if (cells%cells_exist) then

          cell_loop_typeradius: do c = 1, cells%ncells

             ! below, all in the ith cell
             ! loop over all particles i,j in this cell c
             ! cellhead(c) is the first i particle
             ! the first (distinct) j is the next of i

             i = cells%cellhead(c)
             do while (i > 0)
                itype = cells%particle%type(i)
                xtmp = cells%particle%x(1,i)
                ytmp = cells%particle%x(2,i)
                ztmp = cells%particle%x(3,i)
                j = cells%next(i)
                do while (j > 0)
                   jtype = cells%particle%type(j)
                   delx = cells%particle%x(1,j) - xtmp
                   dely = cells%particle%x(2,j) - ytmp
                   delz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delx, dely, delz)
                   rsq = delx*delx + dely*dely + delz*delz
                   rcutsq = rsqmat(itype,jtype)
                   if (rsq < rcutsq) then
                      nc(i) = nc(i) + 1
                      nc(j) = nc(j) + 1
                   end if
                   j = cells%next(j)
                end do
                i = cells%next(i)
             end do

             ! loop over distinct cells
             ! should I make this part of previous loop?

             i = cells%cellhead(c)
             do while (i > 0)
                itype = cells%particle%type(i)
                xtmp = cells%particle%x(1,i)
                ytmp = cells%particle%x(2,i)
                ztmp = cells%particle%x(3,i)
                do d = 1, cells%nstencil
                   neighcell = cells%stencil_neighbor(c,d)
                   if (neighcell == 0) cycle
                   j = cells%cellhead(neighcell)
                   do while (j > 0)
                      jtype = cells%particle%type(j)
                      delx = cells%particle%x(1,j) - xtmp
                      dely = cells%particle%x(2,j) - ytmp
                      delz = cells%particle%x(3,j) - ztmp
                      call cells%domain%minimum_image(delx, dely, delz)
                      rsq = delx*delx + dely*dely + delz*delz
                      rcutsq = rsqmat(itype,jtype)
                      if (rsq < rcutsq) then
                         nc(i) = nc(i) + 1
                         nc(j) = nc(j) + 1
                      end if
                      j = cells%next(j)
                   end do
                end do
                i = cells%next(i)
             end do
                   
          end do cell_loop_typeradius

       else ! not using cells, so loop over all pairs

          do i = 1, nparticles-1
             itype = cells%particle%type(i)
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)
             do j = i+1, nparticles
                jtype = cells%particle%type(j)
                delx = cells%particle%x(1,j) - xtmp
                dely = cells%particle%x(2,j) - ytmp
                delz = cells%particle%x(3,j) - ztmp
                call cells%domain%minimum_image(delx, dely, delz)
                rsq = delx*delx + dely*dely + delz*delz
                rcutsq = rsqmat(itype,jtype)
                if (rsq < rcutsq) then
                   nc(i) = nc(i) + 1
                   nc(j) = nc(j) + 1
                end if
             end do
          end do

       end if if_cells_typeradius

    else ! have to use per-particle cutoff based on delta

       if_cells_perparticle: if (cells%cells_exist) then

          cell_loop_parparticle: do c = 1, cells%ncells

             i = cells%cellhead(c)
             do while (i > 0)
                xtmp = cells%particle%x(1,i)
                ytmp = cells%particle%x(2,i)
                ztmp = cells%particle%x(3,i)
                j = cells%next(i)
                do while (j > 0)
                   delx = cells%particle%x(1,j) - xtmp
                   dely = cells%particle%x(2,j) - ytmp
                   delz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delx, dely, delz)
                   rsq = delx*delx + dely*dely + delz*delz
                   rcutsq = cells%particle%radius(i)
                   rcutsq = rcutsq + cells%particle%radius(j)
                   rcutsq = rcutsq + r
                   rcutsq = rcutsq*rcutsq
                   if (rsq < rcutsq) then
                      nc(i) = nc(i) + 1
                      nc(j) = nc(j) + 1
                   end if
                   j = cells%next(j)
                end do
                i = cells%next(i)
             end do

             ! loop over distinct cells
             ! should I make this part of previous loop?

             i = cells%cellhead(c)
             do while (i > 0)
                xtmp = cells%particle%x(1,i)
                ytmp = cells%particle%x(2,i)
                ztmp = cells%particle%x(3,i)
                do d = 1, cells%nstencil
                   neighcell = cells%stencil_neighbor(c,d)
                   if (neighcell == 0) cycle
                   j = cells%cellhead(neighcell)
                   do while (j > 0)
                      delx = cells%particle%x(1,j) - xtmp
                      dely = cells%particle%x(2,j) - ytmp
                      delz = cells%particle%x(3,j) - ztmp
                      call cells%domain%minimum_image(delx, dely, delz)
                      rsq = delx*delx + dely*dely + delz*delz
                      rcutsq = cells%particle%radius(i)
                      rcutsq = rcutsq + cells%particle%radius(j)
                      rcutsq = rcutsq + r
                      rcutsq = rcutsq*rcutsq
                      if (rsq < rcutsq) then
                         nc(i) = nc(i) + 1
                         nc(j) = nc(j) + 1
                      end if
                      j = cells%next(j)
                   end do
                end do
                i = cells%next(i)
             end do
                   
          end do cell_loop_parparticle

       else ! not using cells, so loop over all pairs

          do i = 1, nparticles-1
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)
             do j = i+1, nparticles
                delx = cells%particle%x(1,j) - xtmp
                dely = cells%particle%x(2,j) - ytmp
                delz = cells%particle%x(3,j) - ztmp
                call cells%domain%minimum_image(delx, dely, delz)
                rsq = delx*delx + dely*dely + delz*delz
                rcutsq = cells%particle%radius(i)
                rcutsq = rcutsq + cells%particle%radius(j)
                rcutsq = rcutsq + r
                rcutsq = rcutsq*rcutsq
                if (rsq < rcutsq) then
                   nc(i) = nc(i) + 1
                   nc(j) = nc(j) + 1
                end if
             end do
          end do

       end if if_cells_perparticle

    end if option_switch

  end subroutine eval_nc

  ! ----------------------------------------------------------------------
  ! Calculate contact numbers, with/without cells in
  ! rectangular/triclinic box in periodic/nonperiodic boundaries
  ! Can apply a cutoff based on a_i+a_j+delta (TRUE) or a uniform
  ! cutoff (FALSE) in the yes_delta argument.  The distance r that follows
  ! is either the delta or the uniform cutoff.
  ! Get the number of bins from the dimensions of the ncdist array
  ! Assuming that r>0 if yes_delta is false
  ! And that r+a_i+a_j > 0 for all i,j if yes_delta is true
  ! ----------------------------------------------------------------------

  subroutine eval_nc_parttype(cells, group, agroup, bgroup, parttypea, parttypeb, r, nc)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    real, intent(in) :: r
    integer, dimension(:), intent(inout) :: nc
    integer, intent(in) :: agroup,bgroup
    integer :: c,d,i,j,neighcell,parttypea,parttypeb,nparticles,ntypes
    logical :: typeradius_flag, yes_delta, groups_identical
    real :: rcutsq,rsq,delx,dely,delz,xtmp,ytmp,ztmp
    integer :: n,groupbit_a,groupbit_b,othergroup,na,nb
    integer, dimension(:), pointer :: mask => null()

    na = group%count(agroup)
    nb = group%count(bgroup)
    if (na < 1) then
       write(*, '(a)') 'group A in pair distribution function has bad count'
       stop
    end if
    if (nb < 1) then
       write(*, '(a)') 'group B in pair distribution function has bad count'
       stop
    end if

    groups_identical = .false.
    if (agroup == bgroup) groups_identical = .true.

    groupbit_a = group%bitmask(agroup)
    groupbit_b = group%bitmask(bgroup)
    mask => group%mask
    rcutsq = r*r
    nc(:) = 0

    if (.not. cells%domain%particles_inside()) then
       call cells%domain%remap()
    end if

    ! keeping logical evaluation out of inner loop
    ! first case: uniform cutoff

    if_cells: if (cells%cells_exist) then

       cell_loop: do c = 1, cells%ncells
          ! below, all in the ith cell
          ! loop over all particles i,j in this cell c
          ! cellhead(c) is the first i particle
          ! the first (distinct) j is the next of i

          i = cells%cellhead(c)
          self_loop: do while (i > 0)

             if (cells%particle%type(i) /= parttypea) then
                i = cells%next(i)
                cycle
             end if
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)
             j = cells%next(i)
             do while (j > 0)
                if (cells%particle%type(j) /= parttypeb) then
                   j = cells%next(j)
                   cycle
                end if
                delx = cells%particle%x(1,j) - xtmp
                dely = cells%particle%x(2,j) - ytmp
                delz = cells%particle%x(3,j) - ztmp
                call cells%domain%minimum_image(delx, dely, delz)
                rsq = delx*delx + dely*dely + delz*delz
                !print *, 'final case ', rsq, rcutsq
                if (rsq < rcutsq) then
                    !print *, 'final case ', rsq, rcutsq, nc(i), i, cells%particle%type(i)
                    nc(i) = nc(i) + 1
                end if
                j = cells%next(j)
              end do
              i = cells%next(i)
          end do self_loop
          ! loop over distinct cells
          ! should I make this part of previous loop?

          i = cells%cellhead(c)
          distinct_loop: do while (i > 0)

             if (cells%particle%type(i) /= parttypea) then
                i = cells%next(i)
                cycle
             end if
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)
             do d = 1, cells%nstencil
                neighcell = cells%stencil_neighbor(c,d)
                if (neighcell == 0) cycle
                j = cells%cellhead(neighcell)
                do while (j > 0)
                   if (cells%particle%type(j) /= parttypeb) then
                      j = cells%next(j)
                      cycle
                   end if
                   delx = cells%particle%x(1,j) - xtmp
                   dely = cells%particle%x(2,j) - ytmp
                   delz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delx, dely, delz)
                   rsq = delx*delx + dely*dely + delz*delz
                   if (rsq < rcutsq) then
                       !print *, 'final case ', rsq, rcutsq, nc(i), i, cells%particle%type(i)
                       nc(i) = nc(i) + 1
                   end if
                   j = cells%next(j)
                end do
             end do
             i = cells%next(i)
          end do distinct_loop
       end do cell_loop
    end if if_cells
    !print*, nc

  end subroutine eval_nc_parttype

  ! ----------------------------------------------------------------------
  ! Get the mean position
  ! Used for subtracting mean drift from displacement
  ! ----------------------------------------------------------------------

  subroutine get_mean_pos(part, mean_pos)

    type(particle_t), intent(in) :: part
    real, dimension(3), intent(out) :: mean_pos
    integer :: nparticles, i

    nparticles = part%nparticles

    mean_pos = 0
    do i = 1, nparticles
       mean_pos = mean_pos + part%x(:,i)
    end do
    mean_pos = mean_pos / nparticles

  end subroutine get_mean_pos

  ! ----------------------------------------------------------------------
  ! Get symmetric tensor of self-displacement variance
  ! Used for getting MSD
  ! IMPORTANT: assumes unwrapped coords
  ! ----------------------------------------------------------------------

  subroutine get_displ_var(from, to, group, var, subtr_mean, igroup)

    type(particle_t), intent(in) :: from,to
    real, dimension(6), intent(out) :: var
    type(group_t), intent(in), target :: group
    logical, intent(in) :: subtr_mean
    integer, intent(in), optional :: igroup
    integer :: nparticles,i,groupsize,jgroup,groupbit
    real, dimension(3) :: mean_pos_from, mean_pos_to, mean_disp
    real :: delx,dely,delz
    integer, dimension(:), pointer :: mask => null()

    mask => group%mask

    ! check for consistency

    nparticles = from%nparticles
    if (to%nparticles /= nparticles) then
       write(*, '(a)') 'cannot get displ. var. between snaps with different nparticles'
       stop
    end if

    ! check if group exists (if present)
    ! exit with error if not
    ! if no group specified, group is all (group 1)

    if (.not. present(igroup)) then
       jgroup = 1
    else
       jgroup = igroup
    end if
    if (jgroup > group%ngroup) then
       write(*, '(a)') 'group ', jgroup, ' not found'
       stop
    end if
    if (jgroup < 1) then
       write(*, '(a)') 'group ', jgroup, ' not found'
       stop
    end if
    groupbit = group%bitmask(jgroup)

    ! loop over particles
    ! subtract mean motion if desired

    call get_mean_pos(from, mean_pos_from)
    call get_mean_pos(to, mean_pos_to)
    mean_disp = mean_pos_to - mean_pos_from

    var(:) = 0.0
    do i = 1, nparticles
       if (iand(mask(i),groupbit) /= 0) then
          delx = to%x(1,i) - from%x(1,i)
          dely = to%x(2,i) - from%x(2,i)
          delz = to%x(3,i) - from%x(3,i)
          if (subtr_mean) then
             delx = delx - mean_disp(1)
             dely = dely - mean_disp(2)
             delz = delz - mean_disp(3)
          end if
          var(1) = var(1) + delx*delx
          var(2) = var(2) + dely*dely
          var(3) = var(3) + delz*delz
          var(4) = var(4) + dely*delz
          var(5) = var(5) + delx*delz
          var(6) = var(6) + delx*dely
       end if
    end do

    ! normalize by size of group

    groupsize = group%count(jgroup)

    if (groupsize > 0) var(:) = var(:) / groupsize

  end subroutine get_displ_var

  ! ----------------------------------------------------------------------
  ! Get symmetric tensor of self-displacement 4th moment
  ! Used for getting non-gaussian parameter
  ! IMPORTANT: assumes unwrapped coords
  ! ----------------------------------------------------------------------

  subroutine get_displ_var_sq(from, to, group, var_sq, subtr_mean, igroup)

    type(particle_t), intent(in) :: from,to
    real, dimension(6), intent(out) :: var_sq
    type(group_t), intent(in), target :: group
    logical, intent(in) :: subtr_mean
    integer, intent(in), optional :: igroup
    integer :: nparticles,i,groupsize,jgroup,groupbit
    real, dimension(3) :: mean_pos_from, mean_pos_to, mean_disp
    real :: delx,dely,delz
    integer, dimension(:), pointer :: mask => null()

    mask => group%mask

    ! check for consistency

    nparticles = from%nparticles
    if (to%nparticles /= nparticles) then
       write(*, '(a)') 'cannot get displ. var. between snaps with different nparticles'
       stop
    end if

    ! check if group exists (if present)
    ! exit with error if not
    ! if no group specified, group is all (group 1)

    if (.not. present(igroup)) then
       jgroup = 1
    else
       jgroup = igroup
    end if
    if (jgroup > group%ngroup) then
       write(*, '(a)') 'group ', jgroup, ' not found'
       stop
    end if
    if (jgroup < 1) then
       write(*, '(a)') 'group ', jgroup, ' not found'
       stop
    end if
    groupbit = group%bitmask(jgroup)

    ! loop over particles
    ! subtract mean motion if desired

    call get_mean_pos(from, mean_pos_from)
    call get_mean_pos(to, mean_pos_to)
    mean_disp = mean_pos_to - mean_pos_from

    var_sq(:) = 0.0
    do i = 1, nparticles
       if (iand(mask(i),groupbit) /= 0) then
          delx = to%x(1,i) - from%x(1,i)
          dely = to%x(2,i) - from%x(2,i)
          delz = to%x(3,i) - from%x(3,i)
          if (subtr_mean) then
             delx = delx - mean_disp(1)
             dely = dely - mean_disp(2)
             delz = delz - mean_disp(3)
          end if
          var_sq(1) = var_sq(1) + delx*delx*delx*delx
          var_sq(2) = var_sq(2) + dely*dely*dely*dely
          var_sq(3) = var_sq(3) + delz*delz*delz*delz
          var_sq(4) = var_sq(4) + dely*delz*dely*delz
          var_sq(5) = var_sq(5) + delx*delz*delx*delz
          var_sq(6) = var_sq(6) + delx*dely*delx*dely
       end if
    end do

    ! normalize by size of group

    groupsize = group%count(jgroup)

    if (groupsize > 0) var_sq(:) = var_sq(:) / groupsize

  end subroutine get_displ_var_sq

  ! ----------------------------------------------------------------------
  ! Get symmetric tensor of self-displacement variance
  ! Used for getting cage-relative(cr) MSD
  ! IMPORTANT: assumes unwrapped coords
  ! ----------------------------------------------------------------------

  subroutine get_displ_var_cr(from, to, group, cells, yes_delta, r, var_cr, igroup)

    type(particle_t), intent(in) :: from,to
    real, dimension(6), intent(out) :: var_cr
    type(group_t), intent(in), target :: group
    integer, intent(in), optional :: igroup
    integer :: groupsize,jgroup,groupbit
    real, dimension(3) :: mean_pos_from_neighbor, mean_pos_to_neighbor
    real :: delx,dely,delz,delxx,delyy,delzz
    integer, dimension(:), pointer :: mask => null()

    type(cells_t), intent(inout) :: cells
    logical, intent(in) :: yes_delta
    real, intent(in) :: r
    integer :: c,d,i,j,neighcell,itype,jtype,nparticles,ntypes,nc
    logical :: typeradius_flag
    real :: rcutsq,rsq,xtmp,ytmp,ztmp
    real, dimension(:,:), allocatable :: rsqmat

    typeradius_flag = cells%particle%typeradius_flag
    ntypes = cells%particle%ntypes

    ! default rcutsq for yes_delta is false
    rcutsq = r*r

    if (yes_delta .and. typeradius_flag) then
       allocate(rsqmat(ntypes,ntypes))
       do i = 1, ntypes
          do j = 1, ntypes
             rsqmat(i,j) = cells%particle%typeradius(i) &
                         + cells%particle%typeradius(j) &
                         + r
          end do
       end do
       rsqmat(:,:) = rsqmat(:,:)*rsqmat(:,:)
    end if

    mask => group%mask

    ! check for consistency

    nparticles = from%nparticles
    if (to%nparticles /= nparticles) then
       write(*, '(a)') 'cannot get displ. var. between snaps with different nparticles'
       stop
    end if

    ! check if group exists (if present)
    ! exit with error if not
    ! if no group specified, group is all (group 1)

    if (.not. present(igroup)) then
       jgroup = 1
    else
       jgroup = igroup
    end if
    if (jgroup > group%ngroup) then
       write(*, '(a)') 'group ', jgroup, ' not found'
       stop
    end if
    if (jgroup < 1) then
       write(*, '(a)') 'group ', jgroup, ' not found'
       stop
    end if
    groupbit = group%bitmask(jgroup)

    ! loop over particles
    ! subtract nearest-neighbor motion

    var_cr(:) = 0.0
    do i = 1, nparticles
       if (iand(mask(i),groupbit) /= 0) then
          delx = to%x(1,i) - from%x(1,i)
          dely = to%x(2,i) - from%x(2,i)
          delz = to%x(3,i) - from%x(3,i)

          xtmp = to%x(1,i)
          ytmp = to%x(2,i)
          ztmp = to%x(3,i)

          nc = 0
          mean_pos_to_neighbor(:) = 0.0
          mean_pos_from_neighbor(:) = 0.0

          ! keeping logical evaluation out of inner loop
          ! first case: uniform cutoff

          option_switch: if (.not. yes_delta) then

             if_cells_uniform: if (cells%cells_exist) then

                cell_loop_uniform: do c = 1, cells%ncells

                   j = cells%cellhead(c)
                   do while (j > 0)
                      if (j==i) cycle
                      delxx = cells%particle%x(1,j) - xtmp
                      delyy = cells%particle%x(2,j) - ytmp
                      delzz = cells%particle%x(3,j) - ztmp
                      call cells%domain%minimum_image(delxx, delyy, delzz)
                      rsq = delxx*delxx + delyy*delyy + delzz*delzz
                      if (rsq < rcutsq) then
                         mean_pos_to_neighbor(:) = mean_pos_to_neighbor(:) + to%x(:,j)
                         mean_pos_from_neighbor(:) = mean_pos_from_neighbor(:) + from%x(:,j)
                         nc = nc + 1
                      end if
                      j = cells%next(j)
                   end do

                   ! loop over distinct cells
                   ! should I make this part of previous loop?

                   j = cells%cellhead(c)
                   do d = 1, cells%nstencil
                      neighcell = cells%stencil_neighbor(c,d)
                      if (neighcell == 0) cycle
                      j = cells%cellhead(neighcell)
                      do while (j > 0)
                         if (j==i) cycle
                         delxx = cells%particle%x(1,j) - xtmp
                         delyy = cells%particle%x(2,j) - ytmp
                         delzz = cells%particle%x(3,j) - ztmp
                         call cells%domain%minimum_image(delxx, delyy, delzz)
                         rsq = delxx*delxx + delyy*delyy + delzz*delzz
                         if (rsq < rcutsq) then
                            mean_pos_to_neighbor(:) = mean_pos_to_neighbor(:) + to%x(:,j)
                            mean_pos_from_neighbor(:) = mean_pos_from_neighbor(:) + from%x(:,j)
                            nc = nc + 1
                         end if
                         j = cells%next(j)
                      end do
                   end do
                         
                end do cell_loop_uniform

             else ! not using cells, so loop over all pairs

                do j = 1, nparticles
                   if (j==i) cycle
                   delxx = cells%particle%x(1,j) - xtmp
                   delyy = cells%particle%x(2,j) - ytmp
                   delzz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delxx, delyy, delzz)
                   rsq = delxx*delxx + delyy*delyy + delzz*delzz
                   if (rsq < rcutsq) then
                      mean_pos_to_neighbor(:) = mean_pos_to_neighbor(:) + to%x(:,j)
                      mean_pos_from_neighbor(:) = mean_pos_from_neighbor(:) + from%x(:,j)
                      nc = nc + 1
                   end if
                end do

             end if if_cells_uniform

          else if (typeradius_flag) then ! can use matrix cutoffs

             if_cells_typeradius: if (cells%cells_exist) then

                cell_loop_typeradius: do c = 1, cells%ncells

                   ! below, all in the ith cell
                   ! loop over all particles i,j in this cell c
                   ! cellhead(c) is the first i particle
                   ! the first (distinct) j is the next of i

                   j = cells%cellhead(c)
                   itype = cells%particle%type(i)
                   do while (j > 0)
                      if (j==i) cycle
                      jtype = cells%particle%type(j)
                      delxx = cells%particle%x(1,j) - xtmp
                      delyy = cells%particle%x(2,j) - ytmp
                      delzz = cells%particle%x(3,j) - ztmp
                      call cells%domain%minimum_image(delxx, delyy, delzz)
                      rsq = delxx*delxx + delyy*delyy + delzz*delzz
                      rcutsq = rsqmat(itype,jtype)
                      if (rsq < rcutsq) then
                         mean_pos_to_neighbor(:) = mean_pos_to_neighbor(:) + to%x(:,j)
                         mean_pos_from_neighbor(:) = mean_pos_from_neighbor(:) + from%x(:,j)
                         nc = nc + 1
                      end if
                      j = cells%next(j)
                   end do

                   ! loop over distinct cells
                   ! should I make this part of previous loop?

                   j = cells%cellhead(c)
                   itype = cells%particle%type(i)
                   do d = 1, cells%nstencil
                      neighcell = cells%stencil_neighbor(c,d)
                      if (neighcell == 0) cycle
                      j = cells%cellhead(neighcell)
                      do while (j > 0)
                         if (j==i) cycle
                         jtype = cells%particle%type(j)
                         delxx = cells%particle%x(1,j) - xtmp
                         delyy = cells%particle%x(2,j) - ytmp
                         delzz = cells%particle%x(3,j) - ztmp
                         call cells%domain%minimum_image(delxx, delyy, delzz)
                         rsq = delxx*delxx + delyy*delyy + delzz*delzz
                         rcutsq = rsqmat(itype,jtype)
                         if (rsq < rcutsq) then
                            mean_pos_to_neighbor(:) = mean_pos_to_neighbor(:) + to%x(:,j)
                            mean_pos_from_neighbor(:) = mean_pos_from_neighbor(:) + from%x(:,j)
                            nc = nc + 1
                         end if
                         j = cells%next(j)
                      end do
                   end do
                         
                end do cell_loop_typeradius

             else ! not using cells, so loop over all pairs

                itype = cells%particle%type(i)
                do j = 1, nparticles
                   if (j==i) cycle
                   jtype = cells%particle%type(j)
                   delxx = cells%particle%x(1,j) - xtmp
                   delyy = cells%particle%x(2,j) - ytmp
                   delzz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delxx, delyy, delzz)
                   rsq = delxx*delxx + delyy*delyy + delzz*delzz
                   rcutsq = rsqmat(itype,jtype)
                   if (rsq < rcutsq) then
                      mean_pos_to_neighbor(:) = mean_pos_to_neighbor(:) + to%x(:,j)
                      mean_pos_from_neighbor(:) = mean_pos_from_neighbor(:) + from%x(:,j)
                      nc = nc + 1
                   end if
                end do

             end if if_cells_typeradius

          else ! have to use per-particle cutoff based on delta

             if_cells_perparticle: if (cells%cells_exist) then

                cell_loop_parparticle: do c = 1, cells%ncells

                   j = cells%cellhead(c)
                   do while (j > 0)
                      if (j==i) cycle
                      delxx = cells%particle%x(1,j) - xtmp
                      delyy = cells%particle%x(2,j) - ytmp
                      delzz = cells%particle%x(3,j) - ztmp
                      call cells%domain%minimum_image(delxx, delyy, delzz)
                      rsq = delxx*delxx + delyy*delyy + delzz*delzz
                      rcutsq = cells%particle%radius(i)
                      rcutsq = rcutsq + cells%particle%radius(j)
                      rcutsq = rcutsq + r
                      rcutsq = rcutsq*rcutsq
                      if (rsq < rcutsq) then
                         mean_pos_to_neighbor(:) = mean_pos_to_neighbor(:) + to%x(:,j)
                         mean_pos_from_neighbor(:) = mean_pos_from_neighbor(:) + from%x(:,j)
                         nc = nc + 1
                      end if
                      j = cells%next(j)
                   end do

                   ! loop over distinct cells
                   ! should I make this part of previous loop?

                   j = cells%cellhead(c)
                   do d = 1, cells%nstencil
                      neighcell = cells%stencil_neighbor(c,d)
                      if (neighcell == 0) cycle
                      j = cells%cellhead(neighcell)
                      do while (j > 0)
                         if (j==i) cycle
                         delxx = cells%particle%x(1,j) - xtmp
                         delyy = cells%particle%x(2,j) - ytmp
                         delzz = cells%particle%x(3,j) - ztmp
                         call cells%domain%minimum_image(delxx, delyy, delzz)
                         rsq = delxx*delxx + delyy*delyy + delzz*delzz
                         rcutsq = cells%particle%radius(i)
                         rcutsq = rcutsq + cells%particle%radius(j)
                         rcutsq = rcutsq + r
                         rcutsq = rcutsq*rcutsq
                         if (rsq < rcutsq) then
                            mean_pos_to_neighbor(:) = mean_pos_to_neighbor(:) + to%x(:,j)
                            mean_pos_from_neighbor(:) = mean_pos_from_neighbor(:) + from%x(:,j)
                            nc = nc + 1
                         end if
                         j = cells%next(j)
                      end do
                   end do
                         
                end do cell_loop_parparticle

             else ! not using cells, so loop over all pairs

                do j = 1, nparticles
                   if (j==i) cycle
                   delxx = cells%particle%x(1,j) - xtmp
                   delyy = cells%particle%x(2,j) - ytmp
                   delzz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delxx, delyy, delzz)
                   rsq = delxx*delxx + delyy*delyy + delzz*delzz
                   rcutsq = cells%particle%radius(i)
                   rcutsq = rcutsq + cells%particle%radius(j)
                   rcutsq = rcutsq + r
                   rcutsq = rcutsq*rcutsq
                   if (rsq < rcutsq) then
                      mean_pos_to_neighbor(:) = mean_pos_to_neighbor(:) + to%x(:,j)
                      mean_pos_from_neighbor(:) = mean_pos_from_neighbor(:) + from%x(:,j)
                      nc = nc + 1
                   end if
                end do

             end if if_cells_perparticle

          end if option_switch

          delx = delx-(mean_pos_from_neighbor(1)-mean_pos_to_neighbor(1))/nc
          dely = dely-(mean_pos_from_neighbor(2)-mean_pos_to_neighbor(2))/nc
          delz = delz-(mean_pos_from_neighbor(3)-mean_pos_to_neighbor(3))/nc

          var_cr(1) = var_cr(1) + delx*delx
          var_cr(2) = var_cr(2) + dely*dely
          var_cr(3) = var_cr(3) + delz*delz
          var_cr(4) = var_cr(4) + dely*delz
          var_cr(5) = var_cr(5) + delx*delz
          var_cr(6) = var_cr(6) + delx*dely
       end if
    end do

    ! normalize by size of group

    groupsize = group%count(jgroup)

    if (groupsize > 0) var_cr(:) = var_cr(:) / groupsize

  end subroutine get_displ_var_cr


! ----------------------------------------------------------------------
! Calculate bound pairs and the state space (0 = unbound, 1 = bound) according to the set cutoff distance, at a given snapshot
! ----------------------------------------------------------------------

  subroutine eval_binding(cells, group, agroup, bgroup, cut_dist, state_space)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    real, intent(in) :: cut_dist
    integer, intent(in) :: agroup,bgroup
    integer, dimension(:), intent(inout) :: state_space
    real, dimension(:,:), allocatable :: rsqmat
    real :: rsq,rcutsq,rminsq,normfac,rgrid,inv_rgrid,delx,dely,delz,xtmp,ytmp,ztmp
    real :: normdum_a,normdum_b,normdum_c
    integer :: i,j,c,d,n,groupbit_a,groupbit_b,othergroup,neighcell,na,nb
    integer :: itype,jtype,count_bind,nparticles,ntypes
    integer, dimension(:), pointer :: mask => null()
    logical :: groups_identical

    mask => group%mask

    nparticles = cells%particle%nparticles
    ntypes = cells%particle%ntypes

  ! If group ids are not equal, cannot have a particle in both groups.

    na = group%count(agroup)
    nb = group%count(bgroup)
    if (na < 1) then
       write(*, '(a)') 'group A in binding calculations has bad count'
       stop
    end if
    if (nb < 1) then
       write(*, '(a)') 'group B in binding calculations has bad count'
       stop
    end if

    groups_identical = .false.
    if (agroup == bgroup) groups_identical = .true.

    groupbit_a = group%bitmask(agroup)
    groupbit_b = group%bitmask(bgroup)

    allocate(rsqmat(ntypes,ntypes))
    do i = 1, ntypes
       do j = 1, ntypes
          rsqmat(i,j) = ( cells%particle%typeradius(i) &
                      + cells%particle%typeradius(j)) * (1 + cut_dist)
       end do
    end do
    rsqmat(:,:) = rsqmat(:,:)*rsqmat(:,:)

    state_space(:) = 0
    count_bind = 1

  ! check if particle coordinates need wrapping, and wrap if they do
  ! can do this here, because no reason to be in lamda coords

    if (.not. cells%domain%particles_inside()) then
       call cells%domain%remap()
    end if

    if_cells: if (cells%cells_exist) then
       cell_loop_typeradius: do c = 1, cells%ncells
      ! below, all in the ith cell
      ! loop over all particles i,j in this cell c
      ! cellhead(c) is the first i particle
      ! the first (distinct) j is the next of i
          i = cells%cellhead(c)
          self_loop: do while (i > 0)
             itype = cells%particle%type(i)
             if (iand(mask(i),groupbit_a) /= 0) then
                othergroup = groupbit_b
             else if (iand(mask(i),groupbit_b) /= 0) then
                othergroup = groupbit_a
             else
                i = cells%next(i)
                cycle
             end if
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)
             j = cells%next(i)
             do while (j > 0)
                jtype = cells%particle%type(j)
                if (iand(mask(j),othergroup) == 0) then
                   j = cells%next(j)
                   cycle
                end if
                delx = cells%particle%x(1,j) - xtmp
                dely = cells%particle%x(2,j) - ytmp
                delz = cells%particle%x(3,j) - ztmp
                call cells%domain%minimum_image(delx, dely, delz)
                rsq = delx*delx + dely*dely + delz*delz
                rcutsq = rsqmat(itype,jtype)
                if (rsq < rcutsq) then
                   state_space(i) = 1
                   state_space(j) = 1
                   count_bind = count_bind + 1
                end if
                j = cells%next(j)
             end do
             i = cells%next(i)
          end do self_loop

           ! loop over distinct cells
           ! should I make this part of previous loop?

          i = cells%cellhead(c)
          distinct_loop: do while (i > 0)
             itype = cells%particle%type(i)
             if (iand(mask(i),groupbit_a) /= 0) then
                othergroup = groupbit_b
             else if (iand(mask(i),groupbit_b) /= 0) then
                othergroup = groupbit_a
             else
                i = cells%next(i)
                cycle
             end if
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)
             do d = 1, cells%nstencil
                neighcell = cells%stencil_neighbor(c,d)
                if (neighcell == 0) cycle
                j = cells%cellhead(neighcell)
                do while (j > 0)
                   jtype = cells%particle%type(j)
                   if (iand(mask(j),othergroup) == 0) then
                      j = cells%next(j)
                      cycle
                   end if
                   delx = cells%particle%x(1,j) - xtmp
                   dely = cells%particle%x(2,j) - ytmp
                   delz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delx, dely, delz)
                   rsq = delx*delx + dely*dely + delz*delz
                   rcutsq = rsqmat(itype,jtype)
                   if (rsq < rcutsq) then
                      state_space(i) = 1
                      state_space(j) = 1
                      count_bind = count_bind + 1
                   end if
                   j = cells%next(j)
                end do
             end do
             i = cells%next(i)
          end do distinct_loop
       end do cell_loop_typeradius

    else ! not using cells, so loop over all pairs

       do i = 1, nparticles-1
          itype = cells%particle%type(i)
          if (iand(mask(i),groupbit_a) /= 0) then
             othergroup = groupbit_b
          else if (iand(mask(i),groupbit_b) /= 0) then
             othergroup = groupbit_a
          else
             cycle
          end if
          xtmp = cells%particle%x(1,i)
          ytmp = cells%particle%x(2,i)
          ztmp = cells%particle%x(3,i)
          do j = i+1, nparticles
             jtype = cells%particle%type(j)
             if (iand(mask(j),othergroup) == 0) cycle
             delx = cells%particle%x(1,j) - xtmp
             dely = cells%particle%x(2,j) - ytmp
             delz = cells%particle%x(3,j) - ztmp
             call cells%domain%minimum_image(delx, dely, delz)
             rsq = delx*delx + dely*dely + delz*delz
             rcutsq = rsqmat(itype,jtype)
             if (rsq < rcutsq) then
                state_space(i) = 1
                state_space(j) = 1
                count_bind = count_bind + 1
             end if
          end do
       end do

    end if if_cells

  end subroutine eval_binding


! ----------------------------------------------------------------------
! Calculate bound pairs and track their separation over time
! ----------------------------------------------------------------------

  subroutine eval_pairdisttime(cells, group, agroup, bgroup, cut_dist, pair_sep_time, pair_sep_id, num_to_calc, timesnap)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    real, intent(in) :: cut_dist
    integer, intent(in) :: agroup,bgroup,num_to_calc, timesnap
    real, dimension(:,:), intent(inout) :: pair_sep_time, pair_sep_id
    real :: rsq,rcutsq,rminsq,normfac,rgrid,inv_rgrid,delx,dely,delz,xtmp,ytmp,ztmp
    real :: normdum_a,normdum_b,normdum_c
    integer :: i,j,c,d,n,groupbit_a,groupbit_b,othergroup,neighcell,na,nb,vali,valj
    integer :: itype,jtype,count_bind,nparticles,ntypes
    integer, dimension(:), pointer :: mask => null()
    real, dimension(:,:), allocatable :: rsqmat
    logical :: groups_identical

    mask => group%mask

    nparticles = cells%particle%nparticles
    ntypes = cells%particle%ntypes

  ! If group ids are not equal, cannot have a particle in both groups.

    na = group%count(agroup)
    nb = group%count(bgroup)
    if (na < 1) then
       write(*, '(a)') 'group A in binding calculations has bad count'
       stop
    end if
    if (nb < 1) then
       write(*, '(a)') 'group B in binding calculations has bad count'
       stop
    end if

    groups_identical = .false.
    if (agroup == bgroup) groups_identical = .true.

    groupbit_a = group%bitmask(agroup)
    groupbit_b = group%bitmask(bgroup)

    allocate(rsqmat(ntypes,ntypes))
    do i = 1, ntypes
       do j = 1, ntypes
          rsqmat(i,j) = ( cells%particle%typeradius(i) &
                      + cells%particle%typeradius(j)) * (1 + cut_dist)
       end do
    end do
    rsqmat(:,:) = rsqmat(:,:)*rsqmat(:,:)
    count_bind = 1

  ! check if particle coordinates need wrapping, and wrap if they do
  ! can do this here, because no reason to be in lamda coords

    if (.not. cells%domain%particles_inside()) then
       call cells%domain%remap()
    end if
    if (pair_sep_id(1,1) == 0) then
       if_cells: if (cells%cells_exist) then
          cell_loop_typeradius: do c = 1, cells%ncells
         ! below, all in the ith cell
         ! loop over all particles i,j in this cell c
         ! cellhead(c) is the first i particle
         ! the first (distinct) j is the next of i
             i = cells%cellhead(c)
             self_loop: do while (i > 0)
                itype = cells%particle%type(i)
                if (iand(mask(i),groupbit_a) /= 0) then
                   othergroup = groupbit_b
                else if (iand(mask(i),groupbit_b) /= 0) then
                   othergroup = groupbit_a
                else
                   i = cells%next(i)
                   cycle
                end if
                xtmp = cells%particle%x(1,i)
                ytmp = cells%particle%x(2,i)
                ztmp = cells%particle%x(3,i)
                j = cells%next(i)
                do while (j > 0)
                   jtype = cells%particle%type(j)
                   if (iand(mask(j),othergroup) == 0) then
                      j = cells%next(j)
                      cycle
                   end if
                   delx = cells%particle%x(1,j) - xtmp
                   dely = cells%particle%x(2,j) - ytmp
                   delz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delx, dely, delz)
                   rsq = delx*delx + dely*dely + delz*delz
                   rcutsq = rsqmat(itype,jtype)
                   if (rsq < rcutsq) then
                      pair_sep_time(timesnap, count_bind) = sqrt(rsq)
                      pair_sep_id(1, count_bind) = i
                      pair_sep_id(2, count_bind) = j
                      count_bind = count_bind + 1
                      if (count_bind > num_to_calc) exit if_cells
                   end if
                   j = cells%next(j)
                end do
                i = cells%next(i)
             end do self_loop
 
              ! loop over distinct cells
              ! should I make this part of previous loop?

             i = cells%cellhead(c)
             distinct_loop: do while (i > 0)
                itype = cells%particle%type(i)
                if (iand(mask(i),groupbit_a) /= 0) then
                   othergroup = groupbit_b
                else if (iand(mask(i),groupbit_b) /= 0) then
                   othergroup = groupbit_a
                else
                   i = cells%next(i)
                   cycle
                end if
                xtmp = cells%particle%x(1,i)
                ytmp = cells%particle%x(2,i)
                ztmp = cells%particle%x(3,i)
                do d = 1, cells%nstencil
                   neighcell = cells%stencil_neighbor(c,d)
                   if (neighcell == 0) cycle
                   j = cells%cellhead(neighcell)
                   do while (j > 0)
                      jtype = cells%particle%type(j)
                      if (iand(mask(j),othergroup) == 0) then
                         j = cells%next(j)
                         cycle
                      end if
                      delx = cells%particle%x(1,j) - xtmp
                      dely = cells%particle%x(2,j) - ytmp
                      delz = cells%particle%x(3,j) - ztmp
                      call cells%domain%minimum_image(delx, dely, delz)
                      rsq = delx*delx + dely*dely + delz*delz
                      rcutsq = rsqmat(itype,jtype)
                      if (rsq < rcutsq) then
                         pair_sep_time(timesnap, count_bind) = sqrt(rsq)
                         pair_sep_id(1, count_bind) = i
                         pair_sep_id(2, count_bind) = j
                         count_bind = count_bind + 1
                         if (count_bind > num_to_calc) exit if_cells
                      end if
                      j = cells%next(j)
                   end do
                end do
                i = cells%next(i)
             end do distinct_loop
          end do cell_loop_typeradius

       else ! not using cells, so loop over all pairs

          do i = 1, nparticles-1
             itype = cells%particle%type(i)
             if (iand(mask(i),groupbit_a) /= 0) then
                othergroup = groupbit_b
             else if (iand(mask(i),groupbit_b) /= 0) then
                othergroup = groupbit_a
             else
                cycle
             end if
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)
             do j = i+1, nparticles
                jtype = cells%particle%type(j)
                if (iand(mask(j),othergroup) == 0) cycle
                delx = cells%particle%x(1,j) - xtmp
                dely = cells%particle%x(2,j) - ytmp
                delz = cells%particle%x(3,j) - ztmp
                call cells%domain%minimum_image(delx, dely, delz)
                rsq = delx*delx + dely*dely + delz*delz
                rcutsq = rsqmat(itype,jtype)
                if (rsq < rcutsq) then
                   pair_sep_time(timesnap, count_bind) = sqrt(rsq)
                   pair_sep_id(1, count_bind) = i
                   pair_sep_id(2, count_bind) = j
                   count_bind = count_bind + 1
                   if (count_bind > num_to_calc) exit if_cells
                end if
             end do
          end do
 
       end if if_cells
    else
       do i = 1, num_to_calc
          vali = pair_sep_id(1,i)
          valj = pair_sep_id(2,i)
          xtmp = cells%particle%x(1,vali)
          ytmp = cells%particle%x(2,vali)
          ztmp = cells%particle%x(3,vali)
          delx = cells%particle%x(1,valj) - xtmp
          dely = cells%particle%x(2,valj) - ytmp
          delz = cells%particle%x(3,valj) - ztmp
          call cells%domain%minimum_image(delx, dely, delz)
          rsq = delx*delx + dely*dely + delz*delz
          pair_sep_time(timesnap, i) = sqrt(rsq)
       end do
    end if
  end subroutine eval_pairdisttime

! ----------------------------------------------------------------------
! Calculate bound pairs and track their separation over time
! ----------------------------------------------------------------------

  subroutine eval_pairdisttime_search(cells, group, agroup, bgroup, cut_dist, timesnap, pair_sep_time, num_to_calc)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    real, intent(in) :: cut_dist
    integer, intent(in) :: agroup,bgroup,timesnap,num_to_calc
    real, dimension(:,:), intent(inout) :: pair_sep_time
    real :: rsq,rcutsq,rminsq,normfac,rgrid,inv_rgrid,delx,dely,delz,xtmp,ytmp,ztmp
    real :: normdum_a,normdum_b,normdum_c
    integer :: i,j,c,d,n,groupbit_a,groupbit_b,othergroup,neighcell,na,nb,vali,valj
    integer :: itype,jtype,count_bind,nparticles,ntypes
    integer, dimension(:), pointer :: mask => null()
    real, dimension(:,:), allocatable :: rsqmat
    logical :: groups_identical

    mask => group%mask

    nparticles = cells%particle%nparticles
    ntypes = cells%particle%ntypes

  ! If group ids are not equal, cannot have a particle in both groups.

    na = group%count(agroup)
    nb = group%count(bgroup)
    if (na < 1) then
       write(*, '(a)') 'group A in binding calculations has bad count'
       stop
    end if
    if (nb < 1) then
       write(*, '(a)') 'group B in binding calculations has bad count'
       stop
    end if

    groups_identical = .false.
    if (agroup == bgroup) groups_identical = .true.

    groupbit_a = group%bitmask(agroup)
    groupbit_b = group%bitmask(bgroup)

    allocate(rsqmat(ntypes,ntypes))
    do i = 1, ntypes
       do j = 1, ntypes
          rsqmat(i,j) = ( cells%particle%typeradius(i) &
                      + cells%particle%typeradius(j)) * (1 + cut_dist)
       end do
    end do
    rsqmat(:,:) = rsqmat(:,:)*rsqmat(:,:)
    count_bind = 1

  ! check if particle coordinates need wrapping, and wrap if they do
  ! can do this here, because no reason to be in lamda coords

    if (.not. cells%domain%particles_inside()) then
       call cells%domain%remap()
    end if
    if (0 == 0) then
       if_cells: if (cells%cells_exist) then
          cell_loop_typeradius: do c = 1, cells%ncells
         ! below, all in the ith cell
         ! loop over all particles i,j in this cell c
         ! cellhead(c) is the first i particle
         ! the first (distinct) j is the next of i
             i = cells%cellhead(c)
             self_loop: do while (i > 0)
                itype = cells%particle%type(i)
                if (iand(mask(i),groupbit_a) /= 0) then
                   othergroup = groupbit_b
                else if (iand(mask(i),groupbit_b) /= 0) then
                   othergroup = groupbit_a
                else
                   i = cells%next(i)
                   cycle
                end if
                if (i > num_to_calc) cycle
                xtmp = cells%particle%x(1,i)
                ytmp = cells%particle%x(2,i)
                ztmp = cells%particle%x(3,i)
                j = cells%next(i)
                do while (j > 0)
                   jtype = cells%particle%type(j)
                   if (iand(mask(j),othergroup) == 0) then
                      j = cells%next(j)
                      cycle
                   end if
                   delx = cells%particle%x(1,j) - xtmp
                   dely = cells%particle%x(2,j) - ytmp
                   delz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delx, dely, delz)
                   rsq = delx*delx + dely*dely + delz*delz
                   rcutsq = rsqmat(itype,jtype)
                   if (rsq < rcutsq) then
                      pair_sep_time(timesnap, i) = 1
                      if (j <= num_to_calc) then
                         pair_sep_time(timesnap, j) = 1
                      end if
                   end if
                   j = cells%next(j)
                end do
                i = cells%next(i)
             end do self_loop
 
              ! loop over distinct cells
              ! should I make this part of previous loop?

             i = cells%cellhead(c)
             distinct_loop: do while (i > 0)
                itype = cells%particle%type(i)
                if (iand(mask(i),groupbit_a) /= 0) then
                   othergroup = groupbit_b
                else if (iand(mask(i),groupbit_b) /= 0) then
                   othergroup = groupbit_a
                else
                   i = cells%next(i)
                   cycle
                end if
                if (i > num_to_calc) cycle
                xtmp = cells%particle%x(1,i)
                ytmp = cells%particle%x(2,i)
                ztmp = cells%particle%x(3,i)
                do d = 1, cells%nstencil
                   neighcell = cells%stencil_neighbor(c,d)
                   if (neighcell == 0) cycle
                   j = cells%cellhead(neighcell)
                   do while (j > 0)
                      jtype = cells%particle%type(j)
                      if (iand(mask(j),othergroup) == 0) then
                         j = cells%next(j)
                         cycle
                      end if
                      delx = cells%particle%x(1,j) - xtmp
                      dely = cells%particle%x(2,j) - ytmp
                      delz = cells%particle%x(3,j) - ztmp
                      call cells%domain%minimum_image(delx, dely, delz)
                      rsq = delx*delx + dely*dely + delz*delz
                      rcutsq = rsqmat(itype,jtype)
                      if (rsq < rcutsq) then
                         pair_sep_time(timesnap, i) = 1
                         if (j <= num_to_calc) then
                            pair_sep_time(timesnap, j) = 1
                         end if
                      end if
                      j = cells%next(j)
                   end do
                end do
                i = cells%next(i)
             end do distinct_loop
          end do cell_loop_typeradius

       else ! not using cells, so loop over all pairs

          do i = 1, nparticles-1
             itype = cells%particle%type(i)
             if (iand(mask(i),groupbit_a) /= 0) then
                othergroup = groupbit_b
             else if (iand(mask(i),groupbit_b) /= 0) then
                othergroup = groupbit_a
             else
                cycle
             end if
             if (i > num_to_calc) cycle
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)
             do j = i+1, nparticles
                jtype = cells%particle%type(j)
                if (iand(mask(j),othergroup) == 0) cycle
                delx = cells%particle%x(1,j) - xtmp
                dely = cells%particle%x(2,j) - ytmp
                delz = cells%particle%x(3,j) - ztmp
                call cells%domain%minimum_image(delx, dely, delz)
                rsq = delx*delx + dely*dely + delz*delz
                rcutsq = rsqmat(itype,jtype)
                if (rsq < rcutsq) then
                   pair_sep_time(timesnap, i) = 1
                   if (j <= num_to_calc) then
                      pair_sep_time(timesnap, j) = 1
                   end if
                end if
             end do
          end do
 
       end if if_cells
    end if
  end subroutine eval_pairdisttime_search


! ----------------------------------------------------------------------
! Calculate bound pairs and track their separation over time
! ----------------------------------------------------------------------

  subroutine eval_bindsearch_traj(pair_sep_time, timesnap, cut_dist, bindtimes_hist, searchtimes_hist, nsnaps)

    real, intent(in) :: cut_dist
    integer, intent(in) :: nsnaps, timesnap
    real, dimension(:), intent(inout) :: bindtimes_hist, searchtimes_hist
    real, dimension(:,:), intent(in) :: pair_sep_time
    real :: rsq,rcutsq,rminsq,normfac,rgrid,inv_rgrid,delx,dely,delz,xtmp,ytmp,ztmp
    real :: normdum_a,normdum_b,normdum_c
    integer :: i,j,c,d,n,groupbit_a,groupbit_b,othergroup,neighcell,na,nb,vali,valj, prev_bound_yn, curr_dur, curr_bound
    integer :: itype,jtype,count_bind,nparticles,ntypes
    real, dimension(:,:), allocatable :: rsqmat

    prev_bound_yn = 1
    curr_bound = 1
    curr_dur = 1
    do i = 1, nsnaps
       ! calculate if pair is currently bound or not
       if (pair_sep_time(i,timesnap) <= cut_dist) then
          curr_bound = 1
       else
          curr_bound = 0
       end if
       ! if current & previous state are the same, increment current event duration
       if (curr_bound == prev_bound_yn) then
          curr_dur = curr_dur + 1
       ! if current & previous state are different, add duration...
       else
          ! if bound, add to histogram
          if (prev_bound_yn == 1) then
             bindtimes_hist(curr_dur) = bindtimes_hist(curr_dur) + 1
             curr_dur = 1
          ! if unbound, add to histogram
          else
             searchtimes_hist(curr_dur) = searchtimes_hist(curr_dur) + 1
             curr_dur = 1
          end if
       end if
       prev_bound_yn = curr_bound
    end do

  end subroutine eval_bindsearch_traj

! ----------------------------------------------------------------------
! Calculate bound pairs and track their separation over time
! ----------------------------------------------------------------------

  subroutine eval_bindsearch_traj_strictloose(pair_sep_time, pair, cut_dist, cut_dist_big, bindtimes_hist, searchtimes_hist, nsnaps)
!       call eval_bindsearch_traj_strictloose(pair_sep_time, i,    rsep,     rsep_big,     bindtimes_hist, searchtimes_hist, nsnaps)

    real, intent(in) :: cut_dist, cut_dist_big
    integer, intent(in) :: nsnaps, pair
    real, dimension(:), intent(inout) :: bindtimes_hist, searchtimes_hist
    real, dimension(:,:), intent(in) :: pair_sep_time
    real :: rsq,rcutsq,rminsq,normfac,rgrid,inv_rgrid,delx,dely,delz,xtmp,ytmp,ztmp
    real :: normdum_a,normdum_b,normdum_c, cut_dist_temp
    integer :: i,j,c,d,n,groupbit_a,groupbit_b,othergroup,neighcell,na,nb,vali,valj, prev_bound_yn, curr_dur, curr_bound
    integer :: itype,jtype,count_bind,nparticles,ntypes
    real, dimension(:,:), allocatable :: rsqmat

    prev_bound_yn = 1
    curr_bound = 1
    curr_dur = 1
    cut_dist_temp = cut_dist_big
    do i = 1, nsnaps
       if (prev_bound_yn == 1) then
          cut_dist_temp = cut_dist_big
       else
          cut_dist_temp = cut_dist
       end if
       !print *, "pair ",pair," time ",i," cut_dist_temp ",cut_dist_temp
       ! calculate if pair is currently bound or not
       if (pair_sep_time(i,pair) <= cut_dist_temp) then
          curr_bound = 1
       else
          curr_bound = 0
       end if
       ! if current & previous state are the same, increment current event duration
       if (curr_bound == prev_bound_yn) then
          curr_dur = curr_dur + 1
       ! if current & previous state are different, add duration...
       else
          ! if bound, add to histogram
          if (prev_bound_yn == 1) then
             bindtimes_hist(curr_dur) = bindtimes_hist(curr_dur) + 1
             curr_dur = 1
          ! if unbound, add to histogram
          else
             searchtimes_hist(curr_dur) = searchtimes_hist(curr_dur) + 1
             curr_dur = 1
          end if
       end if
       prev_bound_yn = curr_bound
    end do

  end subroutine eval_bindsearch_traj_strictloose

! ----------------------------------------------------------------------
! Calculate number of bound "B" according to the state space, (0 = unbound, 1 = bound), and adds to existing num_bound
! ----------------------------------------------------------------------

  subroutine add_numbound(cells, group, bgroup, state_space, num_bound)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    integer, intent(in) :: bgroup
    integer, dimension(:), intent(inout) :: state_space
    integer, intent(inout) :: num_bound
    integer :: i,groupbit_b
    integer, dimension(:), pointer :: mask => null()

    mask => group%mask

    groupbit_b = group%bitmask(bgroup)
    do i = 1, size(state_space,1)
       if (iand(mask(i),groupbit_b) /= 0) then
          if (state_space(i) == 1) then
             num_bound = num_bound + 1
          end if
       end if
    end do

  end subroutine add_numbound


! ----------------------------------------------------------------------
! Calculate number of bound "B" according to the state space in time, (0 = unbound, 1 = bound), and adds to existing num_bound_time
! ----------------------------------------------------------------------

  subroutine add_numbound_time(cells, group, bgroup, state_space_time, num_bound_time)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    integer, intent(in) :: bgroup
    integer, dimension(:), intent(inout) :: state_space_time
    integer, intent(inout) :: num_bound_time
    integer :: i,groupbit_b
    integer, dimension(:), pointer :: mask => null()

    mask => group%mask

    groupbit_b = group%bitmask(bgroup)
    do i = 1, size(state_space_time,1)
       if (iand(mask(i),groupbit_b) /= 0) then
          if (state_space_time(i) == 1) then
             num_bound_time = num_bound_time + 1
          end if
       end if
    end do

  end subroutine add_numbound_time


! ----------------------------------------------------------------------
! Adds the values for bound pairs at the initial timestep to the temp_bound array with a duration of tinc
! ----------------------------------------------------------------------

  subroutine add_tempbound(cells, group, bgroup, &
                            temp_bound_alt, state_space, tinc)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    integer(int64), intent(in) :: tinc
    integer, dimension(:), intent(inout) :: state_space
    integer, intent(in) :: bgroup
    integer, dimension(:,:), intent(inout) :: temp_bound_alt
    integer :: i,groupbit_b
    integer, dimension(:), pointer :: mask => null()

    mask => group%mask
    groupbit_b = group%bitmask(bgroup)

    do i = 1, size(state_space,1)
       if (iand(mask(i),groupbit_b) /= 0) then
          if (state_space(i) == 1) then
             temp_bound_alt(i,1) = 1
             temp_bound_alt(i,2) = tinc
          end if
       end if
    end do


  end subroutine add_tempbound


! ----------------------------------------------------------------------
! Adds the values for unbound B at the initial timestep to the temp_unbound array with a duration of tinc
! ----------------------------------------------------------------------

  subroutine add_state(cells, group, bgroup, state_space, temp_unbound, tinc)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    integer(int64), intent(in) :: tinc
    integer, intent(in) :: bgroup
    integer, dimension(:), intent(inout) :: state_space
    integer, dimension(:,:), intent(inout) :: temp_unbound
    integer :: i,groupbit_b
    integer, dimension(:), pointer :: mask => null()

    mask => group%mask

    groupbit_b = group%bitmask(bgroup)
    do i = 1, size(state_space,1)
       if (iand(mask(i),groupbit_b) /= 0) then
          if (state_space(i) == 0) then
             temp_unbound(i,1) = 1
             temp_unbound(i,2) = tinc
          end if
       end if
    end do

  end subroutine add_state


! ----------------------------------------------------------------------
! Adds duration to pairs that are still bound. Eliminates pairs from temp_bound which have become unbound & adds final duration to histogram of binding times
! ----------------------------------------------------------------------

  subroutine elim_tempbound(tinc, &
                            state_space_time, cut_time, &
                            cells, group, bgroup, state_space, &
                            temp_bound_alt, bindtimes_hist_alt)

    integer, dimension(:,:), intent(inout) :: temp_bound_alt
    integer, dimension(:), intent(inout) :: state_space_time
    integer, dimension(:), intent(inout) :: state_space, bindtimes_hist_alt
    real, intent(in) :: cut_time
    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    integer, dimension(1,2) :: bound_pair_temp, temp_zero
    integer(int64), intent(in) :: tinc
    integer, intent(in) :: bgroup
    integer :: i,j,groupbit_b, bindtimeindex, stillbound_yn
    integer, dimension(:), pointer :: mask => null()
    logical, dimension(2) :: array_match_one, array_match_two, array_match_three
    temp_zero(:,:) = 0

    mask => group%mask

    groupbit_b = group%bitmask(bgroup)
    do i = 1, size(state_space,1) ! for all particles
       if (iand(mask(i),groupbit_b) /= 0) then ! if in group b
          ! store particle as bound for at least cut_time
          if (temp_bound_alt(i,2)>=(cut_time*tinc)) then
             state_space_time(i) = 1
          else
             state_space_time(i) = 0
          end if
          ! check if "unbound" in state_space and was previously bound
          if ((state_space(i) == 0) .and. (temp_bound_alt(i,2)/=0)) then
             bindtimeindex = temp_bound_alt(i,2)/tinc
             bindtimes_hist_alt(bindtimeindex) = bindtimes_hist_alt(bindtimeindex) + 1 ! adding last duration of binding event to histogram of times
             temp_bound_alt(i,:) = 0 ! setting that B to bound
          end if
       end if
    end do

  end subroutine elim_tempbound


! ----------------------------------------------------------------------
! Final check to add the duration to the histogram for binding times
! ----------------------------------------------------------------------

  subroutine final_check_bound(tinc, &
                            state_space_time, cut_time, &
                            cells, group, bgroup, state_space, temp_bound_alt, &
                            bindtimes_hist_alt)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    integer, dimension(:,:), intent(inout) :: temp_bound_alt
    integer, dimension(:), intent(inout) :: state_space_time, state_space
    real, intent(in) :: cut_time
    integer, intent(in) :: bgroup
    integer, dimension(1,2) :: bound_pair_temp, temp_zero
    integer(int64), intent(in) :: tinc
    integer, dimension(:), intent(inout) :: bindtimes_hist_alt
    integer :: i, j, bindtimeindex, stillbound_yn, groupbit_b
    logical, dimension(2) :: array_match_one, array_match_two, array_match_three
    integer, dimension(:), pointer :: mask => null()
    temp_zero(:,:) = 0

    mask => group%mask

    groupbit_b = group%bitmask(bgroup)
    do i = 1, size(state_space,1)
       if (iand(mask(i),groupbit_b) /= 0) then
          if (temp_bound_alt(i,1)==1) then
             bindtimeindex = temp_bound_alt(i,2)/tinc
             bindtimes_hist_alt(bindtimeindex) = bindtimes_hist_alt(bindtimeindex) + 1 ! adding last duration of bound event to histogram of times
          end if
       end if
    end do

  end subroutine final_check_bound

! ----------------------------------------------------------------------
! Eliminates B from temp_unbound which have become bound, adds final duration to histogram of unbinding times
! ----------------------------------------------------------------------

  subroutine final_check_search(cells, group, bgroup, state_space, &
                temp_unbound, tinc, unbindtimes_hist)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    integer(int64), intent(in) :: tinc
    integer, intent(in) :: bgroup
    integer, dimension(:), intent(inout) :: state_space, unbindtimes_hist
    integer, dimension(:,:), intent(inout) :: temp_unbound
    integer :: i,groupbit_b, bindtimeindex
    integer, dimension(:), pointer :: mask => null()

    mask => group%mask

    groupbit_b = group%bitmask(bgroup)
    do i = 1, size(state_space,1)
       ! print *, "tempunbound ",temp_unbound(i,:)
       if (iand(mask(i),groupbit_b) /= 0) then
          if (temp_unbound(i,1)==1) then
             bindtimeindex = temp_unbound(i,2)/tinc
             unbindtimes_hist(bindtimeindex) = unbindtimes_hist(bindtimeindex) + 1 ! adding last duration of unbound event to histogram of times
          end if
       end if
    end do

  end subroutine final_check_search

! ----------------------------------------------------------------------
! For all (A,B) in bound_pairs, add new binding events to temp_bound
! ----------------------------------------------------------------------

  subroutine check_tempbound(tinc, cells, group, bgroup, &
                                state_space, temp_bound_alt)

    integer, dimension(:,:), intent(inout) :: temp_bound_alt
    integer, dimension(:), intent(inout) :: state_space
    integer, dimension(1,2) :: bound_pair_temp, temp_zero
    integer, dimension(1,3) :: temp_bound_clear
    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    integer, intent(in) :: bgroup
    integer(int64), intent(in) :: tinc
    integer :: i, j, bindtimeindex, yn_exist, groupbit_b
    logical, dimension(2) :: array_match_one, array_match_two, array_match_three
    integer, dimension(:), pointer :: mask => null()
    temp_zero(:,:) = 0

    mask => group%mask

    groupbit_b = group%bitmask(bgroup)
    do i = 1, size(state_space,1)
       if (iand(mask(i),groupbit_b) /= 0) then ! if molecule is a B
          if (state_space(i) == 1)  then ! if B is bound in state space
             if (temp_bound_alt(i,1) == 1) then ! if B was bound the last step
                temp_bound_alt(i,2) = temp_bound_alt(i,2) + tinc
             else
                temp_bound_alt(i,1) = 1
                temp_bound_alt(i,2) = tinc
             end if
          end if
       end if
    end do

  end subroutine check_tempbound


! ----------------------------------------------------------------------
! Eliminates B from temp_unbound which have become bound, adds final duration to histogram of unbinding times
! ----------------------------------------------------------------------

  subroutine elim_tempunbound(cells, group, bgroup, state_space, &
                temp_unbound, tinc, unbindtimes_hist)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    integer(int64), intent(in) :: tinc
    integer, intent(in) :: bgroup
    integer, dimension(:), intent(inout) :: state_space, unbindtimes_hist
    integer, dimension(:,:), intent(inout) :: temp_unbound
    integer :: i,groupbit_b, bindtimeindex
    integer, dimension(:), pointer :: mask => null()

    mask => group%mask

    groupbit_b = group%bitmask(bgroup)
    do i = 1, size(state_space,1) ! for all particles
       if (iand(mask(i),groupbit_b) /= 0) then ! if in group b
          ! check if "bound" in state_space and was previously unbound
          if ((state_space(i) == 1) .and. (temp_unbound(i,2)/=0)) then
             bindtimeindex = temp_unbound(i,2)/tinc
             unbindtimes_hist(bindtimeindex) = unbindtimes_hist(bindtimeindex) + 1 ! adding last duration of unbound event to histogram of times
             temp_unbound(i,:) = 0 ! setting that B to bound
          end if
       end if
    end do

  end subroutine elim_tempunbound


! ----------------------------------------------------------------------
! For all B=0 in state_space, check if already exists in temp_unbound and either add B w/ tinc or add tinc to existing duration
! ----------------------------------------------------------------------

  subroutine check_tempunbound(cells, group, bgroup, state_space, temp_unbound, tinc)

    type(cells_t), intent(inout) :: cells
    type(group_t), intent(in), target :: group
    integer, intent(in) :: bgroup
    integer(int64), intent(in) :: tinc
    integer, dimension(:), intent(inout) :: state_space
    integer, dimension(:,:), intent(inout) :: temp_unbound
    integer :: i,groupbit_b, bindtimeindex
    integer, dimension(:), pointer :: mask => null()

    mask => group%mask

    groupbit_b = group%bitmask(bgroup)
    do i = 1, size(state_space,1)
       if (iand(mask(i),groupbit_b) /= 0) then ! if molecule is a B
          if (state_space(i) == 0)  then ! if B is unbound in state space
             if (temp_unbound(i,1) == 1) then ! if B was unbound the last step
                temp_unbound(i,2) = temp_unbound(i,2) + tinc
             else
                temp_unbound(i,1) = 1
                temp_unbound(i,2) = tinc
             end if
          end if
       end if
    end do

  end subroutine check_tempunbound

    ! ----------------------------------------------------------------------
  ! Assign particles with the unique cluster ID that they are a part of at the initial timestep
  ! ----------------------------------------------------------------------

  subroutine eval_cluster_id_time(cells, yes_delta, r, cluster_id_f)

    type(cells_t), intent(in) :: cells
    integer, dimension(:,:), intent(inout) :: cluster_id_f
    logical, intent(in) :: yes_delta
    real, intent(in) :: r
    integer :: c,d,i,j,idi,idj,itype,jtype,neighcell
    integer :: ncells,nparticles,ntypes,cut_type
    real :: rsq,rcutsq,xtmp,ytmp,ztmp,delx,dely,delz,ri,rj
    logical :: done,typeradius_flag
    real, dimension(:,:), allocatable :: rsqmat
    integer :: o, min_val, max_val, count
    integer, dimension(:), allocatable :: val, unique
    integer, dimension(:), allocatable :: final, cluster_id

    ntypes = cells%particle%ntypes
    ncells = cells%ncells
    nparticles = cells%particle%nparticles
    typeradius_flag = cells%particle%typeradius_flag
    allocate(cluster_id(nparticles))
    allocate(val(nparticles))
    allocate(unique(nparticles))
    cluster_id(:) = 0
    o = 0

    if (size(cluster_id) .ne. nparticles) then
       write(*, '(a)') 'cluster_id size must be nparticles'
       stop
    end if

    if (.not. cells%cells_exist) then
       write(*, '(a)') 'eval_cluster_id_time currently requires cells'
       stop
    end if

    ! select rcutsq method by what's loaded

    rcutsq = r*r
    cut_type = 1
    if (yes_delta) then
       if (typeradius_flag) then
          cut_type = 2
          allocate(rsqmat(ntypes,ntypes))
          do i = 1, ntypes
             do j = 1, ntypes
                rsqmat(i,j) = cells%particle%typeradius(i) &
                            + cells%particle%typeradius(j) &
                            + r
             end do
          end do
          rsqmat(:,:) = rsqmat(:,:)*rsqmat(:,:)
       else
          cut_type = 3
          write(*, '(a)') 'currently not supporting clusters and yes_delta w/o typeradius_flag'
          stop
       end if
    end if

    ! all particles start with own cluster id

    do i = 1, nparticles
       cluster_id(i) = i
    end do

    done = .false.
    progress_loop: do while (.not. done)
       done = .true.
       cell_c: do c = 1, ncells

          ! self loop over all particles in cell

          i = cells%cellhead(c)
          self_i: do while (i > 0)

             itype = cells%particle%type(i)
             idi = cluster_id(i)
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)

             j = cells%next(i)
             self_j: do while (j > 0)

                idj = cluster_id(j)

                ! if particles already in same cluster, cycle

                if (idi .eq. idj) then
                   j = cells%next(j)
                   cycle
                end if

                ! check if contacting

                delx = cells%particle%x(1,j) - xtmp
                dely = cells%particle%x(2,j) - ytmp
                delz = cells%particle%x(3,j) - ztmp
                call cells%domain%minimum_image(delx, dely, delz)
                rsq = delx*delx + dely*dely + delz*delz

                select case (cut_type)
                case(1)
                   continue
                case(2)
                   jtype = cells%particle%type(j)
                   rcutsq = rsqmat(itype,jtype)
                case(3)
                   ri = cells%particle%radius(i)
                   rj = cells%particle%radius(j)
                   rcutsq = (ri + rj + r)
                   rcutsq = rcutsq * rcutsq
                end select

                ! if outside of range, go to next j

                if (rsq >= rcutsq) then
                   j = cells%next(j)
                   cycle
                end if

                ! pick new cluster ids from minimum id in pair
                ! flag that we're not done yet
                ! need a cycle without updates to be done

                cluster_id(i) = min(idi,idj)
                cluster_id(j) = cluster_id(i)
                j = cells%next(j)
                done = .false.

             end do self_j
             i = cells%next(i)
          end do self_i

          ! distinct neighbor loop

          i = cells%cellhead(c)
          distinct_i: do while (i > 0)

             itype = cells%particle%type(i)
             idi = cluster_id(i)
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)

             cell_d: do d = 1, cells%nstencil

                neighcell = cells%stencil_neighbor(c,d)
                if (neighcell == 0) cycle
                j = cells%cellhead(neighcell)

                distinct_j: do while (j > 0)

                   idj = cluster_id(j)

                   ! if particles already in same cluster, cycle

                   if (idi .eq. idj) then
                      j = cells%next(j)
                      cycle
                   end if

                   ! check if contacting

                   delx = cells%particle%x(1,j) - xtmp
                   dely = cells%particle%x(2,j) - ytmp
                   delz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delx, dely, delz)
                   rsq = delx*delx + dely*dely + delz*delz

                   select case (cut_type)
                   case(1)
                      continue
                   case(2)
                      jtype = cells%particle%type(j)
                      rcutsq = rsqmat(itype,jtype)
                   case(3)
                      ri = cells%particle%radius(i)
                      rj = cells%particle%radius(j)
                      rcutsq = (ri + rj + r)
                      rcutsq = rcutsq * rcutsq
                   end select

                   ! if outside of range, go to next j

                   if (rsq >= rcutsq) then
                      j = cells%next(j)
                      cycle
                   end if

                   ! pick new cluster ids
                   ! flag that we're not done yet

                   cluster_id(i) = min(idi,idj)
                   cluster_id(j) = cluster_id(i)
                   j = cells%next(j)
                   done = .false.

                end do distinct_j
             end do cell_d
             i = cells%next(i)
          end do distinct_i
       end do cell_c
       if (done) exit
    end do progress_loop

    ! create list of unique cluster ID's in terms of lowest particle ID
    min_val = minval(cluster_id)-1
    max_val = maxval(cluster_id)
    do while (min_val<max_val)
        o = o+1
        min_val = minval(cluster_id, mask=cluster_id>min_val)
        unique(o) = min_val
    end do
    allocate(final(o), source=unique(1:o))
    
    do i = 1, o                                ! for all invidual clusters
       count = 1
       do j = 1, nparticles                    ! for each particle
          if (cluster_id(j) == final(i)) then  ! check if that particle is in the given cluster in i
             cluster_id_f(i,count) = j
             count = count + 1
             !cluster_id_f(j) = i               ! if so, change the value in cluster_id_f to be unique cluster id i
          end if
       end do
    end do

  end subroutine eval_cluster_id_time


! ----------------------------------------------------------------------
! Evaluates if any particles in each unique cluster at time=0 are still
! part of clusters at the current time=t. If so, it adds a unit of time to
! the cluster_duration, and if not, it outputs the current cluster lifetime
! to cluster_decaytime.
! ----------------------------------------------------------------------

  subroutine eval_cluster_breakup(cells, yes_delta, r, &
                cluster_id_init, cluster_duration, cluster_decaytime)

    type(cells_t), intent(in) :: cells
    integer, dimension(:,:), intent(inout) :: cluster_id_init
    integer, dimension(:), intent(inout) :: cluster_duration, cluster_decaytime
    logical, intent(in) :: yes_delta
    real, intent(in) :: r
    integer :: c,d,i,j,idi,idj,itype,jtype,neighcell,cluster_i_count,count_yn
    integer :: ncells,nparticles,ntypes,cut_type
    real :: rsq,rcutsq,xtmp,ytmp,ztmp,delx,dely,delz,ri,rj
    logical :: done,typeradius_flag
    real, dimension(:,:), allocatable :: rsqmat
    integer :: o,min_val,max_val
    integer, dimension(:), allocatable :: val,unique
    integer, dimension(:), allocatable :: final,part_still_clust
    ntypes = cells%particle%ntypes
    ncells = cells%ncells
    nparticles = cells%particle%nparticles
    typeradius_flag = cells%particle%typeradius_flag

    ntypes = cells%particle%ntypes
    ncells = cells%ncells
    nparticles = cells%particle%nparticles
    typeradius_flag = cells%particle%typeradius_flag
    allocate(val(nparticles))
    allocate(unique(nparticles))

    if (.not. cells%cells_exist) then
       write(*, '(a)') 'eval_cluster_id_time currently requires cells'
       stop
    end if

    ! select rcutsq method by what's loaded
    rcutsq = r*r
    cut_type = 1
    if (yes_delta) then
       if (typeradius_flag) then
          cut_type = 2
          allocate(rsqmat(ntypes,ntypes))
          do i = 1, ntypes
             do j = 1, ntypes
                rsqmat(i,j) = cells%particle%typeradius(i) &
                            + cells%particle%typeradius(j) + r
             end do
          end do
          rsqmat(:,:) = rsqmat(:,:)*rsqmat(:,:)
       else
          cut_type = 3
          write(*, '(a)') 'currently not supporting clusters and yes_delta w/o typeradius_flag'
          stop
       end if
    end if

    ! go through each row of cluster_id_init (unique clusters) that is not entirely=0, & check if all their initial particles are still within definition of a cluster
    do i = 1, nparticles ! for each row/cluster in cluster_id_init

       cluster_i_count = count(cluster_id_init(i,:)/=0) ! number of particles in given cluster i
       o = 0
       if (cluster_i_count /= 0) then                ! if particles exist in cluster i
          ! create list of unique particle IDs in given cluster
          min_val = minval(cluster_id_init(i,:))-1
          max_val = maxval(cluster_id_init(i,:))
          do while (min_val<max_val)
             o = o+1
             min_val = minval(cluster_id_init(i,:), mask=cluster_id_init(i,:)>min_val)
             if (min_val /= 0) then
                 unique(o) = min_val
             else
                 o = o-1
             end if
          end do
          allocate(final(o), source=unique(1:o)) ! unique particle IDs in given cluster
          !print*,final
          !print*,o
          allocate(part_still_clust(o))       ! start out with array [<final>] of 0
          part_still_clust(:) = 0
          ! check all particles in <final> is within cutoff distance of another particle in <final>, if so set the IDs from <final> to 1.
          !print*,'o',o
          do j = 1, o-1
             do d = j+1, o
                itype = cells%particle%type(final(j))
                xtmp = cells%particle%x(1,final(j))
                ytmp = cells%particle%x(2,final(j))
                ztmp = cells%particle%x(3,final(j))
                ! check if contacting
                delx = cells%particle%x(1,final(d)) - xtmp
                dely = cells%particle%x(2,final(d)) - ytmp
                delz = cells%particle%x(3,final(d)) - ztmp
                call cells%domain%minimum_image(delx, dely, delz)
                rsq = delx*delx + dely*dely + delz*delz
                ri = cells%particle%radius(final(j))
                rj = cells%particle%radius(final(d))
                rcutsq = (1 + r)
                rcutsq = rcutsq * rcutsq
                ! if particles are within contact separation, set the IDs from <final> to 1.
                !print*,rsq,rcutsq
                if (rsq < rcutsq) then
                   !print*,'partstillclust'
                   part_still_clust(j) = 1
                   part_still_clust(d) = 1
                end if
             end do
          end do
          ! if any particles in original cluster are still bound (any terms in <part_still_clust> are nonzero), add time unit to cluster_duration(i)
          count_yn = count(part_still_clust/=0) ! number of particles in cluster that are still bound
          if (count_yn /= 0) then
             cluster_duration(i) = cluster_duration(i) + 1 ! add to cluster_duration(i)
             !print*,cluster_duration(i)
          else ! if not (all terms in <part_still_clust> are zero), output current cluster lifetime to cluster_decaytime(i)
             cluster_decaytime(i) = cluster_duration(i)
             !cluster_duration(i) = 0
          end if
          deallocate(part_still_clust)       ! start out with array [<final>] of 0
          deallocate(final) ! unique particle IDs in given cluster
       end if
    end do

  end subroutine eval_cluster_breakup


! ----------------------------------------------------------------------
! At final timestep, outputs all of cluster_duration into cluster_decaytime
! ----------------------------------------------------------------------

  subroutine eval_final_cluster(cells, cluster_duration, cluster_decaytime)

    type(cells_t), intent(in) :: cells
    integer, dimension(:), intent(inout) :: cluster_duration, cluster_decaytime
    integer :: c,d,i,j,nparticles

    nparticles = cells%particle%nparticles

    do i = 1, nparticles ! for each row/cluster in cluster_id_init
       if (cluster_duration(i) /= 0) then
          cluster_decaytime(i) = cluster_duration(i)
       end if
    end do

  end subroutine eval_final_cluster


! ----------------------------------------------------------------------
! Calculates the mean cluster decay time for each size cluster size (at initial timestep)
! ----------------------------------------------------------------------

  subroutine eval_clustersize_decay(nparticles, cluster_count, cluster_decaytime, &
                clustersize_decay)

    integer, dimension(:), intent(inout) :: cluster_count, cluster_decaytime
    real, dimension(:), intent(inout) :: clustersize_decay
    integer, intent(in) :: nparticles
    integer :: c,d,i,j,o

    ! add average duration to the cluster
    do j = 1, nparticles ! for each cluster
       if (cluster_decaytime(j) /= 0) then
          c = count(cluster_count .eq. cluster_count(j)) ! find number of occurences of current cluster size
          ! add the average cluster_decaytime for that cluster to overall value for cluster size
          clustersize_decay(cluster_count(j)) = clustersize_decay(cluster_count(j)) + cluster_decaytime(j)/c
       end if
    end do

  end subroutine eval_clustersize_decay


  ! ----------------------------------------------------------------------
  ! Assign particles the lowest id of all the particles in the same
  ! cluster.  E.g. if particles 2--10 are all in the same cluster,
  ! they would all get the number 2.
  ! ----------------------------------------------------------------------

  subroutine eval_cluster_id(cells, yes_delta, r, cluster_id)

    type(cells_t), intent(in) :: cells
    integer, dimension(:), intent(inout) :: cluster_id
    logical, intent(in) :: yes_delta
    real, intent(in) :: r
    integer :: c,d,i,j,idi,idj,itype,jtype,neighcell
    integer :: ncells,nparticles,ntypes,cut_type
    real :: rsq,rcutsq,xtmp,ytmp,ztmp,delx,dely,delz,ri,rj
    logical :: done,typeradius_flag
    real, dimension(:,:), allocatable :: rsqmat

    ntypes = cells%particle%ntypes
    ncells = cells%ncells
    nparticles = cells%particle%nparticles
    typeradius_flag = cells%particle%typeradius_flag

    if (size(cluster_id) .ne. nparticles) then
       write(*, '(a)') 'cluster_id size must be nparticles'
       stop
    end if

    if (.not. cells%cells_exist) then
       write(*, '(a)') 'eval_cluster_id currently requires cells'
       stop
    end if

    ! select rcutsq method by what's loaded

    rcutsq = r*r
    cut_type = 1
    if (yes_delta) then
       if (typeradius_flag) then
          cut_type = 2
          allocate(rsqmat(ntypes,ntypes))
          do i = 1, ntypes
             do j = 1, ntypes
                rsqmat(i,j) = cells%particle%typeradius(i) &
                            + cells%particle%typeradius(j) &
                            + r
             end do
          end do
          rsqmat(:,:) = rsqmat(:,:)*rsqmat(:,:)
       else
          cut_type = 3
          write(*, '(a)') 'currently not supporting clusters and yes_delta w/o typeradius_flag'
          stop
       end if
    end if

    ! all particles start with own cluster id

    do i = 1, nparticles
       cluster_id(i) = i
    end do

    done = .false.
    progress_loop: do while (.not. done)
       done = .true.
       cell_c: do c = 1, ncells

          ! self loop over all particles in cell

          i = cells%cellhead(c)
          self_i: do while (i > 0)

             itype = cells%particle%type(i)
             idi = cluster_id(i)
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)

             j = cells%next(i)
             self_j: do while (j > 0)

                idj = cluster_id(j)

                ! if particles already in same cluster, cycle

                if (idi .eq. idj) then
                   j = cells%next(j)
                   cycle
                end if

                ! check if contacting

                delx = cells%particle%x(1,j) - xtmp
                dely = cells%particle%x(2,j) - ytmp
                delz = cells%particle%x(3,j) - ztmp
                call cells%domain%minimum_image(delx, dely, delz)
                rsq = delx*delx + dely*dely + delz*delz

                select case (cut_type)
                case(1)
                   continue
                case(2)
                   jtype = cells%particle%type(j)
                   rcutsq = rsqmat(itype,jtype)
                case(3)
                   ri = cells%particle%radius(i)
                   rj = cells%particle%radius(j)
                   rcutsq = (ri + rj + r)
                   rcutsq = rcutsq * rcutsq
                end select

                ! if outside of range, go to next j

                if (rsq >= rcutsq) then
                   j = cells%next(j)
                   cycle
                end if

                ! pick new cluster ids from minimum id in pair
                ! flag that we're not done yet
                ! need a cycle without updates to be done

                cluster_id(i) = min(idi,idj)
                cluster_id(j) = cluster_id(i)
                j = cells%next(j)
                done = .false.

             end do self_j
             i = cells%next(i)
          end do self_i

          ! distinct neighbor loop

          i = cells%cellhead(c)
          distinct_i: do while (i > 0)

             itype = cells%particle%type(i)
             idi = cluster_id(i)
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)

             cell_d: do d = 1, cells%nstencil

                neighcell = cells%stencil_neighbor(c,d)
                if (neighcell == 0) cycle
                j = cells%cellhead(neighcell)

                distinct_j: do while (j > 0)

                   idj = cluster_id(j)

                   ! if particles already in same cluster, cycle

                   if (idi .eq. idj) then
                      j = cells%next(j)
                      cycle
                   end if

                   ! check if contacting

                   delx = cells%particle%x(1,j) - xtmp
                   dely = cells%particle%x(2,j) - ytmp
                   delz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delx, dely, delz)
                   rsq = delx*delx + dely*dely + delz*delz

                   select case (cut_type)
                   case(1)
                      continue
                   case(2)
                      jtype = cells%particle%type(j)
                      rcutsq = rsqmat(itype,jtype)
                   case(3)
                      ri = cells%particle%radius(i)
                      rj = cells%particle%radius(j)
                      rcutsq = (ri + rj + r)
                      rcutsq = rcutsq * rcutsq
                   end select

                   ! if outside of range, go to next j

                   if (rsq >= rcutsq) then
                      j = cells%next(j)
                      cycle
                   end if

                   ! pick new cluster ids
                   ! flag that we're not done yet

                   cluster_id(i) = min(idi,idj)
                   cluster_id(j) = cluster_id(i)
                   j = cells%next(j)
                   done = .false.

                end do distinct_j
             end do cell_d
             i = cells%next(i)
          end do distinct_i
       end do cell_c
       if (done) exit
    end do progress_loop

  end subroutine eval_cluster_id


  ! ----------------------------------------------------------------------
  ! Assign particles the lowest id of all the particles in the same
  ! cluster.  E.g. if particles 2--10 are all in the same cluster,
  ! they would all get the number 2.
  ! ----------------------------------------------------------------------

  subroutine eval_cluster_id_parttype(cells, yes_delta, r, cluster_id, typea, typeb)

    type(cells_t), intent(in) :: cells
    integer, dimension(:), intent(inout) :: cluster_id
    logical, intent(in) :: yes_delta
    real, intent(in) :: r
    integer, intent(in) :: typea,typeb
    integer :: c,d,i,j,idi,idj,itype,jtype,neighcell
    integer :: ncells,nparticles,ntypes,cut_type
    real :: rsq,rcutsq,xtmp,ytmp,ztmp,delx,dely,delz,ri,rj
    logical :: done,typeradius_flag
    real, dimension(:,:), allocatable :: rsqmat

    ntypes = cells%particle%ntypes
    ncells = cells%ncells
    nparticles = cells%particle%nparticles
    typeradius_flag = cells%particle%typeradius_flag

    if (size(cluster_id) .ne. nparticles) then
       write(*, '(a)') 'cluster_id size must be nparticles'
       stop
    end if

    if (.not. cells%cells_exist) then
       write(*, '(a)') 'eval_cluster_id currently requires cells'
       stop
    end if

    ! select rcutsq method by what's loaded

    rcutsq = r*r
    cut_type = 1
    if (yes_delta) then
       if (typeradius_flag) then
          cut_type = 2
          allocate(rsqmat(ntypes,ntypes))
          do i = 1, ntypes
             do j = 1, ntypes
                rsqmat(i,j) = cells%particle%typeradius(i) &
                            + cells%particle%typeradius(j) &
                            + r
             end do
          end do
          rsqmat(:,:) = rsqmat(:,:)*rsqmat(:,:)
       else
          cut_type = 3
          write(*, '(a)') 'currently not supporting clusters and yes_delta w/o typeradius_flag'
          stop
       end if
    end if

    ! all particles start with own cluster id

    do i = 1, nparticles
       cluster_id(i) = i
    end do

    done = .false.
    progress_loop: do while (.not. done)
       done = .true.
       cell_c: do c = 1, ncells

          ! self loop over all particles in cell

          i = cells%cellhead(c)
          self_i: do while (i > 0)

             itype = cells%particle%type(i)
             ! if particle i not of the two types to be evaluated, continue
             if ((itype .ne. typea) .and. (itype .ne. typeb)) then
                i = cells%next(i)
                cycle
             end if
             idi = cluster_id(i)
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)

             j = cells%next(i)
             self_j: do while (j > 0)

                jtype = cells%particle%type(j)
                ! if particle j not of the two types to be evaluated, continue
                if ((jtype .ne. typea) .and. (jtype .ne. typeb)) then
                   j = cells%next(j)
                   cycle
                end if
                idj = cluster_id(j)
                ! if particles already in same cluster, cycle

                if (idi .eq. idj) then
                   j = cells%next(j)
                   cycle
                end if

                ! check if contacting

                delx = cells%particle%x(1,j) - xtmp
                dely = cells%particle%x(2,j) - ytmp
                delz = cells%particle%x(3,j) - ztmp
                call cells%domain%minimum_image(delx, dely, delz)
                rsq = delx*delx + dely*dely + delz*delz

                select case (cut_type)
                case(1)
                   continue
                case(2)
                   jtype = cells%particle%type(j)
                   rcutsq = rsqmat(itype,jtype)
                case(3)
                   ri = cells%particle%radius(i)
                   rj = cells%particle%radius(j)
                   rcutsq = (ri + rj + r)
                   rcutsq = rcutsq * rcutsq
                end select

                ! if outside of range, go to next j

                if (rsq >= rcutsq) then
                   j = cells%next(j)
                   cycle
                end if

                ! pick new cluster ids from minimum id in pair
                ! flag that we're not done yet
                ! need a cycle without updates to be done

                cluster_id(i) = min(idi,idj)
                cluster_id(j) = cluster_id(i)
                j = cells%next(j)
                done = .false.

             end do self_j
             i = cells%next(i)
          end do self_i

          ! distinct neighbor loop

          i = cells%cellhead(c)
          distinct_i: do while (i > 0)

             itype = cells%particle%type(i)
             ! if particle i not of the two types to be evaluated, continue
             if ((itype .ne. typea) .and. (itype .ne. typeb)) then
                i = cells%next(i)
                cycle
             end if
             itype = cells%particle%type(i)
             idi = cluster_id(i)
             xtmp = cells%particle%x(1,i)
             ytmp = cells%particle%x(2,i)
             ztmp = cells%particle%x(3,i)

             cell_d: do d = 1, cells%nstencil

                neighcell = cells%stencil_neighbor(c,d)
                if (neighcell == 0) cycle
                j = cells%cellhead(neighcell)

                distinct_j: do while (j > 0)

                   jtype = cells%particle%type(j)
                   ! if particle j not of the two types to be evaluated, continue
                   if ((jtype .ne. typea) .and. (jtype .ne. typeb)) then
                      j = cells%next(j)
                      cycle
                   end if
                   idj = cluster_id(j)

                   ! if particles already in same cluster, cycle

                   if (idi .eq. idj) then
                      j = cells%next(j)
                      cycle
                   end if

                   ! check if contacting

                   delx = cells%particle%x(1,j) - xtmp
                   dely = cells%particle%x(2,j) - ytmp
                   delz = cells%particle%x(3,j) - ztmp
                   call cells%domain%minimum_image(delx, dely, delz)
                   rsq = delx*delx + dely*dely + delz*delz

                   select case (cut_type)
                   case(1)
                      continue
                   case(2)
                      jtype = cells%particle%type(j)
                      rcutsq = rsqmat(itype,jtype)
                   case(3)
                      ri = cells%particle%radius(i)
                      rj = cells%particle%radius(j)
                      rcutsq = (ri + rj + r)
                      rcutsq = rcutsq * rcutsq
                   end select

                   ! if outside of range, go to next j

                   if (rsq >= rcutsq) then
                      j = cells%next(j)
                      cycle
                   end if

                   ! pick new cluster ids
                   ! flag that we're not done yet

                   cluster_id(i) = min(idi,idj)
                   cluster_id(j) = cluster_id(i)
                   j = cells%next(j)
                   done = .false.

                end do distinct_j
             end do cell_d
             i = cells%next(i)
          end do distinct_i
       end do cell_c
       if (done) exit
    end do progress_loop

  end subroutine eval_cluster_id_parttype

  ! ----------------------------------------------------------------------
  ! Add a 3D array g to a version of itself g',
  ! where g'(i,j,k) = g(xmax+1-i,ymax+1-j,zmax+1-k),
  ! and divide by 2, to make a version that is symmetric.
  ! Avoids allocating a copy of the original array.
  ! ----------------------------------------------------------------------

  subroutine symmetrize_3d(g)

    real, dimension(:,:,:), intent(inout) :: g
    integer :: i,j,k,imax,jmax,kmax,imid,jmid,kmid

    imax = size(g,1)
    jmax = size(g,2)
    kmax = size(g,3)
    imid = imax / 2
    jmid = jmax / 2
    kmid = kmax / 2

    ! Divide everything by 2

    g(:,:,:) = 0.5 * g(:,:,:)

    ! Go through elements on the side of the array with small
    ! k-indices, adding elements on the opposide side of the
    ! the array, through its center.  E.g. (imax,jmax,kmax)
    ! added to (1,1,1).  Only go up to kmid = kmax / 2, in
    ! order to avoid double-counting.  If kmax is odd, say
    ! 5, only go up to the rounded-down half, say 2.  This
    ! leaves the middle plane, which must be dealt with below
    ! in a similar fashion.

    do k = 1, kmid
       do j = 1, jmax
          do i = 1, imax
             g(i,j,k) = g(i,j,k) + g(imax+1-i,jmax+1-j,kmax+1-k)
          end do
       end do
    end do

    ! Copy the small-k elements to the large-k elements.

    do k = 1, kmid
        do j = 1, jmax
          do i = 1, imax
             g(imax+1-i,jmax+1-j,kmax+1-k) = g(i,j,k)
          end do
       end do
    end do

    ! If kmax is even, we're done.
    ! If kmax is odd, fix the skipped middle plane.

    if (mod(kmax,2) == 0) return

    do j = 1, jmid
       do i = 1, imax
          g(i,j,kmid+1) = g(i,j,kmid+1) + g(imax+1-i,jmax+1-j,kmid+1)
       end do
    end do

    do j = 1, jmid
       do i = 1, imax
          g(imax+1-i,jmax+1-j,kmid+1) = g(i,j,kmid+1)
       end do
    end do

    ! If jmax is even, we're done.
    ! If jmax is odd, fix the skipped middle line.

    if (mod(jmax,2) == 0) return

    do i = 1, imid
       g(i,jmid+1,kmid+1) = g(i,jmid+1,kmid+1) + g(imax+1-i,jmid+1,kmid+1)
    end do

    do i = 1, imid
       g(imax+1-i,jmid+1,kmid+1) = g(i,jmid+1,kmid+1)
    end do

    ! If imax is even, we're done.
    ! If imax is odd, fix the center point.

    if (mod(imax,2) == 0) return

    g(imid+1,jmid+1,kmid+1) = 2 * g(imid+1,jmid+1,kmid+1)

  end subroutine symmetrize_3d


end module mod_stat_func
