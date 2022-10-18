This folder contains the following Fortran scripts, for analysis of simulation runs:

	'Makefile' -- makefile to compile all necessary scripts on compute cluster or local desktop

	PROGRAMS:

		'bin_msd_tri.f90' - calculation of mean squared displacement as a function of lag time

		'cluster_dist_tri.f90' -- calculation of probability distribution of cluster sizes

		'cluster_dist_tri_time.f90' -- calculation of cluster dissociation times

		'kd_state_percbound.f90' -- calculation of fraction of bound protein and apparent binding affinity

		'nc_dist_tri.f90' -- calculation of coordination number distributions

		'pair_dist_func.f90' -- calculation of radial distribution function, g(r)

		'pair_distance_time.f90' -- calculation of the separation distance between particle pairs, starting as bound


	MODULES:

		'mod_cells.f90' -- creates cells within periodic domain to speed up calculations

		'mod_constants.f90' -- definitions of commonly-used fixed constants, e.g. pi

		'mod_data.f90' -- storage and reading of particle trajectory data

		'mod_domain.f90' -- establishing dimensions and values of periodic boundary domain

		'mod_fftw3.f90' -- FFT inclusion

		'mod_group.f90' -- create groups of particles for evaluation of interaction potentials, etc.

		'mod_particle.f90' -- particle type, with attributes of size, position, etc.

		'mod_stat_func.f90' -- statistical functions for calculation, called in all programs

