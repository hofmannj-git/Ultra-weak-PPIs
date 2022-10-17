# fc must be ifort or gfortran
# mode must be opt or debug
# loc must be mac_air, stampede, or desktop

fc = gfortran
mode = opt
loc = sherlock

# Check for legit Fortran compiler

ifeq ($(fc), ifort)
else
ifeq ($(fc), gfortran)
else
$(error Bad Fortran compiler variable)
endif
endif

# Set paths based on location variable set above
# Produce error if unknown location

ifeq ($(loc), sherlock)
hdf5_lib_dir = /home/users/hofmannj/bin/hdf5/lib
hdf5_include_dir = /home/users/hofmannj/bin/hdf5/include
fftw3_lib_dir = /home/users/hofmannj/bin/fftw/lib
fftw3_include_dir = /home/users/hofmannj/bin/fftw/include
else
ifeq ($(loc), desktop)
hdf5_lib_dir = /usr/local/opt/hdf5/lib
hdf5_include_dir = /usr/local/opt/hdf5/include
fftw3_lib_dir = /usr/local/opt/fftw/lib
fftw3_include_dir = /usr/local/opt/fftw/include
#omp_lib_dir = /usr/lib/gcc/x86_64-unknown-linux-gnu/4.9.2/finclude
#omp_include_dir = /usr/lib/gcc/x86_64-unknown-linux-gnu/4.9.2/finclude
else
$(error Bad location (loc variable))
endif
endif

# Check for legit copiler option

ifeq ($(mode), opt)
else
ifeq ($(mode), debug)
else
$(error Bad compiler mode)
endif
endif

# No machine-specific variables below

programs = bin_msd_tri cluster_dist_tri cluster_dist_tri_time nc_dist_tri pair_dist_func pair_distance_time bindtime_pairtraj kd_state_percbound
modules = mod_particle mod_fftw3 mod_constants mod_domain mod_data mod_group mod_cells mod_stat_func mod_statistics

programs_obj = $(programs:=.o)
modules_obj = $(modules:=.o)

ifort_opt = -O3 -ip -ipo -parallel -stand f08
ifort_debug = -g -check all -fpe0 -traceback -debug extended -stand f08
ifort_compromise = -O3 -g -traceback
gfortran_debug = -O0 -g -fbounds-check -fbacktrace -std=f2008
gfortran_opt = -O2 -std=f2008
fc_flags = $($(fc)_$(mode))

fc_flags += -I$(hdf5_include_dir)
fc_flags += -I$(fftw3_include_dir)
ld_flags = -L$(hdf5_lib_dir) -lhdf5_fortran -lhdf5
ld_flags += -L$(fftw3_lib_dir) -lfftw3_omp -lfftw3 -lm
ld_flags += -fopenmp

all: $(programs)

$(programs_obj): $(modules_obj)

%: %.o
	$(fc) $(fc_flags) -o $@ $^ $(modules_obj) $(ld_flags)

%.o: %.f90
	$(fc) $(fc_flags) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod

veryclean: clean
	rm -f *~ $(programs)
