# ======================================================================
# Makefile for SA_MESH

# ======================================================================
# Make targets (defined below).
.PHONY: default gfortran ifort mingw_static mpi_gcc mpi_intel symbols debug netcdf all clean veryclean
default: all

# ======================================================================
# Options.
#   If the variable is not defined, an assumed value is assigned.
#   Options can be overwritten by user overrides or scripts.

# DIST: Compiler family/distribution.
#   - Blank/undefined to use GNU/GCC compiler (default).
#   - 'intel' to use Intel compiler.
#   - 'mingw' to use GNU/GCC compiler; overrides 'rm' with 'del' for MS-Windows/MS-DOS environment.

# MPI: Parallel/serial compilation.
#   - Blank/undefined to compile in serial (default).
#   - 'ompi' to compile using OMPI compiler.

# LSS: Land surface scheme (LSS).
#   - Blank/undefined to include default versions of CLASS+SVS (default).

# ROUTE: Routing scheme.
#   - Blank/undefined to include default versions of WF_ROUTE,SA_RTE+RTE (default).

# DEBUG: Debugging flags and options.
#   - Blank/undefined to compile with 'o2' optimization (default).
#   - 'yes' to include debug options and disable compiler optimization.

# ======================================================================
# File/object names.
include makefile.def

# ======================================================================
# Pre-configured targets: Compiler options.
ifeq ($(filter gfortran, $(MAKECMDGOALS)), gfortran)
  DIST ?=
  MPI ?=
else ifeq ($(filter ifort, $(MAKECMDGOALS)), ifort)
  DIST ?= intel
  MPI ?=
else ifeq ($(filter mingw_static, $(MAKECMDGOALS)), mingw_static)
  DIST ?= mingw
  MPI ?=
else ifeq ($(filter mpi_gcc, $(MAKECMDGOALS)), mpi_gcc)
  DIST ?=
  MPI ?= ompi
else ifeq ($(filter mpi_intel, $(MAKECMDGOALS)), mpi_intel)
  DIST ?= intel
  MPI ?= ompi
endif

# ======================================================================
# Pre-configured targets: Debug symbols.
ifeq ($(filter debug, $(MAKECMDGOALS)), debug)
  SYMBOLS ?= yes
  DEBUG ?= yes
else ifeq ($(filter symbols, $(MAKECMDGOALS)), symbols)
  SYMBOLS ?= yes
endif

# ======================================================================
# Pre-configured targets: Double precision (where supported).
ifeq ($(filter double, $(MAKECMDGOALS)), double)
  DOUBLE ?= yes
endif

# Summary.
ifdef SUMMARY
  $(info DIST    = $(DIST))
  $(info MPI     = $(MPI))
  $(info SYMBOLS = $(SYMBOLS))
  $(info DEBUG   = $(DEBUG))
  $(info DOUBLE  = $(DOUBLE))
endif

# ======================================================================
# Pre-configured targets: netCDF library.
# This target will call 'nf-config' via the active shell.
# However, if the netCDF library is installed,
# 'nf-config' should be installed as well.
ifeq ($(filter netcdf, $(MAKECMDGOALS)), netcdf)
  ifeq (, $(shell which nf-config))
    $(error The 'netcdf' target is specified but 'nf-config' cannot be found)
  else
    LIBNCO = $(shell nf-config --fflags) -DNETCDF
  endif
  ifeq (, $(shell which nc-config))
    $(error The 'netcdf' target is specified but 'nc-config' cannot be found)
  else
    LIBNCL = $(shell nf-config --flibs)
  endif
  ifdef SUMMARY
    $(info LIBNCO = $(LIBNCO))
    $(info LIBNCL = $(LIBNCL))
  endif
endif

# ======================================================================
# Targets.
gfortran: all
ifort: all
mingw_static: all
mpi_gcc: all
mpi_intel: all
symbols: all
debug: all
double: all
netcdf: all

# ======================================================================
# Compiler overrides.
ifdef FC_OVERRIDE
  FC = FC_OVERRIDE
endif
ifdef CC_OVERRIDE
  CC = CC_OVERRIDE
endif

# ======================================================================
# Compiler check (if not the 'clean' or 'veryclean' targets).
# Minimum requirement.
# Intel 16+ (15+):
#   - Intel 15 introduces full Fortran 2003 support but is untested.
#   - Intel 14 will compile but may stall during run-time. This is
#       presumed to be the result of partial Fortran 2003 support.
# GNU/gcc 5+:
#   - GNU/gcc 4 does not implement the necessary Fortran 2003 features.
ifeq ($(filter clean veryclean, $(MAKECMDGOALS)), )
  ifeq ($(DIST), intel)
    ifneq (, $(shell which ifx))
    else ifneq (, $(shell which icc))
      ifeq ($(shell test $$(icc -dumpversion | cut -d '.' -f 1) -lt 16; echo $$?), 0)
        $(error The code requires Intel compiler version 16 or higher)
      endif
    else
      $(error The 'intel' compiler cannot be found)
    endif
  else
    ifneq (, $(shell which gcc))
      ifeq ($(shell test $$(gcc -dumpversion | cut -d '.' -f 1) -lt 5; echo $$?), 0)
        $(error The code requires GNU/gcc and GNU/gfortran version 5 or higher)
      endif
    else
      $(error The 'gnu' compiler cannot be found)
    endif
  endif
endif

# ======================================================================
# Compiler and options.
ifeq ($(DIST), intel)
  ifneq (, $(shell which ifx))
    FC = ifx
    CC = icx
  else
    FC = ifort
    CC = icc
  endif
  LFLAG = -c -g -O0 -traceback
  FFLAG = -fpp -check bounds -fpe0 -fp-model source
  CFLAG =
  ifeq ($(shell test $$($(CC) -dumpversion | cut -d '.' -f 1) -lt 17; echo $$?), 0)
    CFLAG += -no-multibyte-chars
  endif
else
  FC = gfortran
  CC = gcc
  LFLAG = -c -g -fbacktrace
  FFLAG = -cpp -ffree-form -ffree-line-length-none -fcray-pointer -fbounds-check -ffpe-trap=invalid,zero,overflow -Wconversion -Wsurprising -Wintrinsic-shadow -Wtarget-lifetime
  ifeq ($(shell test $$($(CC) -dumpversion | cut -d '.' -f 1) -gt 5; echo $$?), 0)
    FFLAG += -Winteger-division
  endif
  CFLAG =
endif

# Override debugging options if 'DEBUG' not enabled.
ifndef DEBUG
  ifeq ($(DIST), intel)
    FFLAG = -fpp -fp-model precise
  else
    FFLAG = -cpp -ffree-form -ffree-line-length-none -fcray-pointer
  endif
  ifndef SYMBOLS
    LFLAG = -c -O2
  endif
  CLEANUP = @$(MAKE) -s clean DIST=$(DIST)
endif

# Override compile options if 'DOUBLE' enabled.
ifdef DOUBLE
  ifeq ($(DIST), intel)
    FFLAG += -r8
  else
    FFLAG += -fdefault-double-8 -fdefault-real-8
  endif
endif

# Output: sa_mesh.
ifneq ($(MPI), ompi)
  OUT ?= sa_mesh
else
  OUT ?= mpi_sa_mesh
endif

# If MPI is enabled, switch to OMPI compiler and rename output.
# Otherwise add MPI stub to 'OBJECTS'.
ifeq ($(MPI), ompi)
  ifneq (, $(shell which mpiifx))
    FC = mpiifx
  else ifneq (, $(shell which mpif90))
    FC = mpif90
  else ifneq (, $(shell which mpifort))
    FC = mpifort
  else
    $(error The 'mpi' Fortran compiler cannot be found)
  endif
  ifneq (, $(shell which mpiicx))
    CC = mpiicx
  else ifneq (, $(shell which mpicc))
    CC = mpicc
  else
    $(error The 'mpi' C compiler cannot be found)
  endif
else
  OBJECTS := mpi_stub.o $(OBJECTS)
endif

# Override 'rm' with 'del' and add static option for MinGW.
ifeq ($(DIST), mingw)
  BIN_DEL = del
  LLINK = -static
  FC = gfortran
  CC = gcc
  OUT = sa_mesh_static
else
  BIN_DEL = rm
endif

# Summary.
ifdef SUMMARY
  $(info FC         = $(FC))
  $(info CC         = $(CC))
  $(info LFLAG      = $(LFLAG))
  $(info FFLAG      = $(FFLAG))
  $(info CFLAG      = $(CFLAG))
  $(info OUT        = $(OUT))
  $(info CLEANUP    = $(CLEANUP))
  $(info BIN_DEL    = $(BIN_DEL))
endif

# ======================================================================
# General rules.
%.o: %.f
	$(FC) $(LFLAG) $<
%.o: %.F90
	$(FC) $(LFLAG) $(FFLAG) $(INC_DIRS) $(DFLAG) $(LIBNCO) $<
%.o: %.f90
	$(FC) $(LFLAG) $(FFLAG) $(LIBNCO) $<
%.o: %.for
	$(FC) $(LFLAG) $<
%.o: %.c
	$(CC) $(LFLAG) $(CFLAG) $(INC_DIRS) $<

# ======================================================================
# Special rules for SPS/rpnphy (including SVS).
SPS_BASE_DIR=..
SPS_DIR=$(SPS_BASE_DIR)/sps
SPS_BUILD_DIR=$(SPS_BASE_DIR)/sps_build
LIBSPS = \
  -L$(SPS_BUILD_DIR)/src/rpnphy/rpnphy -lrpnphy \
  -L$(SPS_BUILD_DIR)/src/modelutils/modelutils -lmodelutils -lmodelutils_tmg_stubs \
  -L$(SPS_BUILD_DIR)/src/rmn -lrmn \
  -L$(SPS_BUILD_DIR)/src/rmn/App/src/lib -lApp \
  -L$(SPS_BUILD_DIR)/src/tdpack -ltdpack \
  -L$(SPS_BUILD_DIR)/src/rpncomm/src -lrpncomm
  ifeq ($(DIST), intel)
    LIBSPS += -liomp5 -lpthread
  else
    LIBSPS += -lgomp
  endif
runsvs_mesh.o: runsvs_mesh.F90
	$(FC) $(LFLAG) $(shell $(SPS_BUILD_DIR)/rpnphy-config --fflags) \
	-I$(SPS_DIR)/src/modelutils/include \
	-I$(SPS_DIR)/src/rpnphy/src/utils \
	-I$(SPS_DIR)/src/rpnphy/src/base \
	-I$(SPS_DIR)/src/rpnphy/src/surface \
	-I$(SPS_BUILD_DIR)/src/rpnphy/rpnphy/modules \
	-I$(SPS_BUILD_DIR)/src/modelutils/modelutils/modules \
	-I$(SPS_BUILD_DIR)/src/tdpack/include $<

# ======================================================================
# Make target: all
# Deletes object and modules files unless 'DEBUG' has a value.
all: ${OBJECTS}
	$(FC) $(OBJECTS) -o $(OUT) $(LLINK) $(LIBNCL) $(LIBSPS)
	$(CLEANUP)

# ======================================================================
# Make target: clean
# Remove object and module files.
clean:
	-$(BIN_DEL) *.mod *.o

# ======================================================================
# Make target: veryclean
# Remove object and module files, and output file.
veryclean: clean
	-$(BIN_DEL) $(OUT)
