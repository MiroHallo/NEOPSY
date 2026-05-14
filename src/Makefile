# ######################################################################
# Makefile for the Multizonal Transdimensional Inversion (NEOPSY)
#
# Purpose: Compiles Fortran90 modules and main program for evaluation
#          of Multizonal Transdimensional Bayesian Inversion.
#
# Requirements: Geopsy 3.X, PT library, Fortran90 and MPI compilers
#
# Author:  Miroslav HALLO
# Date:    2026/05
# License: GNU GPL v3
# ######################################################################
#
# Usage:
#    > make all        # Compile all binary files
#    > make clean      # Remove all modules and objects
#    > make distclean  # Remove all modules, objects, and binary files
#
# Example:
#    > make clean MODE=fast all   # Fresh compilation with -O3 optim.
#
# ######################################################################

# Compiler without MPI
F90 = gfortran

# Compiler with MPI
MF90 = mpif90

# Path to the Geopsy 3.X library
GEOPSY_PATH=/opt/geopsy/lib

# Path to the Parallel Tempering (PT) library
PT_PATH=./pt/lib

# Compilation mode (debug/release/fast)
MODE ?= release

# Flags
BASE_FLAGS = -cpp
ifeq ($(MODE), debug)
    FFLAGS = $(BASE_FLAGS) -O0 -g -Wall -fcheck=all
else ifeq ($(MODE), release)
    FFLAGS = $(BASE_FLAGS) -O2
else ifeq ($(MODE), fast)
    FFLAGS = $(BASE_FLAGS) -O3 -march=native
endif

# Program names
TARGETS = ../npdc ../tiser ../timpi ../tires

# Object lists
OBJ = neoLib.o neoInput.o neoData.o io.o
TISER_OBJ = $(OBJ) neoForward.o neoinv.o $(PT_PATH)/pt.serial.o
TIMPI_OBJ = $(OBJ) neoForward.o neoinv.o $(PT_PATH)/pt.mpi.o
TIRES_OBJ = $(OBJ) nr.o sh.o qwl.o neores.o

# ----------------------------------------------------------------------

# Compile final programs
all: $(TARGETS)
	@echo "SUCCESS"

../npdc: npdc.o
	$(F90) $(FFLAGS) -o $@ npdc.o -L$(GEOPSY_PATH) -Wl,-rpath,"$(GEOPSY_PATH)" -lQGpCoreWave

../tiser: $(TISER_OBJ)
	$(F90) $(FFLAGS) -o $@ $(TISER_OBJ)

../timpi: $(TIMPI_OBJ)
	$(MF90) $(FFLAGS) -o $@ $(TIMPI_OBJ)

../tires: $(TIRES_OBJ)
	$(F90) $(FFLAGS) -o $@ $(TIRES_OBJ)

# ----------------------------------------------------------------------

# Global rule for Fortran 90 objects
%.o: %.f90
	$(F90) $(FFLAGS) -c $< -o $@

# Create Fortran 77 object
nr.o: nr.for
	$(F90) $(FFLAGS) -ffixed-form -c $< -o $@

# Create global module (only create the .mod file)
global_module: neoinv.f90
	@sed -n '/MODULE GlobalBox/,/END MODULE/p' neoinv.f90 > globalbox_tmp.f90
	$(F90) $(FFLAGS) -c globalbox_tmp.f90 -o globalbox_tmp.o
	@rm -f globalbox_tmp.f90 globalbox_tmp.o

# Create neoForward.o object (with global_module)
neoForward.o: neoForward.f90 global_module neoInput.o neoData.o neoLib.o io.o
	$(F90) $(FFLAGS) -c $< -o $@

# Create neoinv.o object (after neoForward.o)
neoinv.o: neoinv.f90 neoForward.o neoInput.o neoData.o neoLib.o
	$(F90) $(FFLAGS) -c $< -o $@

# Dependencies
neoInput.o: neoLib.o
neoData.o: neoInput.o
neores.o: sh.o qwl.o neoLib.o neoInput.o neoData.o io.o

# ----------------------------------------------------------------------
# Clean after compilation
clean:
	rm -f *.o *.mod

# ----------------------------------------------------------------------
# Distribution clean
distclean:
	rm -f *.o *.mod $(TARGETS)

# ----------------------------------------------------------------------
.PHONY: all clean distclean