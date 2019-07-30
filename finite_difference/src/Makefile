#------------------------------------------------------------------------------
# BSD 2-Clause License
# 
# Copyright (c) 2018-2019, Science and Technology Facilities Council.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#------------------------------------------------------------------------------
# Authors: A. R. Porter, STFC Daresbury Laboratory
#          J. Henrichs, Australian Bureau of Meteorology

# Makefile for the dl_esm_inf library, intended to be driven by the top-level
# Makefile in dl_esm_inf/finite_difference.
#
# This Makefile picks up the compiler to use plus any flags from
# environment variables. e.g. to use gfortran:
# 
# export F90=mpif90
# export F90FLAGS="-O3"
# export AR=ar
#
# or, to build the serial version:
#
# export F90=gfortran

MPI ?= no
OCL ?= no

ifeq ($(MPI), yes)
	PARALLEL_MOD=parallel_utils_mod
	F90 ?= mpif90
	API_LIB ?= lib_dm_fd.a
else
	PARALLEL_MOD=parallel_utils_stub_mod
	F90 ?= gfortran
	API_LIB ?= lib_fd.a
endif

ifeq ($(OCL), yes)
	FORTCL_TARGET=
else
	FORTCL_TARGET=stub
endif

BASE_SOURCES = argument_mod.f90 \
          global_parameters_mod.f90 \
          decomposition_mod.f90 \
          field_mod.f90 \
          gocean_mod.F90 \
          grid_mod.f90 \
          halo_mod.f90 \
          kernel_mod.f90 \
          kind_params_mod.f90 \
          $(PARALLEL_MOD).f90 \
          parallel_comms_mod.f90 \
          parallel_mod.f90 \
          region_mod.f90 \
          tile_mod.f90

# Get the absolute path to the makefile
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
# Then get the absolute path of the makefile directory,
# which is the source directory.
ROOT  := $(dir $(mkfile_path))

# Convert all source files to use the absolute path.
# This way the makefile can be used in any build directory
SOURCES = $(BASE_SOURCES:%=$(ROOT)/%)

# Create a list of all object files from the F90/f90 source files
MODULES := $(BASE_SOURCES:%.f90=%.o)
MODULES := $(MODULES:%.F90=%.o)

# Location of the FortCL library (git sub-module)
FORTCL_DIR = $(ROOT)../../external/FortCL/src

.PHONY: install clean all distclean

all: ${API_LIB}

# Create the archive. Note that we include the FortCL objects in this
# archive but take care to exclude clfortran.o (as that only has
# interface definitions, not implementation). We could use filter-out
# for this but that is specific to Gnu Make.
${API_LIB}: fortcl ${MODULES}
	${AR} ${ARFLAGS} ${API_LIB} ${FORTCL_DIR}/*.o ${MODULES}
	${AR} d ${API_LIB} clfortran.o

fortcl:
	${MAKE} -C ${FORTCL_DIR} F90=${F90} ${FORTCL_TARGET}

install:
	${MAKE} -C .. install

clean:
	${MAKE} -C ${FORTCL_DIR} clean
	rm -f *.o *.mod *.a

distclean: clean
	${MAKE} -C ${FORTCL_DIR} allclean

# Inter-module dependencies, alphabetical order
argument_mod.o: global_parameters_mod.o

decomposition_mod.o: region_mod.o

field_mod.o: kind_params_mod.o tile_mod.o grid_mod.o \
             halo_mod.o region_mod.o gocean_mod.o

gocean_mod.o: parallel_mod.o

grid_mod.o: kind_params_mod.o region_mod.o gocean_mod.o parallel_mod.o

halo_mod.o: region_mod.o

kernel_mod.o: argument_mod.o

$(PARALLEL_MOD).o: kind_params_mod.o

parallel_mod.o: decomposition_mod.o parallel_comms_mod.o

parallel_comms_mod.o: $(PARALLEL_MOD).o

tile_mod.o: region_mod.o

%.o: $(ROOT)/%.[Ff]90
	$(F90) $(F90FLAGS) -I${FORTCL_DIR} -c $<

$(PARALLEL_MOD).o:	$(ROOT)/parallel/$(PARALLEL_MOD).f90
	$(F90) $(F90FLAGS) -c $< -o $@
