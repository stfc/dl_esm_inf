#------------------------------------------------------------------------------
# BSD 2-Clause License
# 
# Copyright (c) 2018, Science and Technology Facilities Council.
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
# Author: A. R. Porter, STFC Daresbury Laboratory

# Makefile for the example mini-app demonstrating the use of the dl_esm_inf
# library.
#
# This Makefile picks up the compiler to use plus any associated flags from
# environment variables. e.g. to use gfortran:
# 
# export F90=mpif90
# export F90FLAGS="-O3"
# export AR=ar
# make
#
# or, to build the serial version:
#
# export F90=gfortran
# make serial-model

all: model

# Target to build version with MPI support
.PHONY: model
model: inf_lib model.exe

# Target to build serial version
.PHONY: serial-model
serial-model: inf_lib_serial model.exe

# Run tests (check both serial and parallel versions are valgrind clean)
# This is a very basic target and assumes that valgrind and 'mpirun'
# are available.
# TODO use a proper testing framework!
test: serial-test

serial-test: serial-model
	${MAKE} clean
	${MAKE} serial-model
	valgrind ./model.exe

parallel-test: model
	${MAKE} clean
	${MAKE} model
	mpirun -np 2 valgrind ./model.exe

INF_DIR=..
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/lib_fd.a

.PHONY: inf_lib
inf_lib:
	${MAKE} -C ${INF_DIR} F90="${F90}" dm_fd_lib

.PHONY: inf_lib_serial
inf_lib_serial:
	${MAKE} -C ${INF_DIR} F90="${F90}" fd_lib

# Normal targets
model.exe: model.o
	${F90} -o $@ ${LDFLAGS} model.o ${INF_LIB}

%.o: %.f90
	$(F90) $(F90FLAGS) ${OMPFLAGS} -I${INF_INC} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) ${OMPFLAGS} -I${INF_INC} -c $<

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~ *.exe

allclean: clean
	${MAKE} -C ${INF_DIR} distclean
