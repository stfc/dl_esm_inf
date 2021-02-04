#------------------------------------------------------------------------------
# BSD 2-Clause License
# 
# Copyright (c) 2021, Science and Technology Facilities Council.
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
# Author: S. Siso, STFC Daresbury Laboratory

# Makefile for dist_mem tests for the dl_esm_inf finite difference library.
#
# This Makefile picks up the compiler to use plus any associated flags from
# environment variables. e.g. to use gfortran:
# 
# export F90=mpif90
# export MPIRUN=mpirun
# export F90FLAGS="-O3"
# export AR=ar
# make
#
all: tests

# Location of the dl_esm_inf library
INF_DIR=../..
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/lib_fd.a

.PHONY: test_device_io inf_lib clean

# Compile and run tests. This is a very basic target.
tests: test_device_io

test_device_io: inf_lib test_device_io.exe
	./test_device_io.exe

inf_lib:
	${MAKE} -C ${INF_DIR} F90="${F90}" fd_lib

# Normal targets
test_device_io.exe: inf_lib test_device_io.o
	${F90} -o $@ ${LDFLAGS} test_device_io.o ${INF_LIB}

%.o: %.[fF]90
	$(F90) $(F90FLAGS) -I${INF_INC} -c $<

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~ *.exe gnu_opt_report.txt *.optrpt
