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
# Author: A. R. Porter, STFC Daresbury Laboratory

# Top-level Makefile for the dl_esm_inf library.
# Default target is to build the library without MPI. Use the 'dm_fd_lib'
# target to compile with MPI support or 'cl_fd_lib' for OpenCL support.
#
# This Makefile picks up the compiler to use along with any flags from
# environment variables. e.g. to use gfortran:
# 
# export F90=mpif90
# export F90FLAGS="-O3"
# export AR=ar
#
# or, to build the serial version:
#
# export F90=gfortran

API_LIB ?= lib_fd.a
F90 ?= gfortran

.PHONY: doc clean install fd_lib dm_fd_lib cl_fd_lib all

all: fd_lib

# Build the library without MPI or OpenCL
fd_lib:
	${MAKE} -C src API_LIB="${API_LIB}" MPI=no OCL=no

# Build the library with MPI but no OpenCL
dm_fd_lib:
	${MAKE} -C src API_LIB="${API_LIB}" MPI=yes OCL=no

# OpenCL without MPI
cl_fd_lib:
	${MAKE} -C src API_LIB="${API_LIB}" MPI=no OCL=yes

install: fd_lib lib include
	cp src/${API_LIB} ./lib/.
	cp src/*.[Mm][Oo][Dd] ./include/.

lib:
	mkdir lib

include:
	mkdir include

clean:
	${MAKE} -C src clean

allclean: clean
	rm -f lib/${API_LIB}
	rm -f include/*.[Mm][Oo][Dd]

distclean: allclean

doc:
	doxygen src/doxyfile
