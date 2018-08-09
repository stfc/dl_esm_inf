# The name of the library (archive) this Makefile creates.
API_LIB = dl_esm_inf_fd.a
# FortCL is a git submodule of the dl_esm_inf repository.
FORTCL_DIR = ../external/FortCL

.PHONY: all fd_lib fd_lib_cl install clean allclean distclean fortcl fortcl_stub

all: fd_lib

# Default is to build without OpenCL support
fd_lib: fortcl_stub
	${MAKE} -C src API_LIB="${API_LIB}"

# Target requiring OpenCL installation
fd_lib_cl: fortcl
	${MAKE} -C src API_LIB="${API_LIB}"

install: fd_lib lib include
	cp src/${API_LIB} ./lib/.
	cp src/*.[Mm][Oo][Dd] ./include/.

fortcl_stub:
	${MAKE} -C ${FORTCL_DIR} stub

fortcl:
	${MAKE} -C ${FORTCL_DIR}

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

