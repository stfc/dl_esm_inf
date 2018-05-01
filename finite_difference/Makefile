API_LIB = dl_esm_inf_fd.a

all: fd_lib

fd_lib:
	ln -sf src/parallel/parallel_stub_mod.f90 src/parallel_mod.f90
	${MAKE} -C src API_LIB="${API_LIB}"

dm_fd_lib:
	ln -sf src/parallel/parallel_mod.f90 src/parallel_mod.f90
	${MAKE} -C src API_LIB="${API_LIB}" dm_fd_lib

install: fd_lib lib include
	cp src/${API_LIB} ./lib/.
	cp src/*.[Mm][Oo][Dd] ./include/.

lib:
	mkdir lib

include:
	mkdir include

clean:
	rm src/parallel_mod.f90
	${MAKE} -C src clean

distclean: clean
	rm -f lib/${API_LIB}
	rm -f include/*.[Mm][Oo][Dd]
