all: model

.PHONY: model
model: inf_lib model.exe

INF_DIR=..
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/dl_esm_inf_fd.a

.PHONY: inf_lib
inf_lib:
	${MAKE} -C ${INF_DIR} F90=${F90} dm_fd_lib

.PHONY: inf_lib_serial
inf_lib_serial:
	${MAKE} -C ${INF_DIR} F90=${F90} fd_lib

# Normal targets
model.exe: inf_lib model.o
	${F90} -o $@ ${LDFLAGS} model.o ${INF_LIB}

%.o: %.f90
	$(F90) $(F90FLAGS) ${OMPFLAGS} -I${INF_INC} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) ${OMPFLAGS} -I${INF_INC} -c $<

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~ *.exe
