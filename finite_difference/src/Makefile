FORTCL_DIR = ../../external/FortCL/src

MODULES = argument_mod.o \
          kind_params_mod.o \
          global_parameters_mod.o \
          gocean_mod.o \
          grid_mod.o \
          kernel_mod.o \
          halo_mod.o \
          field_mod.o \
          region_mod.o

all: ${API_LIB}

# Create the archive. Note that we include the FortCL objects in this
# archive but take care to exclude clfortran.o (as that only has
# interface definitions, not implementation). We could use filter-out
# for this but that is specific to Gnu Make.
${API_LIB}: ${MODULES}
	${AR} ${ARFLAGS} ${API_LIB} ${FORTCL_DIR}/*.o ${MODULES}
	${AR} d ${API_LIB} clfortran.o

clean:
	rm -f *.o *.mod *.a

# Inter-module dependencies, alphabetical order
argument_mod.o: global_parameters_mod.o

field_mod.o: kind_params_mod.o grid_mod.o \
             halo_mod.o region_mod.o gocean_mod.o

grid_mod.o: kind_params_mod.o region_mod.o gocean_mod.o

halo_mod.o: region_mod.o

%.o: %.[fF]90
	$(F90) $(F90FLAGS) -I${FORTCL_DIR} -c $<
