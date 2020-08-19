# dl_esm_inf
Daresbury Laboratory Earth-System Modelling Infrastructure library.

A small library to aid the creation of earth-system models. Currently
supports two-dimensional, finite-difference models in Fortran.

The first version of this library was developed to support 2D finite-
difference shallow-water models in the GOcean Project.

## Building ##

The ``finite_difference`` directory contains a Makefile which picks up
the compiler and associated flags from environment variables, e.g.:

  export F90=mpif90
  export F90FLAGS=-O2
  export AR=ar

The ``compiler_setup`` directory contains ``gnu.sh`` which may be used
to set-up these environment variables in order to use the Gnu compiler:

  . compiler_setup/gnu.sh

The `fd_lib` target builds the serial version of the library and the
`dm_fd_lib` target builds the distributed-memory (MPI) version. For
the latter you will of course need MPI installed.

## Example ##

The `finite_difference/example` directory contains a very simple example
of how to construct a model grid and associate fields with it.

## Runtime configuration ##

The following environment variables can be set up to tune runtime options:

- `DL_ESM_ALIGNMENT`: Positive integer that specifies a value by which the grid
contiguous dimension length should be divisible. If that is not the case for
the selected problem dimensions, it will add padding elements in that dimension
until the condition holds true. If not specified it defaults to 1.

## Documentation ##

The documentation is in the top-level `doc` directory and can be built
using Sphinx. Therefore, if you wish to build it you will need a
working Python installation with sphinx and sphinxfortran installed:

    pip install sphinx sphinx-fortran

``make latexpdf`` will build the PDF version of the documentation
(using latex and pdflatex) while ``make html`` builds the html
version.
