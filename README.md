# dl_esm_inf
Daresbury Laboratory Earth-System Modelling Infrastructure library.

A small library to aid the creation of earth-system models. Currently
supports two-dimensional, finite-difference models in Fortran.

The first version of this library was developed to support 2D finite-
difference shallow-water models in the GOcean Project.

## Building ##

The `finite_difference` directory contains a Makefile which picks up
the compiler and associated flags from environment variables, e.g.:

  export F90=mpif90
  export F90FLAGS=-O2
  export AR=ar

The `fd_lib` target builds the serial version of the library and the
`dm_fd_lib` target builds the distributed-memory (MPI) version. For
the latter you will of course need MPI installed.

## Example ##

The `finite_difference/example` directory contains a very simple example
of how to construct a model grid and associate fields with it.

## Documentation ##

The documentation is in the top-level `doc` directory and can be built
using Sphinx. Therefore, if you wish to build it you will need a
working Python installation with sphinx and sphinxfortran installed:

    pip install sphinx sphinx-fortran

``make latexpdf`` will build the PDF version of the documentation
(using latex and pdflatex) while ``make html`` builds the html
version.

## Distributed-memory (MPI) support ##

TBD

## OpenCL Support ##

dl_esm_inf incorporates the FortCL library (github.com/stfc/FortCL) as
a git sub-module. Obtaining the source of this submodule requires that
you either supply the ``--recursive`` flag to ``git clone`` or that
you subsequently do ``git submodule init`` and ``git submodule
update``. By default a stub version of this library is compiled so
that there are no dependencies on an OpenCL installation.  The
addition of full OpenCL support to the dl_esm_inf library is the
subject of Issue #10 (github.com/stfc/dl_esm_inf/issues/10) and is
currently a work in progress.
