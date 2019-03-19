# dl_esm_inf
Daresbury Laboratory Earth-System Modelling Infrastructure library.

A small library to aid the creation of earth-system models. Currently
supports two-dimensional, finite-difference models in Fortran.

The first version of this library was developed to support 2D finite-
difference shallow-water models in the GOcean Project.

## Example ##

The `finite_difference/example` directory contains a very simple example
of how to construct a model grid and associate fields with it.

## Concepts ##

dl_esm_inf provides certain types of object from which an earth-system
model may then be constructed

### Grid ###

The grid of points making up the simulation domain.

### Field ###

A field is a representation of some physical quantity (e.g. vorticity)
at points on the grid.

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
