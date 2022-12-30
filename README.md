# Overview

The `linalg` library package is a set of selected functions and
subroutines from [BLAS](https://netlib.org/blas/) and
[LAPACK](https://netlib.org/lapack/) as they are used by and bundled
with the [LAMMPS MD code](https://www.lammps.org/).  While the original
BLAS/LAPACK reference implementation from the Netlib servers are written
in a portable Fortran dialect, this library uses the [f2c
tool](https://netlib.org/f2c/) tools to convert them to C++ while
maintaining the expected binary interface from the most common Fortran
compilers on Linux, macOS, Windows platforms.  This allows to use it
as a drop-in replacement in case there is no compatible BLAS/LAPACK
binary library package or no Fortran compiler anavailable.  It can be
compiled with the same C++ compiler (minimum C++-11 standard) as the
LAMMPS code itself.

The purpose of this package is to provide the source for the translated
C++ files bundled with LAMMPS and easily allow to create updated
versions for new LAPACK releases or add new functionality by updating or
adding to the Fortran sources and re-running the translation to C++.

Most of the original Fortran code can be translated without change.
Those are placed into the `fortran` folder.  A few other Fortran files
need some small modifications to make them compatible with the `f2c`
translator.  In addition a few BLAS/LAPACK support functions
(e.g. DLAMCH or XERBLA) have been directly implemented in C++.  Finally,
a few functions from the `f2c` runtime library are required as well.
Those have been adapted to C++ and for typical platforms this library is
supposed to be used on.

# Version

The included Fortran files correspond to LAPACK version 3.11.0 released
on November 11th, 2022.

