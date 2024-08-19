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
versions for new LAPACK releases or add new functionality if required by
added LAMMPS functionality through updating or adding to the Fortran
sources and re-running the translation to C++.

Most of the original Fortran code can be translated without change.
Those are placed into the `fortran` folder.  A few other Fortran files
need some small modifications to make them compatible with the `f2c`
translator.  In addition a few BLAS/LAPACK support functions
(e.g. DLAMCH or XERBLA) have been directly implemented in C++.  Finally,
a few functions from the `f2c` runtime library are required as well.
Those have been adapted to C++ and for typical platforms this library is
supposed to be used on.

# Version

The included Fortran files correspond to LAPACK version 3.12.0 released
on November 24th, 2023.

# License

Since the bulk of the code is an automated translation of the BLAS/LAPACK
sources, the same licensing terms apply.  Additional or replacement code
is copyright (c) 2022, 2023 by Axel Kohlmeyer `<akohlmey@gmail.com>`.

# Installation

This package uses CMake to set up the translation, and compilation plus
a few tests.  CMake version 3.16 or later is required. Also required are
f2c, the GNU C++ compiler, clang-format, and sed.  This has only been
tested on a Linux machine.

Running the command:

```
cmake -S . -B build
```

Will set up a build folder and

```
cmake --build build
```

will do the translation and also try to compile the generated C++ files into a library
as well as create a compressed tar file with the translated sources in the build folder.

If a suitable Fortran compiler was found, also a couple of tests for BLAS functions
are configured and can be run with ctest.

```
ctest --test-dir build
```
