The C++ files in this folder are either direct C++ implementations of
their Fortran equivalents using the C++ runtime, or they are adapted
copies of functions from the libf2c runtime.  The runtime functions
needed to be renamed to avoid conflics with libgfortran which uses
some of the same function names.  Also the header file `f2c.h` was
renamed to `lmp_f2c.h` to avoid conflicts.

The Fortran files in this folder are modified from their
original versions, so that f2c can correctly translate them.
