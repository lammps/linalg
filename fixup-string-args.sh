#!/bin/sh -vx
# this command replaces typical long string with 
# a single character same as the LAPACK subroutine requires.
# In Fortran there is no difference, but in C++ it can matter.

if [ $# -lt 1 ] ; then \
        echo "usage: $0 <folder with fortran sources>"; \
        exit 1;\
fi

exec sed -i                            \
    -e "s/'All'/'A'/g"                 \
    -e "s/'ALL'/'A'/g"                 \
    -e "s/'Right'/'R'/g"               \
    -e "s/'RIGHT'/'R'/g"               \
    -e "s/'Left'/'L'/g"                \
    -e "s/'Lower'/'L'/g"               \
    -e "s/'Upper'/'U'/g"               \
    -e "s/'No transpose'/'N'/g"        \
    -e "s/'No Transpose'/'N'/g"        \
    -e "s/'NO TRANSPOSE'/'T'/g"        \
    -e "s/'Unit'/'U'/g"                \
    -e "s/'UNIT'/'U'/g"                \
    -e "s/'Non-unit'/'N'/g"            \
    -e "s/'NON-UNIT'/'N'/g"            \
    -e "s/'Conjugate transpose'/'C'/g" \
    -e "s/'Conjugate'/'C'/g"           \
    -e "s/'Full'/'F'/g"                \
    -e "s/'Transpose'/'T'/g"           \
    -e "s/'Backward'/'B'/g"            \
    -e "s/'Forward'/'F'/g"             \
    -e "s/'Columnwise'/'C'/g"          \
    -e "s/'Rowwise'/'R'/g"             \
    "$1"/*.f
