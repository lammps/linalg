# Translate BLAS/LAPACK subset to C++. Compile and setup tests.

cmake_minimum_required(VERSION 3.16..3.25)

project(linalg VERSION 1.0
  DESCRIPTION "BLAS/LAPACK functions translated to C++ as required by the LAMMPS MD package"
  LANGUAGES CXX)

# We require as a minimum the C++11 standard
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
endif()
if(CMAKE_CXX_STANDARD LESS 11)
  message(FATAL_ERROR "C++ standard must be set to at least 11")
endif()

# Set RelWithDebInfo as default build type.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)

########################################################################
# Check for tools
# f2c is required for Fortran to C++ translation
find_program(F2C_COMMAND NAMES f2c f2c.exe REQUIRED DOC "The f2c translater program")
# sed is required to do some automated edits of the C++ code
find_program(SED_COMMAND NAMES sed sed.exe REQUIRED DOC "The sed stream editor program")
# clang-format is required to regularize the resulting files
find_program(CLANG_FORMAT NAMES clang-format REQUIRED DOC "clang-format executable")

# We use a special compiler trick to remove comments
if(NOT (CMAKE_CXX_COMPILER_ID STREQUAL "GNU"))
  message(FATAL_ERROR "Must use GNU C++ compiler")
endif()

########################################################################
# collect Fortran source files
file(GLOB STANDARD_FORTRAN_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/fortran/[^.]*.f)
file(GLOB CUSTOM_FORTRAN_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/static/[^.]*.f)

# skip Fortran files for which we have alternatives
list(FILTER STANDARD_FORTRAN_SOURCES
  EXCLUDE REGEX "(dgetrf2|disnan|dlaisnan|dlamch|dpotrf2|lsame|xerbla|zlarft)\\.f\$")

# collect custom C++ source files
file(GLOB CUSTOM_CXX_SOURCES CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/static/[^.]*.cpp)

list(LENGTH STANDARD_FORTRAN_SOURCES NUM_STANDARD_FORTRAN)
list(LENGTH CUSTOM_FORTRAN_SOURCES NUM_CUSTOM_FORTRAN)
list(LENGTH CUSTOM_CXX_SOURCES NUM_CUSTOM_CXX)

# create list of translated C++ files
list(APPEND TRANSLATED_CXX_FILES ${STANDARD_FORTRAN_SOURCES})
list(APPEND TRANSLATED_CXX_FILES ${CUSTOM_FORTRAN_SOURCES})
list(APPEND TRANSLATED_CXX_FILES ${CUSTOM_CXX_SOURCES})
list(TRANSFORM TRANSLATED_CXX_FILES REPLACE "\\.f\$" ".cpp")
list(TRANSFORM TRANSLATED_CXX_FILES REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/(fortran|static)" "${CMAKE_CURRENT_BINARY_DIR}")

# copy custom C++ files and header
foreach(_OLD ${CUSTOM_CXX_SOURCES})
  string(REGEX REPLACE ".*static.([^.]*)\\.cpp\$" "${CMAKE_BINARY_DIR}/\\1.cpp" _NEW ${_OLD})
  configure_file(${_OLD} ${_NEW} COPYONLY)
endforeach()
configure_file(${CMAKE_SOURCE_DIR}/static/lmp_f2c.h ${CMAKE_BINARY_DIR}/lmp_f2c.h COPYONLY)

# translate Fortran files to C++
# call f2c on the Fortran file then use GCC to remove comments and clang-format to renormalize the sources
# then use sed to fix a few issues:
# - put a (char *) cast in front of all string constants
# - fix up the extern "C" statements that got messed up with the previous fix
# - change including f2c.h to lmp_f2c.h to avoid conflicts with a system f2c
# - rename calls to libf2c functions to include lmp_ so they use the provided custom functions, e.g. s_cat becomes c_lmp_cat
foreach(_F ${STANDARD_FORTRAN_SOURCES} ${CUSTOM_FORTRAN_SOURCES})
  string(REGEX REPLACE ".*(fortran|static).([^.]*)\\.f\$" "${CMAKE_CURRENT_BINARY_DIR}/\\2.cpp" _CXX ${_F})
  add_custom_command(OUTPUT ${_CXX} VERBATIM
    COMMAND ${F2C_COMMAND} -C++ -a -f < ${_F} | 
            ${CMAKE_CXX_COMPILER} -fpreprocessed -dD -P -E - | 
            ${CLANG_FORMAT} --style=file:${CMAKE_SOURCE_DIR}/.clang-format > ${_CXX}
    COMMAND ${SED_COMMAND} -i -e s,\\\(\"\[^\"\]\\\+\"\\\),\(char\ *\)\\1,g
                              -e s,^extern.*\"C\",extern\ \"C\",
                              -e s,^\#include.*\"f2c.h\",\#include\ \"lmp_f2c.h\", ${_CXX}
    COMMAND ${SED_COMMAND} -i -e s,s_\\\(cat\\\|cmp\\\|copy\\\)\(,s_lmp_\\1\(,g
                              -e s,d_\\\(sign\\\|cnjg\\\|imag\\\|lg10\\\)\(,d_lmp_\\1\(,g
                              -e s,z_\\\(abs\\\|div\\\)\(,z_lmp_\\1\(,g
                              -e s,i_\\\(len\\\|nint\\\|dnnt\\\)\(,i_lmp_\\1\(,g
                              -e s,pow_\\\(dd\\\|di\\\|ii\\\)\(,pow_lmp_\\1\(,g ${_CXX}
                              MAIN_DEPENDENCY ${_F})
endforeach()
  
add_library(linalg STATIC ${TRANSLATED_CXX_FILES} ${CMAKE_BINARY_DIR}/lmp_f2c.h)

add_custom_target(tar ALL BYPRODUCTS ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.gz
  DEPENDS ${TRANSLATED_CXX_FILES} COMMAND
  ${CMAKE_COMMAND} -E tar "czf" ${CMAKE_PROJECT_NAME}-${CMAKE_PROJECT_VERSION}.tar.gz
                        ${CMAKE_BINARY_DIR}/lmp_f2c.h ${TRANSLATED_CXX_FILES})

message("\n${NUM_STANDARD_FORTRAN} standard Fortran files, "
    "${NUM_CUSTOM_FORTRAN} custom Fortran files, "
    "${NUM_CUSTOM_CXX} custom C++ files\n")

# set up tests

include(CheckLanguage)
check_language(Fortran)
if(NOT CMAKE_Fortran_COMPILER)
  unset(ENABLE_TESTING FORCE)
  message(STATUS "Could not find working Fortran compiler. Testing not available")
else()
  option(ENABLE_TESTING "Enable or disable testing" ON)
endif()

if(ENABLE_TESTING)
  # valgrind support
  find_program(VALGRIND_BINARY NAMES valgrind)
  if(VALGRIND_BINARY)
    message(STATUS "Found valgrind at ${VALGRIND_BINARY}")
    set(MEMORYCHECK_COMMAND "${VALGRIND_BINARY}"
      CACHE FILEPATH "Memory Check Command")
    set(MEMORYCHECK_COMMAND_OPTIONS "--leak-check=full --show-leak-kinds=all --track-origins=yes"
      CACHE STRING "Memory Check Command Options")
  endif()
  # code coverage
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    option(ENABLE_COVERAGE "Enable collecting code coverage data" OFF)
    mark_as_advanced(ENABLE_COVERAGE)
    if(ENABLE_COVERAGE)
      target_compile_options(linalg PUBLIC --coverage)
      target_link_options(linalg PUBLIC --coverage)
    endif()
  endif()

  include(CTest)
  enable_testing()
  enable_language(Fortran)

  add_executable(test_blas_l3_double tests/dblat3.f)
  target_link_libraries(test_blas_l3_double PRIVATE linalg)
  add_test(NAME BLAS-L3-DOUBLE COMMAND
    ${CMAKE_COMMAND}
    -DTESTCMD=$<TARGET_FILE:test_blas_l3_double> -DTESTINP=${CMAKE_SOURCE_DIR}/tests/dblat3.in
    -DTESTOUT=${CMAKE_BINARY_DIR}/test_blas_l3_double.out -P ${CMAKE_SOURCE_DIR}/tests/runtest.cmake
  )
  set_tests_properties(BLAS-L3-DOUBLE PROPERTIES
    FAIL_REGULAR_EXPRESSION "SUSPECT")

  add_executable(test_blas_zgemm tests/zblat3.f)
  target_link_libraries(test_blas_zgemm PRIVATE linalg)
  add_test(NAME BLAS-ZGEMM COMMAND
    ${CMAKE_COMMAND}
    -DTESTCMD=$<TARGET_FILE:test_blas_zgemm> -DTESTINP=${CMAKE_SOURCE_DIR}/tests/zblat3.in
    -DTESTOUT=${CMAKE_BINARY_DIR}/test_blas_zgemm.out -P ${CMAKE_SOURCE_DIR}/tests/runtest.cmake
  )
  set_tests_properties(BLAS-ZGEMM PROPERTIES
    FAIL_REGULAR_EXPRESSION "SUSPECT")
endif()
