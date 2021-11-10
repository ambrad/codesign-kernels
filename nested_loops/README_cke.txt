Details for C++/Kokkos/EKAT implementations:

GNU for CPU:

1. Clone, build, and install EKAT as follows:

1a. Clone: git clone git@github.com:E3SM-Project/EKAT.git

1b. Configure:

    ekatsrc= # path to EKAT repo
    ekatinstall=ekat-install # or some other path
    rm -rf CMake*
    cmake \
        -D CMAKE_BUILD_TYPE:STRING=DEBUG                  \
        -D CMAKE_CXX_COMPILER:STRING=mpicxx               \
        -D CMAKE_Fortran_COMPILER:STRING=mpifort          \
        -D CMAKE_INSTALL_PREFIX:PATH=$ekatinstall         \
        -D EKAT_ENABLE_TESTS:BOOL=ON                      \
        -D EKAT_TEST_DOUBLE_PRECISION:BOOL=ON             \
        -D EKAT_TEST_SINGLE_PRECISION:BOOL=ON             \
        -D EKAT_TEST_MAX_THREADS:STRING=2                 \
        $ekatsrc

1c. Build and install:

    make install

2. In codesign-kernels/nested_loops, make a make.inc file like this one:

    $ cat make.inc
    EKAT = # path to EKAT install directory

Then

    make gnu-cpu-cke

GNU for V100:

1b.

2. make gnu-gpu-cke
