Details for C++/Kokkos/EKAT implementations:

GNU for CPU:

1. Clone, build, and install EKAT as follows:

1a. Clone:
    git clone git@github.com:E3SM-Project/EKAT.git
    cd EKAT
    git submodule update --init --recursive


1b. Configure:

    ekatsrc= # path to EKAT repo
    ekatinstall=ekat-install # or some other path
    rm -rf CMake*
    cmake \
        -D CMAKE_BUILD_TYPE:STRING=RELEASE                \
        -D CMAKE_CXX_COMPILER:STRING=mpicxx               \
        -D CMAKE_Fortran_COMPILER:STRING=mpifort          \
        -D CMAKE_INSTALL_PREFIX:PATH=$ekatinstall         \
        -D EKAT_ENABLE_TESTS:BOOL=ON                      \
        -D EKAT_DISABLE_TPL_WARNINGS:BOOL=ON              \
        -D EKAT_TEST_DOUBLE_PRECISION:BOOL=ON             \
        -D EKAT_TEST_SINGLE_PRECISION:BOOL=ON             \
        -D EKAT_TEST_MAX_THREADS:STRING=2                 \
        $ekatsrc

1c. Build and install:

    make install
    ctest # verify the build

2. In codesign-kernels/nested_loops, make a make.inc file like this one:

    $ cat make.inc
    EKAT = # path to EKAT install directory

Then

    make gnu-cpu-cke

GNU for Weaver V100:

1b. EKAT config on Weaver:

    ekatsrc= # path to EKAT repo
    ekatinstall=ekat-install
    export OMPI_CXX=${ekatsrc}/extern/kokkos/bin/nvcc_wrapper
    rm -rf CMakeFiles
    rm -f  CMakeCache.txt
    cmake \
        -C ${ekatsrc}/cmake/machine-files/weaver.cmake  \
        -D CMAKE_BUILD_TYPE:STRING=RELEASE              \
        -D CMAKE_CXX_COMPILER:STRING=mpicxx             \
        -D CMAKE_Fortran_COMPILER:STRING=mpifort        \
        -D EKAT_DISABLE_TPL_WARNINGS:BOOL=ON            \
        -D EKAT_DISABLE_TPL_WARNINGS:BOOL=ON            \
        -D EKAT_ENABLE_TESTS:BOOL=ON                    \
        -D EKAT_TEST_SINGLE_PRECISION:BOOL=ON           \
        -D EKAT_TEST_DOUBLE_PRECISION:BOOL=ON           \
        -D CMAKE_INSTALL_PREFIX:PATH=$ekatinstall       \
        ${ekatsrc}

2. make.inc:

    EKAT = # path to EKAT install directory
    CKE_PACK_SIZE = 1

Then

    make gnu-weaver-cke

GNU for Summit V100:

1b. Configure:

    ekatsrc= #...
    ekatinstall=ekat-install
    export OMPI_CXX=${ekatsrc}/extern/kokkos/bin/nvcc_wrapper
    rm -rf CMakeFiles
    rm -f  CMakeCache.txt
    cmake \
        -C ${ekatsrc}/cmake/machine-files/weaver.cmake  \
        -D CMAKE_BUILD_TYPE:STRING=RELEASE              \
        -D CMAKE_CXX_COMPILER:STRING=mpicxx             \
        -D CMAKE_Fortran_COMPILER:STRING=mpifort        \
        -D EKAT_DISABLE_TPL_WARNINGS:BOOL=ON            \
        -D EKAT_DISABLE_TPL_WARNINGS:BOOL=ON            \
        -D EKAT_ENABLE_TESTS:BOOL=ON                    \
        -D EKAT_TEST_SINGLE_PRECISION:BOOL=ON           \
        -D EKAT_TEST_DOUBLE_PRECISION:BOOL=ON           \
        -D CMAKE_INSTALL_PREFIX:PATH=$ekatinstall       \
        ${ekatsrc}

2. make.inc:

    EKAT = # ...
    CKE_PACK_SIZE = 1

Then

    make gnu-summit-cke
