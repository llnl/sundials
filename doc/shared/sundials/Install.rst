..
   Programmer(s): Daniel R. Reynolds @ UMBC
   ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _Installation:

Installing SUNDIALS
===================

In this chapter we discuss two ways for building and installing SUNDIALS from
source. The first is with the `Spack <https://spack.io/>`__ HPC package manager
and the second is with `CMake <https://cmake.org/>`__.

.. _Installation.Spack:

Installing with Spack
---------------------

Spack is a package management tool that provides a simple spec syntax to
configure and install software on a wide variety of platforms and environments.
See the `Getting Started
<https://spack.readthedocs.io/en/latest/getting_started.html>`__ section in the
Spack documentation for more information on installing Spack.

Once Spack is setup on your system, the default SUNDIALS configuration can be
install with the command

.. code-block:: bash

   spack install sundials

Additional options can be enabled through Spack package variants. For
information on the available variants visit the `SUNDIALS Spack package
<https://packages.spack.io/package.html?name=sundials>`__ web page or use the
command

.. code-block:: bash

   spack info sundials

.. _Installation.CMake:

Installing with CMake
---------------------

CMake provides a platform-independent build system capable of generating Unix
and Linux Makefiles, as well as KDevelop, Visual Studio, and (Apple) XCode
project files from the same configuration file. A GUI front end is also
available allowing for an interactive build and installation process.

At a minimum, building SUNDIALS requires CMake version 3.18.0 or higher and a
working C compiler. If a compatible version of CMake is not already installed on
your system, source files or pre-built binary files can be obtained from the
`CMake Download website <https://cmake.org/download/>`__.

When building with CMake, you will need to obtain the SUNDIALS source code. You
can get the source files by either cloning the `SUNDIALS GitHub repository
<https://github.com/LLNL/sundials>`__ with the command

.. code-block:: bash

   git clone https://github.com/LLNL/sundials

or by downloading release compressed archives (``.tar.gz`` files) from the
`SUNDIALS download website
<https://computing.llnl.gov/projects/sundials/sundials-software>`__. The
compressed archives allow for downloading the entire SUNDIALS suite or
individual packages. The name of the distribution archive is of the form
``SOLVER-a.b.c.tar.gz``, where ``SOLVER`` is one of: ``sundials``, ``cvode``,
``cvodes``, ``arkode``, ``ida``, ``idas``, or ``kinsol``, and ``a.b.c``
represents the version number of the SUNDIALS suite or of the individual
package. After downloading the relevant archives, uncompress and expand the
sources. For example, if you downloaded ``sundials-x.y.z.tar.gz``, running the
command

.. code-block:: bash

   tar -zxf sundials-x.y.z.tar.gz

will extract the source files under the ``sundials-x.y.z`` directory.

In the installation steps below we will refer to the following directories:

* ``SOLVER_DIR`` is the ``sundials`` directory created when cloning from GitHub
  or the ``SOLVER-a.b.c`` directory created after uncompressing the release
  archive.

* ``BUILD_DIR`` is the (temporary) directory under which SUNDIALS is built.
  In-source builds are prohibited; the build directory ``BUILD_DIR`` can **not**
  be the same as ``SOLVER_DIR`` and such an attempt will lead to an error. This
  prevents "polluting" the source tree, simplifies building with different
  configurations and/or options, and makes it easy to clean-up all traces of the
  build by simply removing the build directory.

* ``INSTALL_DIR`` is the directory under which the SUNDIALS exported header
  files and libraries will be installed. The installation directory
  ``INSTALL_DIR`` can not be the same as the ``SOLVER_DIR`` directory.
  Typically, header files are exported under a directory ``INSTALL_DIR/include``
  while libraries are typically installed under ``INSTALL_DIR/lib`` or
  ``INSTALL_LIB/lib64``, with ``INSTALL_DIR`` specified at configuration time.

.. _Installation.CMake.Unix:

Linux/Unix systems
^^^^^^^^^^^^^^^^^^

CMake can be used from the command line with the ``cmake`` command, or from
graphical interfaces with the ``ccmake`` or ``cmake-gui`` commands. Below we
present the installation steps using the command line interface.

Using CMake from the command line is simply a matter of generating the build
files for the desired configuration, building, and installing. For example, the
following commands will build and install the default configuration:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR
   cd BUILD_DIR
   make
   make install

The default configuration will install static and shared libraries for all
SUNDIALS packages and install the associated example codes. Additional features
can be enabled by specifying more options in the configuration step. For
example, to enable MPI add ``-D ENABLE_MPI=ON`` to the ``cmake`` command above:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_MPI=ON

See section :numref:`Installation.Options` below for a complete list of SUNDIALS
configuration options and additional configuration examples.

.. _Installation.CMake.Windows:

Windows Systems
^^^^^^^^^^^^^^^

CMake can also be used to build SUNDIALS on Windows. To build SUNDIALS for use
with Visual Studio the following steps should be performed:

#. Create a separate ``BUILD_DIR``

#. Open a Visual Studio Command Prompt and cd to ``BUILD_DIR``

#. Run ``cmake-gui ../SOLVER_DIR``

   a. Hit Configure

   b. Check/Uncheck solvers to be built

   c. Change ``CMAKE_INSTALL_PREFIX`` to ``INSTALL_DIR``

   d. Set other options as desired (see section :numref:`Installation.Options`)

   e. Hit Generate

#. Back in the VS Command Window:

   a. Run ``msbuild ALL_BUILD.vcxproj``

   b. Run ``msbuild INSTALL.vcxproj``

The resulting libraries will be in the ``INSTALL_DIR``.

The SUNDIALS project can also now be opened in Visual Studio.  Double click on
the ``ALL_BUILD.vcxproj`` file to open the project.  Build the whole *solution*
to create the SUNDIALS libraries.  To use the SUNDIALS libraries in your own
projects, you must set the include directories for your project, add the
SUNDIALS libraries to your project solution, and set the SUNDIALS libraries as
dependencies for your project.

.. _Installation.CMake.HPC:

HPC Clusters
^^^^^^^^^^^^

This section is a guide for installing SUNDIALS on specific HPC clusters.  In
general, the procedure is the same as described previously in
:numref:`Installation.CMake.Unix` for Unix/Linux machines. The main differences
are in the modules and environment variables that are specific to different HPC
clusters. We aim to keep this section as up to date as possible, but it may lag
the latest software updates to each cluster.

Frontier
""""""""

`Frontier <https://www.olcf.ornl.gov/frontier/>`__ is an Exascale supercomputer
at the Oak Ridge Leadership Computing Facility. If you are new to this system,
then we recommend that you review the `Frontier user guide
<https://docs.olcf.ornl.gov/systems/frontier_user_guide.html>`__.

**A Standard Installation**

Load the modules and set the environment variables needed to build SUNDIALS.
This configuration enables both MPI and HIP support for distributed and GPU
parallelism. It uses the HIP compiler for C and C++ and the Cray Fortran
compiler. Other configurations are possible.

.. code-block:: bash

   # required dependencies
   module load PrgEnv-cray-amd/8.5.0
   module load craype-accel-amd-gfx90a
   module load rocm/5.3.0
   module load cmake/3.23.2

   # GPU-aware MPI
   export MPICH_GPU_SUPPORT_ENABLED=1

   # compiler environment hints
   export CC=$(which hipcc)
   export CXX=$(which hipcc)
   export FC=$(which ftn)
   export CFLAGS="-I${ROCM_PATH}/include"
   export CXXFLAGS="-I${ROCM_PATH}/include -Wno-pass-failed"
   export LDFLAGS="-L${ROCM_PATH}/lib -lamdhip64 ${PE_MPICH_GTL_DIR_amd_gfx90a} -lmpi_gtl_hsa"

Now we can build SUNDIALS. In general, this is the same procedure described in
the previous sections. The following command builds and installs SUNDIALS with
MPI, HIP, and the Fortran interface enabled, where ``<account>`` is your
allocation account on Frontier:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D AMDGPU_TARGETS=gfx90a \
     -D ENABLE_HIP=ON \
     -D ENABLE_MPI=ON \
     -D BUILD_FORTRAN_MODULE_INTERFACE=ON
   cd BUILD_DIR
   make -j8 install
   # Need an allocation to run the tests:
   salloc -A <account> -t 10 -N 1 -p batch
   make test
   make test_install_all

.. _Installation.Options:

Configuration options
---------------------

All available SUNDIALS CMake options are described in the sections below. The
default values for some options (e.g., compiler flags and installation paths)
are for a Linux system and are provided as illustration only.

.. note::

   When using a CMake graphical interface (``ccmake`` or ``cmake-gui``),
   multiple configuration passes are performed before generating the build
   files. For options where the default value depends on the value of another
   option, the initial value is set on the first configuration pass and is not
   updated automatically if the related option value is changed in subsequent
   passes. For example, the default value of :cmakeop:`EXAMPLES_INSTALL_PATH` is
   ``CMAKE_INSTALL_PREFIX/examples``; if the value of
   :cmakeop:`CMAKE_INSTALL_PREFIX` is updated, then
   :cmakeop:`EXAMPLES_INSTALL_PATH` will also need to be updated as its value
   was set using the :cmakeop:`CMAKE_INSTALL_PREFIX` default.

.. _Installation.Options.BuildType:

Build Type
^^^^^^^^^^

The build type determines the level of compiler optimization, if debug
information is included, and if additional error checking code is generated. The
provided build types are:

* ``Debug`` -- no optimization flags, debugging information included, additional
  error checking enabled

* ``Release`` -- high optimization flags, no debugging information, no
  additional error checks

* ``RelWithDebInfo`` -- high optimization flags, debugging information included,
  no additional error checks

* ``MinSizeRel`` -- minimize size flags, no debugging information, no additional
  error checks

Each build type has a corresponding option for the set of compiler flags that
are appended to the user-specified compiler flags. See section
:numref:`Installation.Options.Compilers` for more information.

.. cmakeoption:: CMAKE_BUILD_TYPE

   Choose the type of build for single-configuration generators (e.g., Makefiles
   or Ninja).

   Default: ``RelWithDebInfo``

.. cmakeoption:: CMAKE_CONFIGURATION_TYPES

   Specifies the build types for multi-config generators (e.g. Visual Studio,
   Xcode, or Ninja Multi-Config) as a semicolon-separated list.

   Default: ``Debug``, ``Release``, ``RelWithDebInfo``, and ``MinSizeRel``

.. _Installation.Options.Compilers:

Compilers and Compiler Flags
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Building SUNDIALS requires a C compiler that supports at least a subset of the
C99 standard (specifically those features implemented by Visual Studio 2015).

Additional SUNDIALS features that interface with external C++ libraries or GPU
programming models require a C++ compiler (e.g., CUDA, HIP, SYCL, Ginkgo,
Trilinos, etc.). The C++ standard required depends on the particular library or
programming model used and is noted with the relevant options below. The C++
convenience classes provided by SUNDIALS require C++14 or newer. C++
applications that require an earlier C++ standard should use the SUNDIALS C API.

When enabling the SUNDIALS Fortran interfaces, a Fortran compiler that supports
the Fortran 2003 or newer standard is required in order to utilize the
``ISO_C_BINDING`` module.

C Compiler
""""""""""

.. cmakeoption:: CMAKE_C_COMPILER

   The full path to the C compiler

   Default: CMake will attempt to automatically locate a C compiler on the
   system (e.g., from the ``CC`` environment variable or common installation
   paths).

.. cmakeoption:: CMAKE_C_FLAGS

   User-specified flags for the C compiler. The value of this option should be a
   string with flags separated by spaces.

   Default: Initialized by the ``CFLAGS`` environment variable.

.. cmakeoption:: CMAKE_C_FLAGS_DEBUG

   C compiler flags appended when the :cmakeop:`CMAKE_BUILD_TYPE` is ``Debug``

   Default: ``-g``

.. cmakeoption:: CMAKE_C_FLAGS_RELEASE

   C compiler flags appended when the :cmakeop:`CMAKE_BUILD_TYPE` is ``Release``

   Default: ``-O3 -DNDEBUG``

.. cmakeoption:: CMAKE_C_FLAGS_RELWITHDEBINFO

   C compiler flags appended when the :cmakeop:`CMAKE_BUILD_TYPE` is
   ``RelWithDebInfo``

   Default: ``-O2 -g -DNDEBUG``

.. cmakeoption:: CMAKE_C_FLAGS_MINSIZEREL

   C compiler flags appended when the :cmakeop:`CMAKE_BUILD_TYPE` is
   ``MinSizeRel``

   Default: ``-Os -DNDEBUG``

.. cmakeoption:: CMAKE_C_STANDARD

   The C standard used when building SUNDIALS C source files.

   Default: ``99``

   Options: ``99``, ``11``, ``17``, or ``23``

.. cmakeoption:: CMAKE_C_EXTENSIONS

   Enable compiler specific C extensions.

   Default: ``ON``

C++ Compiler
""""""""""""

.. cmakeoption:: CMAKE_CXX_COMPILER

   The full path to the C++ compiler

   Default: CMake will attempt to automatically locate a C++ compiler on the
   system (e.g., from the ``CXX`` environment variable or common installation
   paths).

.. cmakeoption:: CMAKE_CXX_FLAGS

   User-specified flags for the C++ compiler. The value of this option should be
   a string with flags separated by spaces.

   Default: Initialized by the ``CXXFLAGS`` environment variable.

.. cmakeoption:: CMAKE_CXX_FLAGS_DEBUG

   C++ compiler flags appended when the :cmakeop:`CMAKE_BUILD_TYPE` is ``Debug``

   Default: ``-g``

.. cmakeoption:: CMAKE_CXX_FLAGS_RELEASE

   C++ compiler flags appended when the :cmakeop:`CMAKE_BUILD_TYPE` is
   ``Release``

   Default: ``-O3 -DNDEBUG``

.. cmakeoption:: CMAKE_CXX_FLAGS_RELWITHDEBINFO

   C++ compiler flags appended when the :cmakeop:`CMAKE_BUILD_TYPE` is
   ``RelWithDebInfo``

   Default: ``-O2 -g -DNDEBUG``

.. cmakeoption:: CMAKE_CXX_FLAGS_MINSIZEREL

   C++ compiler flags appended when the :cmakeop:`CMAKE_BUILD_TYPE` is
   ``MinSizeRel``

   Default: ``-Os -DNDEBUG``

.. cmakeoption:: CMAKE_CXX_STANDARD

   The C++ standard used when building SUNDIALS C++ source files.

   Default: ``14`` or ``17`` if :cmakeop:`ENABLE_GINKGO` or
   :cmakeop:`ENABLE_SYCL` are ``ON``

   Options: ``14``, ``17``, ``20``, or ``23``

.. cmakeoption:: CMAKE_CXX_EXTENSIONS

   Enable compiler specific C++ extensions.

   Default: ``ON``

Fortran Compiler
""""""""""""""""

.. cmakeoption:: CMAKE_Fortran_COMPILER

   The full path to the Fortran compiler

   Default: CMake will attempt to automatically locate a Fortran compiler on the
   system (e.g., from the ``FC`` environment variable or common installation
   paths).

.. cmakeoption:: CMAKE_Fortran_FLAGS

   User-specified flags for the Fortran compiler. The value of this option
   should be a string with flags separated by spaces.

   Default: Initialized by the ``FFLAGS`` environment variable.

.. cmakeoption:: CMAKE_Fortran_FLAGS_DEBUG

   Fortran compiler flags appended when the :cmakeop:`CMAKE_BUILD_TYPE` is
   ``Debug``

   Default: ``-g``

.. cmakeoption:: CMAKE_Fortran_FLAGS_RELEASE

   Fortran compiler flags appended when the :cmakeop:`CMAKE_BUILD_TYPE` is
   ``Release``

   Default: ``-O3``

.. cmakeoption:: CMAKE_Fortran_FLAGS_RELWITHDEBINFO

   Fortran compiler flags appended when the :cmakeop:`CMAKE_BUILD_TYPE` is
   ``RelWithDebInfo``

   Default: ``-O2 -g``

.. cmakeoption:: CMAKE_Fortran_FLAGS_MINSIZEREL

   Fortran compiler flags appended when the :cmakeop:`CMAKE_BUILD_TYPE` is
   ``MinSizeRel``

   Default: ``-Os``

.. _Installation.Options.InstallLocation:

Install Location
^^^^^^^^^^^^^^^^

Use the following options to set where the SUNDIALS headers, library, and CMake
configuration files will be installed.

.. cmakeoption:: CMAKE_INSTALL_PREFIX

   Install path prefix (``INSTALL_DIR``), prepended onto install directories

   Default: ``/usr/local``

   .. note::

      The user must have write access to the location specified through this
      option. Exported SUNDIALS header files and libraries will be installed
      under subdirectories ``include`` and :cmakeop:`CMAKE_INSTALL_LIBDIR` of
      :cmakeop:`CMAKE_INSTALL_PREFIX`, respectively.

.. cmakeoption:: CMAKE_INSTALL_LIBDIR

   The directory under :cmakeop:`CMAKE_INSTALL_PREFIX` where libraries will be
   installed

   Default: Set based on the system as ``lib``, ``lib64``, or
   ``lib/<multiarch-tuple>``

.. cmakeoption:: SUNDIALS_INSTALL_CMAKEDIR

   The directory under :cmakeop:`CMAKE_INSTALL_PREFIX` where the SUNDIALS CMake
   package configuration files will be installed (see section
   :numref:`Installation.CMakeConfigFile` for more information)

   Default: ``CMAKE_INSTALL_LIBDIR/cmake/sundials``

.. _Installation.Options.LibraryTypes:

Shared and Static Libraries
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the following options to set which types of libraries will be installed. By
default both static and shared libraries are installed.

.. cmakeoption:: BUILD_SHARED_LIBS

   Build shared libraries

   Default: ``ON``

.. cmakeoption:: BUILD_STATIC_LIBS

   Build static libraries

   Default: ``ON``

.. _Installation.Options.IndexSize:

Index Size
^^^^^^^^^^

.. cmakeoption:: SUNDIALS_INDEX_SIZE

   The integer size (in bits) used for indices in SUNDIALS (e.g., for vector and
   matrix entries), options are: ``32`` or ``64``

   Default: ``64``

   .. note::

      The build system tries to find an integer type of the appropriate
      size. Candidate 64-bit integer types are (in order of preference):
      ``int64_t``, ``__int64``, ``long long``, and ``long``. Candidate 32-bit
      integers are (in order of preference): ``int32_t``, ``int``, and ``long``.
      The advanced option, :cmakeop:`SUNDIALS_INDEX_TYPE` can be used to provide
      a type not listed here.

.. cmakeoption:: SUNDIALS_INDEX_TYPE

   The integer type used for SUNDIALS indices. The type size must match the size
   provided in the :cmakeop:`SUNDIALS_INDEX_SIZE` option.

   Default: Automatically determined based on :cmakeop:`SUNDIALS_INDEX_SIZE`

   .. versionchanged:: 3.2.0

      In prior versions, this option could be set to ``INT64_T`` to use 64-bit
      integers or ``INT32_T`` to use 32-bit integers. These special values are
      deprecated and a user will only need to use the
      :cmakeop:`SUNDIALS_INDEX_SIZE` option in most cases.


.. .. cmakeoption:: SUNDIALS_COUNTER_TYPE

..    The integer type used for SUNDIALS counters.

..    Default: `int64_t`

..    .. versionadded:: 7.3.0


.. _Installation.Options.Precision:

Precision
^^^^^^^^^

.. cmakeoption:: SUNDIALS_PRECISION

   The floating-point precision used in SUNDIALS packages and class
   implementations, options are: ``single``, ``double``, or ``extended``

   Default: ``double``

.. _Installation.Options.MathLibrary:

Math Library
^^^^^^^^^^^^

.. cmakeoption:: SUNDIALS_MATH_LIBRARY

   The standard C math library (e.g., ``libm``) to link with.

   Default: ``-lm`` on Unix systems, none otherwise

.. _Installation.Options.Packages:

SUNDIALS Packages
^^^^^^^^^^^^^^^^^

The following options can be used to enable/disable particular SUNDIALS
packages.

.. cmakeoption:: BUILD_ARKODE

   Build the ARKODE library

   Default: ``ON``

.. cmakeoption:: BUILD_CVODE

   Build the CVODE library

   Default: ``ON``

.. cmakeoption:: BUILD_CVODES

   Build the CVODES library

   Default: ``ON``

.. cmakeoption:: BUILD_IDA

   Build the IDA library

   Default: ``ON``

.. cmakeoption:: BUILD_IDAS

   Build the IDAS library

   Default: ``ON``

.. cmakeoption:: BUILD_KINSOL

   Build the KINSOL library

   Default: ``ON``

.. _Installation.Options.Examples:

Example Programs
^^^^^^^^^^^^^^^^

.. cmakeoption:: EXAMPLES_ENABLE_C

   Build the SUNDIALS C examples

   Default: ``ON``

.. cmakeoption:: EXAMPLES_ENABLE_CXX

   Build the SUNDIALS C++ examples

   Default: ``OFF``

.. cmakeoption:: EXAMPLES_ENABLE_CUDA

   Build the SUNDIALS CUDA examples

   Default: ``ON`` when :cmakeop:`ENABLE_CUDA` is ``ON``, otherwise ``OFF``

.. cmakeoption:: EXAMPLES_ENABLE_F2003

   Build the SUNDIALS Fortran 2003 examples

   Default: ``ON`` when :cmakeop:`BUILD_FORTRAN_MODULE_INTERFACE` is ``ON``,
   otherwise ``OFF``

.. cmakeoption:: EXAMPLES_INSTALL

   Install example program source files and sample output files. See
   :cmakeop:`EXAMPLES_INSTALL_PATH` for the install location.

   A ``CMakeLists.txt`` file to build the examples will be automatically
   generated and installed with the source files. If building on a Unix-like
   system, a ``Makefile`` for compiling the installed example programs will be
   also generated and installed.

   Default: ``ON``

.. cmakeoption:: EXAMPLES_INSTALL_PATH

   Full path to where example source and output files will be installed

   Default: ``CMAKE_INSTALL_PREFIX/examples``

.. _Installation.Options.Fortran:

Fortran Interfaces
^^^^^^^^^^^^^^^^^^

.. cmakeoption:: BUILD_FORTRAN_MODULE_INTERFACE

   Build the SUNDIALS Fortran 2003 interface

   Default: ``OFF``

   .. note::

      The Fortran interface are only compatible with double precision (i.e.,
      :cmakeop:`SUNDIALS_PRECISION` must be ``double``).

   .. warning::

      There is a known issue with MSYS/gfortran and SUNDIALS shared libraries
      that causes linking the Fortran interfaces to fail when building
      SUNDIALS. For now the work around is to only build with static libraries
      when using MSYS with gfortran on Windows.

.. _Installation.Options.ErrorChecking:

Error Checking
^^^^^^^^^^^^^^

For more information on error handling in SUNDIALS, see
:ref:`SUNDIALS.Errors`.

.. cmakeoption:: SUNDIALS_ENABLE_ERROR_CHECKS

   Build SUNDIALS with more extensive checks for unrecoverable errors.

   Default: ``ON`` when :cmakeop:`CMAKE_BUILD_TYPE` is ``Debug``, otherwise
   ``OFF``

   .. warning::

      Error checks will impact performance, but can be helpful for debugging.

.. _Installation.Options.Logging:

Logging
^^^^^^^

For more information on logging in SUNDIALS, see :ref:`SUNDIALS.Logging`.

.. cmakeoption:: SUNDIALS_LOGGING_LEVEL

   The maximum logging level. The options are:

   * ``0`` -- no logging
   * ``1`` -- log errors
   * ``2`` -- log errors + warnings
   * ``3`` -- log errors + warnings + informational output
   * ``4`` -- log errors + warnings + informational output + debug output
   * ``5`` -- log all of the above and even more (e.g. vector valued variables may be logged)

   Default: ``2``

   .. warning::

      Logging will impact performance, but can be helpful for debugging or
      understanding algorithm performance. The higher the logging level, the
      more output that may be logged, and the more performance may degrade.

   .. versionchanged:: 7.0.0

      Enabling MPI in SUNDIALS enables MPI-aware logging.

.. _Installation.Options.Monitoring:

Monitoring
^^^^^^^^^^

.. cmakeoption:: SUNDIALS_BUILD_WITH_MONITORING

   Build SUNDIALS with capabilities for fine-grained monitoring of solver
   progress and statistics. This is primarily useful for debugging.

   Default: ``OFF``

   .. warning::

      Building with monitoring may result in minor performance degradation even
      if monitoring is not utilized.

.. _Installation.Options.Profiling:

Profiling
^^^^^^^^^

For more information on profiling in SUNDIALS, see :ref:`SUNDIALS.Profiling`.

.. cmakeoption:: SUNDIALS_BUILD_WITH_PROFILING

   Build SUNDIALS with capabilities for fine-grained profiling. This requires
   POSIX timers, the Windows ``profileapi.h`` timers, or enabling Caliper with
   :cmakeop:`ENABLE_CALIPER`.

   Default: ``OFF``

   .. warning::

      Profiling will impact performance, and should be enabled judiciously.

.. _Installation.Options.Adiak:

Building with Adiak
^^^^^^^^^^^^^^^^^^^

`Adiak <http://software.llnl.gov/Adiak/>`__ is a library for recording meta-data
about HPC simulations. Adiak is developed by Lawrence Livermore National
Laboratory and can be obtained from the `Adiak GitHub repository
<https://github.com/LLNL/Adiak>`__.

.. cmakeoption:: ENABLE_ADIAK

   Enable Adiak support

   Default: ``OFF``

.. cmakeoption:: adiak_DIR

   Path to the root of an Adiak installation

   Default: None

.. _Installation.Options.Caliper:

Building with Caliper
^^^^^^^^^^^^^^^^^^^^^

`Caliper <https://software.llnl.gov/Caliper/>`__ is a performance analysis
library providing a code instrumentation and performance measurement framework
for HPC applications. Caliper is developed by Lawrence Livermore National
Laboratory and can be obtained from the `Caliper GitHub repository
<https://github.com/LLNL/Caliper>`__.

When profiling and Caliper are both enabled, SUNDIALS will utilize Caliper for
performance profiling.

To enable Caliper support, set the :cmakeop:`ENABLE_CALIPER` to ``ON`` and set
:cmakeop:`CALIPER_DIR` to the root path of the Caliper installation. For
example, the following command will configure SUNDIALS with profiling and
Caliper support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D SUNDIALS_BUILD_WITH_PROFILING=ON \
     -D ENABLE_CALIPER=ON \
     -D CALIPER_DIR=/path/to/caliper/installation

.. cmakeoption:: ENABLE_CALIPER

   Enable Caliper support

   Default: ``OFF``

 .. note::

      Using Caliper requires setting :cmakeop:`SUNDIALS_BUILD_WITH_PROFILING` to
      ``ON``.

.. cmakeoption:: CALIPER_DIR

   Path to the root of a Caliper installation

   Default: None

.. _Installation.Options.CUDA:

Building with CUDA
^^^^^^^^^^^^^^^^^^

The NVIDIA `CUDA Toolkit <https://developer.nvidia.com/cuda-toolkit>`__ provides
a development environment for GPU-accelerated computing with NVIDIA GPUs. The
CUDA Toolkit and compatible NVIDIA drivers are available from the `NVIDIA
developer website <https://developer.nvidia.com/cuda-downloads>`__. SUNDIALS has
been tested with the CUDA toolkit versions 10, 11, and 12.

When CUDA support is enabled, the :ref:`CUDA NVector <NVectors.CUDA>`, the
:ref:`cuSPARSE SUNMatrix <SUNMatrix.cuSparse>`, and the :ref:`cuSPARSE batched
QR SUNLinearSolver <SUNLinSol.cuSolverSp>` will be built (see sections
:numref:`Installation.LibrariesAndHeaders.Vector.CUDA`,
:numref:`Installation.LibrariesAndHeaders.Matrix.cuSPARSE`, and
:numref:`Installation.LibrariesAndHeaders.LinearSolver.cuSPARSE`, respectively,
for the corresponding header files and libraries). For more information on using
SUNDIALS with GPUs, see :ref:`SUNDIALS.GPU`.

To enable CUDA support, set :cmakeop:`ENABLE_CUDA` to ``ON``. If CUDA is
installed in a nonstandard location, you may need to set
:cmakeop:`CUDA_TOOLKIT_ROOT_DIR` to your CUDA Toolkit installation path. You
will also need to set :cmakeop:`CMAKE_CUDA_ARCHITECTURES` to the CUDA
architecture for your system. For example, the following command will configure
SUNDIALS with CUDA support for a system with an Ampere GPU:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_CUDA=ON \
     -D CMAKE_CUDA_ARCHITECTURES="80"

.. cmakeoption:: ENABLE_CUDA

   Enable CUDA support

   Default: ``OFF``

.. cmakeoption:: CUDA_TOOLKIT_ROOT_DIR

   Path to the CUDA Toolkit installation

   Default: CMake will attempt to automatically locate an installed CUDA Toolkit

.. cmakeoption:: CMAKE_CUDA_ARCHITECTURES

   Specifies the CUDA architecture to compile for i.e., ``60`` for Pascal,
   ``70`` for Volta, ``80`` for Ampere, ``90`` for Hopper, etc. See the `GPU
   compute capability tables <https://developer.nvidia.com/cuda-gpus>`__ on the
   NVIDIA webpage and the `GPU Compilation
   <https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/index.html#gpu-compilation>`__
   section of the CUDA documentation for more information.

   Default: Determined automatically by CMake. Users are encouraged to override
   this value with the architecture for their system as the default varies
   across compilers and compiler versions.

   .. versionchanged:: 7.2.0

      In prior versions :cmakeop:`CMAKE_CUDA_ARCHITECTURES` defaulted to ``70``.

.. _Installation.Options.Ginkgo:

Building with Ginkgo
^^^^^^^^^^^^^^^^^^^^

`Ginkgo <https://ginkgo-project.github.io/>`__ is a high-performance linear
algebra library with a focus on solving sparse linear systems. It is implemented
using modern C++ (you will need at least a C++17 compliant compiler to build
it), with GPU kernels implemented in CUDA (for NVIDIA devices), HIP (for AMD
devices), and SYCL/DPC++ (for Intel devices and other supported
hardware). Ginkgo can be obtained from the `Ginkgo GitHub repository
<https://github.com/ginkgo-project/ginkgo>`__. SUNDIALS requires using Ginkgo
version 1.9.0 or newer and is regularly tested with the latest versions of
Ginkgo, specifically versions 1.9.0 and 1.10.0.

When Ginkgo support is enabled, the :ref:`Ginkgo SUNMatrix <SUNMatrix.Ginkgo>`
and :ref:`SUNLinearSolver <SUNLinSol.Ginkgo>` as well as the :ref:`Ginkgo
Batch SUNMatrix <SUNMatrix.GinkgoBatch>` and :ref:`SUNLinearSolver
<SUNLinSol.GinkgoBatch>` header files will be installed (see sections
:numref:`Installation.LibrariesAndHeaders.Matrix.Ginkgo` and
:numref:`Installation.LibrariesAndHeaders.LinearSolver.Ginkgo`, respectively,
for the corresponding header files). For more information on using SUNDIALS with
GPUs, see :ref:`SUNDIALS.GPU`.

To enable Ginkgo support, set :cmakeop:`ENABLE_GINKGO` to ``ON`` and set
:cmakeop:`Ginkgo_DIR` to the root path of the Ginkgo installation. Additionally,
set :cmakeop:`SUNDIALS_GINKGO_BACKENDS` to a semicolon-separated list of Ginkgo
target architectures/executors. For example, the following command will
configure SUNDIALS with Ginkgo support using the reference, OpenMP, and CUDA
(targeting Ampere GPUs) backends:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_GINKGO=ON \
     -D Ginkgo_DIR=/path/to/ginkgo/installation \
     -D SUNDIALS_GINKGO_BACKENDS="REF;OMP;CUDA" \
     -D ENABLE_CUDA=ON \
     -D CMAKE_CUDA_ARCHITECTURES="80" \
     -D ENABLE_OPENMP=ON

.. note::

   The SUNDIALS interfaces to Ginkgo are not compatible with extended precision
   (i.e., when :cmakeop:`SUNDIALS_PRECISION` is set to ``extended``).

.. cmakeoption:: ENABLE_GINKGO

   Enable Ginkgo support

   Default: ``OFF``

.. cmakeoption:: Ginkgo_DIR

   Path to the Ginkgo installation

   Default: None

.. cmakeoption:: SUNDIALS_GINKGO_BACKENDS

   Semi-colon separated list of Ginkgo target architectures/executors to build
   for. Options currently supported are ``REF`` (the Ginkgo reference executor),
   ``OMP`` (OpenMP), ``CUDA``, ``HIP``, and ``SYCL``.

   Default: ``"REF;OMP"``

   .. versionchanged:: 7.1.0

      The ``DPCPP`` option was changed to ``SYCL`` to align with Ginkgo's naming
      convention.

.. _Installation.Options.HIP:

Building with HIP
^^^^^^^^^^^^^^^^^

The `Heterogeneous-compute Interface for Portability (HIP)
<https://rocm.docs.amd.com/projects/HIP/en/latest/>`__ allows developers to
create portable applications for AMD and NVIDIA GPUs. HIP can be obtained from
the `HIP GitHub repository
<https://github.com/ROCm-Developer-Tools/HIP>`__. SUNDIALS has been tested with
HIP versions between 5.0.0 to 5.4.3.

When HIP support is enabled, the :ref:`HIP NVector <NVectors.HIP>` will be built
(see section :numref:`Installation.LibrariesAndHeaders.Vector.HIP` for the
corresponding header file and library). For more information on using SUNDIALS
with GPUs, see :ref:`SUNDIALS.GPU`.

To enable HIP support, set :cmakeop:`ENABLE_HIP` to ``ON`` and set
:cmakeop:`AMDGPU_TARGETS` to the desired target (e.g., ``gfx705``). In addition,
set :cmakeop:`CMAKE_C_COMPILER` and :cmakeop:`CMAKE_CXX_COMPILER` to a HIP
compatible compiler e.g., ``hipcc``. For example, the following command will
configure SUNDIALS with HIP support for a system with an MI250X GPU:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D CMAKE_C_COMPILER=hipcc \
     -D CMAKE_CXX_COMPILER=hipcc \
     -D ENABLE_HIP=ON \
     -D AMDGPU_TARGETS="gfx90a"

.. cmakeoption:: ENABLE_HIP

   Enable HIP Support

   Default: ``OFF``

.. cmakeoption:: AMDGPU_TARGETS

   Specify which AMD GPUs to target

   Default: None

.. _Installation.Options.hypre:

Building with *hypre*
^^^^^^^^^^^^^^^^^^^^^

`hypre <https://www.llnl.gov/casc/hypre/>`__ is a library of high performance
preconditioners and solvers featuring multigrid methods for the solution of
large, sparse linear systems of equations on massively parallel computers. The
library is developed by Lawrence Livermore National Laboratory and is available
from the `hypre GitHub repository
<https://github.com/hypre-space/hypre>`__. SUNDIALS is regularly tested with the
latest versions of *hypre*, specifically up to version 2.26.0.

When *hypre* support is enabled, the :ref:`ParHyp NVector <NVectors.ParHyp>`
will be built (see section
:numref:`Installation.LibrariesAndHeaders.Vector.ParHyp` for the corresponding
header file and library).

To enable *hypre* support, set :cmakeop:`ENABLE_MPI` to ``ON``, set
:cmakeop:`ENABLE_HYPRE` to ``ON``, and set :cmakeop:`HYPRE_DIR` to the root path
of the *hypre* installation. For example, the following command will configure
SUNDIALS with *hypre* support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_MPI=ON \
     -D ENABLE_HYPRE=ON \
     -D HYPRE_DIR=/path/to/hypre/installation

.. note::

   SUNDIALS must be configured so that :cmakeop:`SUNDIALS_INDEX_SIZE` is
   compatible with ``HYPRE_BigInt`` in the *hypre* installation.

.. cmakeoption:: ENABLE_HYPRE

   Enable *hypre* support

   Default: ``OFF``

.. cmakeoption:: HYPRE_DIR

   Path to the *hypre* installation

   Default: None

.. _Installation.Options.KLU:

Building with KLU
^^^^^^^^^^^^^^^^^

KLU is a software package for the direct solution of sparse nonsymmetric linear
systems of equations that arise in circuit simulation and is part of
`SuiteSparse <https://people.engr.tamu.edu/davis/suitesparse.html>`__, a suite
of sparse matrix software. The library is developed by Texas A&M University and
is available from the `SuiteSparse GitHub repository
<https://github.com/DrTimothyAldenDavis/SuiteSparse>`__. SUNDIALS is regularly
tested with the latest versions of KLU, specifically up to SuiteSparse version
7.7.0.

When KLU support is enabled, the :ref:`KLU SUNLinearSolver <SUNLinSol.KLU>` will
be built (see section
:numref:`Installation.LibrariesAndHeaders.LinearSolver.KLU` for the
corresponding header file and library).

To enable KLU support, set :cmakeop:`ENABLE_KLU` to ``ON``. For SuiteSparse
7.4.0 and newer, set :cmakeop:`KLU_ROOT` to the root of the SuiteSparse
installation. Alternatively, set :cmakeop:`KLU_INCLUDE_DIR` and
:cmakeop:`KLU_LIBRARY_DIR` to the path to the header and library files,
respectively, of the SuiteSparse installation.  For example, the
following command will configure SUNDIALS with KLU support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_KLU=ON \
     -D KLU_ROOT=/path/to/suitesparse/installation

.. cmakeoption:: ENABLE_KLU

   Enable KLU support

   Default: ``OFF``

.. cmakeoption:: KLU_ROOT

   Path to the SuiteSparse installation

   Default: None

.. cmakeoption:: KLU_INCLUDE_DIR

   Path to SuiteSparse header files

   Default: None

.. cmakeoption:: KLU_LIBRARY_DIR

   Path to SuiteSparse installed library files

   Default: None

.. _Installation.Options.Kokkos:

Building with Kokkos
^^^^^^^^^^^^^^^^^^^^

`Kokkos <https://kokkos.github.io/kokkos-core-wiki/>`__ is a modern C++
(requires at least C++14) programming model for witting performance portable
code for multicore CPU and GPU-based systems including NVIDIA, AMD, and Intel
GPUs. Kokkos is developed by Sandia National Laboratory and can be obtained from
the `Kokkos GitHub repository <https://github.com/kokkos/kokkos>`__. The minimum
supported version of Kokkos 3.7.00. SUNDIALS is regularly tested with the latest
versions of Kokkos, specifically up to version 4.3.01.

When Kokkos support is enabled, the :ref:`Kokkos NVector <NVectors.Kokkos>`
header file will be installed (see section
:numref:`Installation.LibrariesAndHeaders.Vector.Kokkos` for the corresponding
header file). For more information on using SUNDIALS with GPUs, see
:ref:`SUNDIALS.GPU`.

To enable Kokkos support, set the :cmakeop:`ENABLE_KOKKOS` to ``ON`` and set
:cmakeop:`Kokkos_DIR` to root path of the Kokkos installation. For example, the
following command will configure SUNDIALS with Kokkos support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_KOKKOS=ON \
     -D Kokkos_DIR=/path/to/kokkos/installation

.. cmakeoption:: ENABLE_KOKKOS

   Enable Kokkos support

   Default: ``OFF``

.. cmakeoption:: Kokkos_DIR

   Path to the Kokkos installation.

   Default: None

.. _Installation.Options.KokkosKernels:

Building with KokkosKernels
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `KokkosKernels <https://github.com/kokkos/kokkos-kernels>`__ library is
built on Kokkos and provides common linear algebra computational kernels.
KokkosKernels is developed by Sandia National Laboratory and can be obtained
from the `KokkosKernels GitHub repository
<https://github.com/kokkos/kokkos-kernels>`__. The minimum supported version of
KokkosKernels 3.7.00. SUNDIALS is regularly tested with the latest versions of
KokkosKernels, specifically up to version 4.3.01.

When KokkosKernels support is enabled, the :ref:`KokkosKernels SUNMatrix
<SUNMatrix.Kokkos>` and :ref:`KokkosKernels SUNLinearSolver <SUNLinSol.Kokkos>`
header files will be installed (see sections
:numref:`Installation.LibrariesAndHeaders.Matrix.KokkosKernels` and
:numref:`Installation.LibrariesAndHeaders.LinearSolver.KokkosKernels`,
respectively, for the corresponding header files). For more information on using
SUNDIALS with GPUs, see :ref:`SUNDIALS.GPU`.

To enable KokkosKernels support, set :cmakeop:`ENABLE_KOKKOS` and
:cmakeop:`ENABLE_KOKKOS_KERNELS` to ``ON`` and set :cmakeop:`Kokkos_DIR` and
:cmakeop:`KokkosKernels_DIR` to the root paths for the Kokkos and KokkosKernels
installations, respectively. For example, the following command will configure
SUNDIALS with Kokkos and KokkosKernels support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_KOKKOS=ON \
     -D Kokkos_DIR=/path/to/kokkos/installation \
     -D ENABLE_KOKKOS_KERNELS=ON \
     -D KokkosKernels_DIR=/path/to/kokkoskernels/installation

.. cmakeoption:: ENABLE_KOKKOS_KERNELS

   Enable KokkosKernels support

   Default: ``OFF``

.. cmakeoption:: KokkosKernels_DIR

   Path to the KokkosKernels installation.

   Default: None

.. _Installation.Options.LAPACK:

Building with LAPACK
^^^^^^^^^^^^^^^^^^^^

The `Linear Algebra PACKage (LAPACK) <https://netlib.org/lapack/>`__ library
interface defines functions for solving systems of linear equations. Several
LAPACK implementations are available e.g., the `Netlib reference implementation
<https://www.netlib.org/lapack/>`__, the `Intel oneAPI Math Kernel Library
<https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html>`__,
or `OpenBLAS <http://www.openmathlib.org/OpenBLAS/>`__ (among others). SUNDIALS
is regularly tested with the latest versions of OpenBLAS, specifically up to
version 0.3.27.

When LAPACK support is enabled, the :ref:`LAPACK banded SUNLinearSolver
<SUNLinSol_LapackBand>` and :ref:`LAPACK dense SUNLinearSolver
<SUNLinSol_LapackDense>` will be built (see sections
:numref:`Installation.LibrariesAndHeaders.LinearSolver.LAPACKBand` and
:numref:`Installation.LibrariesAndHeaders.LinearSolver.LAPACKDense`,
respectively, for the corresponding header files and libraries). Additionally,
the :ref:`Arnoldi iteration SUNDomEigEstimator <SUNDomEigEst.Arnoldi>` will be
build (see :numref:`Installation.LibrariesAndHeaders.DomEigEst.Arnoldi`).

To enable LAPACK support, set :cmakeop:`ENABLE_LAPACK` to ``ON``. CMake will
attempt to find BLAS and LAPACK installations on the system and set the
variables :cmakeop:`BLAS_LIBRARIES`, :cmakeop:`BLAS_LINKER_FLAGS`,
:cmakeop:`LAPACK_LIBRARIES`, and :cmakeop:`LAPACK_LINKER_FLAGS`.  You can set
the :cmakeop:`LAPACK_ROOT` CMake variable to the path of a desired LAPACK
installation, and/or set the option :cmakeop:`BLA_VENDOR` to tell CMake to only
look for LAPACK from a specified vendor (see the `CMake documentation
<https://cmake.org/cmake/help/latest/module/FindLAPACK.html#input-variables>`__).
If necessary, to explicitly override the LAPACK library to build with, manually
set the aforementioned variables to the desired values when configuring the
build. For example, this is sometimes needed when using OpenBLAS:

.. code-block:: bash

  cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_LAPACK=ON \
     -D BLAS_LIBRARIES=/path/to/lapack/installation/lib/libopenblas.so \
     -D LAPACK_LIBRARIES=/path/to/lapack/installation/lib/libopenblas.so

.. note::

   If a working Fortran compiler is not available to infer the name-mangling
   scheme for LAPACK functions, the options :cmakeop:`SUNDIALS_LAPACK_CASE` and
   :cmakeop:`SUNDIALS_LAPACK_UNDERSCORES` *must* be set to bypass the check for
   a Fortran compiler and define the name-mangling scheme. The defaults for
   these options in earlier versions of SUNDIALS were ``lower`` and ``one``,
   respectively.

.. cmakeoption:: ENABLE_LAPACK

   Enable LAPACK support

   Default: ``OFF``

.. cmakeoption:: LAPACK_ROOT

   Path to the LAPACK installation

   Default: None

.. cmakeoption:: BLA_VENDOR

   The LAPACK vendor to search for.

   Default: All vendors

.. cmakeoption:: BLAS_LIBRARIES

   BLAS libraries

   Default: None (CMake will try to find a BLAS installation)

.. cmakeoption:: BLAS_LINKER_FLAGS

   BLAS required linker flags

   Default: None (CMake will try to determine the necessary flags)

.. cmakeoption:: LAPACK_LIBRARIES

   LAPACK libraries

   Default: None (CMake will try to find a LAPACK installation)

.. cmakeoption:: LAPACK_LINKER_FLAGS

   LAPACK required linker flags

   Default: None (CMake will try to determine the necessary flags)

.. cmakeoption:: SUNDIALS_LAPACK_CASE

   Specify the case to use in the Fortran name-mangling scheme,
   options are: ``lower`` or ``upper``

   Default:

   .. note::

      The build system will attempt to infer the Fortran name-mangling scheme
      using the Fortran compiler. This option should only be used if a Fortran
      compiler is not available or to override the inferred or default
      (``lower``) scheme if one can not be determined. If used,
      :cmakeop:`SUNDIALS_LAPACK_UNDERSCORES` must also be set.

.. cmakeoption:: SUNDIALS_LAPACK_UNDERSCORES

   Specify the number of underscores to append in the Fortran
   name-mangling scheme, options are: ``none``, ``one``, or ``two``

   Default:

   .. note::

      The build system will attempt to infer the Fortran name-mangling scheme
      using the Fortran compiler. This option should only be used if a Fortran
      compiler is not available or to override the inferred or default (``one``)
      scheme if one can not be determined. If used,
      :cmakeop:`SUNDIALS_LAPACK_CASE` must also be set.

.. _Installation.Options.MAGMA:

Building with MAGMA
^^^^^^^^^^^^^^^^^^^

The `Matrix Algebra on GPU and Multicore Architectures (MAGMA)
<https://icl.utk.edu/magma/>`__ project provides a dense linear algebra library
similar to LAPACK but targeting heterogeneous architectures. The library is
developed by the University of Tennessee and is available from the `MAGMA GitHub
repository <https://github.com/icl-utk-edu/magma>`__. SUNDIALS is regularly
tested with the latest versions of MAGMA, specifically up to version 2.8.0.

When MAGMA support is enabled, the :ref:`MAGMA dense SUNMatrix
<SUNMatrix.MagmaDense>` and :ref:`MAGMA dense SUNLinearSolver
<SUNLinSol.MagmaDense>` will be built (see sections
:numref:`Installation.LibrariesAndHeaders.Matrix.MAGMADense` and
:numref:`Installation.LibrariesAndHeaders.LinearSolver.MAGMADense`,
respectively, for the corresponding header files and libraries). For more
information on using SUNDIALS with GPUs, see :ref:`SUNDIALS.GPU`.

To enable MAGMA support, set :cmakeop:`ENABLE_MAGMA` to ``ON``,
:cmakeop:`MAGMA_DIR` to the root path of MAGMA installation, and
:cmakeop:`SUNDIALS_MAGMA_BACKENDS` to the desired MAGMA backend to use. For
example, the following command will configure SUNDIALS with MAGMA support with
the CUDA backend (targeting Ampere GPUs):

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_MAGMA=ON \
     -D MAGMA_DIR=/path/to/magma/installation \
     -D SUNDIALS_MAGMA_BACKEND="CUDA" \
     -D ENABLE_CUDA=ON \
     -D CMAKE_CUDA_ARCHITECTURES="80"

.. cmakeoption:: ENABLE_MAGMA

   Enable MAGMA support

   Default: ``OFF``

.. cmakeoption:: MAGMA_DIR

   Path to the MAGMA installation

   Default: None

.. cmakeoption:: SUNDIALS_MAGMA_BACKENDS

   Which MAGMA backend to use under the SUNDIALS MAGMA interface: ``CUDA`` or
   ``HIP``

   Default: ``CUDA``

   .. TODO(DJG): Change this options so it is HIP or SYCL if those options are
      enabled

.. _Installation.Options.MPI:

Building with MPI
^^^^^^^^^^^^^^^^^

The `Message Passing Interface (MPI) <https://www.mpi-forum.org/>`__ is a
standard for communication on parallel computing systems. Several MPI
implementations are available e.g., `OpenMPI <https://www.open-mpi.org/>`__,
`MPICH <https://www.mpich.org/>`__, `MVAPICH
<https://mvapich.cse.ohio-state.edu/>`__, `Cray MPICH
<https://cpe.ext.hpe.com/docs/24.03/mpt/mpich/index.html>`__, `Intel MPI
<https://www.intel.com/content/www/us/en/developer/tools/oneapi/mpi-library.html>`__,
or `IBM Spectrum MPI <https://www.ibm.com/products/spectrum-mpi>`__ (among
others). SUNDIALS is regularly tested with the latest versions of OpenMPI,
specifically up to version 5.0.5.

When MPI support is enabled, the :ref:`parallel NVector <NVectors.NVParallel>`,
:ref:`MPI ManyVector NVector <NVectors.MPIManyVector>`, and :ref:`MPI+X NVector
<NVectors.MPIPlusX>` will be built (see sections
:numref:`Installation.LibrariesAndHeaders.Vector.Parallel`,
:numref:`Installation.LibrariesAndHeaders.Vector.MPIManyVector`, and
:numref:`Installation.LibrariesAndHeaders.Vector.MPIPlusX`,
respectively, for the corresponding header files and libraries).

.. attention::

   .. versionchanged:: 7.0.0

      When MPI is enabled, all SUNDIALS libraries will include MPI symbols and
      applications will need to include the path for MPI headers and link against
      the corresponding MPI library.

To enable MPI support, set :cmakeop:`ENABLE_MPI` to ``ON``. If CMake is unable
to locate an MPI installation, set the relevant ``MPI_<language>_COMPILER``
options to the desired MPI compilers. For example, the following command will
configure SUNDIALS with MPI support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_MPI=ON

.. cmakeoption:: ENABLE_MPI

   Enable MPI support

   Default: ``OFF``

.. cmakeoption:: MPI_C_COMPILER

   The MPI C compiler e.g., ``mpicc``

   Default: CMake will attempt to locate an MPI C compiler

.. cmakeoption:: MPI_CXX_COMPILER

   The MPI C++ compiler e.g., ``mpicxx``

   Default: CMake will attempt to locate an MPI C++ compiler

   .. note::

      This option is only needed if MPI is enabled (:cmakeop:`ENABLE_MPI` is
      ``ON``) and C++ examples are enabled (:cmakeop:`EXAMPLES_ENABLE_CXX` is
      ``ON``). All SUNDIALS solvers can be used from C++ MPI applications by
      without setting any additional configuration options other than
      :cmakeop:`ENABLE_MPI`.

.. cmakeoption:: MPI_Fortran_COMPILER

   The MPI Fortran compiler e.g., ``mpif90``

   Default: CMake will attempt to locate an MPI Fortran compiler

   .. note::

      This option is triggered only needed if MPI is enabled
      (:cmakeop:`ENABLE_MPI` is ``ON``) and the Fortran interfaces are enabled
      (:cmakeop:`BUILD_FORTRAN_MODULE_INTERFACE` is ``ON``).

.. cmakeoption:: MPIEXEC_EXECUTABLE

   Specify the executable for running MPI programs e.g., ``mpiexec``

   Default: CMake will attempt to locate the MPI executable

.. cmakeoption:: MPIEXEC_PREFLAGS

   Specifies flags that come directly after ``MPIEXEC_EXECUTABLE`` and before
   ``MPIEXEC_NUMPROC_FLAG`` and ``MPIEXEC_MAX_NUMPROCS``.

   Default: None

.. cmakeoption:: MPIEXEC_POSTFLAGS

   Specifies flags that come after the executable to run but before any other
   program arguments.

   Default: None

.. _Installation.Options.OneMKL:

Building with oneMKL
^^^^^^^^^^^^^^^^^^^^

The Intel `oneAPI Math Kernel Library (oneMKL)
<https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html>`__
includes CPU and SYCL/DPC++ interfaces for LAPACK dense linear algebra
routines. The SUNDIALS oneMKL interface targets the SYCL/DPC++ routines, to
utilize the CPU routine see section
:numref:`Installation.Options.LAPACK`. SUNDIALS has been tested with oneMKL
version 2021.4.

When oneMKL support is enabled, the :ref:`oneMLK dense SUNMatrix
<SUNMatrix.OneMklDense>` and the :ref:`oneMKL dense SUNLinearSolver
<SUNLinSol.OneMklDense>` will be built (see sections
:numref:`Installation.LibrariesAndHeaders.Matrix.oneMKLDense` and
:numref:`Installation.LibrariesAndHeaders.LinearSolver.oneMKLDense`,
respectively, for the corresponding header files and libraries). For more
information on using SUNDIALS with GPUs, see :ref:`SUNDIALS.GPU`.

To enable the SUNDIALS oneMKL interface set :cmakeop:`ENABLE_ONEMKL` to ``ON``
and :cmakeop:`ONEMKL_DIR` to the root path of oneMKL installation. For example,
the following command will configure SUNDIALS with oneMKL support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_ONEMKL=ON \
     -D ONEMKL_DIR=/path/to/onemkl/installation \

.. cmakeoption:: ENABLE_ONEMKL

   Enable oneMKL support

   Default: ``OFF``

.. cmakeoption:: ONEMKL_DIR

   Path to oneMKL installation.

   Default: None

.. cmakeoption:: SUNDIALS_ONEMKL_USE_GETRF_LOOP

   This advanced debugging option replaces the batched LU factorization with a
   loop over each system in the batch and a non-batched LU factorization.

   Default: ``OFF``

.. cmakeoption:: SUNDIALS_ONEMKL_USE_GETRS_LOOP

   This advanced debugging option replaces the batched LU solve with a loop over
   each system in the batch and a non-batched solve.

   Default: ``OFF``

.. _Installation.Options.OpenMP:

Building with OpenMP
^^^^^^^^^^^^^^^^^^^^

The `OpenMP <https://www.openmp.org/>`__ API defines a directive-based approach
for portable parallel programming across architectures.

When OpenMP support is enabled, the :ref:`OpenMP NVector <NVectors.OpenMP>` will
be built (see section :numref:`Installation.LibrariesAndHeaders.Vector.OpenMP`
for the corresponding header file and library).

To enable OpenMP support, set the :cmakeop:`ENABLE_OPENMP` to ``ON``. For
example, the following command will configure SUNDIALS with OpenMP support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_OPENMP=ON

.. cmakeoption:: ENABLE_OPENMP

   Enable OpenMP support

   Default: ``OFF``

Building with OpenMP Device Offloading
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `OpenMP <https://www.openmp.org/>`__ 4.0 specification added support for
offloading computations to devices (i.e., GPUs). SUNDIALS requires OpenMP 4.5
for GPU offloading support.

When OpenMP offloading support is enabled, the :ref:`OpenMPDEV NVector
<NVectors.OpenMPDEV>` will be built (see section
:numref:`Installation.LibrariesAndHeaders.Vector.OpenMPDEV` for the
corresponding header file and library).

To enable OpenMP device offloading support, set the
:cmakeop:`ENABLE_OPENMP_DEVICE` to ``ON``. For example, the following command
will configure SUNDIALS with OpenMP device offloading support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_OPENMP_DEVICE=ON

.. cmakeoption:: ENABLE_OPENMP_DEVICE

   Enable OpenMP device offloading support

   Default: ``OFF``

.. _Installation.Options.PETSc:

Building with PETSc
^^^^^^^^^^^^^^^^^^^

The `Portable, Extensible Toolkit for Scientific Computation (PETSc)
<https://petsc.org>`__ is a suite of data structures and routines for simulating
applications modeled by partial differential equations. The library is developed
by Argonne National Laboratory and is available from the `PETSc GitLab
repository <https://gitlab.com/petsc/petsc>`__.  SUNDIALS requires PETSc 3.5.0
or newer and is regularly tested with the latest versions of PETSc, specifically
up to version 3.21.4.

When PETSc support is enabled, the :ref:`PETSc NVector <NVectors.NVPETSc>` and
:ref:`PETSc SNES SUNNonlinearSolver <SUNNonlinSol.PetscSNES>` will be built (see
sections :numref:`Installation.LibrariesAndHeaders.Vector.PETSc` and
:numref:`Installation.LibrariesAndHeaders.NonlinearSolver.PETScSNES`,
respectively, for the corresponding header files and libraries).

To enable PETSc support, set :cmakeop:`ENABLE_MPI` to ``ON``, set
:cmakeop:`ENABLE_PETSC` to ``ON``, and set :cmakeop:`PETSC_DIR` to the path of
the PETSc installation. Alternatively, a user can provide a list of include
paths in :cmakeop:`PETSC_INCLUDES` and a list of complete paths to the PETSc
libraries in :cmakeop:`PETSC_LIBRARIES`. For example, the following command will
configure SUNDIALS with PETSc support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_MPI=ON \
     -D ENABLE_PETSC=ON \
     -D PETSC_DIR=/path/to/petsc/installation

.. cmakeoption:: ENABLE_PETSC

   Enable PETSc support

   Default: ``OFF``

.. cmakeoption:: PETSC_DIR

   Path to PETSc installation

   Default: None

.. cmakeoption:: PETSC_LIBRARIES

   Semi-colon separated list of PETSc link libraries. Unless provided by the
   user, this is autopopulated based on the PETSc installation found in
   :cmakeop:`PETSC_DIR`.

   Default: None

.. cmakeoption:: PETSC_INCLUDES

   Semi-colon separated list of PETSc include directories. Unless provided by
   the user, this is autopopulated based on the PETSc installation found in
   :cmakeop:`PETSC_DIR`.

   Default: None

.. _Installation.Options.PThreads:

Building with PThreads
^^^^^^^^^^^^^^^^^^^^^^

POSIX Threads (PThreads) is an API for shared memory programming defined by the
Institute of Electrical and Electronics Engineers (IEEE) standard POSIX.1c.

When PThreads support is enabled, the :ref:`PThreads NVector
<NVectors.PThreads>` will be built (see section
:numref:`Installation.LibrariesAndHeaders.Vector.PThreads` for the corresponding
header file and library).

To enable PThreads support, set :cmakeop:`ENABLE_PTHREAD` to ``ON``. For
example, the following command will configure SUNDIALS with PThreads support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_PTHREAD=ON

.. cmakeoption:: ENABLE_PTHREAD

   Enable PThreads support

   Default: ``OFF``

.. _Installation.Options.RAJA:

Building with RAJA
^^^^^^^^^^^^^^^^^^

`RAJA <https://raja.readthedocs.io/en/develop/>`__ is a performance portability
layer developed by Lawrence Livermore National Laboratory and can be obtained
from the `RAJA GitHub repository <https://github.com/LLNL/RAJA>`__. SUNDIALS is
regularly tested with the latest versions of RAJA, specifically up to version
2024.02.2.

When RAJA support is enabled, the :ref:`RAJA NVector <NVectors.RAJA>` will be
built (see section :numref:`Installation.LibrariesAndHeaders.Vector.RAJA`
for the corresponding header files and libraries).

To enable RAJA support, set :cmakeop:`ENABLE_RAJA` to ``ON``, set
:cmakeop:`RAJA_DIR` to the path of the RAJA installation, set
:cmakeop:`SUNDIALS_RAJA_BACKENDS` to the desired backend (``CUDA``, ``HIP``, or
``SYCL``), and set :cmakeop:`ENABLE_CUDA`, :cmakeop:`ENABLE_HIP`, or
:cmakeop:`ENABLE_SYCL` to ``ON`` depending on the selected backend. For
example, the following command will configure SUNDIALS with RAJA support using
the CUDA backend (targeting Ampere GPUs):

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_RAJA=ON \
     -D RAJA_DIR=/path/to/raja/installation \
     -D SUNDIALS_RAJA_BACKENDS="CUDA" \
     -D ENABLE_CUDA=ON \
     -D CMAKE_CUDA_ARCHITECTURES="80"

.. cmakeoption:: ENABLE_RAJA

   Enable RAJA support

   Default: ``OFF``

.. cmakeoption:: RAJA_DIR

   Path to the RAJA installation

   Default: None

.. cmakeoption:: SUNDIALS_RAJA_BACKENDS

   If building SUNDIALS with RAJA support, this sets the RAJA backend to target.
   Values supported are ``CUDA``, ``HIP``, or ``SYCL``.

   Default: ``CUDA``

.. _Installation.Options.SuperLU_DIST:

Building with SuperLU_DIST
^^^^^^^^^^^^^^^^^^^^^^^^^^

`SuperLU_DIST <https://portal.nersc.gov/project/sparse/superlu/>`__ is a general
purpose library for the direct solution of large, sparse, nonsymmetric systems
of linear equations in a distributed memory setting. The library is developed by
Lawrence Berkeley National Laboratory and is available from the `SuperLU_DIST
GitHub repository <https://github.com/xiaoyeli/superlu_dist>`__. SuperLU_DIST
version 7.0.0 or newer is required. SUNDIALS is regularly tested with the latest
versions of SuperLU_DIST, specifically up to version 8.2.1.

When SuperLU_DIST support is enabled, the :ref:`SuperLU_DIST (SLUNRloc)
SUNMatrix <SUNMatrix.SLUNRloc>` and :ref:`SuperLU_DIST SUNLinearSolver
<SUNLinSol.SuperLUDIST>` will be built (see sections
:numref:`Installation.LibrariesAndHeaders.Matrix.SuperLU_DIST` and
:numref:`Installation.LibrariesAndHeaders.LinearSolver.SuperLU_DIST` for the
corresponding header files and libraries).

To enable SuperLU_DIST support, set :cmakeop:`ENABLE_MPI` to ``ON``, set
:cmakeop:`ENABLE_SUPERLUDIST` to ``ON``, and set :cmakeop:`SUPERLUDIST_DIR` to
the path where SuperLU_DIST is installed. If SuperLU_DIST was built with OpenMP
enabled, set :cmakeop:`SUPERLUDIST_OpenMP` and :cmakeop:`ENABLE_OPENMP` to
``ON``. For example, the following command will configure SUNDIALS with
SuperLU_DIST support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_SUPERLUDIST=ON \
     -D SUPERLUDIST_DIR=/path/to/superludist/installation

.. cmakeoption:: ENABLE_SUPERLUDIST

   Enable SuperLU_DIST support

   Default: ``OFF``

.. cmakeoption:: SUPERLUDIST_DIR

   Path to SuperLU_DIST installation.

   Default: None

.. cmakeoption:: SUPERLUDIST_OpenMP

   Enable SUNDIALS support for SuperLU_DIST built with OpenMP

   Default: None

   .. note::

      SuperLU_DIST must be built with OpenMP support for this option to
      function. Additionally the environment variable ``OMP_NUM_THREADS`` must
      be set to the desired number of threads.

.. cmakeoption:: SUPERLUDIST_INCLUDE_DIRS

   List of include paths for SuperLU_DIST (under a typical SuperLU_DIST
   install, this is typically the SuperLU_DIST ``SRC`` directory)

   Default: None

   .. note::

      This is an advanced option. Prefer to use :cmakeop:`SUPERLUDIST_DIR`.

.. cmakeoption:: SUPERLUDIST_LIBRARIES

   Semi-colon separated list of libraries needed for SuperLU_DIST

   Default: None

   .. note::

      This is an advanced option. Prefer to use :cmakeop:`SUPERLUDIST_DIR`.

.. cmakeoption:: SUPERLUDIST_INCLUDE_DIR

   Path to SuperLU_DIST header files (under a typical SuperLU_DIST
   install, this is typically the SuperLU_DIST ``SRC`` directory)

   Default: None

   .. note::

      This is an advanced option. This option is deprecated. Use
      :cmakeop:`SUPERLUDIST_INCLUDE_DIRS`.

.. cmakeoption:: SUPERLUDIST_LIBRARY_DIR

   Path to SuperLU_DIST installed library files

   Default: None

   .. note::

      This option is deprecated. Use :cmakeop:`SUPERLUDIST_DIR`.

.. _Installation.Options.SuperLU_MT:

Building with SuperLU_MT
^^^^^^^^^^^^^^^^^^^^^^^^

`SuperLU_MT <https://portal.nersc.gov/project/sparse/superlu/>`__ is a general
purpose library for the direct solution of large, sparse, nonsymmetric systems
of linear equations on shared memory parallel machines. The library is developed
by Lawrence Berkeley National Laboratory and is available from the `SuperLU_MT
GitHub repository <https://github.com/xiaoyeli/superlu_mt>`__. SUNDIALS is
regularly tested with the latest versions of SuperLU_MT, specifically up to
version 4.0.1.

When SuperLU_MT support is enabled, the :ref:`SuperLU_MT SUNLinearSolver
<SUNLinSol.SuperLUMT>` will be built (see section
:numref:`Installation.LibrariesAndHeaders.LinearSolver.SuperLU_MT` for the
corresponding header file and library).

To enable SuperLU_MT support, set :cmakeop:`ENABLE_SUPERLUMT` to ``ON``, set
:cmakeop:`SUPERLUMT_INCLUDE_DIR` and :cmakeop:`SUPERLUMT_LIBRARY_DIR` to the
location of the header and library files, respectively, of the SuperLU_MT
installation. Depending on the SuperLU_MT installation, it may also be necessary
to set :cmakeop:`SUPERLUMT_LIBRARIES` to a semi-colon separated list of other
libraries SuperLU_MT depends on. For example, if SuperLU_MT was build with an
external blas library, then include the full path to the blas library in this
list. Additionally, the variable :cmakeop:`SUPERLUMT_THREAD_TYPE` must be set to
either ``Pthread`` or ``OpenMP``.  For example, the following command will
configure SUNDIALS with SuperLU_MT support using PThreads:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_SUPERLUMT=ON \
     -D SUPERLUMT_INCLUDE_DIR=/path/to/superlumt/installation/include/dir \
     -D SUPERLUMT_LIBRARY_DIR=/path/to/superlumt/installation/library/dir \
     -D SUPERLUMT_THREAD_TYPE="Pthread"

.. warning::

   Do not mix thread types when using SUNDIALS packages. For example, if using the
   OpenMP or PThreads NVector then the SuperLU_MT installation should use the same
   threading type.

.. cmakeoption:: ENABLE_SUPERLUMT

   Enable SuperLU_MT support

   Default: ``OFF``

.. cmakeoption:: SUPERLUMT_INCLUDE_DIR

   Path to SuperLU_MT header files (under a typical SuperLU_MT
   install, this is typically the SuperLU_MT ``SRC`` directory)

   Default: None

.. cmakeoption:: SUPERLUMT_LIBRARY_DIR

   Path to SuperLU_MT installed library files

   Default: None

.. cmakeoption:: SUPERLUMT_LIBRARIES

   Semi-colon separated list of libraries needed for SuperLU_MT

   Default: None

.. cmakeoption:: SUPERLUMT_THREAD_TYPE

   Must be set to Pthread or OpenMP, depending on how SuperLU_MT was compiled.

   Default: Pthread

.. _Installation.Options.SYCL:

Building with SYCL
^^^^^^^^^^^^^^^^^^

`SYCL <https://www.khronos.org/sycl/>`__ is an abstraction layer for programming
heterogeneous parallel computing based on C++17.

When SYCL support is enabled, the :ref:`SYCL NVector <NVectors.SYCL>` will
be built (see section :numref:`Installation.LibrariesAndHeaders.Vector.SYCL`
for the corresponding header file and library).

To enable SYCL support, set the :cmakeop:`ENABLE_SYCL` to ``ON``. For example,
the following command will configure SUNDIALS with SYCL support using Intel
compilers:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D CMAKE_C_COMPILER=icx \
     -D CMAKE_CXX_COMPILER=icpx \
     -D CMAKE_CXX_FLAGS="-fsycl" \
     -D ENABLE_SYCL=ON

.. cmakeoption:: ENABLE_SYCL

   Enable SYCL support

   Default: ``OFF``

   .. note::

      Building with SYCL enabled requires a compiler that supports a subset of
      the of SYCL 2020 specification (specifically ``sycl/sycl.hpp`` must be
      available).

      CMake does not currently support autodetection of SYCL compilers and
      :cmakeop:`CMAKE_CXX_COMPILER` must be set to a valid SYCL compiler. At
      present the only supported SYCL compilers are the Intel oneAPI compilers
      i.e., ``dpcpp`` and ``icpx``. When using ``icpx`` the ``-fsycl`` flag and
      any ahead of time compilation flags must be added to
      :cmakeop:`CMAKE_CXX_FLAGS`.

.. cmakeoption:: SUNDIALS_SYCL_2020_UNSUPPORTED

   This advanced option disables the use of *some* features from the SYCL 2020
   standard in SUNDIALS libraries and examples. This can be used to work around
   some cases of incomplete compiler support for SYCL 2020.

   Default: ``OFF``

.. _Installation.Options.Trilinos:

Building with Trilinos
^^^^^^^^^^^^^^^^^^^^^^

`Trilinos <https://trilinos.github.io/>`__ is a collection of C++ libraries of
linear solvers, non-linear solvers, optimization solvers, etc. developed by
Sandia National Laboratory and available from the `Trilinos GitHub repository
<https://github.com/trilinos/Trilinos>`__. SUNDIALS is regularly tested with
the latest versions of Trilinos, specifically up to version 16.0.0.

When Trilinos support is enabled, the :ref:`Trilinos Tpetra NVector
<NVectors.NVTrilinos>` will be built (see section
:numref:`Installation.LibrariesAndHeaders.Vector.Trilinos` for the corresponding
header file and library).

To enable Trilinos support, set the :cmakeop:`ENABLE_TRILINOS` to ``ON`` and set
:cmakeop:`Trilinos_DIR` to root path of the Trilinos installation. For example,
the following command will configure SUNDIALS with Trilinos support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D ENABLE_TRILONOS=ON \
     -D TRILINOS_DIR=/path/to/trilinos/installation

.. cmakeoption:: ENABLE_TRILINOS

   Enable Trilinos support

   Default: ``OFF``

.. cmakeoption:: Trilinos_DIR

   Path to the Trilinos installation

   Default: None

.. _Installation.Options.XBraid:

Building with XBraid
^^^^^^^^^^^^^^^^^^^^

XBraid is parallel-in-time library implementing an optimal-scaling multigrid
reduction in time (MGRIT) solver. The library is developed by Lawrence Livermore
National Laboratory and is available from the `XBraid GitHub repository
<https://github.com/XBraid/xbraid>`__. SUNDIALS is regularly tested with the
latest versions of XBraid, specifically up to version 3.0.0.

To enable XBraid support, set :cmakeop:`ENABLE_MPI` to ``ON``, set
:cmakeop:`ENABLE_XBRAID` to ``ON``, set :cmakeop:`XBRAID_DIR` to the root path
of the XBraid installation. For example, the following command will configure
SUNDIALS with XBraid support:

.. code-block:: bash

   cmake \
     -S SOLVER_DIR \
     -B BUILD_DIR \
     -D CMAKE_INSTALL_PREFIX=INSTALL_DIR \
     -D SUNDIALS_INDEX_SIZE="32" \
     -D ENABLE_MPI=ON \
     -D ENABLE_XBRAID=ON \
     -D XBRAID_DIR=/path/to/xbraid/installation

.. note::

   At this time the XBraid types ``braid_Int`` and ``braid_Real`` are hard-coded
   to ``int`` and ``double`` respectively. As such SUNDIALS must be configured
   with :cmakeop:`SUNDIALS_INDEX_SIZE` set to ``32`` and
   :cmakeop:`SUNDIALS_PRECISION` set to ``double``. Additionally, SUNDIALS must
   be configured with :cmakeop:`ENABLE_MPI` set to ``ON``.

.. cmakeoption:: ENABLE_XBRAID

   Enable or disable the ARKStep + XBraid interface.

   Default: ``OFF``

.. cmakeoption:: XBRAID_DIR

   The root directory of the XBraid installation.

   Default: ``OFF``

.. cmakeoption:: XBRAID_INCLUDES

   Semi-colon separated list of XBraid include directories. Unless provided by
   the user, this is autopopulated based on the XBraid installation found in
   :cmakeop:`XBRAID_DIR`.

   Default: None

.. cmakeoption:: XBRAID_LIBRARIES

   Semi-colon separated list of XBraid link libraries. Unless provided by
   the user, this is autopopulated based on the XBraid installation found in
   :cmakeop:`XBRAID_DIR`.

   Default: None


.. _Installation.Options.xSDK:

Building with xSDK Defaults
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `Extreme-scale Scientific Software Development Kit (xSDK)
<https://xsdk.info>`__ is a community of HPC libraries and applications
developing best practices and standards for scientific software.

.. cmakeoption:: USE_XSDK_DEFAULTS

   Enable xSDK default configuration settings. This sets the default value for
   :cmakeop:`CMAKE_BUILD_TYPE` to ``Debug``, :cmakeop:`SUNDIALS_INDEX_SIZE` to
   ``32``, and :cmakeop:`SUNDIALS_PRECISION` to ``double``.

   Default: ``OFF``

.. _Installation.Options.Addons:

Building with External Addons
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SUNDIALS "addons" are community developed code additions for SUNDIALS that can
be subsumed by the SUNDIALS build system so that they have full access to all
internal SUNDIALS symbols. The intent is for SUNDIALS addons to function as if
they are part of the SUNDIALS library, while allowing them to potentially have
different licenses (although we encourage BSD-3-Clause still), code style
(although we encourage them to follow the SUNDIALS style outlined :ref:`here
<SourceCode>`).

.. warning::

   SUNDIALS addons are not maintained by the SUNDIALS team and may come with
   different licenses. Use them at your own risk.

To build with SUNDIALS addons,

1. Clone/copy the addon(s) into ``SOLVER_DIR/external/``

2. Copy the ``sundials-addon-example`` block in the
   ``SOLVER_DIR/external/CMakeLists.txt``, paste it below the example block, and
   modify the path listed for your own external addon(s).

3. When building SUNDIALS, set the CMake option
   :cmakeop:`SUNDIALS_ENABLE_EXTERNAL_ADDONS` to ``ON``

4. Build SUNDIALS as usual.

.. cmakeoption:: SUNDIALS_ENABLE_EXTERNAL_ADDONS

   Build SUNDIALS with any external addons that you have put in
   ``SOLVER_DIR/external``.

   Default: ``OFF``

.. _Installation.Testing:

Testing the Build and Installation
----------------------------------

If SUNDIALS was configured with any ``EXAMPLES_ENABLE_<language>`` options set
to ``ON``, then a set of regression tests can be run after building with the
command:

.. code-block:: bash

   make test

Additionally, if :cmakeop:`EXAMPLES_INSTALL` is set to ``ON``, then a set of
smoke tests can be run after installing with the command:

.. code-block:: bash

   make test_install

.. _Installation.BuildRunExamples:

Building and Running Examples
-----------------------------

Each of the SUNDIALS solvers is distributed with a set of examples demonstrating
basic usage. To build and install the examples, set at least one of the
``EXAMPLES_ENABLE_<language>`` options to ``ON``, and set
:cmakeop:`EXAMPLES_INSTALL` to ``ON``. Along side the example sources and
outputs, automatically generated ``CMakeLists.txt`` configuration files (and
``Makefile`` files if on Linux/Unix systems) are installed referencing the
*installed* SUNDIALS headers and libraries.

Either the ``CMakeLists.txt`` file or the traditional ``Makefile`` may be used
to build the examples and serve as a template for building user developed
problems. To use the supplied ``Makefile`` simply run ``make`` to compile and
generate the executables. To use CMake from within the installed example
directory, run ``cmake`` (or ``ccmake`` or ``cmake-gui`` to use the GUI)
followed by ``make`` to compile the example code.  Note that if CMake is used,
it will overwrite the traditional ``Makefile`` with a new CMake-generated
``Makefile``.

The resulting output from running the examples can be compared with example
output bundled in the SUNDIALS distribution.

.. note::

   There will potentially be differences in the output due to machine
   architecture, compiler versions, use of third party libraries, etc.

.. _Installation.UsingSUNDIALS:

Using SUNDIALS In Your Project
------------------------------

After installing SUNDIALS, building your application with SUNDIALS involves two
steps: including the right header files and linking to the right libraries.
Depending on what features of SUNDIALS that your application uses, the header
files and libraries needed will vary. For example, if you want to use CVODE for
serial computations you need the following includes:

.. code-block:: c

   #include <cvode/cvode.h>
   #include <nvector/nvector_serial.h>

and must link to ``libsundials_cvode`` and ``libsundials_nvecserial``. If you
wanted to use CVODE with the GMRES linear solver and the CUDA NVector, you need
the following includes:

.. code-block:: c

   #include <cvode/cvode.h>
   #include <nvector/nvector_cuda.h>
   #include <sunlinsol/sunlinsol_spgmr.h>

and must link to ``libsundials_cvode``, ``libsundials_nveccuda``, and
``libsundials_sunlinsolspgmr``.

.. attention::

   .. versionadded:: 7.0.0

      All applications must also link to ``libsundials_core``. For projects
      using SUNDIALS CMake targets (see section
      :numref:`Installation.CMakeConfigFile`), this dependency is automatically
      included.

Refer to section :numref:`Installation.LibrariesAndHeaders` below or the
documentations sections for the individual SUNDIALS packages and modules of
interest for the proper includes and libraries to link against.

.. _Installation.CMakeConfigFile:

CMake Projects
^^^^^^^^^^^^^^

For projects that use CMake, the SUNDIALS `CMake package configuration file
<https://cmake.org/cmake/help/v3.18/manual/cmake-packages.7.html>`__ provides
CMake targets for the consuming project. Use the CMake ``find_package`` command
to search for the configuration file, ``SUNDIALSConfig.cmake``, which is
installed alongside a package version file, ``SUNDIALSConfigVersion.cmake``,
under the ``INSTALL_DIR/SUNDIALS_INSTALL_CMAKEDIR`` directory. The SUNDIALS
CMake targets follow the same naming convention as the generated library
binaries with the ``libsundials_`` prefix replaced by ``SUNDIALS::``. For
example, the exported target for ``libsundials_cvode`` is
``SUNDIALS::cvode``. See section :numref:`Installation.LibrariesAndHeaders` for
a complete list of CMake targets. The CMake code snippit below shows how a
consuming project might leverage the SUNDIALS package configuration file to
build against SUNDIALS in their own CMake project.

.. code-block:: cmake

  project(MyProject)

  # Set the variable SUNDIALS_DIR to the SUNDIALS instdir.
  # When using the cmake CLI command, this can be done like so:
  #   cmake -D SUNDIALS_DIR=/path/to/sundials/installation

  # Find any SUNDIALS version...
  find_package(SUNDIALS REQUIRED)

  # ... or find any version newer than some minimum...
  find_package(SUNDIALS 7.1.0 REQUIRED)

  # ... or find a version in a range
  find_package(SUNDIALS 7.0.0...7.1.0 REQUIRED)

  # To check if specific components are available in the SUNDIALS installation,
  # use the COMPONENTS option followed by the desired target names
  find_package(SUNDIALS REQUIRED COMPONENTS cvode nvecpetsc)

  add_executable(myexec main.c)

  # Link to SUNDIALS libraries through the exported targets.
  # This is just an example, users should link to the targets appropriate
  # for their use case.
  target_link_libraries(myexec PUBLIC SUNDIALS::cvode SUNDIALS::nvecpetsc)

.. note::

   .. versionchanged:: 7.1.0

      A single version provided to ``find_package`` denotes the minimum version
      of SUNDIALS to look for, and any version equal or newer than what is
      specified will match. In prior versions ``SUNDIALSConfig.cmake`` required
      the version found to have the same major version number as the single
      version provided to ``find_package``.

To accommodate installing both static and shared libraries simultaneously,
targets are created with ``_static`` and ``_shared`` suffixes, respectively, and
the un-suffixed target is an alias to the ``_shared`` version. For example,
``SUNDIALS::cvode`` is an alias to ``SUNDIALS::cvode_shared`` in this
case. Projects that wish to use static libraries should use the ``_static``
version of the target when both library types are installed. When only static or
shared libraries (not both) are installed the un-suffixed alias corresponds to
the library type chosen at configuration time (see section
:numref:`Installation.Options.LibraryTypes`).

.. _Installation.LibrariesAndHeaders:

Libraries and Header Files
--------------------------

As noted above, the SUNDIALS the header files and libraries are installed under
the :cmakeop:`CMAKE_INSTALL_PREFIX` path in the ``include`` and
:cmakeop:`CMAKE_INSTALL_LIBDIR` subdirectories, respectively. The public header
files are further organized into subdirectories under the ``include`` directory.
The installed public header files and libraries are listed for reference in the
sections below. Additionally, the exported CMake targets are also listed for
projects using CMake (see section :numref:`Installation.CMakeConfigFile`). The
file extension ``.LIB`` used below is typically ``.so``, ``.dll``, or ``.dylib``
for shared libraries and ``.a`` or ``.lib`` for static libraries.

.. warning::

   SUNDIALS installs some header files to
   ``CMAKE_INSTALL_PREFIX/include/sundials/priv``. All of the header files in
   this directory are private and **should not be included in user code**. The
   private headers are subject to change without any notice and relying on them
   may break your code.

.. _Installation.LibrariesAndHeaders.Core:

SUNDIALS Core
^^^^^^^^^^^^^

The core library contains the shared infrastructure utilized by SUNDIALS
packages. All applications using SUNDIALS must link against the core
library. For codes using the SUNDIALS CMake targets, the core target is
automatically included as needed by other targets.

.. table:: The SUNDIALS core library, header, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_core.LIB``                     |
   +--------------+----------------------------------------------+
   | Headers      | ``sundials/sundials_core.h``                 |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::core``                           |
   +--------------+----------------------------------------------+

The core header file is a convenient way to include all the header files that
make up the SUNDIALS core infrastructure.

.. table:: Header files included by ``sundials_core.h``
   :align: center

   +--------------+-------------------------------------------------+
   | Headers      | ``sundials/sundials_adaptcontroller.h``         |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_adjointstepper.h``          |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_adjointcheckpointscheme.h`` |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_config.h``                  |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_context.h``                 |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_domeigestimator.h``         |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_errors.h``                  |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_iterative.h``               |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_linearsolver.h``            |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_logger.h``                  |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_math.h``                    |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_matrix.h``                  |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_memory.h``                  |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_nonlinearsolver.h``         |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_nvector.h``                 |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_profiler.h``                |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_types.h``                   |
   |              +-------------------------------------------------+
   |              | ``sundials/sundials_version.h``                 |
   +--------------+-------------------------------------------------+

For C++ applications, several convenience classes are provided for interacting
with SUNDIALS objects. These can be accessed by including the C++ core header
file.

.. table:: The SUNDIALS C++ core header file
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sundials/sundials_core.hpp``               |
   +--------------+----------------------------------------------+

Like the C core header file, the C++ core header file is a convenient way to
include all the header files for the core C++ classes.

.. warning::

   Features in the ``sundials::experimental`` namespace are not yet part of the
   public API and are subject to change or removal without notice.

.. table:: Header files included by ``sundials_core.hpp``
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sundials/sundials_context.hpp``            |
   |              +----------------------------------------------+
   |              | ``sundials/sundials_core.h``                 |
   |              +----------------------------------------------+
   |              | ``sundials/sundials_linearsolver.hpp``       |
   |              +----------------------------------------------+
   |              | ``sundials/sundials_matrix.hpp``             |
   |              +----------------------------------------------+
   |              | ``sundials/sundials_memory.hpp``             |
   |              +----------------------------------------------+
   |              | ``sundials/sundials_nonlinearsolver.hpp``    |
   |              +----------------------------------------------+
   |              | ``sundials/sundials_nvector.hpp``            |
   |              +----------------------------------------------+
   |              | ``sundials/sundials_profiler.hpp``           |
   +--------------+----------------------------------------------+

When MPI support is enabled (:cmakeop:`ENABLE_MPI` is ``ON``), the following
header file provides aliases between MPI data types and SUNDIALS types. The
alias ``MPI_SUNREALTYPE`` is one of ``MPI_FLOAT``, ``MPI_DOUBLE``, or
``MPI_LONG_DOUBLE`` depending on the value of :cmakeop:`SUNDIALS_PRECISION`. The
alias ``MPI_SUNINDEXTYPE`` is either ``MPI_INT32_T`` or ``MPI_INT64_T``
depending on the value of :cmakeop:`SUNDIALS_INDEX_SIZE`.

.. table:: Header file defining aliases between SUNDIALS and MPI data types
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sundials/sundials_mpi_types.h``            |
   +--------------+----------------------------------------------+

When XBraid support is enabled (:cmakeop:`ENABLE_XBRAID` is ``ON``), the
following header file defines types and functions for interfacing SUNDIALS with
XBraid.

.. table:: SUNDIALS header for interfacing with XBraid
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sundials/sundials_xbraid.h``               |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Packages:

SUNDIALS Packages
^^^^^^^^^^^^^^^^^

.. _Installation.LibrariesAndHeaders.Packages.CVODE:

CVODE
"""""

To use the :ref:`CVODE <CVODE>` package, include the header file and link to the
library given below.

.. table:: CVODE library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_cvode.LIB``                    |
   +--------------+----------------------------------------------+
   | Headers      | ``cvode/cvode.h``                            |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::cvode``                          |
   +--------------+----------------------------------------------+

The CVODE header file includes the files below which define functions, types,
and constants for the CVODE linear solver interface and using projection methods
with CVODE.

.. table:: Additional header files included by ``cvode.h``
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``cvode/cvode_ls.h``                         |
   |              +----------------------------------------------+
   |              | ``cvode/cvode_proj.h``                       |
   +--------------+----------------------------------------------+

CVODE provides a specialized linear solver module for diagonal linear
systems. Include the header file below to access the related functions.

.. table:: CVODE diagonal linear solver
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``cvode/cvode_diag.h``                       |
   +--------------+----------------------------------------------+

For problems in which the user cannot define a more effective, problem-specific
preconditioner for Krylov iterative linear solvers, CVODE provides banded
(``bandpre``) and band-block-diagonal (``bbdpre``) preconditioner
modules. Include the header files below to access the related functions.

.. table:: CVODE preconditioner modules
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``cvode/cvode_bandpre.h``                    |
   |              +----------------------------------------------+
   |              | ``cvode/cvode_bbdpre.h``                     |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Packages.CVODES:

CVODES
""""""

To use the :ref:`CVODES <CVODES>` package, include the header file and link to
the library given below.

.. warning::

   CVODES is a superset of CVODE and defines the same functions as provided by
   CVODE. As such, applications should not link to both CVODES and CVODE.

.. table:: CVODES library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_cvodes.LIB``                   |
   +--------------+----------------------------------------------+
   | Headers      | ``cvodes/cvodes.h``                          |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::cvodes``                         |
   +--------------+----------------------------------------------+

The CVODES header file includes the files below which define functions, types,
and constants for the CVODES linear solver interface and using projection
methods with CVODES.

.. table:: Additional header files included by ``cvodes.h``
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``cvodes/cvodes_ls.h``                       |
   |              +----------------------------------------------+
   |              | ``cvodes/cvodes_proj.h``                     |
   +--------------+----------------------------------------------+

CVODES provides a specialized linear solver module for diagonal linear
systems. Include the header file below to access the related functions.

.. table:: CVODES diagonal linear solver
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``cvodes/cvodes_diag.h``                     |
   +--------------+----------------------------------------------+

For problems in which the user cannot define a more effective, problem-specific
preconditioner for Krylov iterative linear solvers, CVODES provides banded
(``bandpre``) and band-block-diagonal (``bbdpre``) preconditioner
modules. Include the header files below to access the related functions.

.. table:: CVODES preconditioner modules
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``cvodes/cvodes_bandpre.h``                  |
   |              +----------------------------------------------+
   |              | ``cvodes/cvodes_bbdpre.h``                   |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Packages.ARKODE:

ARKODE
""""""

To use the :ref:`ARKODE <ARKODE>` package, link to the library below and include
the header file for the desired module.

.. table:: ARKODE library, header files, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_arkode.LIB``                   |
   +--------------+----------------------------------------------+
   | Headers      | ``arkode/arkode_arkstep.h``                  |
   |              +----------------------------------------------+
   |              | ``arkode/arkode_erkstep.h``                  |
   |              +----------------------------------------------+
   |              | ``arkode/arkode_forcingstep.h``              |
   |              +----------------------------------------------+
   |              | ``arkode/arkode_lsrkstep.h``                 |
   |              +----------------------------------------------+
   |              | ``arkode/arkode_mristep.h``                  |
   |              +----------------------------------------------+
   |              | ``arkode/arkode_splittingstep.h``            |
   |              +----------------------------------------------+
   |              | ``arkode/arkode_sprkstep.h``                 |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::arkode``                         |
   +--------------+----------------------------------------------+

The ARKODE module header files include the header file for the shared ARKODE
interface functions, constants, and types (``arkode.h``). As appropriate, the
module header files also include the ARKODE linear solver interface as well as
the header files defining method coefficients.

.. table:: Additional header files included by ``arkode_*step.h`` header files
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``arkode/arkode.h``                          |
   |              +----------------------------------------------+
   |              | ``arkode/arkode_butcher.h``                  |
   |              +----------------------------------------------+
   |              | ``arkode/arkode_butcher_dirk.h``             |
   |              +----------------------------------------------+
   |              | ``arkode/arkode_butcher_erk.h``              |
   |              +----------------------------------------------+
   |              | ``arkode/arkode_ls.h``                       |
   |              +----------------------------------------------+
   |              | ``arkode/arkode_sprk.h``                     |
   +--------------+----------------------------------------------+

For problems in which the user cannot define a more effective, problem-specific
preconditioner for Krylov iterative linear solvers, ARKODE provides banded
(``bandpre``) and band-block-diagonal (``bbdpre``) preconditioner
modules. Include the header files below to access the related functions.

.. table:: ARKODE preconditioner modules
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``arkode/arkode_bandpre.h``                  |
   |              +----------------------------------------------+
   |              | ``arkode/arkode_bbdpre.h``                   |
   +--------------+----------------------------------------------+

When XBraid support is enabled (:cmakeop:`ENABLE_XBRAID` is ``ON``), include the
ARKODE-XBraid interface header file and link to the interface library given
below to use ARKODE and XBraid together.

.. table:: ARKODE library, header, and CMake target for interfacing with XBraid
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_arkode_xbraid.LIB``            |
   +--------------+----------------------------------------------+
   | Headers      | ``arkode/arkode_xbraid.h``                   |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::arkode_xbraid``                  |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Packages.IDA:

IDA
"""

To use the :ref:`IDA <IDA>` package, include the header file and link to the
library given below.

.. table:: IDA library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_ida.LIB``                      |
   +--------------+----------------------------------------------+
   | Headers      | ``ida/ida.h``                                |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::ida``                            |
   +--------------+----------------------------------------------+

The IDA header file includes the header file below which defines functions,
types, and constants for the IDA linear solver interface.

.. table:: Additional header files included by ``ida.h``
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``ida/ida_ls.h``                             |
   +--------------+----------------------------------------------+

For problems in which the user cannot define a more effective, problem-specific
preconditioner for Krylov iterative linear solvers, IDA provides a
band-block-diagonal (``bbdpre``) preconditioner module. Include the header
file below to access the related functions.

.. table:: IDA preconditioner modules
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``ida/ida_bbdpre.h``                         |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Packages.IDAS:

IDAS
""""

To use the :ref:`IDAS <IDAS>` package, include the header file and link to the
library given below.

.. warning::

   IDAS is a superset of IDA and defines the same functions as provided by
   IDA. As such, applications should not link to both IDAS and IDA.


.. table:: IDAS library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_idas.LIB``                     |
   +--------------+----------------------------------------------+
   | Headers      | ``idas/idas.h``                              |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::idas``                           |
   +--------------+----------------------------------------------+

The IDAS header file includes the header file below which defines functions,
types, and constants for the IDAS linear solver interface.

.. table:: Additional header files included by ``idas.h``
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``idas/idas_ls.h``                           |
   +--------------+----------------------------------------------+

For problems in which the user cannot define a more effective, problem-specific
preconditioner for Krylov iterative linear solvers, IDAS provides a
band-block-diagonal (``bbdpre``) preconditioner module. Include the header
file below to access the related functions.

.. table:: IDAS preconditioner modules
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``idas/idas_bbdpre.h``                       |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Packages.KINSOL:

KINSOL
""""""

To use the :ref:`KINSOL <KINSOL>` package, include the header file and link to
the library given below.

.. table:: KINSOL library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_kinsol.LIB``                   |
   +--------------+----------------------------------------------+
   | Headers      | ``kinsol/kinsol.h``                          |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::kinsol``                         |
   +--------------+----------------------------------------------+

The KINSOL header file includes the header file below which defines functions,
types, and constants for the KINSOL linear solver interface.

.. table:: Additional header files included by ``kinsol.h``
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``kinsol/kinsol_ls.h``                       |
   +--------------+----------------------------------------------+

For problems in which the user cannot define a more effective, problem-specific
preconditioner for Krylov iterative linear solvers, KINSOL provides a
band-block-diagonal (``bbdpre``) preconditioner module. Include the header
file below to access the related functions.

.. table:: KINSOL preconditioner modules
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``kinsol/kinsol_bbdpre.h``                   |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector:

Vectors
^^^^^^^

.. _Installation.LibrariesAndHeaders.Vector.Serial:

Serial
""""""

To use the :ref:`serial NVector <NVectors.NVSerial>`, include the header file
and link to the library given below.

When using SUNDIALS time integration packages or the KINSOL package, the serial
NVector is bundled with the package library and it is not necessary to link to
the library below when using those packages.

.. table:: The serial NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvecserial.LIB``               |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_serial.h``                 |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvecserial``                     |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.ManyVector:

ManyVector
""""""""""

To use the :ref:`ManyVector NVector <NVectors.ManyVector>`, include the header
file and link to the library given below.

.. table:: The ManyVector NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvecmanyvector.LIB``           |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_manyvector.h``             |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvecmanyvector``                 |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.Parallel:

Parallel (MPI)
""""""""""""""

To use the :ref:`parallel (MPI) NVector <NVectors.NVParallel>`, include the
header file and link to the library given below.

.. table:: The parallel (MPI) NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvecparallel.LIB``             |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_parallel.h``               |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvecparallel``                   |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.MPIManyVector:

MPI ManyVector
""""""""""""""

To use the :ref:`MPI ManyVector NVector <NVectors.MPIManyVector>`, include the
header file and link to the library given below.

.. table:: The MPI ManyVector NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvecmpimanyvector.LIB``        |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_mpimanyvector.h``          |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvecmpimanyvector``              |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.MPIPlusX:

MPI+X
"""""

To use the :ref:`MPI+X NVector <NVectors.MPIPlusX>`, include the header
file and link to the library given below.

.. table:: The MPI+X NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvecmpiplusx.LIB``             |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_mpiplusx.h``               |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvecmpiplusx``                   |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.OpenMP:

OpenMP
""""""

To use the :ref:`OpenMP NVector <NVectors.OpenMP>`, include the header file and
link to the library given below.

.. table:: The OpenMP NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvecopenmp.LIB``               |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_openmp.h``                 |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvecopenmp``                     |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.OpenMPDEV:

OpenMPDEV
"""""""""

To use the :ref:`OpenMP device offload NVector <NVectors.OpenMPDEV>`, include
the header file and link to the library given below.

.. table:: The OpenMP device offload NVector library, header file, and CMake
           target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvecopenmpdev.LIB``            |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_openmpdev.h``              |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvecopenmpdev``                  |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.PThreads:

PThreads
""""""""

To use the :ref:`POSIX Threads NVector <NVectors.PThreads>`, include the header
file and link to the library given below.

.. table:: The POSIX Threads NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvecpthreads.LIB``             |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_pthreads.h``               |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvecpthreads``                   |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.ParHyp:

*hypre* (ParHyp)
""""""""""""""""

To use the :ref:`hypre (ParHyp) NVector <NVectors.ParHyp>`, include the header
file and link to the library given below.

.. table:: The *hypre* (ParHyp) NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvecparhyp.LIB``               |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_parhyp.h``                 |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvecparhyp``                     |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.PETSc:

PETSc
"""""

To use the :ref:`PETSc NVector <NVectors.ParHyp>`, include the header file and
link to the library given below.

.. table:: The PETSc NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvecpetsc.LIB``                |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_petsc.h``                  |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvecpetsc``                      |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.CUDA:

CUDA
""""

To use the :ref:`CUDA NVector <NVectors.CUDA>`, include the header
file and link to the library given below.

.. table:: The CUDA NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nveccuda.LIB``                 |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_cuda.h``                   |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nveccuda``                       |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.HIP:

HIP
"""

To use the :ref:`HIP NVector <NVectors.HIP>`, include the header
file and link to the library given below.

.. table:: The HIP NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvechip.LIB``                  |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_hip.h``                    |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvechip``                        |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.RAJA:

RAJA
""""

To use the :ref:`RAJA NVector <NVectors.RAJA>`, include the header
file and link to the library given below for the desired backend.

.. table:: The RAJA NVector libraries, header file, and CMake targets
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nveccudaraja.LIB``             |
   |              +----------------------------------------------+
   |              | ``libsundials_nvechipraja.LIB``              |
   |              +----------------------------------------------+
   |              | ``libsundials_nvecsyclraja.LIB``             |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_raja.h``                   |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nveccudaraja``                   |
   |              +----------------------------------------------+
   |              | ``SUNDIALS::nvechipraja``                    |
   |              +----------------------------------------------+
   |              | ``SUNDIALS::nvecsyclraja``                   |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.SYCL:

SYCL
""""

To use the :ref:`SYCL NVector <NVectors.SYCL>`, include the header
file and link to the library given below.

.. table:: The SYCL NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvecsycl.LIB``                 |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_sycl.h``                   |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvecsycl``                       |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.Trilinos:

Trilinos (Tpetra)
"""""""""""""""""

To use the :ref:`Trilinos (Tpetra) NVector <NVectors.NVTrilinos>`, include the
header file and link to the library given below.

.. table:: The Trilinos (Tpetra) NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_nvectrilinos.LIB``             |
   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_trilinos.h``               |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nvectrilinos``                   |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Vector.Kokkos:

Kokkos
""""""

To use the :ref:`Kokkos NVector <NVectors.Kokkos>`, include the header
file and link to the library given below.

.. table:: The Kokkos NVector library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``nvector/nvector_kokkos.hpp``               |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::nveckokkos``                     |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Matrix:

Matrices
^^^^^^^^

.. _Installation.LibrariesAndHeaders.Matrix.Band:

Banded
""""""

To use the :ref:`banded SUNMatrix <SUNMatrix.Band>`, include the header file and
link to the library given below.

When using SUNDIALS time integration packages or the KINSOL package, the banded
SUNMatrix is bundled with the package library and it is not necessary to link to
the library below when using those packages.

.. table:: The banded SUNMatrix library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunmatrixband.LIB``            |
   +--------------+----------------------------------------------+
   | Headers      | ``sunmatrix/sunmatrix_band.h``               |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunmatrixband``                  |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Matrix.cuSPARSE:

cuSPARSE
""""""""

To use the :ref:`cuSPARSE SUNMatrix <SUNMatrix.cuSPARSE>`, include the header
file and link to the library given below.

.. table:: The cuSPARSE SUNMatrix library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunmatrixcusparse.LIB``        |
   +--------------+----------------------------------------------+
   | Headers      | ``sunmatrix/sunmatrix_cusparse.h``           |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunmatrixcusparse``              |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Matrix.Dense:

Dense
"""""

To use the :ref:`dense SUNMatrix <SUNMatrix.Dense>`, include the header file and
link to the library given below.

When using SUNDIALS time integration packages or the KINSOL package, the dense
SUNMatrix is bundled with the package library and it is not necessary to link to
the library below when using those packages.

.. table:: The dense SUNMatrix library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunmatrixdense.LIB``           |
   +--------------+----------------------------------------------+
   | Headers      | ``sunmatrix/sunmatrix_dense.h``              |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunmatrixdense``                 |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Matrix.Ginkgo:

Ginkgo
""""""

To use the :ref:`Ginkgo SUNMatrix <SUNMatrix.Ginkgo>` or :ref:`Ginkgo Batch
SUNMatrix <SUNMatrix.GinkgoBatch>`, include the corresponding header file given
below.

.. table:: The Ginkgo and Ginkgo Batch SUNMatrix header files and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sunmatrix/sunmatrix_ginkgo.hpp``           |
   |              +----------------------------------------------+
   |              | ``sunmatrix/sunmatrix_ginkgobatch.hpp``      |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunmatrixginkgo``                |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Matrix.KokkosKernels:

KokkosKernels Dense
"""""""""""""""""""

To use the :ref:`KokkosKernels dense SUNMatrix <SUNMatrix.Kokkos>`, include the
header file given below.

.. table:: The dense KokkosKernels SUNMatrix library, header file, and CMake
           target
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sunmatrix/sunmatrix_kokkosdense.hpp``      |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunmatrixkokkosdense``           |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Matrix.MAGMADense:

MAGMA Dense
"""""""""""

To use the :ref:`MAGMA dense SUNMatrix <SUNMatrix.MAGMADense>`, include the
header file and link to the library given below.

.. table:: The dense MAGMA SUNMatrix library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunmatrixmagmadense.LIB``      |
   +--------------+----------------------------------------------+
   | Headers      | ``sunmatrix/sunmatrix_magmadense.h``         |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunmatrixmagmadense``            |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Matrix.oneMKLDense:

oneMKL Dense
""""""""""""

To use the :ref:`oneMKL dense SUNMatrix <SUNMatrix.OneMklDense>`, include the
header file and link to the library given below.

.. table:: The dense oneMKL SUNMatrix library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunmatrixonemkldense.LIB``     |
   +--------------+----------------------------------------------+
   | Headers      | ``sunmatrix/sunmatrix_onemkldense.h``        |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunmatrixonemkldense``           |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Matrix.Sparse:

Sparse
""""""

To use the :ref:`sparse SUNMatrix <SUNMatrix.Sparse>`, include the header file
and link to the library given below.

When using SUNDIALS time integration packages or the KINSOL package, the sparse
SUNMatrix is bundled with the package library and it is not necessary to link to
the library below when using those packages.

.. table:: The sparse SUNMatrix library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunmatrixsparse.LIB``          |
   +--------------+----------------------------------------------+
   | Headers      | ``sunmatrix/sunmatrix_sparse.h``             |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunmatrixsparse``                |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.Matrix.SuperLU_DIST:

SuperLU_DIST (SLUNRloc)
"""""""""""""""""""""""

To use the :ref:`SuperLU_DIST (SLUNRloc) SUNMatrix <SUNMatrix.SLUNRloc>`,
include the header file and link to the library given below.

.. table:: The SuperLU_DIST (SLUNRloc) SUNMatrix library, header file, and CMake
           target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunmatrixslunrloc.LIB``        |
   +--------------+----------------------------------------------+
   | Headers      | ``sunmatrix/sunmatrix_slunrloc.h``           |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunmatrixslunrloc``              |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver:

Linear Solvers
^^^^^^^^^^^^^^

.. _Installation.LibrariesAndHeaders.LinearSolver.Band:

Banded
""""""

To use the :ref:`banded SUNLinearSolver <SUNLinSol_Band>`, include the header
file and link to the library given below.

When using SUNDIALS time integration packages or the KINSOL package, the banded
SUNLinearSolver is bundled with the package library and it is not necessary to
link to the library below when using those packages.

.. table:: The banded SUNLinearSolver library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsolband.LIB``            |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_band.h``               |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolband``                  |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.cuSPARSE:

cuSPARSE Batched QR
"""""""""""""""""""

To use the :ref:`cuSPARSE batched QR SUNLinearSolver <SUNLinSol.cuSolverSp>`,
include the header file and link to the library given below.

.. table:: The cuSPARSE batched QR SUNLinearSolver library, header file, and
           CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsolcusolversp.LIB``      |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_cusolversp_batchqr.h`` |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolcusolversp``            |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.Dense:

Dense
"""""

To use the :ref:`dense SUNLinearSolver <SUNLinSol_Dense>`, include the header
file and link to the library given below.

When using SUNDIALS time integration packages or the KINSOL package, the dense
SUNLinearSolver is bundled with the package library and it is not necessary to
link to the library below when using those packages.

.. table:: The dense SUNLinearSolver library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsoldense.LIB``           |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_dense.h``              |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsoldense``                 |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.Ginkgo:

Ginkgo
""""""

To use the :ref:`Ginkgo SUNLinearSolver <SUNLinSol.Ginkgo>` or :ref:`Ginkgo
Batch SUNLinearSolver <SUNLinSol.Ginkgo>`, include the corresponding header file
given below.

.. table:: The Ginkgo and Ginkgo Batch SUNLinearSolver header files and CMake
           target
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_ginkgo.hpp``           |
   |              +----------------------------------------------+
   |              | ``sunlinsol/sunlinsol_ginkgobatch.hpp``      |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolginkgo``                |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.KLU:

KLU
"""

To use the :ref:`KLU SUNLinearSolver <SUNLinSol.KLU>`, include the header file
and link to the library given below.

.. table:: The KLU SUNLinearSolver library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsolklu.LIB``             |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_klu.h``                |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolklu``                   |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.KokkosKernels:

KokkosKernels Dense
"""""""""""""""""""

To use the :ref:`KokkosKernels dense SUNLinearSolver <SUNLinSol.Kokkos>`, include the
header file given below.

.. table:: The KokkosKernels dense SUNLinearSolver header file and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_kokkosdense.hpp``      |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolkokkosdense``           |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.LAPACKBand:

LAPACK Banded
"""""""""""""

To use the :ref:`LAPACK banded SUNLinearSolver <SUNLinSol_LapackBand>`, include
the header file and link to the library given below.

.. table:: The LAPACK banded SUNLinearSolver library, header file, and CMake
           target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsollapackband.LIB``      |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_lapackband.h``         |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsollapackband``            |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.LAPACKDense:

LAPACK Dense
""""""""""""

To use the :ref:`LAPACK dense SUNLinearSolver <SUNLinSol_LapackDense>`, include
the header file and link to the library given below.

.. table:: The LAPACK dense SUNLinearSolver library, header file, and CMake
           target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsollapackdense.LIB``     |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_lapackdense.h``        |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsollapackdense``           |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.MAGMADense:

MAGMA Dense
"""""""""""

To use the :ref:`MAGMA dense SUNLinearSolver <SUNLinSol.MAGMADense>`, include
the header file and link to the library given below.

.. table:: The MAGMA dense SUNLinearSolver library, header file, and CMake
           target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsolmagmadense.LIB``      |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_magmadense.h``         |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolmagmadense``            |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.oneMKLDense:

oneMKL Dense
""""""""""""

To use the :ref:`oneMKL dense SUNLinearSolver <SUNLinSol.OneMklDense>`, include
the header file and link to the library given below.

.. table:: The oneMKL dense SUNLinearSolver library, header file, and CMake
           target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsolonemkldense.LIB``     |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_onemkldense.h``        |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolonemkldense``           |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.PCG:

Preconditioned Conjugate Gradient (PCG)
"""""""""""""""""""""""""""""""""""""""

To use the :ref:`PCG SUNLinearSolver <SUNLinSol.PCG>`, include the header file
and link to the library given below.

When using SUNDIALS time integration packages or the KINSOL package, the PCG
SUNLinearSolver is bundled with the package library and it is not necessary to
link to the library below when using those packages.

.. table:: The PCG SUNLinearSolver library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsolpcg.LIB``             |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_pcg.h``                |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolpcg``                   |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.SPBCGS:

Scaled, Preconditioned Bi-Conjugate Gradient, Stabilized (SPBCGS)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To use the :ref:`SPBCGS SUNLinearSolver <SUNLinSol.SPBCGS>`, include the header
file and link to the library given below.

When using SUNDIALS time integration packages or the KINSOL package, the SPBCGS
SUNLinearSolver is bundled with the package library and it is not necessary to
link to the library below when using those packages.

.. table:: The SPBCGS SUNLinearSolver library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsolspbcgs.LIB``          |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_spbcgs.h``             |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolspbcgs``                |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.SPFGMR:

Scaled, Preconditioned, Flexible, Generalized Minimum Residual (SPFGMR)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To use the :ref:`SPFGMR SUNLinearSolver <SUNLinSol.SPFGMR>`, include the header
file and link to the library given below.

When using SUNDIALS time integration packages or the KINSOL package, the SPFGMR
SUNLinearSolver is bundled with the package library and it is not necessary to
link to the library below when using those packages.

.. table:: The SPFGMR SUNLinearSolver library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsolspfgmr.LIB``          |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_spfgmr.h``             |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolspfgmr``                |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.SPGMR:

Scaled, Preconditioned, Generalized Minimum Residual (SPGMR)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To use the :ref:`SPGMR SUNLinearSolver <SUNLinSol.SPGMR>`, include the header
file and link to the library given below.

When using SUNDIALS time integration packages or the KINSOL package, the SPGMR
SUNLinearSolver is bundled with the package library and it is not necessary to
link to the library below when using those packages.

.. table:: The SPGMR SUNLinearSolver library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsolspgmr.LIB``           |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_spgmr.h``              |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolspgmr``                 |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.SPTFQMR:

Scaled, Preconditioned, Transpose-Free Quasi-Minimum Residual (SPTFQMR)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To use the :ref:`SPTFQMR SUNLinearSolver <SUNLinSol.SPTFQMR>`, include the
header file and link to the library given below.

When using SUNDIALS time integration packages or the KINSOL package, the SPTFQMR
SUNLinearSolver is bundled with the package library and it is not necessary to
link to the library below when using those packages.

.. table:: The SPTFQMR SUNLinearSolver library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsolsptfqmr.LIB``         |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_sptfqmr.h``            |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolsptfqmr``               |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.SuperLU_DIST:

SuperLU_DIST
""""""""""""

To use the :ref:`SuperLU_DIST SUNLinearSolver <SUNLinSol.SuperLUDIST>`, include
the header file and link to the library given below.

.. table:: The SuperLU_DIST SUNLinearSolver library, header file, and CMake
           target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsolsuperludist.LIB``     |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_superludist.h``        |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolsuperludist``           |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.LinearSolver.SuperLU_MT:

SuperLU_MT
""""""""""

To use the :ref:`SuperLU_MT SUNLinearSolver <SUNLinSol.SuperLUMT>`, include
the header file and link to the library given below.

.. table:: The SuperLU_MT SUNLinearSolver library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunlinsolsuperlumt.LIB``       |
   +--------------+----------------------------------------------+
   | Headers      | ``sunlinsol/sunlinsol_superlumt.h``          |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunlinsolsuperlumt``             |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.NonlinearSolver:

Nonlinear Solvers
^^^^^^^^^^^^^^^^^

.. _Installation.LibrariesAndHeaders.NonlinearSolver.Newton:

Newton
""""""

To use the :ref:`Newton SUNNonlinearSolver <SUNNonlinSol.Newton>`, include the
header file and link to the library given below.

When using SUNDIALS time integration packages, the Newton SUNNonlinearSolver is
bundled with the package library and it is not necessary to link to the library
below when using those packages.

.. table:: The Newton SUNNonlinearSolver library, header file, and CMake target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunnonlinsolnewton.LIB``       |
   +--------------+----------------------------------------------+
   | Headers      | ``sunnonlinsol/sunnonlinsol_newton.h``       |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunnonlinsolnewton``             |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.NonlinearSolver.FixedPoint:

Fixed-point
"""""""""""

To use the :ref:`fixed-point SUNNonlinearSolver <SUNNonlinSol.FixedPoint>`,
include the header file and link to the library given below.

When using SUNDIALS time integration packages, the fixed-point
SUNNonlinearSolver is bundled with the package library and it is not necessary
to link to the library below when using those packages.

.. table:: The Fixed-point SUNNonlinearSolver library, header file, and CMake
           target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunnonlinsolfixedpoint.LIB``   |
   +--------------+----------------------------------------------+
   | Headers      | ``sunnonlinsol/sunnonlinsol_fixedpoint.h``   |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunnonlinsolfixedpoint``         |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.NonlinearSolver.PETScSNES:

PETSc SNES
""""""""""

To use the :ref:`PETSc SNES SUNNonlinearSolver <SUNNonlinSol.PetscSNES>`, include
the header file and link to the library given below.

.. table:: The PETSc SNES SUNNonlinearSolver library, header file, and CMake
           target
   :align: center

   +--------------+----------------------------------------------+
   | Libraries    | ``libsundials_sunnonlinsolpetscsnes.LIB``    |
   +--------------+----------------------------------------------+
   | Headers      | ``sunnonlinsol/sunnonlinsol_petscsnes.h``    |
   +--------------+----------------------------------------------+
   | CMake target | ``SUNDIALS::sunnonlinsolpetscsnes``          |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.MemoryHelper:

Memory Helpers
^^^^^^^^^^^^^^

System
""""""

When using SUNDIALS time integration packages or the KINSOL package, the system
SUNMemoryHelper is bundled with the package library and it is not necessary to
link to the library below when using those packages.

.. table:: SUNDIALS system memory helper header file
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sunmemory/sunmemory_system.h``             |
   +--------------+----------------------------------------------+

CUDA
""""

To use the :ref:`CUDA SUNMemoryHelper <SUNMemory.CUDA>`, include the header file
given below when using a CUDA-enabled NVector or SUNMatrix.

.. table:: SUNDIALS CUDA memory helper header file
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sunmemory/sunmemory_cuda.h``               |
   +--------------+----------------------------------------------+

HIP
"""

To use the :ref:`HIP SUNMemoryHelper <SUNMemory.HIP>`, include the header file
given below when using a HIP-enabled NVector or SUNMatrix.

.. table:: SUNDIALS HIP memory helper header file
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sunmemory/sunmemory_hip.h``                |
   +--------------+----------------------------------------------+

SYCL
""""

To use the :ref:`SYCL SUNMemoryHelper <SUNMemory.SYCL>`, include the header file
given below when using a SYCL-enabled NVector or SUNMatrix.

.. table:: SUNDIALS SYCL memory helper header file
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sunmemory/sunmemory_sycl.h``               |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.ExecutionPolicies:

Execution Policies
^^^^^^^^^^^^^^^^^^

CUDA
""""

When using a CUDA-enabled NVector or SUNMatrix, include the header file below to
access the CUDA execution policy C++ classes.

.. table:: SUNDIALS CUDA execution policies header file
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sundials/sundials_cuda_policies.hpp``      |
   +--------------+----------------------------------------------+

HIP
"""

When using a HIP-enabled NVector or SUNMatrix, include the header file below to
access the HIP execution policy C++ classes.

.. table:: SUNDIALS HIP execution policies header file
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sundials/sundials_hip_policies.hpp``       |
   +--------------+----------------------------------------------+

SYCL
""""

When using a SYCL-enabled NVector or SUNMatrix, include the header file below to
access the SYCL execution policy C++ classes.

.. table:: SUNDIALS SYCL execution policies header file
   :align: center

   +--------------+----------------------------------------------+
   | Headers      | ``sundials/sundials_sycl_policies.hpp``      |
   +--------------+----------------------------------------------+

.. _Installation.LibrariesAndHeaders.AdjointSensitivityCheckpointing:

Adjoint Sensitivity Checkpointing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _Installation.LibrariesAndHeaders.ASACheckpointing.Fixed:

Fixed ASA Checkpointing
"""""""""""""""""""""""

For fixed-interval adjoint checkpointing, include the header file below:

.. table:: SUNDIALS fixed adjoint checkpointing header files
   :align: center

   +--------------+---------------------------------------------------------------------+
   | Headers      | ``sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h``   |
   +--------------+---------------------------------------------------------------------+

.. _Installation.LibrariesAndHeaders.DomEigEst:

Dominant Eigenvalue Estimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _Installation.LibrariesAndHeaders.DomEigEst.Power:

Power Iteration
"""""""""""""""

To use the :ref:`Power iteration SUNDomEigEstimator <SUNDomEigEst.Power>`,
include the header file and link to the library given below.

.. table:: The SUNDIALS Power iteration SUNDomEigEstimator library, header file,
           and CMake target
   :align: center

   +--------------+---------------------------------------+
   | Libraries    | ``libsundials_sundomeigestpower.LIB`` |
   +--------------+---------------------------------------+
   | Headers      | ``sundomeigest/sundomeigest_power.h`` |
   +--------------+---------------------------------------+
   | CMake target | ``SUNDIALS::sundomeigestpower``       |
   +--------------+---------------------------------------+

.. _Installation.LibrariesAndHeaders.DomEigEst.Arnoldi:

Arnoldi Iteration
"""""""""""""""""

To use the :ref:`Arnoldi iteration SUNDomEigEstimator <SUNDomEigEst.Arnoldi>`,
include the header file and link to the library given below.

.. table:: The SUNDIALS Arnoldi iteration SUNDomEigEstimator library, header
           file, and CMake target
   :align: center

   +--------------+-----------------------------------------+
   | Libraries    | ``libsundials_sundomeigestarnoldi.LIB`` |
   +--------------+-----------------------------------------+
   | Headers      | ``sundomeigest/sundomeigest_arnoldi.h`` |
   +--------------+-----------------------------------------+
   | CMake target | ``SUNDIALS::sundomeigestarnoldi``       |
   +--------------+-----------------------------------------+
