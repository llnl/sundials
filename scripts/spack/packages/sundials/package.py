# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

import os
import re
import sys

from llnl.util import tty

from spack.package import *


class Sundials(CachedCMakePackage, CudaPackage, ROCmPackage):
    """SUNDIALS (SUite of Nonlinear and DIfferential/ALgebraic equation
    Solvers)"""

    homepage = "https://computing.llnl.gov/projects/sundials"
    url = "https://github.com/LLNL/sundials/releases/download/v2.7.0/sundials-2.7.0.tar.gz"
    git = "https://github.com/llnl/sundials.git"
    tags = ["radiuss", "e4s"]
    test_requires_compiler = True

    maintainers = ["balos1", "cswoodward", "gardner48"]

    # ==========================================================================
    # Versions
    # ==========================================================================
    version("develop", branch="develop")
    version("6.5.1", sha256="4252303805171e4dbdd19a01e52c1dcfe0dafc599c3cfedb0a5c2ffb045a8a75")
    version("6.5.0", sha256="4e0b998dff292a2617e179609b539b511eb80836f5faacf800e688a886288502")
    version("6.4.1", sha256="7bf10a8d2920591af3fba2db92548e91ad60eb7241ab23350a9b1bc51e05e8d0")
    version("6.4.0", sha256="0aff803a12c6d298d05b56839197dd09858631864017e255ed89e28b49b652f1")
    version("6.3.0", sha256="89a22bea820ff250aa7239f634ab07fa34efe1d2dcfde29cc8d3af11455ba2a7")
    version("6.2.0", sha256="195d5593772fc483f63f08794d79e4bab30c2ec58e6ce4b0fb6bcc0e0c48f31d")
    version("6.1.1", sha256="cfaf637b792c330396a25ef787eb59d58726c35918ebbc08e33466e45d50470c")
    version("6.1.0", sha256="eea49f52140640e54931c779e73aece65f34efa996a26b2263db6a1e27d0901c")
    version("6.0.0", sha256="c7178e54df20a9363ae3e5ac5b3ee9db756a4ddd4b8fff045127e93b73b151f4")
    version("5.8.0", sha256="d4ed403351f72434d347df592da6c91a69452071860525385b3339c824e8a213")
    version("5.7.0", sha256="48da7baa8152ddb22aed1b02d82d1dbb4fbfea22acf67634011aa0303a100a43")
    version("5.6.1", sha256="16b77999ec7e7f2157aa1d04ca1de4a2371ca8150e056d24951d0c58966f2a83")
    version("5.6.0", sha256="95e4201912e150f29c6f6f7625de763385e2073dae7f929c4a544561ea29915d")
    version("5.5.0", sha256="2a755e89aab96d2ff096a4e30bf00bb162e80be20e9e99f424dccfb249098237")
    version("5.4.0", sha256="04d8a2ebe02cdaeef5a9e22ff7e3146bb563d8400f65772b6c7af80001413ffa")
    version("5.3.0", sha256="88dff7e11a366853d8afd5de05bf197a8129a804d9d4461fb64297f1ef89bca7")
    version("5.2.0", sha256="95f058acce5bd66e654de65acdbb1c9f44c90cf1b4e28f8d933cdb4415ebba3e")
    version("5.1.0", sha256="fb22d14fad42203809dc46d046b001149ec4e901b23882bd4a80619157fd9b21")
    version("5.0.0", sha256="345141ec01c641d0bdfb3476c478b7e74fd6a7192a478a27cafe75d9da2d7dd3")
    version("4.1.0", sha256="280de1c27b2360170a6f46cb3799b2aee9dff3bddbafc8b08c291a47ab258aa5")
    version("4.0.1", sha256="29e409c8620e803990edbda1ebf49e03a38c08b9187b90658d86bddae913aed4")
    version("3.2.1", sha256="47d94d977ab2382cdcdd02f72a25ebd4ba8ca2634bbb2f191fe1636e71c86808")
    version("3.2.0", sha256="d2b690afecadf8b5a048bb27ab341de591d714605b98d3518985dfc2250e93f9")
    version("3.1.2", sha256="a8985bb1e851d90e24260450667b134bc13d71f5c6effc9e1d7183bd874fe116")
    version("3.1.1", sha256="a24d643d31ed1f31a25b102a1e1759508ce84b1e4739425ad0e18106ab471a24")
    version("3.1.0", sha256="18d52f8f329626f77b99b8bf91e05b7d16b49fde2483d3a0ea55496ce4cdd43a")
    version("3.0.0", sha256="28b8e07eecfdef66e2c0d0ea0cb1b91af6e4e94d71008abfe80c27bf39f63fde")
    version("2.7.0", sha256="d39fcac7175d701398e4eb209f7e92a5b30a78358d4a0c0fcc23db23c11ba104")
    version("2.6.2", sha256="d8ed0151509dd2b0f317b318a4175f8b95a174340fc3080b8c20617da8aa4d2f")

    # ==========================================================================
    # Variants
    # ==========================================================================

    # SUNDIALS solvers
    sun_solvers = ["CVODE", "CVODES", "ARKODE", "IDA", "IDAS", "KINSOL"]

    for pkg in sun_solvers:
        variant(pkg, default=True, description="Enable %s solver" % pkg)

    # Language standards
    variant(
        "cstd", default="99", description="C language standard", values=("90", "99", "11", "17")
    )

    variant("cxxstd", default="14", description="C++ language standard", values=("14", "17"))

    # Logging (default=0 when "@6.2.0:6.7.0", default=2 when "@7.0.0:")
    variant(
        "logging-level",
        default="2",
        description="logging level\n 0 = no logging,\n 1 = errors,\n "
        "2 = errors + warnings,\n 3 = errors + "
        "warnings + info,\n 4 = errors + warnings + info + debugging, "
        "\n 5 = all of the above + more",
        values=("0", "1", "2", "3", "4", "5"),
        multi=False,
        when="@6.2.0:",
    )

    # MPI logging (option removed in 7.0)
    variant(
        "logging-mpi",
        default="OFF",
        description="enable MPI support in the logger",
        when="@6.2.0:6.7.0",
    )

    # Real type
    variant(
        "precision",
        default="double",
        description="real type precision",
        values=("single", "double", "extended"),
        multi=False,
    )

    # Index type
    variant("int64", default=False, description="Use 64bit integers for indices")

    # Parallelism
    variant("mpi", default=True, description="Enable MPI parallel vector")
    variant("openmp", default=False, description="Enable OpenMP parallel vector")
    variant("pthread", default=False, description="Enable Pthreads parallel vector")
    variant("raja", default=False, when="@3.0.0:", description="Enable RAJA vector")
    variant("sycl", default=False, when="@5.7.0:", description="Enable SYCL vector")

    # External libraries
    variant("adiak", default=False, when="@6.6.0:", description="Enable Adiak interfaces")
    variant(
        "caliper",
        default=False,
        when="@6.0.0: +profiling",
        description="Enable Caliper instrumentation/profiling",
    )
    variant("ginkgo", default=False, when="@6.4.0:", description="Enable Ginkgo interfaces")
    variant("hypre", default=False, when="@2.7.0:", description="Enable Hypre MPI parallel vector")
    variant("kokkos", default=False, when="@6.4.0:", description="Enable Kokkos vector")
    variant(
        "kokkos-kernels",
        default=False,
        when="@6.4.0:",
        description="Enable KokkosKernels based matrix and linear solver",
    )
    variant("klu", default=False, description="Enable KLU sparse, direct solver")
    variant("lapack", default=False, description="Enable LAPACK direct solvers")
    variant("petsc", default=False, when="@2.7.0:", description="Enable PETSc interfaces")
    variant("magma", default=False, when="@5.7.0:", description="Enable MAGMA interface")
    variant("superlu-mt", default=False, description="Enable SuperLU_MT sparse, direct solver")
    variant(
        "superlu-dist",
        default=False,
        when="@5.0.0:",
        description="Enable SuperLU_DIST sparse, direct solver",
    )
    variant("trilinos", default=False, when="@5.0.0:", description="Enable Trilinos interfaces")

    # Library type
    variant("shared", default=True, description="Build shared libraries")
    variant("static", default=True, description="Build static libraries")

    # Fortran interfaces
    variant("fcmix", default=False, description="Enable Fortran 77 interface")
    variant("f2003", default=False, description="Enable Fortran 2003 interface")

    # Examples
    variant("examples", default=True, description="Enable examples")
    variant("examples-install", default=True, description="Install examples")

    # Generic (std-c) math libraries (UNIX only)
    variant(
        "generic-math",
        default=True,
        description="Use generic (std-c) math libraries on unix systems",
    )

    # Monitoring
    variant(
        "monitoring",
        default=False,
        when="@5.5.0:",
        description="Build with simulation monitoring capabilities",
    )

    # Profiling
    variant(
        "profiling", default=False, when="@6.0.0:", description="Build with profiling capabilities"
    )

    # Scheduler
    variant(
        "scheduler",
        default="slurm",
        description="Specify which scheduler the system runs on",
        values=("flux", "lsf", "slurm"),
    )

    # Benchmarking
    variant("benchmarks", default=False, description="Build benchmark programs")

    # Profiling examples
    variant(
        "profile-examples",
        default=False,
        when="+adiak +caliper",
        description="Build examples with profiling capabilities",
    )

    # Caliper Directory
    variant("caliper-dir", default="none", description="Specify where to place Caliper files")

    # ==========================================================================
    # Dependencies
    # ==========================================================================

    # Build dependencies
    depends_on("cmake@3.12:", when="~cuda", type="build")
    depends_on("cmake@3.18:", when="+cuda", type="build")

    # MPI related dependencies
    depends_on("mpi", when="+mpi")
    depends_on("mpi", when="+hypre")
    depends_on("mpi", when="+petsc")
    depends_on("mpi", when="+superlu-dist")

    # Other parallelism dependencies
    depends_on("raja", when="+raja")
    depends_on("raja+cuda", when="+raja +cuda")
    depends_on("raja+rocm", when="+raja +rocm")

    # External libraries
    depends_on("adiak", when="+adiak")
    depends_on("caliper", when="+caliper")
    depends_on("ginkgo@1.9.0:", when="+ginkgo")
    depends_on("kokkos", when="+kokkos")
    depends_on("kokkos-kernels", when="+kokkos-kernels")
    for cuda_arch in CudaPackage.cuda_arch_values:
        depends_on(
            "kokkos+cuda+cuda_lambda+cuda_constexpr cuda_arch=%s" % cuda_arch,
            when="+kokkos +cuda cuda_arch=%s" % cuda_arch,
        )
        depends_on(
            "kokkos-kernels+cuda cuda_arch=%s" % cuda_arch,
            when="+kokkos-kernels +cuda cuda_arch=%s" % cuda_arch,
        )
    for rocm_arch in ROCmPackage.amdgpu_targets:
        depends_on(
            "kokkos+rocm amdgpu_target=%s" % rocm_arch,
            when="+kokkos +rocm amdgpu_target=%s" % rocm_arch,
        )
    depends_on("lapack", when="+lapack")
    depends_on("hypre+mpi@2.22.1:", when="@5.7.1: +hypre")
    depends_on("hypre+mpi@:2.22.0", when="@:5.7.0 +hypre")
    depends_on("magma", when="+magma")
    depends_on("petsc+mpi", when="+petsc")
    depends_on("suite-sparse", when="+klu")
    depends_on("superlu-dist@7.0.0:", when="@6.4.0: +superlu-dist")
    depends_on("superlu-dist@6.3.0:", when="@5.5.0:6.3.0 +superlu-dist")
    depends_on("superlu-dist@6.1.1:", when="@:5.4.0 +superlu-dist")
    depends_on("superlu-mt+blas", when="+superlu-mt")
    depends_on("trilinos+tpetra@:13", when="@:7.1.1 +trilinos")
    depends_on("trilinos+tpetra@14:", when="@7.2.0: +trilinos")

    # Require that external libraries built with the same precision
    depends_on("petsc~double~complex", when="+petsc precision=single")
    depends_on("petsc+double~complex", when="+petsc precision=double")

    # Require that external libraries built with the same index type
    with when("+int64"):
        depends_on("hypre+mpi+int64", when="+hypre +int64")
        depends_on("petsc+int64", when="+petsc +int64")
        depends_on("superlu-dist+int64", when="+superlu-dist +int64")

    with when("~int64"):
        depends_on("hypre+mpi~int64", when="+hypre ~int64")
        depends_on("petsc~int64", when="+petsc ~int64")
        depends_on("superlu-dist~int64", when="+superlu-dist ~int64")

    # ==========================================================================
    # Conflicts
    # ==========================================================================

    conflicts("+cuda", when="@:2.7.0")
    conflicts("+f2003", when="@:4.1.0")
    conflicts("~int64", when="@:2.7.0")
    conflicts("+rocm", when="@:5.6.0")
    conflicts("~openmp", when="^superlu-dist+openmp")

    # External libraries incompatible with 64-bit indices
    conflicts("+lapack", when="@3.0.0: +int64")
    conflicts("+hypre", when="+hypre@:2.6.1a +int64")

    # External libraries incompatible with single precision
    with when("precision=single"):
        conflicts("+hypre", when="+hypre@:2.12.0")
        conflicts("+klu")
        conflicts("+superlu-dist")

    # External libraries incompatible with extended (quad) precision
    with when("precision=extended"):
        conflicts("+hypre", when="+hypre@:2.12.0")
        conflicts("+klu")
        conflicts("+lapack")
        conflicts("+superlu-dist")
        conflicts("+superlu-mt")

    # SuperLU_MT interface requires lapack for external blas (before v3.0.0)
    conflicts("+superlu-mt", when="@:2.7.0 ~lapack")

    # rocm+examples and cstd do not work together in 6.0.0
    conflicts("+rocm+examples", when="@6.0.0")

    # Ginkgo requires at least C++17 starting with 1.9.0
    conflicts("cxxstd=14", when="+ginkgo@1.9.0:")

    # SYCL requires at least C++17
    conflicts("cxxstd=14", when="+sycl")

    # ==========================================================================
    # Patches
    # ==========================================================================

    # remove OpenMP header file and function from hypre vector test code
    patch("test_nvector_parhyp.patch", when="@2.7.0:3.0.0")
    patch("FindPackageMultipass.cmake.patch", when="@5.0.0")
    patch("5.5.0-xsdk-patches.patch", when="@5.5.0")
    patch("0001-add-missing-README-to-examples-cvode-hip.patch", when="@5.6.0:5.7.0")
    # remove sundials_nvecopenmp target from ARKODE SuperLU_DIST example
    patch("remove-links-to-OpenMP-vector.patch", when="@5.5.0:5.7.0")
    # fix issues with exported PETSc target(s) in SUNDIALSConfig.cmake
    patch("sundials-v5.8.0.patch", when="@5.8.0")
    # https://github.com/spack/spack/issues/29526
    patch("nvector-pic.patch", when="@6.1.0:6.2.0 +rocm")

    # ==========================================================================
    # Post Install Actions
    # ==========================================================================

    @run_after("install")
    def post_install(self):
        """Run after install to fix install name of dynamic libraries
        on Darwin to have full path and install the LICENSE file."""
        spec = self.spec
        prefix = self.spec.prefix

        if sys.platform == "darwin":
            fix_darwin_install_name(prefix.lib)

        if spec.satisfies("@:3.0.0"):
            install("LICENSE", prefix)

    @run_after("install")
    def filter_compilers(self):
        """Run after install to tell the example program Makefiles
        to use the compilers that Spack built the package with.

        If this isn't done, they'll have CC, CPP, and F77 set to
        Spack's generic cc and f77. We want them to be bound to
        whatever compiler they were built with."""

        spec = self.spec

        kwargs = {"ignore_absent": True, "backup": False, "string": True}
        dirname = os.path.join(self.prefix, "examples")

        cc_files = [
            "arkode/C_openmp/Makefile",
            "arkode/C_parallel/Makefile",
            "arkode/C_parhyp/Makefile",
            "arkode/C_petsc/Makefile",
            "arkode/C_serial/Makefile",
            "cvode/C_openmp/Makefile",
            "cvode/parallel/Makefile",
            "cvode/parhyp/Makefile",
            "cvode/petsc/Makefile",
            "cvode/serial/Makefile",
            "cvodes/C_openmp/Makefile",
            "cvodes/parallel/Makefile",
            "cvodes/serial/Makefile",
            "ida/C_openmp/Makefile",
            "ida/parallel/Makefile",
            "ida/petsc/Makefile",
            "ida/serial/Makefile",
            "idas/C_openmp/Makefile",
            "idas/parallel/Makefile",
            "idas/serial/Makefile",
            "kinsol/C_openmp/Makefile",
            "kinsol/parallel/Makefile",
            "kinsol/serial/Makefile",
            "nvector/C_openmp/Makefile",
            "nvector/parallel/Makefile",
            "nvector/parhyp/Makefile",
            "nvector/petsc/Makefile",
            "nvector/pthreads/Makefile",
            "nvector/serial/Makefile",
            "sunlinsol/band/Makefile",
            "sunlinsol/dense/Makefile",
            "sunlinsol/klu/Makefile",
            "sunlinsol/lapackband/Makefile",
            "sunlinsol/lapackdense/Makefile",
            "sunlinsol/pcg/parallel/Makefile",
            "sunlinsol/pcg/serial/Makefile",
            "sunlinsol/spbcgs/parallel/Makefile",
            "sunlinsol/spbcgs/serial/Makefile",
            "sunlinsol/spfgmr/parallel/Makefile",
            "sunlinsol/spfgmr/serial/Makefile",
            "sunlinsol/spgmr/parallel/Makefile",
            "sunlinsol/spgmr/serial/Makefile",
            "sunlinsol/sptfqmr/parallel/Makefile",
            "sunlinsol/sptfqmr/serial/Makefile",
            "sunlinsol/superlumt/Makefile",
            "sunlinsol/superludist/Makefile",
            "sunmatrix/band/Makefile",
            "sunmatrix/dense/Makefile",
            "sunmatrix/sparse/Makefile",
        ]

        cxx_files = [
            "arkode/CXX_parallel/Makefile",
            "arkode/CXX_serial/Makefile" "cvode/cuda/Makefile",
            "cvode/raja/Makefile",
            "nvector/cuda/Makefile",
            "nvector/raja/Makefile",
        ]

        f77_files = [
            "arkode/F77_parallel/Makefile",
            "arkode/F77_serial/Makefile",
            "cvode/fcmix_parallel/Makefile",
            "cvode/fcmix_serial/Makefile",
            "ida/fcmix_openmp/Makefile",
            "ida/fcmix_parallel/Makefile",
            "ida/fcmix_pthreads/Makefile",
            "ida/fcmix_serial/Makefile",
            "kinsol/fcmix_parallel/Makefile",
            "kinsol/fcmix_serial/Makefile",
        ]

        f90_files = ["arkode/F90_parallel/Makefile", "arkode/F90_serial/Makefile"]

        f2003_files = [
            "arkode/F2003_serial/Makefile",
            "cvode/F2003_serial/Makefile",
            "cvodes/F2003_serial/Makefike",
            "ida/F2003_serial/Makefile",
            "idas/F2003_serial/Makefile",
            "kinsol/F2003_serial/Makefile",
        ]

        for filename in cc_files:
            filter_file(
                os.environ["CC"], self.compiler.cc, os.path.join(dirname, filename), **kwargs
            )

        for filename in cc_files:
            filter_file(r"^CPP\s*=.*", self.compiler.cc, os.path.join(dirname, filename), **kwargs)

        for filename in cxx_files:
            filter_file(
                os.environ["CXX"], self.compiler.cxx, os.path.join(dirname, filename), **kwargs
            )

        for filename in cxx_files:
            filter_file(r"^CPP\s*=.*", self.compiler.cc, os.path.join(dirname, filename), **kwargs)

        if ("+fcmix" in spec) and ("+examples" in spec):
            for filename in f77_files:
                filter_file(
                    os.environ["F77"], self.compiler.f77, os.path.join(dirname, filename), **kwargs
                )

        if ("+fcmix" in spec) and ("+examples" in spec):
            for filename in f90_files:
                filter_file(
                    os.environ["FC"], self.compiler.fc, os.path.join(dirname, filename), **kwargs
                )

        if ("+f2003" in spec) and ("+examples" in spec):
            for filename in f2003_files:
                filter_file(
                    os.environ["FC"], self.compiler.fc, os.path.join(dirname, filename), **kwargs
                )

    @property
    def headers(self):
        """Export the headers and defines of SUNDIALS.
        Sample usage: spec['sundials'].headers.cpp_flags
        """
        # SUNDIALS headers are inside subdirectories, so we use a fake header
        # in the include directory.
        hdr = find(self.prefix.include.nvector, "nvector_serial.h", recursive=False)
        return HeaderList(join_path(self.spec.prefix.include, "fake.h")) if hdr else None

    @property
    def libs(self):
        """Export the libraries of SUNDIALS.
        Sample usage: spec['sundials'].libs.ld_flags
                      spec['sundials:arkode,cvode'].libs.ld_flags
        """
        query_parameters = self.spec.last_query.extra_parameters
        if not query_parameters:
            sun_libs = "libsundials_*[!0-9]"
            # Q: should the result be ordered by dependency?
        else:
            sun_libs = ["libsundials_" + p for p in query_parameters]
        is_shared = "+shared" in self.spec

        libs = find_libraries(sun_libs, root=self.prefix, shared=is_shared, recursive=True)

        return libs or None  # Raise an error if no libs are found

    @run_after("install")
    @on_package_attributes(run_tests=True)
    def test_install(self):
        """Perform make test_install."""
        with working_dir(self.build_directory):
            make("test_install")

    @property
    def _smoke_tests(self):
        # smoke_tests tuple: exe, args, purpose, use cmake (true/false)
        smoke_tests = []
        if "+CVODE" in self.spec:
            smoke_tests.append(("cvode/serial/cvAdvDiff_bnd", [], "Test CVODE", True))

        if "+cuda" in self.spec:
            if "+CVODE" in self.spec:
                smoke_tests.append(
                    ("cvode/cuda/cvAdvDiff_kry_cuda", [], "Test CVODE with CUDA", True)
                )

        if "+hip" in self.spec:
            if "+CVODE" in self.spec:
                smoke_tests.append(
                    ("cvode/hip/cvAdvDiff_kry_hip", [], "Test CVODE with HIP", True)
                )

        if "+sycl" in self.spec:
            if "+CVODE" in self.spec:
                smoke_tests.append(
                    ("cvode/sycl/cvAdvDiff_kry_sycl", [], "Test CVODE with SYCL", True)
                )

        return smoke_tests

    @property
    def _smoke_tests_path(self):
        # examples/smoke-tests are cached for testing
        return self.prefix.examples

    # TODO: Replace this method and its 'get' use for cmake path with
    #   join_path(self.spec['cmake'].prefix.bin, 'cmake') once stand-alone
    #   tests can access build dependencies through self.spec['cmake'].
    def cmake_bin(self, set=True):
        """(Hack) Set/get cmake dependency path."""
        filepath = join_path(self.install_test_root, "cmake_bin_path.txt")
        if set:
            with open(filepath, "w") as out_file:
                cmake_bin = join_path(self.spec["cmake"].prefix.bin, "cmake")
                out_file.write("{0}\n".format(cmake_bin))
        elif os.path.isfile(filepath):
            with open(filepath, "r") as in_file:
                return in_file.read().strip()

    @run_after("install")
    def setup_smoke_tests(self):
        install_tree(self._smoke_tests_path, join_path(self.install_test_root, "testing"))
        self.cmake_bin(set=True)

    def build_smoke_tests(self):
        cmake_bin = self.cmake_bin(set=False)

        if not cmake_bin:
            tty.msg("Skipping sundials test: cmake_bin_path.txt not found")
            return

        for smoke_test in self._smoke_tests:
            work_dir = join_path(self._smoke_tests_path, os.path.dirname(smoke_test[0]))
            with working_dir(work_dir):
                if smoke_test[3]:  # use cmake
                    self.run_test(exe=cmake_bin, options=["."])
                self.run_test(exe="make")

    def run_smoke_tests(self):
        for smoke_test in self._smoke_tests:
            self.run_test(
                exe=join_path(self._smoke_tests_path, smoke_test[0]),
                options=smoke_test[1],
                status=[0],
                installed=True,
                skip_missing=True,
                purpose=smoke_test[2],
            )

    def clean_smoke_tests(self):
        for smoke_test in self._smoke_tests:
            work_dir = join_path(self._smoke_tests_path, os.path.dirname(smoke_test[0]))
            with working_dir(work_dir):
                self.run_test(exe="make", options=["clean"])

    def test(self):
        self.build_smoke_tests()
        self.run_smoke_tests()
        self.clean_smoke_tests()
        return
        """Run the smoke tests."""
        if "+examples" not in self.spec:
            print("Smoke tests were skipped: install with examples enabled")
        return

        self.run_test(
            "examples/nvector/serial/test_nvector_serial",
            options=["10", "0"],
            work_dir=self._extra_tests_path,
        )
        if "+cuda" in self.spec:
            self.run_test("examples/cvode/cuda/cvAdvDiff_ky_cuda", work_dir=self._extra_tests_path)
        if "+rocm" in self.spec:
            self.run_test("examples/cvode/hip/cvAdvDiff_kry_hip", work_dir=self._extra_tests_path)
        if "+sycl" in self.spec:
            self.run_test(
                "examples/cvode/CXX_sycl/cvAdvDiff_kry_sycl", work_dir=self._extra_tests_path
            )
        return

    # ==========================================================================
    # SUNDIALS Settings
    # ==========================================================================

    def cmake_args(self):
        spec = self.spec

        options = []

        return options

    def _from_variant_helper(self, cmake_var, variant):
        if variant is None:
            variant = cmake_var.lower()
        if variant not in self.variants:
            raise KeyError('"{0}" is not a variant of "{1}"'.format(variant, self.name))
        if variant not in self.spec.variants:
            return ""
        value = self.spec.variants[variant].value
        if isinstance(value, (tuple, list)):
            # Sort multi-valued variants for reproducibility
            value = sorted(value)
        return value

    def cache_string_from_variant(self, cmake_var, variant):
        value = self._from_variant_helper(cmake_var, variant)
        return cmake_cache_string(cmake_var, value)

    def cache_option_from_variant(self, cmake_var, variant):
        value = self._from_variant_helper(cmake_var, variant)
        return cmake_cache_option(cmake_var, value)

    def initconfig_compiler_entries(self):
        entries = []

        entries.extend(
            [
                # compilers
                cmake_cache_path("CMAKE_C_COMPILER", self.compiler.cc),
                cmake_cache_path("CMAKE_CXX_COMPILER", self.compiler.cxx),
                cmake_cache_path("CMAKE_Fortran_COMPILER", self.compiler.fc),
                # language standard
                self.cache_string_from_variant("CMAKE_C_STANDARD", "cstd"),
                self.cache_string_from_variant("CMAKE_CXX_STANDARD", "cxxstd"),
            ]
        )

        return entries

    def initconfig_mpi_entries(self):
        spec = self.spec
        entries = []

        if "+mpi" in spec:
            entries.extend(
                [
                    self.cache_option_from_variant("MPI_ENABLE", "mpi"),
                    cmake_cache_path("MPI_MPICC", spec["mpi"].mpicc),
                    cmake_cache_path("MPI_MPICXX", spec["mpi"].mpicxx),
                    cmake_cache_path("MPI_MPIF77", spec["mpi"].mpif77),
                    cmake_cache_path("MPI_MPIF90", spec["mpi"].mpifc),
                ]
            )
            if "scheduler=flux" in spec:
                entries.append(cmake_cache_string("SUNDIALS_SCHEDULER_COMMAND", "flux run"))
            if "scheduler=slurm" in spec:
                entries.append(cmake_cache_string("SUNDIALS_SCHEDULER_COMMAND", "srun"))
            if "scheduler=lsf" in spec:
                entries.append(cmake_cache_string("SUNDIALS_SCHEDULER_COMMAND", "jsrun"))

        return entries

    def initconfig_hardware_entries(self):
        spec = self.spec
        entries = []

        if "+cuda" in spec:
            entries.append(self.cache_option_from_variant("CUDA_ENABLE", "cuda"))
            if not spec.satisfies("cuda_arch=none"):
                cuda_arch = spec.variants["cuda_arch"].value
                entries.append(
                    cmake_cache_string("CMAKE_CUDA_ARCHITECTURES", "{0}".format(cuda_arch[0]))
                )

        if "+rocm" in spec:
            entries.extend(
                [
                    self.cache_option_from_variant("ENABLE_HIP", "rocm"),
                    cmake_cache_path("HIP_PATH", spec["hip"].prefix),
                    cmake_cache_path("HIP_DIR", spec["hip"].prefix.cmake),
                    cmake_cache_path("HIP_CLANG_INCLUDE_PATH", spec["llvm-amdgpu"].prefix.include),
                    cmake_cache_path("ROCM_PATH", spec["hsa-rocr-dev"].prefix),
                    cmake_cache_string(
                        "AMDGPU_TARGETS", ";".join(spec.variants["amdgpu_target"].value)
                    ),
                ]
            )
        return entries

    def initconfig_package_entries(self):
        spec = self.spec
        entries = []

        # SUNDIALS solvers
        for pkg in self.sun_solvers:
            entries.append(self.cache_option_from_variant("BUILD_" + pkg, pkg))

        entries.extend(
            [
                # Precision
                self.cache_string_from_variant("SUNDIALS_PRECISION", "precision"),
                # Fortran interface
                self.cache_option_from_variant("F77_INTERFACE_ENABLE", "fcmix"),
                self.cache_option_from_variant("F2003_INTERFACE_ENABLE", "f2003"),
                # library type
                self.cache_option_from_variant("BUILD_SHARED_LIBS", "shared"),
                self.cache_option_from_variant("BUILD_STATIC_LIBS", "static"),
                # Generic (std-c) math libraries
                self.cache_option_from_variant("USE_GENERIC_MATH", "generic-math"),
                # Logging
                self.cache_string_from_variant("SUNDIALS_LOGGING_LEVEL", "logging-level"),
                # Monitoring
                self.cache_option_from_variant("SUNDIALS_BUILD_WITH_MONITORING", "monitoring"),
                # Profiling
                self.cache_option_from_variant("SUNDIALS_BUILD_WITH_PROFILING", "profiling"),
                self.cache_option_from_variant("ENABLE_CALIPER", "caliper"),
                self.cache_option_from_variant("ENABLE_ADIAK", "adiak"),
                # Benchmarking
                self.cache_option_from_variant("BUILD_BENCHMARKS", "benchmarks"),
                # Profile examples
                self.cache_option_from_variant(
                    "SUNDIALS_TEST_ENABLE_PROFILING", "profile-examples"
                ),
                self.cache_option_from_variant(
                    "SUNDIALS_TEST_ENABLE_DEV_TESTS", "profile-examples"
                ),
                cmake_cache_string("SPACK_VERSION", ".".join(map(str, spack.spack_version_info))),
            ]
        )

        # index type (v3.0.0 or later)
        if spec.satisfies("@3:"):
            intsize = "64" if "+int64" in spec else "32"
            entries.extend(
                [
                    cmake_cache_string("SUNDIALS_INDEX_SIZE", intsize),
                    cmake_cache_string("SUNDIALS_INDEX_TYPE", "int{}_t".format(intsize)),
                ]
            )

        # TPLs
        entries.extend(
            [
                self.cache_option_from_variant("ENABLE_GINKGO", "ginkgo"),
                self.cache_option_from_variant("ENABLE_KOKKOS_KERNELS", "kokkos-kernels"),
                self.cache_option_from_variant("ENABLE_KOKKOS", "kokkos"),
                self.cache_option_from_variant("ENABLE_SYCL", "sycl"),
                self.cache_option_from_variant("EXAMPLES_INSTALL", "examples-install"),
                self.cache_option_from_variant("HYPRE_ENABLE", "hypre"),
                self.cache_option_from_variant("KLU_ENABLE", "klu"),
                self.cache_option_from_variant("LAPACK_ENABLE", "lapack"),
                self.cache_option_from_variant("OPENMP_ENABLE", "openmp"),
                self.cache_option_from_variant("PETSC_ENABLE", "petsc"),
                self.cache_option_from_variant("PTHREAD_ENABLE", "pthread"),
                self.cache_option_from_variant("RAJA_ENABLE", "raja"),
                self.cache_option_from_variant("SUPERLUDIST_ENABLE", "superlu-dist"),
                self.cache_option_from_variant("SUPERLUMT_ENABLE", "superlu-mt"),
                self.cache_option_from_variant("Trilinos_ENABLE", "trilinos"),
            ]
        )

        # Building with Adiak
        if "+adiak" in spec:
            entries.append(
                cmake_cache_path("adiak_DIR", spec["adiak"].prefix.lib.cmake + "/adiak")
            )

        # Building with Caliper
        if "+caliper" in spec:
            entries.append(cmake_cache_path("CALIPER_DIR", spec["caliper"].prefix))
            if "+adiak" in spec["caliper"]:
                entries.append(
                    cmake_cache_path("adiak_DIR", spec["adiak"].prefix.lib.cmake + "/adiak")
                )

            if not "caliper-dir=none" in spec:
                entries.append(
                    self.cache_string_from_variant(
                        "SUNDIALS_TEST_CALIPER_OUTPUT_DIR", "caliper-dir"
                    )
                )
                entries.append(
                    self.cache_string_from_variant(
                        "SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR", "caliper-dir"
                    )
                )

        # Building with Ginkgo
        if "+ginkgo" in spec:
            gko_backends = ["REF"]
            if "+openmp" in spec["ginkgo"] and "+openmp" in spec:
                gko_backends.append("OMP")
            if "+cuda" in spec["ginkgo"] and "+cuda" in spec:
                gko_backends.append("CUDA")
            if "+rocm" in spec["ginkgo"] and "+rocm" in spec:
                gko_backends.append("HIP")
            if "+oneapi" in spec["ginkgo"] and "+sycl" in spec:
                gko_backends.append("DPCPP")
            entries.extend(
                [
                    cmake_cache_path("Ginkgo_DIR", spec["ginkgo"].prefix),
                    cmake_cache_string("SUNDIALS_GINKGO_BACKENDS", ";".join(gko_backends)),
                ]
            )

        # Building with Hypre
        if "+hypre" in spec:
            entries.extend(
                [
                    cmake_cache_path("HYPRE_INCLUDE_DIR", spec["hypre"].prefix.include),
                    cmake_cache_path("HYPRE_LIBRARY_DIR", spec["hypre"].prefix.lib),
                ]
            )
            if not spec["hypre"].variants["shared"].value:
                hypre_libs = spec["blas"].libs + spec["lapack"].libs
                entries.extend([cmake_cache_string("HYPRE_LIBRARIES", hypre_libs.joined(";"))])

        # Building with Kokkos and KokkosKernels
        if "+kokkos" in spec:
            entries.extend([cmake_cache_path("Kokkos_DIR", spec["kokkos"].prefix)])
        if "+kokkos-kernels" in spec:
            entries.extend([cmake_cache_path("KokkosKernels_DIR", spec["kokkos-kernels"].prefix)])

        # Building with KLU
        if "+klu" in spec:
            entries.extend(
                [
                    cmake_cache_path("KLU_INCLUDE_DIR", spec["suite-sparse"].prefix.include),
                    cmake_cache_path("KLU_LIBRARY_DIR", spec["suite-sparse"].prefix.lib),
                ]
            )

        # Building with LAPACK
        if "+lapack" in spec:
            lapack_libs = []
            lapack_libs.extend(spec["lapack"].libs)
            lapack_libs.extend(spec["blas"].libs)
            entries.append(cmake_cache_string("LAPACK_LIBRARIES", ";".join(lapack_libs)))

        # Building with MAGMA
        if "+magma" in spec:
            entries.extend(
                [
                    cmake_cache_option("ENABLE_MAGMA", True),
                    cmake_cache_path("MAGMA_DIR", spec["magma"].prefix),
                ]
            )
            if "+cuda" in spec:
                entries.append(cmake_cache_string("SUNDIALS_MAGMA_BACKENDS", "CUDA"))
            if "+rocm" in spec:
                entries.append(cmake_cache_string("SUNDIALS_MAGMA_BACKENDS", "HIP"))

        # Building with PETSc
        if "+petsc" in spec:
            if spec.version >= Version("5"):
                entries.append(cmake_cache_path("PETSC_DIR", spec["petsc"].prefix))
                if "+kokkos" in spec["petsc"]:
                    entries.extend(
                        [
                            cmake_cache_path("Kokkos_DIR", spec["kokkos"].prefix),
                            cmake_cache_path("KokkosKernels_DIR", spec["kokkos-kernels"].prefix),
                        ]
                    )
            else:
                entries.extend(
                    [
                        cmake_cache_path("PETSC_INCLUDE_DIR", spec["petsc"].prefix.include),
                        cmake_cache_path("PETSC_LIBRARY_DIR", spec["petsc"].prefix.lib),
                    ]
                )

        # Building with RAJA
        if "+raja" in spec:
            entries.append(cmake_cache_path("RAJA_DIR", spec["raja"].prefix))
            if "camp" in spec:
                entries.append(
                    cmake_cache_path("camp_DIR", spec["camp"].prefix.lib.cmake + "/camp")
                )
            if "+rocm" in spec:
                entries.append(cmake_cache_string("SUNDIALS_RAJA_BACKENDS", "HIP"))

        # Building with SuperLU_DIST
        if "+superlu-dist" in spec:
            # if spec.satisfies("@6.4.0:"):
            if False:
                entries.extend(
                    [
                        cmake_cache_path("SUPERLUDIST_DIR", spec["superlu-dist"].prefix),
                        cmake_cache_string("SUPERLUDIST_OpenMP", "^superlu-dist+openmp" in spec),
                    ]
                )
            else:
                superludist_libs = []
                superludist_libs.extend(spec["parmetis"].libs)
                superludist_libs.extend(spec["metis"].libs)
                superludist_libs.extend(spec["superlu-dist"].libs)
                entries.extend(
                    [
                        cmake_cache_path(
                            "SUPERLUDIST_INCLUDE_DIR", spec["superlu-dist"].prefix.include
                        ),
                        cmake_cache_path(
                            "SUPERLUDIST_LIBRARY_DIR", spec["superlu-dist"].prefix.lib
                        ),
                        cmake_cache_string("SUPERLUDIST_LIBRARIES", ";".join(superludist_libs)),
                        cmake_cache_string("SUPERLUDIST_OpenMP", "^superlu-dist+openmp" in spec),
                    ]
                )

        # Building with SuperLU_MT
        if "+superlu-mt" in spec:
            if spec.satisfies("@3:"):
                entries.extend(
                    [
                        cmake_cache_string("BLAS_ENABLE", True),
                        cmake_cache_string("BLAS_LIBRARIES", spec["blas"].libs),
                    ]
                )
            entries.extend(
                [
                    cmake_cache_path("SUPERLUMT_INCLUDE_DIR", spec["superlu-mt"].prefix.include),
                    cmake_cache_path("SUPERLUMT_LIBRARY_DIR", spec["superlu-mt"].prefix.lib),
                    cmake_cache_string(
                        "SUPERLUMT_THREAD_TYPE",
                        "OpenMP" if "^superlu-mt+openmp" in spec else "Pthread",
                    ),
                ]
            )

        # Building with Trilinos
        if "+trilinos" in spec:
            entries.append(cmake_cache_path("Trilinos_DIR", spec["trilinos"].prefix))

        # Examples
        if spec.satisfies("@3:"):
            entries.extend(
                [
                    self.cache_option_from_variant("EXAMPLES_ENABLE_C", "examples"),
                    self.cache_option_from_variant("EXAMPLES_ENABLE_CXX", "examples"),
                    cmake_cache_option("EXAMPLES_ENABLE_CUDA", "+examples+cuda" in spec),
                    cmake_cache_option("EXAMPLES_ENABLE_F77", "+examples+fcmix" in spec),
                    cmake_cache_option("EXAMPLES_ENABLE_F90", "+examples+fcmix" in spec),
                    cmake_cache_option("EXAMPLES_ENABLE_F2003", "+examples+f2003" in spec),
                ]
            )
        else:
            entries.extend(
                [
                    self.cache_option_from_variant("EXAMPLES_ENABLE", "examples"),
                    self.cache_option_from_variant("CXX_ENABLE", "examples"),
                    cmake_cache_option("F90_ENABLE", "+examples+fcmix" in spec),
                ]
            )
        return entries
