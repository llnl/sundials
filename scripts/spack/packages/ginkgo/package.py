# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *
import spack.pkg.builtin.ginkgo


class Ginkgo(spack.pkg.builtin.ginkgo.Ginkgo):

    # The version of Spack we are using does not include Ginkgo 1.9.0
    version("1.9.0", commit="20cfd68795f58078898da9890baa311b46845a8b")
