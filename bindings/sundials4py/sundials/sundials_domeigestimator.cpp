/* -----------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2025-2025, Lawrence Livermore National Security,
 * University of Maryland Baltimore County, and the SUNDIALS contributors.
 * Copyright (c) 2013-2025, Lawrence Livermore National Security
 * and Southern Methodist University.
 * Copyright (c) 2002-2013, Lawrence Livermore National Security.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * This file is the entrypoint for the Python binding code for the
 * SUNDIALS N_Vector class. It contains hand-written code for functions
 * that require special treatment, and includes the generated code
 * produced with the generate.py script.
 * -----------------------------------------------------------------*/

#include <cstdlib>
#include <cstring>

#include "sundials4py.hpp"

#include <sundials/sundials_classview.hpp>
#include <sundials/sundials_domeigestimator.hpp>

#include "sundials_domeigestimator_usersupplied.hpp"

namespace nb = nanobind;
using namespace sundials::experimental;

namespace sundials4py {

void bind_sundomeigestimator(nb::module_& m)
{
#include "sundials_domeigestimator_generated.hpp"

  m.def(
    "SUNDomEigEstimator_SetOptions",
    [](SUNDomEigEstimator self, const std::string& id,
       const std::string& file_name, int argc,
       const std::vector<std::string>& args)
    {
      std::vector<char*> argv;
      argv.reserve(args.size());

      for (const auto& arg : args)
      {
        // We need a non-const char*, so we use data() and an explicit cast.
        // This is safe as long as the underlying std::string is not modified.
        argv.push_back(const_cast<char*>(arg.data()));
      }

      return SUNDomEigEstimator_SetOptions(self,
                                           id.empty() ? nullptr : id.c_str(),
                                           file_name.empty() ? nullptr
                                                             : file_name.c_str(),
                                           argc, argv.data());
    },
    nb::arg("self"), nb::arg("id"), nb::arg("file_name"), nb::arg("argc"),
    nb::arg("args"));

  m.def(
    "SUNDomEigEstimator_SetATimes",
    [](SUNDomEigEstimator dee,
       std::function<std::remove_pointer_t<SUNATimesFn>> ATimes) -> SUNErrCode
    {
      if (!dee->python) { dee->python = new SUNDomEigEstimatorFunctionTable; }

      auto fntable = static_cast<SUNDomEigEstimatorFunctionTable*>(dee->python);

      fntable->atimes = nb::cast(ATimes);

      if (ATimes)
      {
        return SUNDomEigEstimator_SetATimes(dee, fntable,
                                            sundomeigestimator_atimes_wrapper);
      }
      else { return SUNDomEigEstimator_SetATimes(dee, fntable, nullptr); }
    },
    nb::arg("DEE"), nb::arg("ATimes").none());
}

} // namespace sundials4py

extern "C" void SUNDomEigEstimatorFunctionTable_Destroy(void* ptr)
{
  delete static_cast<SUNDomEigEstimatorFunctionTable*>(ptr);
}
