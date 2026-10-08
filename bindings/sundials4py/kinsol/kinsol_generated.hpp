// #ifndef _KINSOL_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("KIN_SUCCESS")             = 0;
m.attr("KIN_INITIAL_GUESS_OK")    = 1;
m.attr("KIN_STEP_LT_STPTOL")      = 2;
m.attr("KIN_WARNING")             = 99;
m.attr("KIN_MEM_NULL")            = -1;
m.attr("KIN_ILL_INPUT")           = -2;
m.attr("KIN_NO_MALLOC")           = -3;
m.attr("KIN_MEM_FAIL")            = -4;
m.attr("KIN_LINESEARCH_NONCONV")  = -5;
m.attr("KIN_MAXITER_REACHED")     = -6;
m.attr("KIN_MXNEWT_5X_EXCEEDED")  = -7;
m.attr("KIN_LINESEARCH_BCFAIL")   = -8;
m.attr("KIN_LINSOLV_NO_RECOVERY") = -9;
m.attr("KIN_LINIT_FAIL")          = -10;
m.attr("KIN_LSETUP_FAIL")         = -11;
m.attr("KIN_LSOLVE_FAIL")         = -12;
m.attr("KIN_SYSFUNC_FAIL")        = -13;
m.attr("KIN_FIRST_SYSFUNC_ERR")   = -14;
m.attr("KIN_REPTD_SYSFUNC_ERR")   = -15;
m.attr("KIN_VECTOROP_ERR")        = -16;
m.attr("KIN_CONTEXT_ERR")         = -17;
m.attr("KIN_DAMPING_FN_ERR")      = -18;
m.attr("KIN_DEPTH_FN_ERR")        = -19;
m.attr("KIN_ORTH_MGS")            = 0;
m.attr("KIN_ORTH_ICWY")           = 1;
m.attr("KIN_ORTH_CGS2")           = 2;
m.attr("KIN_ORTH_DCGS2")          = 3;
m.attr("KIN_ETACHOICE1")          = 1;
m.attr("KIN_ETACHOICE2")          = 2;
m.attr("KIN_ETACONSTANT")         = 3;
m.attr("KIN_NONE")                = 0;
m.attr("KIN_LINESEARCH")          = 1;
m.attr("KIN_PICARD")              = 2;
m.attr("KIN_FP")                  = 3;

sundials4py::scoped_def(m, "KINSetUserData", KINSetUserData, nb::arg("kinmem"),
                        nb::arg("user_data"));

sundials4py::scoped_def(m, "KINSetDamping", KINSetDamping, nb::arg("kinmem"),
                        nb::arg("beta"));

sundials4py::scoped_def(m, "KINSetMAA", KINSetMAA, nb::arg("kinmem"),
                        nb::arg("maa"));

sundials4py::scoped_def(m, "KINSetOrthAA", KINSetOrthAA, nb::arg("kinmem"),
                        nb::arg("orthaa"));

sundials4py::scoped_def(m, "KINSetDelayAA", KINSetDelayAA, nb::arg("kinmem"),
                        nb::arg("delay"));

sundials4py::scoped_def(m, "KINSetDampingAA", KINSetDampingAA,
                        nb::arg("kinmem"), nb::arg("beta"));

sundials4py::scoped_def(m, "KINSetReturnNewest", KINSetReturnNewest,
                        nb::arg("kinmem"), nb::arg("ret_newest"));

sundials4py::scoped_def(m, "KINSetNumMaxIters", KINSetNumMaxIters,
                        nb::arg("kinmem"), nb::arg("mxiter"));

sundials4py::scoped_def(m, "KINSetNoInitSetup", KINSetNoInitSetup,
                        nb::arg("kinmem"), nb::arg("noInitSetup"));

sundials4py::scoped_def(m, "KINSetNoResMon", KINSetNoResMon, nb::arg("kinmem"),
                        nb::arg("noNNIResMon"));

sundials4py::scoped_def(m, "KINSetMaxSetupCalls", KINSetMaxSetupCalls,
                        nb::arg("kinmem"), nb::arg("msbset"));

sundials4py::scoped_def(m, "KINSetMaxSubSetupCalls", KINSetMaxSubSetupCalls,
                        nb::arg("kinmem"), nb::arg("msbsetsub"));

sundials4py::scoped_def(m, "KINSetEtaForm", KINSetEtaForm, nb::arg("kinmem"),
                        nb::arg("etachoice"));

sundials4py::scoped_def(m, "KINSetEtaConstValue", KINSetEtaConstValue,
                        nb::arg("kinmem"), nb::arg("eta"));

sundials4py::scoped_def(m, "KINSetEtaParams", KINSetEtaParams,
                        nb::arg("kinmem"), nb::arg("egamma"), nb::arg("ealpha"));

sundials4py::scoped_def(m, "KINSetResMonParams", KINSetResMonParams,
                        nb::arg("kinmem"), nb::arg("omegamin"),
                        nb::arg("omegamax"));

sundials4py::scoped_def(m, "KINSetResMonConstValue", KINSetResMonConstValue,
                        nb::arg("kinmem"), nb::arg("omegaconst"));

sundials4py::scoped_def(m, "KINSetNoMinEps", KINSetNoMinEps, nb::arg("kinmem"),
                        nb::arg("noMinEps"));

sundials4py::scoped_def(m, "KINSetMaxNewtonStep", KINSetMaxNewtonStep,
                        nb::arg("kinmem"), nb::arg("mxnewtstep"));

sundials4py::scoped_def(m, "KINSetMaxBetaFails", KINSetMaxBetaFails,
                        nb::arg("kinmem"), nb::arg("mxnbcf"));

sundials4py::scoped_def(m, "KINSetRelErrFunc", KINSetRelErrFunc,
                        nb::arg("kinmem"), nb::arg("relfunc"));

sundials4py::scoped_def(m, "KINSetFuncNormTol", KINSetFuncNormTol,
                        nb::arg("kinmem"), nb::arg("fnormtol"));

sundials4py::scoped_def(m, "KINSetScaledStepTol", KINSetScaledStepTol,
                        nb::arg("kinmem"), nb::arg("scsteptol"));

sundials4py::scoped_def(m, "KINSetConstraints", KINSetConstraints,
                        nb::arg("kinmem"), nb::arg("constraints"));

sundials4py::scoped_def(
  m, "KINGetNumNonlinSolvIters",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nniters_adapt_modifiable;

      int r = KINGetNumNonlinSolvIters(kinmem, &nniters_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable);
    };

    return KINGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumFuncEvals",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumFuncEvals_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nfevals_adapt_modifiable;

      int r = KINGetNumFuncEvals(kinmem, &nfevals_adapt_modifiable);
      return std::make_tuple(r, nfevals_adapt_modifiable);
    };

    return KINGetNumFuncEvals_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumBetaCondFails",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumBetaCondFails_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nbcfails_adapt_modifiable;

      int r = KINGetNumBetaCondFails(kinmem, &nbcfails_adapt_modifiable);
      return std::make_tuple(r, nbcfails_adapt_modifiable);
    };

    return KINGetNumBetaCondFails_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumBacktrackOps",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumBacktrackOps_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nbacktr_adapt_modifiable;

      int r = KINGetNumBacktrackOps(kinmem, &nbacktr_adapt_modifiable);
      return std::make_tuple(r, nbacktr_adapt_modifiable);
    };

    return KINGetNumBacktrackOps_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetFuncNorm",
  [](void* kinmem) -> std::tuple<int, sunrealtype>
  {
    auto KINGetFuncNorm_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype fnorm_adapt_modifiable;

      int r = KINGetFuncNorm(kinmem, &fnorm_adapt_modifiable);
      return std::make_tuple(r, fnorm_adapt_modifiable);
    };

    return KINGetFuncNorm_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetStepLength",
  [](void* kinmem) -> std::tuple<int, sunrealtype>
  {
    auto KINGetStepLength_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype steplength_adapt_modifiable;

      int r = KINGetStepLength(kinmem, &steplength_adapt_modifiable);
      return std::make_tuple(r, steplength_adapt_modifiable);
    };

    return KINGetStepLength_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(m, "KINPrintAllStats", KINPrintAllStats,
                        nb::arg("kinmem"), nb::arg("outfile"), nb::arg("fmt"));

sundials4py::scoped_def(m, "KINGetReturnFlagName", KINGetReturnFlagName,
                        nb::arg("flag"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _KINLS_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("KINLS_SUCCESS")     = 0;
m.attr("KINLS_MEM_NULL")    = -1;
m.attr("KINLS_LMEM_NULL")   = -2;
m.attr("KINLS_ILL_INPUT")   = -3;
m.attr("KINLS_MEM_FAIL")    = -4;
m.attr("KINLS_PMEM_NULL")   = -5;
m.attr("KINLS_JACFUNC_ERR") = -6;
m.attr("KINLS_SUNMAT_FAIL") = -7;
m.attr("KINLS_SUNLS_FAIL")  = -8;

sundials4py::scoped_def(
  m, "KINSetLinearSolver",
  [](void* kinmem, SUNLinearSolver LS,
     std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto KINSetLinearSolver_adapt_optional_arg_with_default_null =
      [](void* kinmem, SUNLinearSolver LS,
         std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = KINSetLinearSolver(kinmem, LS, A_adapt_default_null);
      return lambda_result;
    };

    return KINSetLinearSolver_adapt_optional_arg_with_default_null(kinmem, LS, A);
  },
  nb::arg("kinmem"), nb::arg("LS"), nb::arg("A").none() = nb::none());

sundials4py::scoped_def(
  m, "KINGetJac",
  [](void* kinmem) -> std::tuple<int, SUNMatrix>
  {
    auto KINGetJac_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, SUNMatrix>
    {
      SUNMatrix J_adapt_modifiable;

      int r = KINGetJac(kinmem, &J_adapt_modifiable);
      return std::make_tuple(r, J_adapt_modifiable);
    };

    return KINGetJac_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"), "nb::rv_policy::reference", nb::rv_policy::reference);

sundials4py::scoped_def(
  m, "KINGetJacNumIters",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetJacNumIters_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nni_J_adapt_modifiable;

      int r = KINGetJacNumIters(kinmem, &nni_J_adapt_modifiable);
      return std::make_tuple(r, nni_J_adapt_modifiable);
    };

    return KINGetJacNumIters_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumJacEvals",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumJacEvals_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long njevals_adapt_modifiable;

      int r = KINGetNumJacEvals(kinmem, &njevals_adapt_modifiable);
      return std::make_tuple(r, njevals_adapt_modifiable);
    };

    return KINGetNumJacEvals_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumLinFuncEvals",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumLinFuncEvals_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nfevals_adapt_modifiable;

      int r = KINGetNumLinFuncEvals(kinmem, &nfevals_adapt_modifiable);
      return std::make_tuple(r, nfevals_adapt_modifiable);
    };

    return KINGetNumLinFuncEvals_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumPrecEvals",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumPrecEvals_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long npevals_adapt_modifiable;

      int r = KINGetNumPrecEvals(kinmem, &npevals_adapt_modifiable);
      return std::make_tuple(r, npevals_adapt_modifiable);
    };

    return KINGetNumPrecEvals_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumPrecSolves",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumPrecSolves_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long npsolves_adapt_modifiable;

      int r = KINGetNumPrecSolves(kinmem, &npsolves_adapt_modifiable);
      return std::make_tuple(r, npsolves_adapt_modifiable);
    };

    return KINGetNumPrecSolves_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumLinIters",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumLinIters_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nliters_adapt_modifiable;

      int r = KINGetNumLinIters(kinmem, &nliters_adapt_modifiable);
      return std::make_tuple(r, nliters_adapt_modifiable);
    };

    return KINGetNumLinIters_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumLinConvFails",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumLinConvFails_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long nlcfails_adapt_modifiable;

      int r = KINGetNumLinConvFails(kinmem, &nlcfails_adapt_modifiable);
      return std::make_tuple(r, nlcfails_adapt_modifiable);
    };

    return KINGetNumLinConvFails_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetNumJtimesEvals",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetNumJtimesEvals_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long njvevals_adapt_modifiable;

      int r = KINGetNumJtimesEvals(kinmem, &njvevals_adapt_modifiable);
      return std::make_tuple(r, njvevals_adapt_modifiable);
    };

    return KINGetNumJtimesEvals_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(
  m, "KINGetLastLinFlag",
  [](void* kinmem) -> std::tuple<int, long>
  {
    auto KINGetLastLinFlag_adapt_modifiable_immutable_to_return =
      [](void* kinmem) -> std::tuple<int, long>
    {
      long flag_adapt_modifiable;

      int r = KINGetLastLinFlag(kinmem, &flag_adapt_modifiable);
      return std::make_tuple(r, flag_adapt_modifiable);
    };

    return KINGetLastLinFlag_adapt_modifiable_immutable_to_return(kinmem);
  },
  nb::arg("kinmem"));

sundials4py::scoped_def(m, "KINGetLinReturnFlagName", KINGetLinReturnFlagName,
                        nb::arg("flag"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
