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

m.def("KINSol", KINSol, nb::arg("kinmem"), nb::arg("uu"), nb::arg("strategy"),
      nb::arg("u_scale"), nb::arg("f_scale"), "Solver function");

m.def("KINSetUserData", KINSetUserData, nb::arg("kinmem"), nb::arg("user_data"));

m.def("KINSetDamping", KINSetDamping, nb::arg("kinmem"), nb::arg("beta"));

m.def("KINSetMAA", KINSetMAA, nb::arg("kinmem"), nb::arg("maa"));

m.def("KINSetOrthAA", KINSetOrthAA, nb::arg("kinmem"), nb::arg("orthaa"));

m.def("KINSetDelayAA", KINSetDelayAA, nb::arg("kinmem"), nb::arg("delay"));

m.def("KINSetDampingAA", KINSetDampingAA, nb::arg("kinmem"), nb::arg("beta"));

m.def("KINSetReturnNewest", KINSetReturnNewest, nb::arg("kinmem"),
      nb::arg("ret_newest"));

m.def("KINSetNumMaxIters", KINSetNumMaxIters, nb::arg("kinmem"),
      nb::arg("mxiter"));

m.def("KINSetNoInitSetup", KINSetNoInitSetup, nb::arg("kinmem"),
      nb::arg("noInitSetup"));

m.def("KINSetNoResMon", KINSetNoResMon, nb::arg("kinmem"),
      nb::arg("noNNIResMon"));

m.def("KINSetMaxSetupCalls", KINSetMaxSetupCalls, nb::arg("kinmem"),
      nb::arg("msbset"));

m.def("KINSetMaxSubSetupCalls", KINSetMaxSubSetupCalls, nb::arg("kinmem"),
      nb::arg("msbsetsub"));

m.def("KINSetEtaForm", KINSetEtaForm, nb::arg("kinmem"), nb::arg("etachoice"));

m.def("KINSetEtaConstValue", KINSetEtaConstValue, nb::arg("kinmem"),
      nb::arg("eta"));

m.def("KINSetEtaParams", KINSetEtaParams, nb::arg("kinmem"), nb::arg("egamma"),
      nb::arg("ealpha"));

m.def("KINSetResMonParams", KINSetResMonParams, nb::arg("kinmem"),
      nb::arg("omegamin"), nb::arg("omegamax"));

m.def("KINSetResMonConstValue", KINSetResMonConstValue, nb::arg("kinmem"),
      nb::arg("omegaconst"));

m.def("KINSetNoMinEps", KINSetNoMinEps, nb::arg("kinmem"), nb::arg("noMinEps"));

m.def("KINSetMaxNewtonStep", KINSetMaxNewtonStep, nb::arg("kinmem"),
      nb::arg("mxnewtstep"));

m.def("KINSetMaxBetaFails", KINSetMaxBetaFails, nb::arg("kinmem"),
      nb::arg("mxnbcf"));

m.def("KINSetRelErrFunc", KINSetRelErrFunc, nb::arg("kinmem"),
      nb::arg("relfunc"));

m.def("KINSetFuncNormTol", KINSetFuncNormTol, nb::arg("kinmem"),
      nb::arg("fnormtol"));

m.def("KINSetScaledStepTol", KINSetScaledStepTol, nb::arg("kinmem"),
      nb::arg("scsteptol"));

m.def("KINSetConstraints", KINSetConstraints, nb::arg("kinmem"),
      nb::arg("constraints"));

m.def(
  "KINGetNumNonlinSolvIters",
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

m.def(
  "KINGetNumFuncEvals",
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

m.def(
  "KINGetNumBetaCondFails",
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

m.def(
  "KINGetNumBacktrackOps",
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

m.def(
  "KINGetFuncNorm",
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

m.def(
  "KINGetStepLength",
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

m.def("KINPrintAllStats", KINPrintAllStats, nb::arg("kinmem"),
      nb::arg("outfile"), nb::arg("fmt"));

m.def("KINGetReturnFlagName", KINGetReturnFlagName, nb::arg("flag"));
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

m.def(
  "KINSetLinearSolver",
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

m.def(
  "KINGetJac",
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

m.def(
  "KINGetJacNumIters",
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

m.def(
  "KINGetNumJacEvals",
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

m.def(
  "KINGetNumLinFuncEvals",
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

m.def(
  "KINGetNumPrecEvals",
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

m.def(
  "KINGetNumPrecSolves",
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

m.def(
  "KINGetNumLinIters",
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

m.def(
  "KINGetNumLinConvFails",
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

m.def(
  "KINGetNumJtimesEvals",
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

m.def(
  "KINGetLastLinFlag",
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

m.def("KINGetLinReturnFlagName", KINGetLinReturnFlagName, nb::arg("flag"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
