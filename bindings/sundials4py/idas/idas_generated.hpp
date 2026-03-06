// #ifndef _IDAS_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("IDA_NORMAL")             = 1;
m.attr("IDA_ONE_STEP")           = 2;
m.attr("IDA_YA_YDP_INIT")        = 1;
m.attr("IDA_Y_INIT")             = 2;
m.attr("IDA_SIMULTANEOUS")       = 1;
m.attr("IDA_STAGGERED")          = 2;
m.attr("IDA_CENTERED")           = 1;
m.attr("IDA_FORWARD")            = 2;
m.attr("IDA_HERMITE")            = 1;
m.attr("IDA_POLYNOMIAL")         = 2;
m.attr("IDA_SUCCESS")            = 0;
m.attr("IDA_TSTOP_RETURN")       = 1;
m.attr("IDA_ROOT_RETURN")        = 2;
m.attr("IDA_WARNING")            = 99;
m.attr("IDA_TOO_MUCH_WORK")      = -1;
m.attr("IDA_TOO_MUCH_ACC")       = -2;
m.attr("IDA_ERR_FAIL")           = -3;
m.attr("IDA_CONV_FAIL")          = -4;
m.attr("IDA_LINIT_FAIL")         = -5;
m.attr("IDA_LSETUP_FAIL")        = -6;
m.attr("IDA_LSOLVE_FAIL")        = -7;
m.attr("IDA_RES_FAIL")           = -8;
m.attr("IDA_REP_RES_ERR")        = -9;
m.attr("IDA_RTFUNC_FAIL")        = -10;
m.attr("IDA_CONSTR_FAIL")        = -11;
m.attr("IDA_FIRST_RES_FAIL")     = -12;
m.attr("IDA_LINESEARCH_FAIL")    = -13;
m.attr("IDA_NO_RECOVERY")        = -14;
m.attr("IDA_NLS_INIT_FAIL")      = -15;
m.attr("IDA_NLS_SETUP_FAIL")     = -16;
m.attr("IDA_NLS_FAIL")           = -17;
m.attr("IDA_MEM_NULL")           = -20;
m.attr("IDA_MEM_FAIL")           = -21;
m.attr("IDA_ILL_INPUT")          = -22;
m.attr("IDA_NO_MALLOC")          = -23;
m.attr("IDA_BAD_EWT")            = -24;
m.attr("IDA_BAD_K")              = -25;
m.attr("IDA_BAD_T")              = -26;
m.attr("IDA_BAD_DKY")            = -27;
m.attr("IDA_VECTOROP_ERR")       = -28;
m.attr("IDA_CONTEXT_ERR")        = -29;
m.attr("IDA_NO_QUAD")            = -30;
m.attr("IDA_QRHS_FAIL")          = -31;
m.attr("IDA_FIRST_QRHS_ERR")     = -32;
m.attr("IDA_REP_QRHS_ERR")       = -33;
m.attr("IDA_NO_SENS")            = -40;
m.attr("IDA_SRES_FAIL")          = -41;
m.attr("IDA_REP_SRES_ERR")       = -42;
m.attr("IDA_BAD_IS")             = -43;
m.attr("IDA_NO_QUADSENS")        = -50;
m.attr("IDA_QSRHS_FAIL")         = -51;
m.attr("IDA_FIRST_QSRHS_ERR")    = -52;
m.attr("IDA_REP_QSRHS_ERR")      = -53;
m.attr("IDA_TOO_CLOSE")          = -60;
m.attr("IDA_UNRECOGNIZED_ERROR") = -99;
m.attr("IDA_NO_ADJ")             = -101;
m.attr("IDA_NO_FWD")             = -102;
m.attr("IDA_NO_BCK")             = -103;
m.attr("IDA_BAD_TB0")            = -104;
m.attr("IDA_REIFWD_FAIL")        = -105;
m.attr("IDA_FWD_FAIL")           = -106;
m.attr("IDA_GETY_BADT")          = -107;

m.def("IDAReInit", IDAReInit, nb::arg("ida_mem"), nb::arg("t0"), nb::arg("yy0"),
      nb::arg("yp0"));

m.def("IDASStolerances", IDASStolerances, nb::arg("ida_mem"), nb::arg("reltol"),
      nb::arg("abstol"));

m.def("IDASVtolerances", IDASVtolerances, nb::arg("ida_mem"), nb::arg("reltol"),
      nb::arg("abstol"));

m.def("IDACalcIC", IDACalcIC, nb::arg("ida_mem"), nb::arg("icopt"),
      nb::arg("tout1"), "Initial condition calculation function");

m.def("IDASetNonlinConvCoefIC", IDASetNonlinConvCoefIC, nb::arg("ida_mem"),
      nb::arg("epiccon"));

m.def("IDASetMaxNumStepsIC", IDASetMaxNumStepsIC, nb::arg("ida_mem"),
      nb::arg("maxnh"));

m.def("IDASetMaxNumJacsIC", IDASetMaxNumJacsIC, nb::arg("ida_mem"),
      nb::arg("maxnj"));

m.def("IDASetMaxNumItersIC", IDASetMaxNumItersIC, nb::arg("ida_mem"),
      nb::arg("maxnit"));

m.def("IDASetLineSearchOffIC", IDASetLineSearchOffIC, nb::arg("ida_mem"),
      nb::arg("lsoff"));

m.def("IDASetStepToleranceIC", IDASetStepToleranceIC, nb::arg("ida_mem"),
      nb::arg("steptol"));

m.def("IDASetMaxBacksIC", IDASetMaxBacksIC, nb::arg("ida_mem"),
      nb::arg("maxbacks"));

m.def("IDASetDeltaCjLSetup", IDASetDeltaCjLSetup, nb::arg("ida_max"),
      nb::arg("dcj"));

m.def("IDASetMaxOrd", IDASetMaxOrd, nb::arg("ida_mem"), nb::arg("maxord"));

m.def("IDASetMaxNumSteps", IDASetMaxNumSteps, nb::arg("ida_mem"),
      nb::arg("mxsteps"));

m.def("IDASetInitStep", IDASetInitStep, nb::arg("ida_mem"), nb::arg("hin"));

m.def("IDASetMaxStep", IDASetMaxStep, nb::arg("ida_mem"), nb::arg("hmax"));

m.def("IDASetMinStep", IDASetMinStep, nb::arg("ida_mem"), nb::arg("hmin"));

m.def("IDASetStopTime", IDASetStopTime, nb::arg("ida_mem"), nb::arg("tstop"));

m.def("IDAClearStopTime", IDAClearStopTime, nb::arg("ida_mem"));

m.def("IDASetMaxErrTestFails", IDASetMaxErrTestFails, nb::arg("ida_mem"),
      nb::arg("maxnef"));

m.def("IDASetSuppressAlg", IDASetSuppressAlg, nb::arg("ida_mem"),
      nb::arg("suppressalg"));

m.def("IDASetId", IDASetId, nb::arg("ida_mem"), nb::arg("id"));

m.def("IDASetConstraints", IDASetConstraints, nb::arg("ida_mem"),
      nb::arg("constraints"));

m.def("IDASetMaxNumConstraintFails", IDASetMaxNumConstraintFails,
      nb::arg("ida_mem"), nb::arg("max_fails"));

m.def("IDASetEtaFixedStepBounds", IDASetEtaFixedStepBounds, nb::arg("ida_mem"),
      nb::arg("eta_min_fx"), nb::arg("eta_max_fx"));

m.def("IDASetEtaMin", IDASetEtaMin, nb::arg("ida_mem"), nb::arg("eta_min"));

m.def("IDASetEtaMax", IDASetEtaMax, nb::arg("ida_mem"), nb::arg("eta_max"));

m.def("IDASetEtaLow", IDASetEtaLow, nb::arg("ida_mem"), nb::arg("eta_low"));

m.def("IDASetEtaMinErrFail", IDASetEtaMinErrFail, nb::arg("ida_mem"),
      nb::arg("eta_min_ef"));

m.def("IDASetEtaConvFail", IDASetEtaConvFail, nb::arg("ida_mem"),
      nb::arg("eta_cf"));

m.def("IDASetMaxConvFails", IDASetMaxConvFails, nb::arg("ida_mem"),
      nb::arg("maxncf"));

m.def("IDASetMaxNonlinIters", IDASetMaxNonlinIters, nb::arg("ida_mem"),
      nb::arg("maxcor"));

m.def("IDASetNonlinConvCoef", IDASetNonlinConvCoef, nb::arg("ida_mem"),
      nb::arg("epcon"));

m.def("IDASetNonlinearSolver", IDASetNonlinearSolver, nb::arg("ida_mem"),
      nb::arg("NLS"));

m.def(
  "IDASetRootDirection",
  [](void* ida_mem, std::vector<int> rootdir_1d) -> int
  {
    auto IDASetRootDirection_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<int> rootdir_1d) -> int
    {
      int* rootdir_1d_ptr = rootdir_1d.empty() ? nullptr : rootdir_1d.data();

      auto lambda_result = IDASetRootDirection(ida_mem, rootdir_1d_ptr);
      return lambda_result;
    };

    return IDASetRootDirection_adapt_arr_ptr_to_std_vector(ida_mem, rootdir_1d);
  },
  nb::arg("ida_mem"), nb::arg("rootdir_1d"));

m.def("IDASetNoInactiveRootWarn", IDASetNoInactiveRootWarn, nb::arg("ida_mem"));

m.def(
  "IDASolve",
  [](void* ida_mem, sunrealtype tout, N_Vector yret, N_Vector ypret,
     int itask) -> std::tuple<int, sunrealtype>
  {
    auto IDASolve_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, sunrealtype tout, N_Vector yret, N_Vector ypret,
         int itask) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDASolve(ida_mem, tout, &tret_adapt_modifiable, yret, ypret, itask);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return IDASolve_adapt_modifiable_immutable_to_return(ida_mem, tout, yret,
                                                         ypret, itask);
  },
  nb::arg("ida_mem"), nb::arg("tout"), nb::arg("yret"), nb::arg("ypret"),
  nb::arg("itask"), "Solver function");

m.def("IDAComputeY", IDAComputeY, nb::arg("ida_mem"), nb::arg("ycor"),
      nb::arg("y"));

m.def("IDAComputeYp", IDAComputeYp, nb::arg("ida_mem"), nb::arg("ycor"),
      nb::arg("yp"));

m.def(
  "IDAComputeYSens",
  [](void* ida_mem, std::vector<N_Vector> ycor_1d,
     std::vector<N_Vector> yyS_1d) -> int
  {
    auto IDAComputeYSens_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> ycor_1d,
         std::vector<N_Vector> yyS_1d) -> int
    {
      N_Vector* ycor_1d_ptr = ycor_1d.empty() ? nullptr : ycor_1d.data();
      N_Vector* yyS_1d_ptr  = yyS_1d.empty() ? nullptr : yyS_1d.data();

      auto lambda_result = IDAComputeYSens(ida_mem, ycor_1d_ptr, yyS_1d_ptr);
      return lambda_result;
    };

    return IDAComputeYSens_adapt_arr_ptr_to_std_vector(ida_mem, ycor_1d, yyS_1d);
  },
  nb::arg("ida_mem"), nb::arg("ycor_1d"), nb::arg("yyS_1d"));

m.def(
  "IDAComputeYpSens",
  [](void* ida_mem, std::vector<N_Vector> ycor_1d,
     std::vector<N_Vector> ypS_1d) -> int
  {
    auto IDAComputeYpSens_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> ycor_1d,
         std::vector<N_Vector> ypS_1d) -> int
    {
      N_Vector* ycor_1d_ptr = ycor_1d.empty() ? nullptr : ycor_1d.data();
      N_Vector* ypS_1d_ptr  = ypS_1d.empty() ? nullptr : ypS_1d.data();

      auto lambda_result = IDAComputeYpSens(ida_mem, ycor_1d_ptr, ypS_1d_ptr);
      return lambda_result;
    };

    return IDAComputeYpSens_adapt_arr_ptr_to_std_vector(ida_mem, ycor_1d, ypS_1d);
  },
  nb::arg("ida_mem"), nb::arg("ycor_1d"), nb::arg("ypS_1d"));

m.def("IDAGetDky", IDAGetDky, nb::arg("ida_mem"), nb::arg("t"), nb::arg("k"),
      nb::arg("dky"));

m.def(
  "IDAGetNumSteps",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumSteps_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nsteps_adapt_modifiable;

      int r = IDAGetNumSteps(ida_mem, &nsteps_adapt_modifiable);
      return std::make_tuple(r, nsteps_adapt_modifiable);
    };

    return IDAGetNumSteps_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumResEvals",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumResEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nrevals_adapt_modifiable;

      int r = IDAGetNumResEvals(ida_mem, &nrevals_adapt_modifiable);
      return std::make_tuple(r, nrevals_adapt_modifiable);
    };

    return IDAGetNumResEvals_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumLinSolvSetups",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nlinsetups_adapt_modifiable;

      int r = IDAGetNumLinSolvSetups(ida_mem, &nlinsetups_adapt_modifiable);
      return std::make_tuple(r, nlinsetups_adapt_modifiable);
    };

    return IDAGetNumLinSolvSetups_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumErrTestFails",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long netfails_adapt_modifiable;

      int r = IDAGetNumErrTestFails(ida_mem, &netfails_adapt_modifiable);
      return std::make_tuple(r, netfails_adapt_modifiable);
    };

    return IDAGetNumErrTestFails_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumBacktrackOps",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumBacktrackOps_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nbacktr_adapt_modifiable;

      int r = IDAGetNumBacktrackOps(ida_mem, &nbacktr_adapt_modifiable);
      return std::make_tuple(r, nbacktr_adapt_modifiable);
    };

    return IDAGetNumBacktrackOps_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def("IDAGetConsistentIC", IDAGetConsistentIC, nb::arg("ida_mem"),
      nb::arg("yy0_mod"), nb::arg("yp0_mod"));

m.def(
  "IDAGetLastOrder",
  [](void* ida_mem) -> std::tuple<int, int>
  {
    auto IDAGetLastOrder_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, int>
    {
      int klast_adapt_modifiable;

      int r = IDAGetLastOrder(ida_mem, &klast_adapt_modifiable);
      return std::make_tuple(r, klast_adapt_modifiable);
    };

    return IDAGetLastOrder_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetCurrentOrder",
  [](void* ida_mem) -> std::tuple<int, int>
  {
    auto IDAGetCurrentOrder_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, int>
    {
      int kcur_adapt_modifiable;

      int r = IDAGetCurrentOrder(ida_mem, &kcur_adapt_modifiable);
      return std::make_tuple(r, kcur_adapt_modifiable);
    };

    return IDAGetCurrentOrder_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetCurrentCj",
  [](void* ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetCurrentCj_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype cj_adapt_modifiable;

      int r = IDAGetCurrentCj(ida_mem, &cj_adapt_modifiable);
      return std::make_tuple(r, cj_adapt_modifiable);
    };

    return IDAGetCurrentCj_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetCurrentY",
  [](void* ida_mem) -> std::tuple<int, N_Vector>
  {
    auto IDAGetCurrentY_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, N_Vector>
    {
      N_Vector ycur_adapt_modifiable;

      int r = IDAGetCurrentY(ida_mem, &ycur_adapt_modifiable);
      return std::make_tuple(r, ycur_adapt_modifiable);
    };

    return IDAGetCurrentY_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

m.def(
  "IDAGetCurrentYp",
  [](void* ida_mem) -> std::tuple<int, N_Vector>
  {
    auto IDAGetCurrentYp_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, N_Vector>
    {
      N_Vector ypcur_adapt_modifiable;

      int r = IDAGetCurrentYp(ida_mem, &ypcur_adapt_modifiable);
      return std::make_tuple(r, ypcur_adapt_modifiable);
    };

    return IDAGetCurrentYp_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

m.def(
  "IDAGetActualInitStep",
  [](void* ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetActualInitStep_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hinused_adapt_modifiable;

      int r = IDAGetActualInitStep(ida_mem, &hinused_adapt_modifiable);
      return std::make_tuple(r, hinused_adapt_modifiable);
    };

    return IDAGetActualInitStep_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetLastStep",
  [](void* ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetLastStep_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hlast_adapt_modifiable;

      int r = IDAGetLastStep(ida_mem, &hlast_adapt_modifiable);
      return std::make_tuple(r, hlast_adapt_modifiable);
    };

    return IDAGetLastStep_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetCurrentStep",
  [](void* ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetCurrentStep_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hcur_adapt_modifiable;

      int r = IDAGetCurrentStep(ida_mem, &hcur_adapt_modifiable);
      return std::make_tuple(r, hcur_adapt_modifiable);
    };

    return IDAGetCurrentStep_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetCurrentTime",
  [](void* ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetCurrentTime_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tcur_adapt_modifiable;

      int r = IDAGetCurrentTime(ida_mem, &tcur_adapt_modifiable);
      return std::make_tuple(r, tcur_adapt_modifiable);
    };

    return IDAGetCurrentTime_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetTolScaleFactor",
  [](void* ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetTolScaleFactor_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tolsfact_adapt_modifiable;

      int r = IDAGetTolScaleFactor(ida_mem, &tolsfact_adapt_modifiable);
      return std::make_tuple(r, tolsfact_adapt_modifiable);
    };

    return IDAGetTolScaleFactor_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def("IDAGetErrWeights", IDAGetErrWeights, nb::arg("ida_mem"),
      nb::arg("eweight"));

m.def("IDAGetEstLocalErrors", IDAGetEstLocalErrors, nb::arg("ida_mem"),
      nb::arg("ele"));

m.def(
  "IDAGetNumConstraintFails",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumConstraintFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long num_fails_out_adapt_modifiable;

      int r = IDAGetNumConstraintFails(ida_mem, &num_fails_out_adapt_modifiable);
      return std::make_tuple(r, num_fails_out_adapt_modifiable);
    };

    return IDAGetNumConstraintFails_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumConstraintCorrections",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumConstraintCorrections_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long num_corrections_out_adapt_modifiable;

      int r =
        IDAGetNumConstraintCorrections(ida_mem,
                                       &num_corrections_out_adapt_modifiable);
      return std::make_tuple(r, num_corrections_out_adapt_modifiable);
    };

    return IDAGetNumConstraintCorrections_adapt_modifiable_immutable_to_return(
      ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumGEvals",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumGEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long ngevals_adapt_modifiable;

      int r = IDAGetNumGEvals(ida_mem, &ngevals_adapt_modifiable);
      return std::make_tuple(r, ngevals_adapt_modifiable);
    };

    return IDAGetNumGEvals_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetRootInfo",
  [](void* ida_mem, std::vector<int> rootsfound_1d) -> int
  {
    auto IDAGetRootInfo_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<int> rootsfound_1d) -> int
    {
      int* rootsfound_1d_ptr = rootsfound_1d.empty() ? nullptr
                                                     : rootsfound_1d.data();

      auto lambda_result = IDAGetRootInfo(ida_mem, rootsfound_1d_ptr);
      return lambda_result;
    };

    return IDAGetRootInfo_adapt_arr_ptr_to_std_vector(ida_mem, rootsfound_1d);
  },
  nb::arg("ida_mem"), nb::arg("rootsfound_1d"));

m.def(
  "IDAGetIntegratorStats",
  [](void* ida_mem) -> std::tuple<int, long, long, long, long, int, int,
                                  sunrealtype, sunrealtype, sunrealtype, sunrealtype>
  {
    auto IDAGetIntegratorStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem)
      -> std::tuple<int, long, long, long, long, int, int, sunrealtype,
                    sunrealtype, sunrealtype, sunrealtype>
    {
      long nsteps_adapt_modifiable;
      long nrevals_adapt_modifiable;
      long nlinsetups_adapt_modifiable;
      long netfails_adapt_modifiable;
      int qlast_adapt_modifiable;
      int qcur_adapt_modifiable;
      sunrealtype hinused_adapt_modifiable;
      sunrealtype hlast_adapt_modifiable;
      sunrealtype hcur_adapt_modifiable;
      sunrealtype tcur_adapt_modifiable;

      int r =
        IDAGetIntegratorStats(ida_mem, &nsteps_adapt_modifiable,
                              &nrevals_adapt_modifiable,
                              &nlinsetups_adapt_modifiable,
                              &netfails_adapt_modifiable,
                              &qlast_adapt_modifiable, &qcur_adapt_modifiable,
                              &hinused_adapt_modifiable, &hlast_adapt_modifiable,
                              &hcur_adapt_modifiable, &tcur_adapt_modifiable);
      return std::make_tuple(r, nsteps_adapt_modifiable, nrevals_adapt_modifiable,
                             nlinsetups_adapt_modifiable,
                             netfails_adapt_modifiable, qlast_adapt_modifiable,
                             qcur_adapt_modifiable, hinused_adapt_modifiable,
                             hlast_adapt_modifiable, hcur_adapt_modifiable,
                             tcur_adapt_modifiable);
    };

    return IDAGetIntegratorStats_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumNonlinSolvIters",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nniters_adapt_modifiable;

      int r = IDAGetNumNonlinSolvIters(ida_mem, &nniters_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable);
    };

    return IDAGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumNonlinSolvConvFails",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nnfails_adapt_modifiable;

      int r = IDAGetNumNonlinSolvConvFails(ida_mem, &nnfails_adapt_modifiable);
      return std::make_tuple(r, nnfails_adapt_modifiable);
    };

    return IDAGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(
      ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNonlinSolvStats",
  [](void* ida_mem) -> std::tuple<int, long, long>
  {
    auto IDAGetNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long, long>
    {
      long nniters_adapt_modifiable;
      long nnfails_adapt_modifiable;

      int r = IDAGetNonlinSolvStats(ida_mem, &nniters_adapt_modifiable,
                                    &nnfails_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable,
                             nnfails_adapt_modifiable);
    };

    return IDAGetNonlinSolvStats_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumStepSolveFails",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumStepSolveFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nncfails_adapt_modifiable;

      int r = IDAGetNumStepSolveFails(ida_mem, &nncfails_adapt_modifiable);
      return std::make_tuple(r, nncfails_adapt_modifiable);
    };

    return IDAGetNumStepSolveFails_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def("IDAPrintAllStats", IDAPrintAllStats, nb::arg("ida_mem"),
      nb::arg("outfile"), nb::arg("fmt"));

m.def("IDAGetReturnFlagName", IDAGetReturnFlagName, nb::arg("flag"));

m.def("IDAQuadReInit", IDAQuadReInit, nb::arg("ida_mem"), nb::arg("yQ0"));

m.def("IDAQuadSStolerances", IDAQuadSStolerances, nb::arg("ida_mem"),
      nb::arg("reltolQ"), nb::arg("abstolQ"));

m.def("IDAQuadSVtolerances", IDAQuadSVtolerances, nb::arg("ida_mem"),
      nb::arg("reltolQ"), nb::arg("abstolQ"));

m.def("IDASetQuadErrCon", IDASetQuadErrCon, nb::arg("ida_mem"),
      nb::arg("errconQ"));

m.def(
  "IDAGetQuad",
  [](void* ida_mem, N_Vector yQout) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetQuad_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, N_Vector yQout) -> std::tuple<int, sunrealtype>
    {
      sunrealtype t_adapt_modifiable;

      int r = IDAGetQuad(ida_mem, &t_adapt_modifiable, yQout);
      return std::make_tuple(r, t_adapt_modifiable);
    };

    return IDAGetQuad_adapt_modifiable_immutable_to_return(ida_mem, yQout);
  },
  nb::arg("ida_mem"), nb::arg("yQout"));

m.def("IDAGetQuadDky", IDAGetQuadDky, nb::arg("ida_mem"), nb::arg("t"),
      nb::arg("k"), nb::arg("dky"));

m.def(
  "IDAGetQuadNumRhsEvals",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nrhsQevals_adapt_modifiable;

      int r = IDAGetQuadNumRhsEvals(ida_mem, &nrhsQevals_adapt_modifiable);
      return std::make_tuple(r, nrhsQevals_adapt_modifiable);
    };

    return IDAGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetQuadNumErrTestFails",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nQetfails_adapt_modifiable;

      int r = IDAGetQuadNumErrTestFails(ida_mem, &nQetfails_adapt_modifiable);
      return std::make_tuple(r, nQetfails_adapt_modifiable);
    };

    return IDAGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def("IDAGetQuadErrWeights", IDAGetQuadErrWeights, nb::arg("ida_mem"),
      nb::arg("eQweight"));

m.def(
  "IDAGetQuadStats",
  [](void* ida_mem) -> std::tuple<int, long, long>
  {
    auto IDAGetQuadStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long, long>
    {
      long nrhsQevals_adapt_modifiable;
      long nQetfails_adapt_modifiable;

      int r = IDAGetQuadStats(ida_mem, &nrhsQevals_adapt_modifiable,
                              &nQetfails_adapt_modifiable);
      return std::make_tuple(r, nrhsQevals_adapt_modifiable,
                             nQetfails_adapt_modifiable);
    };

    return IDAGetQuadStats_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDASensReInit",
  [](void* ida_mem, int ism, std::vector<N_Vector> yS0_1d,
     std::vector<N_Vector> ypS0_1d) -> int
  {
    auto IDASensReInit_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, int ism, std::vector<N_Vector> yS0_1d,
         std::vector<N_Vector> ypS0_1d) -> int
    {
      N_Vector* yS0_1d_ptr  = yS0_1d.empty() ? nullptr : yS0_1d.data();
      N_Vector* ypS0_1d_ptr = ypS0_1d.empty() ? nullptr : ypS0_1d.data();

      auto lambda_result = IDASensReInit(ida_mem, ism, yS0_1d_ptr, ypS0_1d_ptr);
      return lambda_result;
    };

    return IDASensReInit_adapt_arr_ptr_to_std_vector(ida_mem, ism, yS0_1d,
                                                     ypS0_1d);
  },
  nb::arg("ida_mem"), nb::arg("ism"), nb::arg("yS0_1d"), nb::arg("ypS0_1d"));

m.def(
  "IDASensSStolerances",
  [](void* ida_mem, sunrealtype reltolS, sundials4py::Array1d abstolS_1d) -> int
  {
    auto IDASensSStolerances_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype reltolS, sundials4py::Array1d abstolS_1d) -> int
    {
      sunrealtype* abstolS_1d_ptr = abstolS_1d.size() == 0 ? nullptr
                                                           : abstolS_1d.data();

      auto lambda_result = IDASensSStolerances(ida_mem, reltolS, abstolS_1d_ptr);
      return lambda_result;
    };

    return IDASensSStolerances_adapt_arr_ptr_to_std_vector(ida_mem, reltolS,
                                                           abstolS_1d);
  },
  nb::arg("ida_mem"), nb::arg("reltolS"), nb::arg("abstolS_1d"));

m.def(
  "IDASensSVtolerances",
  [](void* ida_mem, sunrealtype reltolS, std::vector<N_Vector> abstolS_1d) -> int
  {
    auto IDASensSVtolerances_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype reltolS,
         std::vector<N_Vector> abstolS_1d) -> int
    {
      N_Vector* abstolS_1d_ptr = abstolS_1d.empty() ? nullptr : abstolS_1d.data();

      auto lambda_result = IDASensSVtolerances(ida_mem, reltolS, abstolS_1d_ptr);
      return lambda_result;
    };

    return IDASensSVtolerances_adapt_arr_ptr_to_std_vector(ida_mem, reltolS,
                                                           abstolS_1d);
  },
  nb::arg("ida_mem"), nb::arg("reltolS"), nb::arg("abstolS_1d"));

m.def("IDASensEEtolerances", IDASensEEtolerances, nb::arg("ida_mem"));

m.def(
  "IDAGetSensConsistentIC",
  [](void* ida_mem, std::vector<N_Vector> yyS0_1d,
     std::vector<N_Vector> ypS0_1d) -> int
  {
    auto IDAGetSensConsistentIC_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> yyS0_1d,
         std::vector<N_Vector> ypS0_1d) -> int
    {
      N_Vector* yyS0_1d_ptr = yyS0_1d.empty() ? nullptr : yyS0_1d.data();
      N_Vector* ypS0_1d_ptr = ypS0_1d.empty() ? nullptr : ypS0_1d.data();

      auto lambda_result = IDAGetSensConsistentIC(ida_mem, yyS0_1d_ptr,
                                                  ypS0_1d_ptr);
      return lambda_result;
    };

    return IDAGetSensConsistentIC_adapt_arr_ptr_to_std_vector(ida_mem, yyS0_1d,
                                                              ypS0_1d);
  },
  nb::arg("ida_mem"), nb::arg("yyS0_1d"), nb::arg("ypS0_1d"));

m.def("IDASetSensDQMethod", IDASetSensDQMethod, nb::arg("ida_mem"),
      nb::arg("DQtype"), nb::arg("DQrhomax"));

m.def("IDASetSensErrCon", IDASetSensErrCon, nb::arg("ida_mem"),
      nb::arg("errconS"));

m.def("IDASetSensMaxNonlinIters", IDASetSensMaxNonlinIters, nb::arg("ida_mem"),
      nb::arg("maxcorS"));

m.def(
  "IDASetSensParams",
  [](void* ida_mem, sundials4py::Array1d p_1d, sundials4py::Array1d pbar_1d,
     std::vector<int> plist_1d) -> int
  {
    auto IDASetSensParams_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sundials4py::Array1d p_1d, sundials4py::Array1d pbar_1d,
         std::vector<int> plist_1d) -> int
    {
      sunrealtype* p_1d_ptr    = p_1d.size() == 0 ? nullptr : p_1d.data();
      sunrealtype* pbar_1d_ptr = pbar_1d.size() == 0 ? nullptr : pbar_1d.data();
      int* plist_1d_ptr        = plist_1d.empty() ? nullptr : plist_1d.data();

      auto lambda_result = IDASetSensParams(ida_mem, p_1d_ptr, pbar_1d_ptr,
                                            plist_1d_ptr);
      return lambda_result;
    };

    return IDASetSensParams_adapt_arr_ptr_to_std_vector(ida_mem, p_1d, pbar_1d,
                                                        plist_1d);
  },
  nb::arg("ida_mem"), nb::arg("p_1d"), nb::arg("pbar_1d"), nb::arg("plist_1d"));

m.def("IDASetNonlinearSolverSensSim", IDASetNonlinearSolverSensSim,
      nb::arg("ida_mem"), nb::arg("NLS"));

m.def("IDASetNonlinearSolverSensStg", IDASetNonlinearSolverSensStg,
      nb::arg("ida_mem"), nb::arg("NLS"));

m.def("IDASensToggleOff", IDASensToggleOff, nb::arg("ida_mem"));

m.def(
  "IDAGetSens",
  [](void* ida_mem, std::vector<N_Vector> yySout_1d) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetSens_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype* tret, std::vector<N_Vector> yySout_1d) -> int
    {
      N_Vector* yySout_1d_ptr = yySout_1d.empty() ? nullptr : yySout_1d.data();

      auto lambda_result = IDAGetSens(ida_mem, tret, yySout_1d_ptr);
      return lambda_result;
    };
    auto IDAGetSens_adapt_modifiable_immutable_to_return =
      [&IDAGetSens_adapt_arr_ptr_to_std_vector](void* ida_mem,
                                                std::vector<N_Vector> yySout_1d)
      -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDAGetSens_adapt_arr_ptr_to_std_vector(ida_mem,
                                                     &tret_adapt_modifiable,
                                                     yySout_1d);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return IDAGetSens_adapt_modifiable_immutable_to_return(ida_mem, yySout_1d);
  },
  nb::arg("ida_mem"), nb::arg("yySout_1d"));

m.def(
  "IDAGetSens1",
  [](void* ida_mem, int is, N_Vector yySret) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetSens1_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int is, N_Vector yySret) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDAGetSens1(ida_mem, &tret_adapt_modifiable, is, yySret);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return IDAGetSens1_adapt_modifiable_immutable_to_return(ida_mem, is, yySret);
  },
  nb::arg("ida_mem"), nb::arg("is_"), nb::arg("yySret"));

m.def(
  "IDAGetSensDky",
  [](void* ida_mem, sunrealtype t, int k, std::vector<N_Vector> dkyS_1d) -> int
  {
    auto IDAGetSensDky_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype t, int k, std::vector<N_Vector> dkyS_1d) -> int
    {
      N_Vector* dkyS_1d_ptr = dkyS_1d.empty() ? nullptr : dkyS_1d.data();

      auto lambda_result = IDAGetSensDky(ida_mem, t, k, dkyS_1d_ptr);
      return lambda_result;
    };

    return IDAGetSensDky_adapt_arr_ptr_to_std_vector(ida_mem, t, k, dkyS_1d);
  },
  nb::arg("ida_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dkyS_1d"));

m.def("IDAGetSensDky1", IDAGetSensDky1, nb::arg("ida_mem"), nb::arg("t"),
      nb::arg("k"), nb::arg("is_"), nb::arg("dkyS"));

m.def(
  "IDAGetSensNumResEvals",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetSensNumResEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nresSevals_adapt_modifiable;

      int r = IDAGetSensNumResEvals(ida_mem, &nresSevals_adapt_modifiable);
      return std::make_tuple(r, nresSevals_adapt_modifiable);
    };

    return IDAGetSensNumResEvals_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumResEvalsSens",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumResEvalsSens_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nresevalsS_adapt_modifiable;

      int r = IDAGetNumResEvalsSens(ida_mem, &nresevalsS_adapt_modifiable);
      return std::make_tuple(r, nresevalsS_adapt_modifiable);
    };

    return IDAGetNumResEvalsSens_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetSensNumErrTestFails",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetSensNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nSetfails_adapt_modifiable;

      int r = IDAGetSensNumErrTestFails(ida_mem, &nSetfails_adapt_modifiable);
      return std::make_tuple(r, nSetfails_adapt_modifiable);
    };

    return IDAGetSensNumErrTestFails_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetSensNumLinSolvSetups",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nlinsetupsS_adapt_modifiable;

      int r = IDAGetSensNumLinSolvSetups(ida_mem, &nlinsetupsS_adapt_modifiable);
      return std::make_tuple(r, nlinsetupsS_adapt_modifiable);
    };

    return IDAGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return(
      ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetSensErrWeights",
  [](void* ida_mem, std::vector<N_Vector> eSweight_1d) -> int
  {
    auto IDAGetSensErrWeights_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> eSweight_1d) -> int
    {
      N_Vector* eSweight_1d_ptr = eSweight_1d.empty() ? nullptr
                                                      : eSweight_1d.data();

      auto lambda_result = IDAGetSensErrWeights(ida_mem, eSweight_1d_ptr);
      return lambda_result;
    };

    return IDAGetSensErrWeights_adapt_arr_ptr_to_std_vector(ida_mem, eSweight_1d);
  },
  nb::arg("ida_mem"), nb::arg("eSweight_1d"));

m.def(
  "IDAGetSensStats",
  [](void* ida_mem) -> std::tuple<int, long, long, long, long>
  {
    auto IDAGetSensStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long, long, long, long>
    {
      long nresSevals_adapt_modifiable;
      long nresevalsS_adapt_modifiable;
      long nSetfails_adapt_modifiable;
      long nlinsetupsS_adapt_modifiable;

      int r = IDAGetSensStats(ida_mem, &nresSevals_adapt_modifiable,
                              &nresevalsS_adapt_modifiable,
                              &nSetfails_adapt_modifiable,
                              &nlinsetupsS_adapt_modifiable);
      return std::make_tuple(r, nresSevals_adapt_modifiable,
                             nresevalsS_adapt_modifiable,
                             nSetfails_adapt_modifiable,
                             nlinsetupsS_adapt_modifiable);
    };

    return IDAGetSensStats_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetSensNumNonlinSolvIters",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nSniters_adapt_modifiable;

      int r = IDAGetSensNumNonlinSolvIters(ida_mem, &nSniters_adapt_modifiable);
      return std::make_tuple(r, nSniters_adapt_modifiable);
    };

    return IDAGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return(
      ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetSensNumNonlinSolvConvFails",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nSnfails_adapt_modifiable;

      int r = IDAGetSensNumNonlinSolvConvFails(ida_mem,
                                               &nSnfails_adapt_modifiable);
      return std::make_tuple(r, nSnfails_adapt_modifiable);
    };

    return IDAGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(
      ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetSensNonlinSolvStats",
  [](void* ida_mem) -> std::tuple<int, long, long>
  {
    auto IDAGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long, long>
    {
      long nSniters_adapt_modifiable;
      long nSnfails_adapt_modifiable;

      int r = IDAGetSensNonlinSolvStats(ida_mem, &nSniters_adapt_modifiable,
                                        &nSnfails_adapt_modifiable);
      return std::make_tuple(r, nSniters_adapt_modifiable,
                             nSnfails_adapt_modifiable);
    };

    return IDAGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumStepSensSolveFails",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nSncfails_adapt_modifiable;

      int r = IDAGetNumStepSensSolveFails(ida_mem, &nSncfails_adapt_modifiable);
      return std::make_tuple(r, nSncfails_adapt_modifiable);
    };

    return IDAGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return(
      ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAQuadSensReInit",
  [](void* ida_mem, std::vector<N_Vector> yQS0_1d) -> int
  {
    auto IDAQuadSensReInit_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> yQS0_1d) -> int
    {
      N_Vector* yQS0_1d_ptr = yQS0_1d.empty() ? nullptr : yQS0_1d.data();

      auto lambda_result = IDAQuadSensReInit(ida_mem, yQS0_1d_ptr);
      return lambda_result;
    };

    return IDAQuadSensReInit_adapt_arr_ptr_to_std_vector(ida_mem, yQS0_1d);
  },
  nb::arg("ida_mem"), nb::arg("yQS0_1d"));

m.def(
  "IDAQuadSensSStolerances",
  [](void* ida_mem, sunrealtype reltolQS, sundials4py::Array1d abstolQS_1d) -> int
  {
    auto IDAQuadSensSStolerances_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype reltolQS,
         sundials4py::Array1d abstolQS_1d) -> int
    {
      sunrealtype* abstolQS_1d_ptr =
        abstolQS_1d.size() == 0 ? nullptr : abstolQS_1d.data();

      auto lambda_result = IDAQuadSensSStolerances(ida_mem, reltolQS,
                                                   abstolQS_1d_ptr);
      return lambda_result;
    };

    return IDAQuadSensSStolerances_adapt_arr_ptr_to_std_vector(ida_mem, reltolQS,
                                                               abstolQS_1d);
  },
  nb::arg("ida_mem"), nb::arg("reltolQS"), nb::arg("abstolQS_1d"));

m.def(
  "IDAQuadSensSVtolerances",
  [](void* ida_mem, sunrealtype reltolQS, std::vector<N_Vector> abstolQS_1d) -> int
  {
    auto IDAQuadSensSVtolerances_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype reltolQS,
         std::vector<N_Vector> abstolQS_1d) -> int
    {
      N_Vector* abstolQS_1d_ptr = abstolQS_1d.empty() ? nullptr
                                                      : abstolQS_1d.data();

      auto lambda_result = IDAQuadSensSVtolerances(ida_mem, reltolQS,
                                                   abstolQS_1d_ptr);
      return lambda_result;
    };

    return IDAQuadSensSVtolerances_adapt_arr_ptr_to_std_vector(ida_mem, reltolQS,
                                                               abstolQS_1d);
  },
  nb::arg("ida_mem"), nb::arg("reltolQS"), nb::arg("abstolQS_1d"));

m.def("IDAQuadSensEEtolerances", IDAQuadSensEEtolerances, nb::arg("ida_mem"));

m.def("IDASetQuadSensErrCon", IDASetQuadSensErrCon, nb::arg("ida_mem"),
      nb::arg("errconQS"));

m.def(
  "IDAGetQuadSens",
  [](void* ida_mem,
     std::vector<N_Vector> yyQSout_1d) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetQuadSens_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype* tret, std::vector<N_Vector> yyQSout_1d) -> int
    {
      N_Vector* yyQSout_1d_ptr = yyQSout_1d.empty() ? nullptr : yyQSout_1d.data();

      auto lambda_result = IDAGetQuadSens(ida_mem, tret, yyQSout_1d_ptr);
      return lambda_result;
    };
    auto IDAGetQuadSens_adapt_modifiable_immutable_to_return =
      [&IDAGetQuadSens_adapt_arr_ptr_to_std_vector](void* ida_mem,
                                                    std::vector<N_Vector> yyQSout_1d)
      -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDAGetQuadSens_adapt_arr_ptr_to_std_vector(ida_mem,
                                                         &tret_adapt_modifiable,
                                                         yyQSout_1d);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return IDAGetQuadSens_adapt_modifiable_immutable_to_return(ida_mem,
                                                               yyQSout_1d);
  },
  nb::arg("ida_mem"), nb::arg("yyQSout_1d"));

m.def(
  "IDAGetQuadSens1",
  [](void* ida_mem, int is, N_Vector yyQSret) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetQuadSens1_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int is, N_Vector yyQSret) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDAGetQuadSens1(ida_mem, &tret_adapt_modifiable, is, yyQSret);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return IDAGetQuadSens1_adapt_modifiable_immutable_to_return(ida_mem, is,
                                                                yyQSret);
  },
  nb::arg("ida_mem"), nb::arg("is_"), nb::arg("yyQSret"));

m.def(
  "IDAGetQuadSensDky",
  [](void* ida_mem, sunrealtype t, int k, std::vector<N_Vector> dkyQS_1d) -> int
  {
    auto IDAGetQuadSensDky_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, sunrealtype t, int k, std::vector<N_Vector> dkyQS_1d) -> int
    {
      N_Vector* dkyQS_1d_ptr = dkyQS_1d.empty() ? nullptr : dkyQS_1d.data();

      auto lambda_result = IDAGetQuadSensDky(ida_mem, t, k, dkyQS_1d_ptr);
      return lambda_result;
    };

    return IDAGetQuadSensDky_adapt_arr_ptr_to_std_vector(ida_mem, t, k, dkyQS_1d);
  },
  nb::arg("ida_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dkyQS_1d"));

m.def("IDAGetQuadSensDky1", IDAGetQuadSensDky1, nb::arg("ida_mem"),
      nb::arg("t"), nb::arg("k"), nb::arg("is_"), nb::arg("dkyQS"));

m.def(
  "IDAGetQuadSensNumRhsEvals",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nrhsQSevals_adapt_modifiable;

      int r = IDAGetQuadSensNumRhsEvals(ida_mem, &nrhsQSevals_adapt_modifiable);
      return std::make_tuple(r, nrhsQSevals_adapt_modifiable);
    };

    return IDAGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetQuadSensNumErrTestFails",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nQSetfails_adapt_modifiable;

      int r = IDAGetQuadSensNumErrTestFails(ida_mem,
                                            &nQSetfails_adapt_modifiable);
      return std::make_tuple(r, nQSetfails_adapt_modifiable);
    };

    return IDAGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return(
      ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetQuadSensErrWeights",
  [](void* ida_mem, std::vector<N_Vector> eQSweight_1d) -> int
  {
    auto IDAGetQuadSensErrWeights_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, std::vector<N_Vector> eQSweight_1d) -> int
    {
      N_Vector* eQSweight_1d_ptr = eQSweight_1d.empty() ? nullptr
                                                        : eQSweight_1d.data();

      auto lambda_result = IDAGetQuadSensErrWeights(ida_mem, eQSweight_1d_ptr);
      return lambda_result;
    };

    return IDAGetQuadSensErrWeights_adapt_arr_ptr_to_std_vector(ida_mem,
                                                                eQSweight_1d);
  },
  nb::arg("ida_mem"), nb::arg("eQSweight_1d"));

m.def(
  "IDAGetQuadSensStats",
  [](void* ida_mem) -> std::tuple<int, long, long>
  {
    auto IDAGetQuadSensStats_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long, long>
    {
      long nrhsQSevals_adapt_modifiable;
      long nQSetfails_adapt_modifiable;

      int r = IDAGetQuadSensStats(ida_mem, &nrhsQSevals_adapt_modifiable,
                                  &nQSetfails_adapt_modifiable);
      return std::make_tuple(r, nrhsQSevals_adapt_modifiable,
                             nQSetfails_adapt_modifiable);
    };

    return IDAGetQuadSensStats_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def("IDAAdjInit", IDAAdjInit, nb::arg("ida_mem"), nb::arg("steps"),
      nb::arg("interp"));

m.def("IDAAdjReInit", IDAAdjReInit, nb::arg("ida_mem"));

m.def(
  "IDACreateB",
  [](void* ida_mem) -> std::tuple<int, int>
  {
    auto IDACreateB_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, int>
    {
      int which_adapt_modifiable;

      int r = IDACreateB(ida_mem, &which_adapt_modifiable);
      return std::make_tuple(r, which_adapt_modifiable);
    };

    return IDACreateB_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def("IDAReInitB", IDAReInitB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("tB0"), nb::arg("yyB0"), nb::arg("ypB0"));

m.def("IDASStolerancesB", IDASStolerancesB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("relTolB"), nb::arg("absTolB"));

m.def("IDASVtolerancesB", IDASVtolerancesB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("relTolB"), nb::arg("absTolB"));

m.def("IDAQuadReInitB", IDAQuadReInitB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("yQB0"));

m.def("IDAQuadSStolerancesB", IDAQuadSStolerancesB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("reltolQB"), nb::arg("abstolQB"));

m.def("IDAQuadSVtolerancesB", IDAQuadSVtolerancesB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("reltolQB"), nb::arg("abstolQB"));

m.def("IDACalcICB", IDACalcICB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("tout1"), nb::arg("yy0"), nb::arg("yp0"));

m.def(
  "IDACalcICBS",
  [](void* ida_mem, int which, sunrealtype tout1, N_Vector yy0, N_Vector yp0,
     std::vector<N_Vector> yyS0_1d, std::vector<N_Vector> ypS0_1d) -> int
  {
    auto IDACalcICBS_adapt_arr_ptr_to_std_vector =
      [](void* ida_mem, int which, sunrealtype tout1, N_Vector yy0, N_Vector yp0,
         std::vector<N_Vector> yyS0_1d, std::vector<N_Vector> ypS0_1d) -> int
    {
      N_Vector* yyS0_1d_ptr = yyS0_1d.empty() ? nullptr : yyS0_1d.data();
      N_Vector* ypS0_1d_ptr = ypS0_1d.empty() ? nullptr : ypS0_1d.data();

      auto lambda_result = IDACalcICBS(ida_mem, which, tout1, yy0, yp0,
                                       yyS0_1d_ptr, ypS0_1d_ptr);
      return lambda_result;
    };

    return IDACalcICBS_adapt_arr_ptr_to_std_vector(ida_mem, which, tout1, yy0,
                                                   yp0, yyS0_1d, ypS0_1d);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("tout1"), nb::arg("yy0"),
  nb::arg("yp0"), nb::arg("yyS0_1d"), nb::arg("ypS0_1d"));

m.def(
  "IDASolveF",
  [](void* ida_mem, sunrealtype tout, N_Vector yret, N_Vector ypret,
     int itask) -> std::tuple<int, sunrealtype, int>
  {
    auto IDASolveF_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, sunrealtype tout, N_Vector yret, N_Vector ypret,
         int itask) -> std::tuple<int, sunrealtype, int>
    {
      sunrealtype tret_adapt_modifiable;
      int ncheckPtr_adapt_modifiable;

      int r = IDASolveF(ida_mem, tout, &tret_adapt_modifiable, yret, ypret,
                        itask, &ncheckPtr_adapt_modifiable);
      return std::make_tuple(r, tret_adapt_modifiable,
                             ncheckPtr_adapt_modifiable);
    };

    return IDASolveF_adapt_modifiable_immutable_to_return(ida_mem, tout, yret,
                                                          ypret, itask);
  },
  nb::arg("ida_mem"), nb::arg("tout"), nb::arg("yret"), nb::arg("ypret"),
  nb::arg("itask"));

m.def("IDASolveB", IDASolveB, nb::arg("ida_mem"), nb::arg("tBout"),
      nb::arg("itaskB"));

m.def("IDAAdjSetNoSensi", IDAAdjSetNoSensi, nb::arg("ida_mem"));

m.def("IDASetMaxOrdB", IDASetMaxOrdB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("maxordB"));

m.def("IDASetMaxNumStepsB", IDASetMaxNumStepsB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("mxstepsB"));

m.def("IDASetInitStepB", IDASetInitStepB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("hinB"));

m.def("IDASetMaxStepB", IDASetMaxStepB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("hmaxB"));

m.def("IDASetSuppressAlgB", IDASetSuppressAlgB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("suppressalgB"));

m.def("IDASetIdB", IDASetIdB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("idB"));

m.def("IDASetConstraintsB", IDASetConstraintsB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("constraintsB"));

m.def("IDASetQuadErrConB", IDASetQuadErrConB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("errconQB"));

m.def("IDASetNonlinearSolverB", IDASetNonlinearSolverB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("NLS"));

m.def(
  "IDAGetB",
  [](void* ida_mem, int which, N_Vector yy,
     N_Vector yp) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetB_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which, N_Vector yy,
         N_Vector yp) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDAGetB(ida_mem, which, &tret_adapt_modifiable, yy, yp);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return IDAGetB_adapt_modifiable_immutable_to_return(ida_mem, which, yy, yp);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("yy"), nb::arg("yp"));

m.def(
  "IDAGetQuadB",
  [](void* ida_mem, int which, N_Vector qB) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetQuadB_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which, N_Vector qB) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = IDAGetQuadB(ida_mem, which, &tret_adapt_modifiable, qB);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return IDAGetQuadB_adapt_modifiable_immutable_to_return(ida_mem, which, qB);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("qB"));

m.def("IDAGetAdjIDABmem", IDAGetAdjIDABmem, nb::arg("ida_mem"), nb::arg("which"));

m.def("IDAGetConsistentICB", IDAGetConsistentICB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("yyB0"), nb::arg("ypB0"));

m.def("IDAGetAdjY", IDAGetAdjY, nb::arg("ida_mem"), nb::arg("t"), nb::arg("yy"),
      nb::arg("yp"));

m.def(
  "IDAGetAdjDataPointHermite",
  [](void* ida_mem, int which, std::optional<N_Vector> yy = std::nullopt,
     std::optional<N_Vector> yd = std::nullopt) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetAdjDataPointHermite_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which, N_Vector yy,
         N_Vector yd) -> std::tuple<int, sunrealtype>
    {
      sunrealtype t_adapt_modifiable;

      int r = IDAGetAdjDataPointHermite(ida_mem, which, &t_adapt_modifiable, yy,
                                        yd);
      return std::make_tuple(r, t_adapt_modifiable);
    };
    auto IDAGetAdjDataPointHermite_adapt_optional_arg_with_default_null =
      [&IDAGetAdjDataPointHermite_adapt_modifiable_immutable_to_return](void* ida_mem,
                                                                        int which,
                                                                        std::optional<N_Vector>
                                                                          yy =
                                                                            std::nullopt,
                                                                        std::optional<N_Vector>
                                                                          yd =
                                                                            std::nullopt)
      -> std::tuple<int, sunrealtype>
    {
      N_Vector yy_adapt_default_null = nullptr;
      if (yy.has_value()) yy_adapt_default_null = yy.value();
      N_Vector yd_adapt_default_null = nullptr;
      if (yd.has_value()) yd_adapt_default_null = yd.value();

      auto lambda_result =
        IDAGetAdjDataPointHermite_adapt_modifiable_immutable_to_return(ida_mem,
                                                                       which,
                                                                       yy_adapt_default_null,
                                                                       yd_adapt_default_null);
      return lambda_result;
    };

    return IDAGetAdjDataPointHermite_adapt_optional_arg_with_default_null(ida_mem,
                                                                          which,
                                                                          yy, yd);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("yy").none() = nb::none(),
  nb::arg("yd").none() = nb::none());

m.def(
  "IDAGetAdjDataPointPolynomial",
  [](void* ida_mem, int which,
     std::optional<N_Vector> y = std::nullopt) -> std::tuple<int, sunrealtype, int>
  {
    auto IDAGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return =
      [](void* ida_mem, int which, N_Vector y) -> std::tuple<int, sunrealtype, int>
    {
      sunrealtype t_adapt_modifiable;
      int order_adapt_modifiable;

      int r = IDAGetAdjDataPointPolynomial(ida_mem, which, &t_adapt_modifiable,
                                           &order_adapt_modifiable, y);
      return std::make_tuple(r, t_adapt_modifiable, order_adapt_modifiable);
    };
    auto IDAGetAdjDataPointPolynomial_adapt_optional_arg_with_default_null =
      [&IDAGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return](void* ida_mem,
                                                                           int which,
                                                                           std::optional<N_Vector>
                                                                             y =
                                                                               std::nullopt)
      -> std::tuple<int, sunrealtype, int>
    {
      N_Vector y_adapt_default_null = nullptr;
      if (y.has_value()) y_adapt_default_null = y.value();

      auto lambda_result =
        IDAGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return(ida_mem,
                                                                          which,
                                                                          y_adapt_default_null);
      return lambda_result;
    };

    return IDAGetAdjDataPointPolynomial_adapt_optional_arg_with_default_null(ida_mem,
                                                                             which,
                                                                             y);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("y").none() = nb::none());
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _IDASLS_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("IDALS_SUCCESS")         = 0;
m.attr("IDALS_MEM_NULL")        = -1;
m.attr("IDALS_LMEM_NULL")       = -2;
m.attr("IDALS_ILL_INPUT")       = -3;
m.attr("IDALS_MEM_FAIL")        = -4;
m.attr("IDALS_PMEM_NULL")       = -5;
m.attr("IDALS_JACFUNC_UNRECVR") = -6;
m.attr("IDALS_JACFUNC_RECVR")   = -7;
m.attr("IDALS_SUNMAT_FAIL")     = -8;
m.attr("IDALS_SUNLS_FAIL")      = -9;
m.attr("IDALS_NO_ADJ")          = -101;
m.attr("IDALS_LMEMB_NULL")      = -102;

m.def(
  "IDASetLinearSolver",
  [](void* ida_mem, SUNLinearSolver LS,
     std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto IDASetLinearSolver_adapt_optional_arg_with_default_null =
      [](void* ida_mem, SUNLinearSolver LS,
         std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = IDASetLinearSolver(ida_mem, LS, A_adapt_default_null);
      return lambda_result;
    };

    return IDASetLinearSolver_adapt_optional_arg_with_default_null(ida_mem, LS,
                                                                   A);
  },
  nb::arg("ida_mem"), nb::arg("LS"), nb::arg("A").none() = nb::none());

m.def("IDASetEpsLin", IDASetEpsLin, nb::arg("ida_mem"), nb::arg("eplifac"));

m.def("IDASetLSNormFactor", IDASetLSNormFactor, nb::arg("ida_mem"),
      nb::arg("nrmfac"));

m.def("IDASetLinearSolutionScaling", IDASetLinearSolutionScaling,
      nb::arg("ida_mem"), nb::arg("onoff"));

m.def("IDASetIncrementFactor", IDASetIncrementFactor, nb::arg("ida_mem"),
      nb::arg("dqincfac"));

m.def(
  "IDAGetJac",
  [](void* ida_mem) -> std::tuple<int, SUNMatrix>
  {
    auto IDAGetJac_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, SUNMatrix>
    {
      SUNMatrix J_adapt_modifiable;

      int r = IDAGetJac(ida_mem, &J_adapt_modifiable);
      return std::make_tuple(r, J_adapt_modifiable);
    };

    return IDAGetJac_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

m.def(
  "IDAGetJacCj",
  [](void* ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetJacCj_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype cj_J_adapt_modifiable;

      int r = IDAGetJacCj(ida_mem, &cj_J_adapt_modifiable);
      return std::make_tuple(r, cj_J_adapt_modifiable);
    };

    return IDAGetJacCj_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetJacTime",
  [](void* ida_mem) -> std::tuple<int, sunrealtype>
  {
    auto IDAGetJacTime_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype t_J_adapt_modifiable;

      int r = IDAGetJacTime(ida_mem, &t_J_adapt_modifiable);
      return std::make_tuple(r, t_J_adapt_modifiable);
    };

    return IDAGetJacTime_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetJacNumSteps",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetJacNumSteps_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nst_J_adapt_modifiable;

      int r = IDAGetJacNumSteps(ida_mem, &nst_J_adapt_modifiable);
      return std::make_tuple(r, nst_J_adapt_modifiable);
    };

    return IDAGetJacNumSteps_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumJacEvals",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumJacEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long njevals_adapt_modifiable;

      int r = IDAGetNumJacEvals(ida_mem, &njevals_adapt_modifiable);
      return std::make_tuple(r, njevals_adapt_modifiable);
    };

    return IDAGetNumJacEvals_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumPrecEvals",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumPrecEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long npevals_adapt_modifiable;

      int r = IDAGetNumPrecEvals(ida_mem, &npevals_adapt_modifiable);
      return std::make_tuple(r, npevals_adapt_modifiable);
    };

    return IDAGetNumPrecEvals_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumPrecSolves",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumPrecSolves_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long npsolves_adapt_modifiable;

      int r = IDAGetNumPrecSolves(ida_mem, &npsolves_adapt_modifiable);
      return std::make_tuple(r, npsolves_adapt_modifiable);
    };

    return IDAGetNumPrecSolves_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumLinIters",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumLinIters_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nliters_adapt_modifiable;

      int r = IDAGetNumLinIters(ida_mem, &nliters_adapt_modifiable);
      return std::make_tuple(r, nliters_adapt_modifiable);
    };

    return IDAGetNumLinIters_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumLinConvFails",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumLinConvFails_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nlcfails_adapt_modifiable;

      int r = IDAGetNumLinConvFails(ida_mem, &nlcfails_adapt_modifiable);
      return std::make_tuple(r, nlcfails_adapt_modifiable);
    };

    return IDAGetNumLinConvFails_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumJTSetupEvals",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumJTSetupEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long njtsetups_adapt_modifiable;

      int r = IDAGetNumJTSetupEvals(ida_mem, &njtsetups_adapt_modifiable);
      return std::make_tuple(r, njtsetups_adapt_modifiable);
    };

    return IDAGetNumJTSetupEvals_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumJtimesEvals",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumJtimesEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long njvevals_adapt_modifiable;

      int r = IDAGetNumJtimesEvals(ida_mem, &njvevals_adapt_modifiable);
      return std::make_tuple(r, njvevals_adapt_modifiable);
    };

    return IDAGetNumJtimesEvals_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetNumLinResEvals",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetNumLinResEvals_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long nrevalsLS_adapt_modifiable;

      int r = IDAGetNumLinResEvals(ida_mem, &nrevalsLS_adapt_modifiable);
      return std::make_tuple(r, nrevalsLS_adapt_modifiable);
    };

    return IDAGetNumLinResEvals_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def(
  "IDAGetLastLinFlag",
  [](void* ida_mem) -> std::tuple<int, long>
  {
    auto IDAGetLastLinFlag_adapt_modifiable_immutable_to_return =
      [](void* ida_mem) -> std::tuple<int, long>
    {
      long flag_adapt_modifiable;

      int r = IDAGetLastLinFlag(ida_mem, &flag_adapt_modifiable);
      return std::make_tuple(r, flag_adapt_modifiable);
    };

    return IDAGetLastLinFlag_adapt_modifiable_immutable_to_return(ida_mem);
  },
  nb::arg("ida_mem"));

m.def("IDAGetLinReturnFlagName", IDAGetLinReturnFlagName, nb::arg("flag"));

m.def(
  "IDASetLinearSolverB",
  [](void* ida_mem, int which, SUNLinearSolver LS,
     std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto IDASetLinearSolverB_adapt_optional_arg_with_default_null =
      [](void* ida_mem, int which, SUNLinearSolver LS,
         std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = IDASetLinearSolverB(ida_mem, which, LS,
                                               A_adapt_default_null);
      return lambda_result;
    };

    return IDASetLinearSolverB_adapt_optional_arg_with_default_null(ida_mem,
                                                                    which, LS, A);
  },
  nb::arg("ida_mem"), nb::arg("which"), nb::arg("LS"),
  nb::arg("A").none() = nb::none());

m.def("IDASetEpsLinB", IDASetEpsLinB, nb::arg("ida_mem"), nb::arg("which"),
      nb::arg("eplifacB"));

m.def("IDASetLSNormFactorB", IDASetLSNormFactorB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("nrmfacB"));

m.def("IDASetLinearSolutionScalingB", IDASetLinearSolutionScalingB,
      nb::arg("ida_mem"), nb::arg("which"), nb::arg("onoffB"));

m.def("IDASetIncrementFactorB", IDASetIncrementFactorB, nb::arg("ida_mem"),
      nb::arg("which"), nb::arg("dqincfacB"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
