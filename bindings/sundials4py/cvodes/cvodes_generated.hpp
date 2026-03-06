// #ifndef _CVODES_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("CV_ADAMS")               = 1;
m.attr("CV_BDF")                 = 2;
m.attr("CV_NORMAL")              = 1;
m.attr("CV_ONE_STEP")            = 2;
m.attr("CV_SIMULTANEOUS")        = 1;
m.attr("CV_STAGGERED")           = 2;
m.attr("CV_STAGGERED1")          = 3;
m.attr("CV_CENTERED")            = 1;
m.attr("CV_FORWARD")             = 2;
m.attr("CV_HERMITE")             = 1;
m.attr("CV_POLYNOMIAL")          = 2;
m.attr("CV_SUCCESS")             = 0;
m.attr("CV_TSTOP_RETURN")        = 1;
m.attr("CV_ROOT_RETURN")         = 2;
m.attr("CV_WARNING")             = 99;
m.attr("CV_TOO_MUCH_WORK")       = -1;
m.attr("CV_TOO_MUCH_ACC")        = -2;
m.attr("CV_ERR_FAILURE")         = -3;
m.attr("CV_CONV_FAILURE")        = -4;
m.attr("CV_LINIT_FAIL")          = -5;
m.attr("CV_LSETUP_FAIL")         = -6;
m.attr("CV_LSOLVE_FAIL")         = -7;
m.attr("CV_RHSFUNC_FAIL")        = -8;
m.attr("CV_FIRST_RHSFUNC_ERR")   = -9;
m.attr("CV_REPTD_RHSFUNC_ERR")   = -10;
m.attr("CV_UNREC_RHSFUNC_ERR")   = -11;
m.attr("CV_RTFUNC_FAIL")         = -12;
m.attr("CV_NLS_INIT_FAIL")       = -13;
m.attr("CV_NLS_SETUP_FAIL")      = -14;
m.attr("CV_CONSTR_FAIL")         = -15;
m.attr("CV_NLS_FAIL")            = -16;
m.attr("CV_MEM_FAIL")            = -20;
m.attr("CV_MEM_NULL")            = -21;
m.attr("CV_ILL_INPUT")           = -22;
m.attr("CV_NO_MALLOC")           = -23;
m.attr("CV_BAD_K")               = -24;
m.attr("CV_BAD_T")               = -25;
m.attr("CV_BAD_DKY")             = -26;
m.attr("CV_TOO_CLOSE")           = -27;
m.attr("CV_VECTOROP_ERR")        = -28;
m.attr("CV_NO_QUAD")             = -30;
m.attr("CV_QRHSFUNC_FAIL")       = -31;
m.attr("CV_FIRST_QRHSFUNC_ERR")  = -32;
m.attr("CV_REPTD_QRHSFUNC_ERR")  = -33;
m.attr("CV_UNREC_QRHSFUNC_ERR")  = -34;
m.attr("CV_NO_SENS")             = -40;
m.attr("CV_SRHSFUNC_FAIL")       = -41;
m.attr("CV_FIRST_SRHSFUNC_ERR")  = -42;
m.attr("CV_REPTD_SRHSFUNC_ERR")  = -43;
m.attr("CV_UNREC_SRHSFUNC_ERR")  = -44;
m.attr("CV_BAD_IS")              = -45;
m.attr("CV_NO_QUADSENS")         = -50;
m.attr("CV_QSRHSFUNC_FAIL")      = -51;
m.attr("CV_FIRST_QSRHSFUNC_ERR") = -52;
m.attr("CV_REPTD_QSRHSFUNC_ERR") = -53;
m.attr("CV_UNREC_QSRHSFUNC_ERR") = -54;
m.attr("CV_CONTEXT_ERR")         = -55;
m.attr("CV_PROJ_MEM_NULL")       = -56;
m.attr("CV_PROJFUNC_FAIL")       = -57;
m.attr("CV_REPTD_PROJFUNC_ERR")  = -58;
m.attr("CV_BAD_TINTERP")         = -59;
m.attr("CV_UNRECOGNIZED_ERR")    = -99;
m.attr("CV_NO_ADJ")              = -101;
m.attr("CV_NO_FWD")              = -102;
m.attr("CV_NO_BCK")              = -103;
m.attr("CV_BAD_TB0")             = -104;
m.attr("CV_REIFWD_FAIL")         = -105;
m.attr("CV_FWD_FAIL")            = -106;
m.attr("CV_GETY_BADT")           = -107;

m.def("CVodeReInit", CVodeReInit, nb::arg("cvode_mem"), nb::arg("t0"),
      nb::arg("y0"));

m.def(
  "CVodeResizeHistory",
  [](void* cvode_mem, sundials4py::Array1d t_hist_1d,
     std::vector<N_Vector> y_hist_1d, std::vector<N_Vector> f_hist_1d,
     int num_y_hist, int num_f_hist) -> int
  {
    auto CVodeResizeHistory_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sundials4py::Array1d t_hist_1d,
         std::vector<N_Vector> y_hist_1d, std::vector<N_Vector> f_hist_1d,
         int num_y_hist, int num_f_hist) -> int
    {
      sunrealtype* t_hist_1d_ptr = t_hist_1d.size() == 0 ? nullptr
                                                         : t_hist_1d.data();
      N_Vector* y_hist_1d_ptr = y_hist_1d.empty() ? nullptr : y_hist_1d.data();
      N_Vector* f_hist_1d_ptr = f_hist_1d.empty() ? nullptr : f_hist_1d.data();

      auto lambda_result = CVodeResizeHistory(cvode_mem, t_hist_1d_ptr,
                                              y_hist_1d_ptr, f_hist_1d_ptr,
                                              num_y_hist, num_f_hist);
      return lambda_result;
    };

    return CVodeResizeHistory_adapt_arr_ptr_to_std_vector(cvode_mem, t_hist_1d,
                                                          y_hist_1d, f_hist_1d,
                                                          num_y_hist, num_f_hist);
  },
  nb::arg("cvode_mem"), nb::arg("t_hist_1d"), nb::arg("y_hist_1d"),
  nb::arg("f_hist_1d"), nb::arg("num_y_hist"), nb::arg("num_f_hist"));

m.def("CVodeSStolerances", CVodeSStolerances, nb::arg("cvode_mem"),
      nb::arg("reltol"), nb::arg("abstol"));

m.def("CVodeSVtolerances", CVodeSVtolerances, nb::arg("cvode_mem"),
      nb::arg("reltol"), nb::arg("abstol"));

m.def("CVodeSetConstraints", CVodeSetConstraints, nb::arg("cvode_mem"),
      nb::arg("constraints"));

m.def("CVodeSetMaxNumConstraintFails", CVodeSetMaxNumConstraintFails,
      nb::arg("cvode_mem"), nb::arg("max_fails"));

m.def("CVodeSetDeltaGammaMaxLSetup", CVodeSetDeltaGammaMaxLSetup,
      nb::arg("cvode_mem"), nb::arg("dgmax_lsetup"));

m.def("CVodeSetInitStep", CVodeSetInitStep, nb::arg("cvode_mem"), nb::arg("hin"));

m.def("CVodeSetLSetupFrequency", CVodeSetLSetupFrequency, nb::arg("cvode_mem"),
      nb::arg("msbp"));

m.def("CVodeSetMaxConvFails", CVodeSetMaxConvFails, nb::arg("cvode_mem"),
      nb::arg("maxncf"));

m.def("CVodeSetMaxErrTestFails", CVodeSetMaxErrTestFails, nb::arg("cvode_mem"),
      nb::arg("maxnef"));

m.def("CVodeSetMaxHnilWarns", CVodeSetMaxHnilWarns, nb::arg("cvode_mem"),
      nb::arg("mxhnil"));

m.def("CVodeSetMaxNonlinIters", CVodeSetMaxNonlinIters, nb::arg("cvode_mem"),
      nb::arg("maxcor"));

m.def("CVodeSetMaxNumSteps", CVodeSetMaxNumSteps, nb::arg("cvode_mem"),
      nb::arg("mxsteps"));

m.def("CVodeSetMaxOrd", CVodeSetMaxOrd, nb::arg("cvode_mem"), nb::arg("maxord"));

m.def("CVodeSetMaxStep", CVodeSetMaxStep, nb::arg("cvode_mem"), nb::arg("hmax"));

m.def("CVodeSetMinStep", CVodeSetMinStep, nb::arg("cvode_mem"), nb::arg("hmin"));

m.def("CVodeSetMonitorFrequency", CVodeSetMonitorFrequency,
      nb::arg("cvode_mem"), nb::arg("nst"));

m.def("CVodeSetNonlinConvCoef", CVodeSetNonlinConvCoef, nb::arg("cvode_mem"),
      nb::arg("nlscoef"));

m.def("CVodeSetNonlinearSolver", CVodeSetNonlinearSolver, nb::arg("cvode_mem"),
      nb::arg("NLS"));

m.def("CVodeSetStabLimDet", CVodeSetStabLimDet, nb::arg("cvode_mem"),
      nb::arg("stldet"));

m.def("CVodeSetStopTime", CVodeSetStopTime, nb::arg("cvode_mem"),
      nb::arg("tstop"));

m.def("CVodeSetInterpolateStopTime", CVodeSetInterpolateStopTime,
      nb::arg("cvode_mem"), nb::arg("interp"));

m.def("CVodeClearStopTime", CVodeClearStopTime, nb::arg("cvode_mem"));

m.def("CVodeSetEtaFixedStepBounds", CVodeSetEtaFixedStepBounds,
      nb::arg("cvode_mem"), nb::arg("eta_min_fx"), nb::arg("eta_max_fx"));

m.def("CVodeSetEtaMaxFirstStep", CVodeSetEtaMaxFirstStep, nb::arg("cvode_mem"),
      nb::arg("eta_max_fs"));

m.def("CVodeSetEtaMaxEarlyStep", CVodeSetEtaMaxEarlyStep, nb::arg("cvode_mem"),
      nb::arg("eta_max_es"));

m.def("CVodeSetNumStepsEtaMaxEarlyStep", CVodeSetNumStepsEtaMaxEarlyStep,
      nb::arg("cvode_mem"), nb::arg("small_nst"));

m.def("CVodeSetEtaMax", CVodeSetEtaMax, nb::arg("cvode_mem"),
      nb::arg("eta_max_gs"));

m.def("CVodeSetEtaMin", CVodeSetEtaMin, nb::arg("cvode_mem"), nb::arg("eta_min"));

m.def("CVodeSetEtaMinErrFail", CVodeSetEtaMinErrFail, nb::arg("cvode_mem"),
      nb::arg("eta_min_ef"));

m.def("CVodeSetEtaMaxErrFail", CVodeSetEtaMaxErrFail, nb::arg("cvode_mem"),
      nb::arg("eta_max_ef"));

m.def("CVodeSetNumFailsEtaMaxErrFail", CVodeSetNumFailsEtaMaxErrFail,
      nb::arg("cvode_mem"), nb::arg("small_nef"));

m.def("CVodeSetEtaConvFail", CVodeSetEtaConvFail, nb::arg("cvode_mem"),
      nb::arg("eta_cf"));

m.def(
  "CVodeSetRootDirection",
  [](void* cvode_mem, std::vector<int> rootdir_1d) -> int
  {
    auto CVodeSetRootDirection_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<int> rootdir_1d) -> int
    {
      int* rootdir_1d_ptr = rootdir_1d.empty() ? nullptr : rootdir_1d.data();

      auto lambda_result = CVodeSetRootDirection(cvode_mem, rootdir_1d_ptr);
      return lambda_result;
    };

    return CVodeSetRootDirection_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                             rootdir_1d);
  },
  nb::arg("cvode_mem"), nb::arg("rootdir_1d"));

m.def("CVodeSetNoInactiveRootWarn", CVodeSetNoInactiveRootWarn,
      nb::arg("cvode_mem"));

m.def(
  "CVode",
  [](void* cvode_mem, sunrealtype tout, N_Vector yout,
     int itask) -> std::tuple<int, sunrealtype>
  {
    auto CVode_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, sunrealtype tout, N_Vector yout,
         int itask) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = CVode(cvode_mem, tout, yout, &tret_adapt_modifiable, itask);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return CVode_adapt_modifiable_immutable_to_return(cvode_mem, tout, yout,
                                                      itask);
  },
  nb::arg("cvode_mem"), nb::arg("tout"), nb::arg("yout"), nb::arg("itask"),
  "Solver function");

m.def("CVodeComputeState", CVodeComputeState, nb::arg("cvode_mem"),
      nb::arg("ycor"), nb::arg("y"));

m.def(
  "CVodeComputeStateSens",
  [](void* cvode_mem, std::vector<N_Vector> yScor_1d,
     std::vector<N_Vector> yS_1d) -> int
  {
    auto CVodeComputeStateSens_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<N_Vector> yScor_1d,
         std::vector<N_Vector> yS_1d) -> int
    {
      N_Vector* yScor_1d_ptr = yScor_1d.empty() ? nullptr : yScor_1d.data();
      N_Vector* yS_1d_ptr    = yS_1d.empty() ? nullptr : yS_1d.data();

      auto lambda_result = CVodeComputeStateSens(cvode_mem, yScor_1d_ptr,
                                                 yS_1d_ptr);
      return lambda_result;
    };

    return CVodeComputeStateSens_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                             yScor_1d, yS_1d);
  },
  nb::arg("cvode_mem"), nb::arg("yScor_1d"), nb::arg("yS_1d"));

m.def("CVodeComputeStateSens1", CVodeComputeStateSens1, nb::arg("cvode_mem"),
      nb::arg("idx"), nb::arg("yScor1"), nb::arg("yS1"));

m.def("CVodeGetDky", CVodeGetDky, nb::arg("cvode_mem"), nb::arg("t"),
      nb::arg("k"), nb::arg("dky"), "Dense output function");

m.def(
  "CVodeGetNumSteps",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumSteps_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nsteps_adapt_modifiable;

      int r = CVodeGetNumSteps(cvode_mem, &nsteps_adapt_modifiable);
      return std::make_tuple(r, nsteps_adapt_modifiable);
    };

    return CVodeGetNumSteps_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumRhsEvals",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nfevals_adapt_modifiable;

      int r = CVodeGetNumRhsEvals(cvode_mem, &nfevals_adapt_modifiable);
      return std::make_tuple(r, nfevals_adapt_modifiable);
    };

    return CVodeGetNumRhsEvals_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumLinSolvSetups",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nlinsetups_adapt_modifiable;

      int r = CVodeGetNumLinSolvSetups(cvode_mem, &nlinsetups_adapt_modifiable);
      return std::make_tuple(r, nlinsetups_adapt_modifiable);
    };

    return CVodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumErrTestFails",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long netfails_adapt_modifiable;

      int r = CVodeGetNumErrTestFails(cvode_mem, &netfails_adapt_modifiable);
      return std::make_tuple(r, netfails_adapt_modifiable);
    };

    return CVodeGetNumErrTestFails_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetLastOrder",
  [](void* cvode_mem) -> std::tuple<int, int>
  {
    auto CVodeGetLastOrder_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, int>
    {
      int qlast_adapt_modifiable;

      int r = CVodeGetLastOrder(cvode_mem, &qlast_adapt_modifiable);
      return std::make_tuple(r, qlast_adapt_modifiable);
    };

    return CVodeGetLastOrder_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetCurrentOrder",
  [](void* cvode_mem) -> std::tuple<int, int>
  {
    auto CVodeGetCurrentOrder_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, int>
    {
      int qcur_adapt_modifiable;

      int r = CVodeGetCurrentOrder(cvode_mem, &qcur_adapt_modifiable);
      return std::make_tuple(r, qcur_adapt_modifiable);
    };

    return CVodeGetCurrentOrder_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetCurrentGamma",
  [](void* cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetCurrentGamma_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype gamma_adapt_modifiable;

      int r = CVodeGetCurrentGamma(cvode_mem, &gamma_adapt_modifiable);
      return std::make_tuple(r, gamma_adapt_modifiable);
    };

    return CVodeGetCurrentGamma_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumStabLimOrderReds",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumStabLimOrderReds_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nslred_adapt_modifiable;

      int r = CVodeGetNumStabLimOrderReds(cvode_mem, &nslred_adapt_modifiable);
      return std::make_tuple(r, nslred_adapt_modifiable);
    };

    return CVodeGetNumStabLimOrderReds_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetActualInitStep",
  [](void* cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetActualInitStep_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hinused_adapt_modifiable;

      int r = CVodeGetActualInitStep(cvode_mem, &hinused_adapt_modifiable);
      return std::make_tuple(r, hinused_adapt_modifiable);
    };

    return CVodeGetActualInitStep_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetLastStep",
  [](void* cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetLastStep_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hlast_adapt_modifiable;

      int r = CVodeGetLastStep(cvode_mem, &hlast_adapt_modifiable);
      return std::make_tuple(r, hlast_adapt_modifiable);
    };

    return CVodeGetLastStep_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetCurrentStep",
  [](void* cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetCurrentStep_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hcur_adapt_modifiable;

      int r = CVodeGetCurrentStep(cvode_mem, &hcur_adapt_modifiable);
      return std::make_tuple(r, hcur_adapt_modifiable);
    };

    return CVodeGetCurrentStep_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetCurrentState",
  [](void* cvode_mem) -> std::tuple<int, N_Vector>
  {
    auto CVodeGetCurrentState_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, N_Vector>
    {
      N_Vector y_adapt_modifiable;

      int r = CVodeGetCurrentState(cvode_mem, &y_adapt_modifiable);
      return std::make_tuple(r, y_adapt_modifiable);
    };

    return CVodeGetCurrentState_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

m.def(
  "CVodeGetCurrentSensSolveIndex",
  [](void* cvode_mem) -> std::tuple<int, int>
  {
    auto CVodeGetCurrentSensSolveIndex_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, int>
    {
      int index_adapt_modifiable;

      int r = CVodeGetCurrentSensSolveIndex(cvode_mem, &index_adapt_modifiable);
      return std::make_tuple(r, index_adapt_modifiable);
    };

    return CVodeGetCurrentSensSolveIndex_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetCurrentTime",
  [](void* cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetCurrentTime_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tcur_adapt_modifiable;

      int r = CVodeGetCurrentTime(cvode_mem, &tcur_adapt_modifiable);
      return std::make_tuple(r, tcur_adapt_modifiable);
    };

    return CVodeGetCurrentTime_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetTolScaleFactor",
  [](void* cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetTolScaleFactor_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tolsfac_adapt_modifiable;

      int r = CVodeGetTolScaleFactor(cvode_mem, &tolsfac_adapt_modifiable);
      return std::make_tuple(r, tolsfac_adapt_modifiable);
    };

    return CVodeGetTolScaleFactor_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def("CVodeGetErrWeights", CVodeGetErrWeights, nb::arg("cvode_mem"),
      nb::arg("eweight"));

m.def("CVodeGetEstLocalErrors", CVodeGetEstLocalErrors, nb::arg("cvode_mem"),
      nb::arg("ele"));

m.def(
  "CVodeGetNumGEvals",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumGEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long ngevals_adapt_modifiable;

      int r = CVodeGetNumGEvals(cvode_mem, &ngevals_adapt_modifiable);
      return std::make_tuple(r, ngevals_adapt_modifiable);
    };

    return CVodeGetNumGEvals_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetRootInfo",
  [](void* cvode_mem, std::vector<int> rootsfound_1d) -> int
  {
    auto CVodeGetRootInfo_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<int> rootsfound_1d) -> int
    {
      int* rootsfound_1d_ptr = rootsfound_1d.empty() ? nullptr
                                                     : rootsfound_1d.data();

      auto lambda_result = CVodeGetRootInfo(cvode_mem, rootsfound_1d_ptr);
      return lambda_result;
    };

    return CVodeGetRootInfo_adapt_arr_ptr_to_std_vector(cvode_mem, rootsfound_1d);
  },
  nb::arg("cvode_mem"), nb::arg("rootsfound_1d"));

m.def(
  "CVodeGetIntegratorStats",
  [](void* cvode_mem) -> std::tuple<int, long, long, long, long, int, int,
                                    sunrealtype, sunrealtype, sunrealtype, sunrealtype>
  {
    auto CVodeGetIntegratorStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem)
      -> std::tuple<int, long, long, long, long, int, int, sunrealtype,
                    sunrealtype, sunrealtype, sunrealtype>
    {
      long nsteps_adapt_modifiable;
      long nfevals_adapt_modifiable;
      long nlinsetups_adapt_modifiable;
      long netfails_adapt_modifiable;
      int qlast_adapt_modifiable;
      int qcur_adapt_modifiable;
      sunrealtype hinused_adapt_modifiable;
      sunrealtype hlast_adapt_modifiable;
      sunrealtype hcur_adapt_modifiable;
      sunrealtype tcur_adapt_modifiable;

      int r =
        CVodeGetIntegratorStats(cvode_mem, &nsteps_adapt_modifiable,
                                &nfevals_adapt_modifiable,
                                &nlinsetups_adapt_modifiable,
                                &netfails_adapt_modifiable,
                                &qlast_adapt_modifiable, &qcur_adapt_modifiable,
                                &hinused_adapt_modifiable,
                                &hlast_adapt_modifiable, &hcur_adapt_modifiable,
                                &tcur_adapt_modifiable);
      return std::make_tuple(r, nsteps_adapt_modifiable, nfevals_adapt_modifiable,
                             nlinsetups_adapt_modifiable,
                             netfails_adapt_modifiable, qlast_adapt_modifiable,
                             qcur_adapt_modifiable, hinused_adapt_modifiable,
                             hlast_adapt_modifiable, hcur_adapt_modifiable,
                             tcur_adapt_modifiable);
    };

    return CVodeGetIntegratorStats_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumNonlinSolvIters",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nniters_adapt_modifiable;

      int r = CVodeGetNumNonlinSolvIters(cvode_mem, &nniters_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable);
    };

    return CVodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumNonlinSolvConvFails",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nnfails_adapt_modifiable;

      int r = CVodeGetNumNonlinSolvConvFails(cvode_mem,
                                             &nnfails_adapt_modifiable);
      return std::make_tuple(r, nnfails_adapt_modifiable);
    };

    return CVodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNonlinSolvStats",
  [](void* cvode_mem) -> std::tuple<int, long, long>
  {
    auto CVodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long, long>
    {
      long nniters_adapt_modifiable;
      long nnfails_adapt_modifiable;

      int r = CVodeGetNonlinSolvStats(cvode_mem, &nniters_adapt_modifiable,
                                      &nnfails_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable,
                             nnfails_adapt_modifiable);
    };

    return CVodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumStepSolveFails",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nncfails_adapt_modifiable;

      int r = CVodeGetNumStepSolveFails(cvode_mem, &nncfails_adapt_modifiable);
      return std::make_tuple(r, nncfails_adapt_modifiable);
    };

    return CVodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumConstraintFails",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumConstraintFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long num_fails_out_adapt_modifiable;

      int r = CVodeGetNumConstraintFails(cvode_mem,
                                         &num_fails_out_adapt_modifiable);
      return std::make_tuple(r, num_fails_out_adapt_modifiable);
    };

    return CVodeGetNumConstraintFails_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumConstraintCorrections",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumConstraintCorrections_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long num_corrections_out_adapt_modifiable;

      int r =
        CVodeGetNumConstraintCorrections(cvode_mem,
                                         &num_corrections_out_adapt_modifiable);
      return std::make_tuple(r, num_corrections_out_adapt_modifiable);
    };

    return CVodeGetNumConstraintCorrections_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def("CVodePrintAllStats", CVodePrintAllStats, nb::arg("cvode_mem"),
      nb::arg("outfile"), nb::arg("fmt"));

m.def("CVodeGetReturnFlagName", CVodeGetReturnFlagName, nb::arg("flag"));

m.def("CVodeQuadReInit", CVodeQuadReInit, nb::arg("cvode_mem"), nb::arg("yQ0"));

m.def("CVodeQuadSStolerances", CVodeQuadSStolerances, nb::arg("cvode_mem"),
      nb::arg("reltolQ"), nb::arg("abstolQ"));

m.def("CVodeQuadSVtolerances", CVodeQuadSVtolerances, nb::arg("cvode_mem"),
      nb::arg("reltolQ"), nb::arg("abstolQ"));

m.def("CVodeSetQuadErrCon", CVodeSetQuadErrCon, nb::arg("cvode_mem"),
      nb::arg("errconQ"));

m.def(
  "CVodeGetQuad",
  [](void* cvode_mem, N_Vector yQout) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetQuad_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, N_Vector yQout) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = CVodeGetQuad(cvode_mem, &tret_adapt_modifiable, yQout);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return CVodeGetQuad_adapt_modifiable_immutable_to_return(cvode_mem, yQout);
  },
  nb::arg("cvode_mem"), nb::arg("yQout"));

m.def("CVodeGetQuadDky", CVodeGetQuadDky, nb::arg("cvode_mem"), nb::arg("t"),
      nb::arg("k"), nb::arg("dky"));

m.def(
  "CVodeGetQuadNumRhsEvals",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nfQevals_adapt_modifiable;

      int r = CVodeGetQuadNumRhsEvals(cvode_mem, &nfQevals_adapt_modifiable);
      return std::make_tuple(r, nfQevals_adapt_modifiable);
    };

    return CVodeGetQuadNumRhsEvals_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetQuadNumErrTestFails",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nQetfails_adapt_modifiable;

      int r = CVodeGetQuadNumErrTestFails(cvode_mem, &nQetfails_adapt_modifiable);
      return std::make_tuple(r, nQetfails_adapt_modifiable);
    };

    return CVodeGetQuadNumErrTestFails_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def("CVodeGetQuadErrWeights", CVodeGetQuadErrWeights, nb::arg("cvode_mem"),
      nb::arg("eQweight"));

m.def(
  "CVodeGetQuadStats",
  [](void* cvode_mem) -> std::tuple<int, long, long>
  {
    auto CVodeGetQuadStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long, long>
    {
      long nfQevals_adapt_modifiable;
      long nQetfails_adapt_modifiable;

      int r = CVodeGetQuadStats(cvode_mem, &nfQevals_adapt_modifiable,
                                &nQetfails_adapt_modifiable);
      return std::make_tuple(r, nfQevals_adapt_modifiable,
                             nQetfails_adapt_modifiable);
    };

    return CVodeGetQuadStats_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeSensReInit",
  [](void* cvode_mem, int ism, std::vector<N_Vector> yS0_1d) -> int
  {
    auto CVodeSensReInit_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, int ism, std::vector<N_Vector> yS0_1d) -> int
    {
      N_Vector* yS0_1d_ptr = yS0_1d.empty() ? nullptr : yS0_1d.data();

      auto lambda_result = CVodeSensReInit(cvode_mem, ism, yS0_1d_ptr);
      return lambda_result;
    };

    return CVodeSensReInit_adapt_arr_ptr_to_std_vector(cvode_mem, ism, yS0_1d);
  },
  nb::arg("cvode_mem"), nb::arg("ism"), nb::arg("yS0_1d"));

m.def(
  "CVodeSensSStolerances",
  [](void* cvode_mem, sunrealtype reltolS, sundials4py::Array1d abstolS_1d) -> int
  {
    auto CVodeSensSStolerances_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype reltolS,
         sundials4py::Array1d abstolS_1d) -> int
    {
      sunrealtype* abstolS_1d_ptr = abstolS_1d.size() == 0 ? nullptr
                                                           : abstolS_1d.data();

      auto lambda_result = CVodeSensSStolerances(cvode_mem, reltolS,
                                                 abstolS_1d_ptr);
      return lambda_result;
    };

    return CVodeSensSStolerances_adapt_arr_ptr_to_std_vector(cvode_mem, reltolS,
                                                             abstolS_1d);
  },
  nb::arg("cvode_mem"), nb::arg("reltolS"), nb::arg("abstolS_1d"));

m.def(
  "CVodeSensSVtolerances",
  [](void* cvode_mem, sunrealtype reltolS, std::vector<N_Vector> abstolS_1d) -> int
  {
    auto CVodeSensSVtolerances_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype reltolS,
         std::vector<N_Vector> abstolS_1d) -> int
    {
      N_Vector* abstolS_1d_ptr = abstolS_1d.empty() ? nullptr : abstolS_1d.data();

      auto lambda_result = CVodeSensSVtolerances(cvode_mem, reltolS,
                                                 abstolS_1d_ptr);
      return lambda_result;
    };

    return CVodeSensSVtolerances_adapt_arr_ptr_to_std_vector(cvode_mem, reltolS,
                                                             abstolS_1d);
  },
  nb::arg("cvode_mem"), nb::arg("reltolS"), nb::arg("abstolS_1d"));

m.def("CVodeSensEEtolerances", CVodeSensEEtolerances, nb::arg("cvode_mem"));

m.def("CVodeSetSensDQMethod", CVodeSetSensDQMethod, nb::arg("cvode_mem"),
      nb::arg("DQtype"), nb::arg("DQrhomax"));

m.def("CVodeSetSensErrCon", CVodeSetSensErrCon, nb::arg("cvode_mem"),
      nb::arg("errconS"));

m.def("CVodeSetSensMaxNonlinIters", CVodeSetSensMaxNonlinIters,
      nb::arg("cvode_mem"), nb::arg("maxcorS"));

m.def(
  "CVodeSetSensParams",
  [](void* cvode_mem, sundials4py::Array1d p_1d, sundials4py::Array1d pbar_1d,
     std::vector<int> plist_1d) -> int
  {
    auto CVodeSetSensParams_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sundials4py::Array1d p_1d,
         sundials4py::Array1d pbar_1d, std::vector<int> plist_1d) -> int
    {
      sunrealtype* p_1d_ptr    = p_1d.size() == 0 ? nullptr : p_1d.data();
      sunrealtype* pbar_1d_ptr = pbar_1d.size() == 0 ? nullptr : pbar_1d.data();
      int* plist_1d_ptr        = plist_1d.empty() ? nullptr : plist_1d.data();

      auto lambda_result = CVodeSetSensParams(cvode_mem, p_1d_ptr, pbar_1d_ptr,
                                              plist_1d_ptr);
      return lambda_result;
    };

    return CVodeSetSensParams_adapt_arr_ptr_to_std_vector(cvode_mem, p_1d,
                                                          pbar_1d, plist_1d);
  },
  nb::arg("cvode_mem"), nb::arg("p_1d"), nb::arg("pbar_1d"), nb::arg("plist_1d"));

m.def("CVodeSetNonlinearSolverSensSim", CVodeSetNonlinearSolverSensSim,
      nb::arg("cvode_mem"), nb::arg("NLS"));

m.def("CVodeSetNonlinearSolverSensStg", CVodeSetNonlinearSolverSensStg,
      nb::arg("cvode_mem"), nb::arg("NLS"));

m.def("CVodeSetNonlinearSolverSensStg1", CVodeSetNonlinearSolverSensStg1,
      nb::arg("cvode_mem"), nb::arg("NLS"));

m.def("CVodeSensToggleOff", CVodeSensToggleOff, nb::arg("cvode_mem"));

m.def(
  "CVodeGetSens",
  [](void* cvode_mem,
     std::vector<N_Vector> ySout_1d) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetSens_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype* tret, std::vector<N_Vector> ySout_1d) -> int
    {
      N_Vector* ySout_1d_ptr = ySout_1d.empty() ? nullptr : ySout_1d.data();

      auto lambda_result = CVodeGetSens(cvode_mem, tret, ySout_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetSens_adapt_modifiable_immutable_to_return =
      [&CVodeGetSens_adapt_arr_ptr_to_std_vector](void* cvode_mem,
                                                  std::vector<N_Vector> ySout_1d)
      -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = CVodeGetSens_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                       &tret_adapt_modifiable,
                                                       ySout_1d);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return CVodeGetSens_adapt_modifiable_immutable_to_return(cvode_mem, ySout_1d);
  },
  nb::arg("cvode_mem"), nb::arg("ySout_1d"));

m.def(
  "CVodeGetSens1",
  [](void* cvode_mem, int is, N_Vector ySout) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetSens1_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int is, N_Vector ySout) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = CVodeGetSens1(cvode_mem, &tret_adapt_modifiable, is, ySout);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return CVodeGetSens1_adapt_modifiable_immutable_to_return(cvode_mem, is,
                                                              ySout);
  },
  nb::arg("cvode_mem"), nb::arg("is_"), nb::arg("ySout"));

m.def(
  "CVodeGetSensDky",
  [](void* cvode_mem, sunrealtype t, int k, std::vector<N_Vector> dkyA_1d) -> int
  {
    auto CVodeGetSensDky_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype t, int k,
         std::vector<N_Vector> dkyA_1d) -> int
    {
      N_Vector* dkyA_1d_ptr = dkyA_1d.empty() ? nullptr : dkyA_1d.data();

      auto lambda_result = CVodeGetSensDky(cvode_mem, t, k, dkyA_1d_ptr);
      return lambda_result;
    };

    return CVodeGetSensDky_adapt_arr_ptr_to_std_vector(cvode_mem, t, k, dkyA_1d);
  },
  nb::arg("cvode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dkyA_1d"));

m.def("CVodeGetSensDky1", CVodeGetSensDky1, nb::arg("cvode_mem"), nb::arg("t"),
      nb::arg("k"), nb::arg("is_"), nb::arg("dky"));

m.def(
  "CVodeGetSensNumRhsEvals",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nfSevals_adapt_modifiable;

      int r = CVodeGetSensNumRhsEvals(cvode_mem, &nfSevals_adapt_modifiable);
      return std::make_tuple(r, nfSevals_adapt_modifiable);
    };

    return CVodeGetSensNumRhsEvals_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumRhsEvalsSens",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumRhsEvalsSens_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nfevalsS_adapt_modifiable;

      int r = CVodeGetNumRhsEvalsSens(cvode_mem, &nfevalsS_adapt_modifiable);
      return std::make_tuple(r, nfevalsS_adapt_modifiable);
    };

    return CVodeGetNumRhsEvalsSens_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetSensNumErrTestFails",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nSetfails_adapt_modifiable;

      int r = CVodeGetSensNumErrTestFails(cvode_mem, &nSetfails_adapt_modifiable);
      return std::make_tuple(r, nSetfails_adapt_modifiable);
    };

    return CVodeGetSensNumErrTestFails_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetSensNumLinSolvSetups",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nlinsetupsS_adapt_modifiable;

      int r = CVodeGetSensNumLinSolvSetups(cvode_mem,
                                           &nlinsetupsS_adapt_modifiable);
      return std::make_tuple(r, nlinsetupsS_adapt_modifiable);
    };

    return CVodeGetSensNumLinSolvSetups_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetSensErrWeights",
  [](void* cvode_mem, std::vector<N_Vector> eSweight_1d) -> int
  {
    auto CVodeGetSensErrWeights_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<N_Vector> eSweight_1d) -> int
    {
      N_Vector* eSweight_1d_ptr = eSweight_1d.empty() ? nullptr
                                                      : eSweight_1d.data();

      auto lambda_result = CVodeGetSensErrWeights(cvode_mem, eSweight_1d_ptr);
      return lambda_result;
    };

    return CVodeGetSensErrWeights_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                              eSweight_1d);
  },
  nb::arg("cvode_mem"), nb::arg("eSweight_1d"));

m.def(
  "CVodeGetSensStats",
  [](void* cvode_mem) -> std::tuple<int, long, long, long, long>
  {
    auto CVodeGetSensStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long, long, long, long>
    {
      long nfSevals_adapt_modifiable;
      long nfevalsS_adapt_modifiable;
      long nSetfails_adapt_modifiable;
      long nlinsetupsS_adapt_modifiable;

      int r = CVodeGetSensStats(cvode_mem, &nfSevals_adapt_modifiable,
                                &nfevalsS_adapt_modifiable,
                                &nSetfails_adapt_modifiable,
                                &nlinsetupsS_adapt_modifiable);
      return std::make_tuple(r, nfSevals_adapt_modifiable,
                             nfevalsS_adapt_modifiable,
                             nSetfails_adapt_modifiable,
                             nlinsetupsS_adapt_modifiable);
    };

    return CVodeGetSensStats_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetSensNumNonlinSolvIters",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nSniters_adapt_modifiable;

      int r = CVodeGetSensNumNonlinSolvIters(cvode_mem,
                                             &nSniters_adapt_modifiable);
      return std::make_tuple(r, nSniters_adapt_modifiable);
    };

    return CVodeGetSensNumNonlinSolvIters_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetSensNumNonlinSolvConvFails",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nSnfails_adapt_modifiable;

      int r = CVodeGetSensNumNonlinSolvConvFails(cvode_mem,
                                                 &nSnfails_adapt_modifiable);
      return std::make_tuple(r, nSnfails_adapt_modifiable);
    };

    return CVodeGetSensNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetSensNonlinSolvStats",
  [](void* cvode_mem) -> std::tuple<int, long, long>
  {
    auto CVodeGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long, long>
    {
      long nSniters_adapt_modifiable;
      long nSnfails_adapt_modifiable;

      int r = CVodeGetSensNonlinSolvStats(cvode_mem, &nSniters_adapt_modifiable,
                                          &nSnfails_adapt_modifiable);
      return std::make_tuple(r, nSniters_adapt_modifiable,
                             nSnfails_adapt_modifiable);
    };

    return CVodeGetSensNonlinSolvStats_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumStepSensSolveFails",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nSncfails_adapt_modifiable;

      int r = CVodeGetNumStepSensSolveFails(cvode_mem,
                                            &nSncfails_adapt_modifiable);
      return std::make_tuple(r, nSncfails_adapt_modifiable);
    };

    return CVodeGetNumStepSensSolveFails_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetStgrSensNumNonlinSolvIters",
  [](void* cvode_mem, std::vector<long> nSTGR1niters_1d) -> int
  {
    auto CVodeGetStgrSensNumNonlinSolvIters_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<long> nSTGR1niters_1d) -> int
    {
      long* nSTGR1niters_1d_ptr =
        nSTGR1niters_1d.empty() ? nullptr : nSTGR1niters_1d.data();

      auto lambda_result =
        CVodeGetStgrSensNumNonlinSolvIters(cvode_mem, nSTGR1niters_1d_ptr);
      return lambda_result;
    };

    return CVodeGetStgrSensNumNonlinSolvIters_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                                          nSTGR1niters_1d);
  },
  nb::arg("cvode_mem"), nb::arg("nSTGR1niters_1d"));

m.def(
  "CVodeGetStgrSensNumNonlinSolvConvFails",
  [](void* cvode_mem, std::vector<long> nSTGR1nfails_1d) -> int
  {
    auto CVodeGetStgrSensNumNonlinSolvConvFails_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<long> nSTGR1nfails_1d) -> int
    {
      long* nSTGR1nfails_1d_ptr =
        nSTGR1nfails_1d.empty() ? nullptr : nSTGR1nfails_1d.data();

      auto lambda_result =
        CVodeGetStgrSensNumNonlinSolvConvFails(cvode_mem, nSTGR1nfails_1d_ptr);
      return lambda_result;
    };

    return CVodeGetStgrSensNumNonlinSolvConvFails_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                                              nSTGR1nfails_1d);
  },
  nb::arg("cvode_mem"), nb::arg("nSTGR1nfails_1d"));

m.def(
  "CVodeGetStgrSensNonlinSolvStats",
  [](void* cvode_mem, std::vector<long> nSTGR1niters_1d,
     std::vector<long> nSTGR1nfails_1d) -> int
  {
    auto CVodeGetStgrSensNonlinSolvStats_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<long> nSTGR1niters_1d,
         std::vector<long> nSTGR1nfails_1d) -> int
    {
      long* nSTGR1niters_1d_ptr =
        nSTGR1niters_1d.empty() ? nullptr : nSTGR1niters_1d.data();
      long* nSTGR1nfails_1d_ptr =
        nSTGR1nfails_1d.empty() ? nullptr : nSTGR1nfails_1d.data();

      auto lambda_result = CVodeGetStgrSensNonlinSolvStats(cvode_mem,
                                                           nSTGR1niters_1d_ptr,
                                                           nSTGR1nfails_1d_ptr);
      return lambda_result;
    };

    return CVodeGetStgrSensNonlinSolvStats_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                                       nSTGR1niters_1d,
                                                                       nSTGR1nfails_1d);
  },
  nb::arg("cvode_mem"), nb::arg("nSTGR1niters_1d"), nb::arg("nSTGR1nfails_1d"));

m.def(
  "CVodeGetNumStepStgrSensSolveFails",
  [](void* cvode_mem, std::vector<long> nSTGR1ncfails_1d) -> int
  {
    auto CVodeGetNumStepStgrSensSolveFails_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<long> nSTGR1ncfails_1d) -> int
    {
      long* nSTGR1ncfails_1d_ptr =
        nSTGR1ncfails_1d.empty() ? nullptr : nSTGR1ncfails_1d.data();

      auto lambda_result =
        CVodeGetNumStepStgrSensSolveFails(cvode_mem, nSTGR1ncfails_1d_ptr);
      return lambda_result;
    };

    return CVodeGetNumStepStgrSensSolveFails_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                                         nSTGR1ncfails_1d);
  },
  nb::arg("cvode_mem"), nb::arg("nSTGR1ncfails_1d"));

m.def(
  "CVodeQuadSensReInit",
  [](void* cvode_mem, std::vector<N_Vector> yQS0_1d) -> int
  {
    auto CVodeQuadSensReInit_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<N_Vector> yQS0_1d) -> int
    {
      N_Vector* yQS0_1d_ptr = yQS0_1d.empty() ? nullptr : yQS0_1d.data();

      auto lambda_result = CVodeQuadSensReInit(cvode_mem, yQS0_1d_ptr);
      return lambda_result;
    };

    return CVodeQuadSensReInit_adapt_arr_ptr_to_std_vector(cvode_mem, yQS0_1d);
  },
  nb::arg("cvode_mem"), nb::arg("yQS0_1d"));

m.def(
  "CVodeQuadSensSStolerances",
  [](void* cvode_mem, sunrealtype reltolQS, sundials4py::Array1d abstolQS_1d) -> int
  {
    auto CVodeQuadSensSStolerances_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype reltolQS,
         sundials4py::Array1d abstolQS_1d) -> int
    {
      sunrealtype* abstolQS_1d_ptr =
        abstolQS_1d.size() == 0 ? nullptr : abstolQS_1d.data();

      auto lambda_result = CVodeQuadSensSStolerances(cvode_mem, reltolQS,
                                                     abstolQS_1d_ptr);
      return lambda_result;
    };

    return CVodeQuadSensSStolerances_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                                 reltolQS,
                                                                 abstolQS_1d);
  },
  nb::arg("cvode_mem"), nb::arg("reltolQS"), nb::arg("abstolQS_1d"));

m.def(
  "CVodeQuadSensSVtolerances",
  [](void* cvode_mem, sunrealtype reltolQS, std::vector<N_Vector> abstolQS_1d) -> int
  {
    auto CVodeQuadSensSVtolerances_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype reltolQS,
         std::vector<N_Vector> abstolQS_1d) -> int
    {
      N_Vector* abstolQS_1d_ptr = abstolQS_1d.empty() ? nullptr
                                                      : abstolQS_1d.data();

      auto lambda_result = CVodeQuadSensSVtolerances(cvode_mem, reltolQS,
                                                     abstolQS_1d_ptr);
      return lambda_result;
    };

    return CVodeQuadSensSVtolerances_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                                 reltolQS,
                                                                 abstolQS_1d);
  },
  nb::arg("cvode_mem"), nb::arg("reltolQS"), nb::arg("abstolQS_1d"));

m.def("CVodeQuadSensEEtolerances", CVodeQuadSensEEtolerances,
      nb::arg("cvode_mem"));

m.def("CVodeSetQuadSensErrCon", CVodeSetQuadSensErrCon, nb::arg("cvode_mem"),
      nb::arg("errconQS"));

m.def(
  "CVodeGetQuadSens",
  [](void* cvode_mem,
     std::vector<N_Vector> yQSout_1d) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetQuadSens_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype* tret, std::vector<N_Vector> yQSout_1d) -> int
    {
      N_Vector* yQSout_1d_ptr = yQSout_1d.empty() ? nullptr : yQSout_1d.data();

      auto lambda_result = CVodeGetQuadSens(cvode_mem, tret, yQSout_1d_ptr);
      return lambda_result;
    };
    auto CVodeGetQuadSens_adapt_modifiable_immutable_to_return =
      [&CVodeGetQuadSens_adapt_arr_ptr_to_std_vector](void* cvode_mem,
                                                      std::vector<N_Vector> yQSout_1d)
      -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = CVodeGetQuadSens_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                           &tret_adapt_modifiable,
                                                           yQSout_1d);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return CVodeGetQuadSens_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                 yQSout_1d);
  },
  nb::arg("cvode_mem"), nb::arg("yQSout_1d"));

m.def(
  "CVodeGetQuadSens1",
  [](void* cvode_mem, int is, N_Vector yQSout) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetQuadSens1_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int is, N_Vector yQSout) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = CVodeGetQuadSens1(cvode_mem, &tret_adapt_modifiable, is, yQSout);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return CVodeGetQuadSens1_adapt_modifiable_immutable_to_return(cvode_mem, is,
                                                                  yQSout);
  },
  nb::arg("cvode_mem"), nb::arg("is_"), nb::arg("yQSout"));

m.def(
  "CVodeGetQuadSensDky",
  [](void* cvode_mem, sunrealtype t, int k,
     std::vector<N_Vector> dkyQS_all_1d) -> int
  {
    auto CVodeGetQuadSensDky_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, sunrealtype t, int k,
         std::vector<N_Vector> dkyQS_all_1d) -> int
    {
      N_Vector* dkyQS_all_1d_ptr = dkyQS_all_1d.empty() ? nullptr
                                                        : dkyQS_all_1d.data();

      auto lambda_result = CVodeGetQuadSensDky(cvode_mem, t, k, dkyQS_all_1d_ptr);
      return lambda_result;
    };

    return CVodeGetQuadSensDky_adapt_arr_ptr_to_std_vector(cvode_mem, t, k,
                                                           dkyQS_all_1d);
  },
  nb::arg("cvode_mem"), nb::arg("t"), nb::arg("k"), nb::arg("dkyQS_all_1d"));

m.def("CVodeGetQuadSensDky1", CVodeGetQuadSensDky1, nb::arg("cvode_mem"),
      nb::arg("t"), nb::arg("k"), nb::arg("is_"), nb::arg("dkyQS"));

m.def(
  "CVodeGetQuadSensNumRhsEvals",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nfQSevals_adapt_modifiable;

      int r = CVodeGetQuadSensNumRhsEvals(cvode_mem, &nfQSevals_adapt_modifiable);
      return std::make_tuple(r, nfQSevals_adapt_modifiable);
    };

    return CVodeGetQuadSensNumRhsEvals_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetQuadSensNumErrTestFails",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nQSetfails_adapt_modifiable;

      int r = CVodeGetQuadSensNumErrTestFails(cvode_mem,
                                              &nQSetfails_adapt_modifiable);
      return std::make_tuple(r, nQSetfails_adapt_modifiable);
    };

    return CVodeGetQuadSensNumErrTestFails_adapt_modifiable_immutable_to_return(
      cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetQuadSensErrWeights",
  [](void* cvode_mem, std::vector<N_Vector> eQSweight_1d) -> int
  {
    auto CVodeGetQuadSensErrWeights_adapt_arr_ptr_to_std_vector =
      [](void* cvode_mem, std::vector<N_Vector> eQSweight_1d) -> int
    {
      N_Vector* eQSweight_1d_ptr = eQSweight_1d.empty() ? nullptr
                                                        : eQSweight_1d.data();

      auto lambda_result = CVodeGetQuadSensErrWeights(cvode_mem,
                                                      eQSweight_1d_ptr);
      return lambda_result;
    };

    return CVodeGetQuadSensErrWeights_adapt_arr_ptr_to_std_vector(cvode_mem,
                                                                  eQSweight_1d);
  },
  nb::arg("cvode_mem"), nb::arg("eQSweight_1d"));

m.def(
  "CVodeGetQuadSensStats",
  [](void* cvode_mem) -> std::tuple<int, long, long>
  {
    auto CVodeGetQuadSensStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long, long>
    {
      long nfQSevals_adapt_modifiable;
      long nQSetfails_adapt_modifiable;

      int r = CVodeGetQuadSensStats(cvode_mem, &nfQSevals_adapt_modifiable,
                                    &nQSetfails_adapt_modifiable);
      return std::make_tuple(r, nfQSevals_adapt_modifiable,
                             nQSetfails_adapt_modifiable);
    };

    return CVodeGetQuadSensStats_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def("CVodeAdjInit", CVodeAdjInit, nb::arg("cvode_mem"), nb::arg("steps"),
      nb::arg("interp"));

m.def("CVodeAdjReInit", CVodeAdjReInit, nb::arg("cvode_mem"));

m.def(
  "CVodeCreateB",
  [](void* cvode_mem, int lmmB) -> std::tuple<int, int>
  {
    auto CVodeCreateB_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int lmmB) -> std::tuple<int, int>
    {
      int which_adapt_modifiable;

      int r = CVodeCreateB(cvode_mem, lmmB, &which_adapt_modifiable);
      return std::make_tuple(r, which_adapt_modifiable);
    };

    return CVodeCreateB_adapt_modifiable_immutable_to_return(cvode_mem, lmmB);
  },
  nb::arg("cvode_mem"), nb::arg("lmmB"));

m.def("CVodeReInitB", CVodeReInitB, nb::arg("cvode_mem"), nb::arg("which"),
      nb::arg("tB0"), nb::arg("yB0"));

m.def("CVodeSStolerancesB", CVodeSStolerancesB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("reltolB"), nb::arg("abstolB"));

m.def("CVodeSVtolerancesB", CVodeSVtolerancesB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("reltolB"), nb::arg("abstolB"));

m.def("CVodeQuadReInitB", CVodeQuadReInitB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("yQB0"));

m.def("CVodeQuadSStolerancesB", CVodeQuadSStolerancesB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("reltolQB"), nb::arg("abstolQB"));

m.def("CVodeQuadSVtolerancesB", CVodeQuadSVtolerancesB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("reltolQB"), nb::arg("abstolQB"));

m.def(
  "CVodeF",
  [](void* cvode_mem, sunrealtype tout, N_Vector yout,
     int itask) -> std::tuple<int, sunrealtype, int>
  {
    auto CVodeF_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, sunrealtype tout, N_Vector yout,
         int itask) -> std::tuple<int, sunrealtype, int>
    {
      sunrealtype tret_adapt_modifiable;
      int ncheckPtr_adapt_modifiable;

      int r = CVodeF(cvode_mem, tout, yout, &tret_adapt_modifiable, itask,
                     &ncheckPtr_adapt_modifiable);
      return std::make_tuple(r, tret_adapt_modifiable,
                             ncheckPtr_adapt_modifiable);
    };

    return CVodeF_adapt_modifiable_immutable_to_return(cvode_mem, tout, yout,
                                                       itask);
  },
  nb::arg("cvode_mem"), nb::arg("tout"), nb::arg("yout"), nb::arg("itask"));

m.def("CVodeB", CVodeB, nb::arg("cvode_mem"), nb::arg("tBout"),
      nb::arg("itaskB"));

m.def("CVodeSetAdjNoSensi", CVodeSetAdjNoSensi, nb::arg("cvode_mem"));

m.def("CVodeSetMaxOrdB", CVodeSetMaxOrdB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("maxordB"));

m.def("CVodeSetMaxNumStepsB", CVodeSetMaxNumStepsB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("mxstepsB"));

m.def("CVodeSetStabLimDetB", CVodeSetStabLimDetB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("stldetB"));

m.def("CVodeSetInitStepB", CVodeSetInitStepB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("hinB"));

m.def("CVodeSetMinStepB", CVodeSetMinStepB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("hminB"));

m.def("CVodeSetMaxStepB", CVodeSetMaxStepB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("hmaxB"));

m.def("CVodeSetConstraintsB", CVodeSetConstraintsB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("constraintsB"));

m.def("CVodeSetQuadErrConB", CVodeSetQuadErrConB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("errconQB"));

m.def("CVodeSetNonlinearSolverB", CVodeSetNonlinearSolverB,
      nb::arg("cvode_mem"), nb::arg("which"), nb::arg("NLS"));

m.def(
  "CVodeGetB",
  [](void* cvode_mem, int which, N_Vector yB) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetB_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int which, N_Vector yB) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tBret_adapt_modifiable;

      int r = CVodeGetB(cvode_mem, which, &tBret_adapt_modifiable, yB);
      return std::make_tuple(r, tBret_adapt_modifiable);
    };

    return CVodeGetB_adapt_modifiable_immutable_to_return(cvode_mem, which, yB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("yB"));

m.def(
  "CVodeGetQuadB",
  [](void* cvode_mem, int which, N_Vector qB) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetQuadB_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int which, N_Vector qB) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tBret_adapt_modifiable;

      int r = CVodeGetQuadB(cvode_mem, which, &tBret_adapt_modifiable, qB);
      return std::make_tuple(r, tBret_adapt_modifiable);
    };

    return CVodeGetQuadB_adapt_modifiable_immutable_to_return(cvode_mem, which,
                                                              qB);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("qB"));

m.def("CVodeGetAdjCVodeBmem", CVodeGetAdjCVodeBmem, nb::arg("cvode_mem"),
      nb::arg("which"));

m.def("CVodeGetAdjY", CVodeGetAdjY, nb::arg("cvode_mem"), nb::arg("t"),
      nb::arg("y"));

m.def(
  "CVodeGetAdjDataPointHermite",
  [](void* cvode_mem, int which, std::optional<N_Vector> y = std::nullopt,
     std::optional<N_Vector> yd = std::nullopt) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetAdjDataPointHermite_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int which, N_Vector y,
         N_Vector yd) -> std::tuple<int, sunrealtype>
    {
      sunrealtype t_adapt_modifiable;

      int r = CVodeGetAdjDataPointHermite(cvode_mem, which, &t_adapt_modifiable,
                                          y, yd);
      return std::make_tuple(r, t_adapt_modifiable);
    };
    auto CVodeGetAdjDataPointHermite_adapt_optional_arg_with_default_null =
      [&CVodeGetAdjDataPointHermite_adapt_modifiable_immutable_to_return](void* cvode_mem,
                                                                          int which,
                                                                          std::optional<N_Vector>
                                                                            y =
                                                                              std::nullopt,
                                                                          std::optional<N_Vector>
                                                                            yd =
                                                                              std::nullopt)
      -> std::tuple<int, sunrealtype>
    {
      N_Vector y_adapt_default_null = nullptr;
      if (y.has_value()) y_adapt_default_null = y.value();
      N_Vector yd_adapt_default_null = nullptr;
      if (yd.has_value()) yd_adapt_default_null = yd.value();

      auto lambda_result =
        CVodeGetAdjDataPointHermite_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                         which,
                                                                         y_adapt_default_null,
                                                                         yd_adapt_default_null);
      return lambda_result;
    };

    return CVodeGetAdjDataPointHermite_adapt_optional_arg_with_default_null(cvode_mem,
                                                                            which,
                                                                            y,
                                                                            yd);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("y").none() = nb::none(),
  nb::arg("yd").none() = nb::none());

m.def(
  "CVodeGetAdjDataPointPolynomial",
  [](void* cvode_mem, int which,
     std::optional<N_Vector> y = std::nullopt) -> std::tuple<int, sunrealtype, int>
  {
    auto CVodeGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem, int which,
         N_Vector y) -> std::tuple<int, sunrealtype, int>
    {
      sunrealtype t_adapt_modifiable;
      int order_adapt_modifiable;

      int r = CVodeGetAdjDataPointPolynomial(cvode_mem, which,
                                             &t_adapt_modifiable,
                                             &order_adapt_modifiable, y);
      return std::make_tuple(r, t_adapt_modifiable, order_adapt_modifiable);
    };
    auto CVodeGetAdjDataPointPolynomial_adapt_optional_arg_with_default_null =
      [&CVodeGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return](void* cvode_mem,
                                                                             int which,
                                                                             std::optional<N_Vector>
                                                                               y =
                                                                                 std::nullopt)
      -> std::tuple<int, sunrealtype, int>
    {
      N_Vector y_adapt_default_null = nullptr;
      if (y.has_value()) y_adapt_default_null = y.value();

      auto lambda_result =
        CVodeGetAdjDataPointPolynomial_adapt_modifiable_immutable_to_return(cvode_mem,
                                                                            which,
                                                                            y_adapt_default_null);
      return lambda_result;
    };

    return CVodeGetAdjDataPointPolynomial_adapt_optional_arg_with_default_null(cvode_mem,
                                                                               which,
                                                                               y);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("y").none() = nb::none());
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _CVSLS_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("CVLS_SUCCESS")         = 0;
m.attr("CVLS_MEM_NULL")        = -1;
m.attr("CVLS_LMEM_NULL")       = -2;
m.attr("CVLS_ILL_INPUT")       = -3;
m.attr("CVLS_MEM_FAIL")        = -4;
m.attr("CVLS_PMEM_NULL")       = -5;
m.attr("CVLS_JACFUNC_UNRECVR") = -6;
m.attr("CVLS_JACFUNC_RECVR")   = -7;
m.attr("CVLS_SUNMAT_FAIL")     = -8;
m.attr("CVLS_SUNLS_FAIL")      = -9;
m.attr("CVLS_NO_ADJ")          = -101;
m.attr("CVLS_LMEMB_NULL")      = -102;

m.def(
  "CVodeSetLinearSolver",
  [](void* cvode_mem, SUNLinearSolver LS,
     std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto CVodeSetLinearSolver_adapt_optional_arg_with_default_null =
      [](void* cvode_mem, SUNLinearSolver LS,
         std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = CVodeSetLinearSolver(cvode_mem, LS,
                                                A_adapt_default_null);
      return lambda_result;
    };

    return CVodeSetLinearSolver_adapt_optional_arg_with_default_null(cvode_mem,
                                                                     LS, A);
  },
  nb::arg("cvode_mem"), nb::arg("LS"), nb::arg("A").none() = nb::none());

m.def("CVodeSetJacEvalFrequency", CVodeSetJacEvalFrequency,
      nb::arg("cvode_mem"), nb::arg("msbj"));

m.def("CVodeSetLinearSolutionScaling", CVodeSetLinearSolutionScaling,
      nb::arg("cvode_mem"), nb::arg("onoff"));

m.def("CVodeSetDeltaGammaMaxBadJac", CVodeSetDeltaGammaMaxBadJac,
      nb::arg("cvode_mem"), nb::arg("dgmax_jbad"));

m.def("CVodeSetEpsLin", CVodeSetEpsLin, nb::arg("cvode_mem"), nb::arg("eplifac"));

m.def("CVodeSetLSNormFactor", CVodeSetLSNormFactor, nb::arg("arkode_mem"),
      nb::arg("nrmfac"));

m.def(
  "CVodeGetJac",
  [](void* cvode_mem) -> std::tuple<int, SUNMatrix>
  {
    auto CVodeGetJac_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, SUNMatrix>
    {
      SUNMatrix J_adapt_modifiable;

      int r = CVodeGetJac(cvode_mem, &J_adapt_modifiable);
      return std::make_tuple(r, J_adapt_modifiable);
    };

    return CVodeGetJac_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

m.def(
  "CVodeGetJacTime",
  [](void* cvode_mem) -> std::tuple<int, sunrealtype>
  {
    auto CVodeGetJacTime_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype t_J_adapt_modifiable;

      int r = CVodeGetJacTime(cvode_mem, &t_J_adapt_modifiable);
      return std::make_tuple(r, t_J_adapt_modifiable);
    };

    return CVodeGetJacTime_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetJacNumSteps",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetJacNumSteps_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nst_J_adapt_modifiable;

      int r = CVodeGetJacNumSteps(cvode_mem, &nst_J_adapt_modifiable);
      return std::make_tuple(r, nst_J_adapt_modifiable);
    };

    return CVodeGetJacNumSteps_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumJacEvals",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumJacEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long njevals_adapt_modifiable;

      int r = CVodeGetNumJacEvals(cvode_mem, &njevals_adapt_modifiable);
      return std::make_tuple(r, njevals_adapt_modifiable);
    };

    return CVodeGetNumJacEvals_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumPrecEvals",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumPrecEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long npevals_adapt_modifiable;

      int r = CVodeGetNumPrecEvals(cvode_mem, &npevals_adapt_modifiable);
      return std::make_tuple(r, npevals_adapt_modifiable);
    };

    return CVodeGetNumPrecEvals_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumPrecSolves",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumPrecSolves_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long npsolves_adapt_modifiable;

      int r = CVodeGetNumPrecSolves(cvode_mem, &npsolves_adapt_modifiable);
      return std::make_tuple(r, npsolves_adapt_modifiable);
    };

    return CVodeGetNumPrecSolves_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumLinIters",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumLinIters_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nliters_adapt_modifiable;

      int r = CVodeGetNumLinIters(cvode_mem, &nliters_adapt_modifiable);
      return std::make_tuple(r, nliters_adapt_modifiable);
    };

    return CVodeGetNumLinIters_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumLinConvFails",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumLinConvFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nlcfails_adapt_modifiable;

      int r = CVodeGetNumLinConvFails(cvode_mem, &nlcfails_adapt_modifiable);
      return std::make_tuple(r, nlcfails_adapt_modifiable);
    };

    return CVodeGetNumLinConvFails_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumJTSetupEvals",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long njtsetups_adapt_modifiable;

      int r = CVodeGetNumJTSetupEvals(cvode_mem, &njtsetups_adapt_modifiable);
      return std::make_tuple(r, njtsetups_adapt_modifiable);
    };

    return CVodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumJtimesEvals",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long njvevals_adapt_modifiable;

      int r = CVodeGetNumJtimesEvals(cvode_mem, &njvevals_adapt_modifiable);
      return std::make_tuple(r, njvevals_adapt_modifiable);
    };

    return CVodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumLinRhsEvals",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nfevalsLS_adapt_modifiable;

      int r = CVodeGetNumLinRhsEvals(cvode_mem, &nfevalsLS_adapt_modifiable);
      return std::make_tuple(r, nfevalsLS_adapt_modifiable);
    };

    return CVodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetLinSolveStats",
  [](void* cvode_mem)
    -> std::tuple<int, long, long, long, long, long, long, long, long>
  {
    auto CVodeGetLinSolveStats_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem)
      -> std::tuple<int, long, long, long, long, long, long, long, long>
    {
      long njevals_adapt_modifiable;
      long nfevalsLS_adapt_modifiable;
      long nliters_adapt_modifiable;
      long nlcfails_adapt_modifiable;
      long npevals_adapt_modifiable;
      long npsolves_adapt_modifiable;
      long njtsetups_adapt_modifiable;
      long njtimes_adapt_modifiable;

      int r = CVodeGetLinSolveStats(cvode_mem, &njevals_adapt_modifiable,
                                    &nfevalsLS_adapt_modifiable,
                                    &nliters_adapt_modifiable,
                                    &nlcfails_adapt_modifiable,
                                    &npevals_adapt_modifiable,
                                    &npsolves_adapt_modifiable,
                                    &njtsetups_adapt_modifiable,
                                    &njtimes_adapt_modifiable);
      return std::make_tuple(r, njevals_adapt_modifiable,
                             nfevalsLS_adapt_modifiable,
                             nliters_adapt_modifiable, nlcfails_adapt_modifiable,
                             npevals_adapt_modifiable, npsolves_adapt_modifiable,
                             njtsetups_adapt_modifiable,
                             njtimes_adapt_modifiable);
    };

    return CVodeGetLinSolveStats_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetLastLinFlag",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetLastLinFlag_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long flag_adapt_modifiable;

      int r = CVodeGetLastLinFlag(cvode_mem, &flag_adapt_modifiable);
      return std::make_tuple(r, flag_adapt_modifiable);
    };

    return CVodeGetLastLinFlag_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def("CVodeGetLinReturnFlagName", CVodeGetLinReturnFlagName, nb::arg("flag"));

m.def(
  "CVodeSetLinearSolverB",
  [](void* cvode_mem, int which, SUNLinearSolver LS,
     std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto CVodeSetLinearSolverB_adapt_optional_arg_with_default_null =
      [](void* cvode_mem, int which, SUNLinearSolver LS,
         std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = CVodeSetLinearSolverB(cvode_mem, which, LS,
                                                 A_adapt_default_null);
      return lambda_result;
    };

    return CVodeSetLinearSolverB_adapt_optional_arg_with_default_null(cvode_mem,
                                                                      which, LS,
                                                                      A);
  },
  nb::arg("cvode_mem"), nb::arg("which"), nb::arg("LS"),
  nb::arg("A").none() = nb::none());

m.def("CVodeSetEpsLinB", CVodeSetEpsLinB, nb::arg("cvode_mem"),
      nb::arg("which"), nb::arg("eplifacB"));

m.def("CVodeSetLSNormFactorB", CVodeSetLSNormFactorB, nb::arg("arkode_mem"),
      nb::arg("which"), nb::arg("nrmfacB"));

m.def("CVodeSetLinearSolutionScalingB", CVodeSetLinearSolutionScalingB,
      nb::arg("cvode_mem"), nb::arg("which"), nb::arg("onoffB"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _CVPROJ_H
//
// #ifdef __cplusplus
// #endif
//

m.def("CVodeSetProjErrEst", CVodeSetProjErrEst, nb::arg("cvode_mem"),
      nb::arg("onoff"));

m.def("CVodeSetProjFrequency", CVodeSetProjFrequency, nb::arg("cvode_mem"),
      nb::arg("proj_freq"));

m.def("CVodeSetMaxNumProjFails", CVodeSetMaxNumProjFails, nb::arg("cvode_mem"),
      nb::arg("max_fails"));

m.def("CVodeSetEpsProj", CVodeSetEpsProj, nb::arg("cvode_mem"), nb::arg("eps"));

m.def("CVodeSetProjFailEta", CVodeSetProjFailEta, nb::arg("cvode_mem"),
      nb::arg("eta"));

m.def(
  "CVodeGetNumProjEvals",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumProjEvals_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nproj_adapt_modifiable;

      int r = CVodeGetNumProjEvals(cvode_mem, &nproj_adapt_modifiable);
      return std::make_tuple(r, nproj_adapt_modifiable);
    };

    return CVodeGetNumProjEvals_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));

m.def(
  "CVodeGetNumProjFails",
  [](void* cvode_mem) -> std::tuple<int, long>
  {
    auto CVodeGetNumProjFails_adapt_modifiable_immutable_to_return =
      [](void* cvode_mem) -> std::tuple<int, long>
    {
      long nprf_adapt_modifiable;

      int r = CVodeGetNumProjFails(cvode_mem, &nprf_adapt_modifiable);
      return std::make_tuple(r, nprf_adapt_modifiable);
    };

    return CVodeGetNumProjFails_adapt_modifiable_immutable_to_return(cvode_mem);
  },
  nb::arg("cvode_mem"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
