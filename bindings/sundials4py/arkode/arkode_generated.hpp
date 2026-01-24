// #ifndef _ARKODE_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("ARK_NORMAL")                 = 1;
m.attr("ARK_ONE_STEP")               = 2;
m.attr("ARK_ADAPT_CUSTOM")           = -1;
m.attr("ARK_ADAPT_PID")              = 0;
m.attr("ARK_ADAPT_PI")               = 1;
m.attr("ARK_ADAPT_I")                = 2;
m.attr("ARK_ADAPT_EXP_GUS")          = 3;
m.attr("ARK_ADAPT_IMP_GUS")          = 4;
m.attr("ARK_ADAPT_IMEX_GUS")         = 5;
m.attr("ARK_FULLRHS_START")          = 0;
m.attr("ARK_FULLRHS_END")            = 1;
m.attr("ARK_FULLRHS_OTHER")          = 2;
m.attr("ARK_INTERP_MAX_DEGREE")      = 5;
m.attr("ARK_INTERP_NONE")            = -1;
m.attr("ARK_INTERP_HERMITE")         = 0;
m.attr("ARK_INTERP_LAGRANGE")        = 1;
m.attr("ARK_SUCCESS")                = 0;
m.attr("ARK_TSTOP_RETURN")           = 1;
m.attr("ARK_ROOT_RETURN")            = 2;
m.attr("ARK_WARNING")                = 99;
m.attr("ARK_TOO_MUCH_WORK")          = -1;
m.attr("ARK_TOO_MUCH_ACC")           = -2;
m.attr("ARK_ERR_FAILURE")            = -3;
m.attr("ARK_CONV_FAILURE")           = -4;
m.attr("ARK_LINIT_FAIL")             = -5;
m.attr("ARK_LSETUP_FAIL")            = -6;
m.attr("ARK_LSOLVE_FAIL")            = -7;
m.attr("ARK_RHSFUNC_FAIL")           = -8;
m.attr("ARK_FIRST_RHSFUNC_ERR")      = -9;
m.attr("ARK_REPTD_RHSFUNC_ERR")      = -10;
m.attr("ARK_UNREC_RHSFUNC_ERR")      = -11;
m.attr("ARK_RTFUNC_FAIL")            = -12;
m.attr("ARK_LFREE_FAIL")             = -13;
m.attr("ARK_MASSINIT_FAIL")          = -14;
m.attr("ARK_MASSSETUP_FAIL")         = -15;
m.attr("ARK_MASSSOLVE_FAIL")         = -16;
m.attr("ARK_MASSFREE_FAIL")          = -17;
m.attr("ARK_MASSMULT_FAIL")          = -18;
m.attr("ARK_CONSTR_FAIL")            = -19;
m.attr("ARK_MEM_FAIL")               = -20;
m.attr("ARK_MEM_NULL")               = -21;
m.attr("ARK_ILL_INPUT")              = -22;
m.attr("ARK_NO_MALLOC")              = -23;
m.attr("ARK_BAD_K")                  = -24;
m.attr("ARK_BAD_T")                  = -25;
m.attr("ARK_BAD_DKY")                = -26;
m.attr("ARK_TOO_CLOSE")              = -27;
m.attr("ARK_VECTOROP_ERR")           = -28;
m.attr("ARK_NLS_INIT_FAIL")          = -29;
m.attr("ARK_NLS_SETUP_FAIL")         = -30;
m.attr("ARK_NLS_SETUP_RECVR")        = -31;
m.attr("ARK_NLS_OP_ERR")             = -32;
m.attr("ARK_INNERSTEP_ATTACH_ERR")   = -33;
m.attr("ARK_INNERSTEP_FAIL")         = -34;
m.attr("ARK_OUTERTOINNER_FAIL")      = -35;
m.attr("ARK_INNERTOOUTER_FAIL")      = -36;
m.attr("ARK_POSTPROCESS_FAIL")       = -37;
m.attr("ARK_POSTPROCESS_STEP_FAIL")  = -37;
m.attr("ARK_POSTPROCESS_STAGE_FAIL") = -38;
m.attr("ARK_USER_PREDICT_FAIL")      = -39;
m.attr("ARK_INTERP_FAIL")            = -40;
m.attr("ARK_INVALID_TABLE")          = -41;
m.attr("ARK_CONTEXT_ERR")            = -42;
m.attr("ARK_RELAX_FAIL")             = -43;
m.attr("ARK_RELAX_MEM_NULL")         = -44;
m.attr("ARK_RELAX_FUNC_FAIL")        = -45;
m.attr("ARK_RELAX_JAC_FAIL")         = -46;
m.attr("ARK_CONTROLLER_ERR")         = -47;
m.attr("ARK_STEPPER_UNSUPPORTED")    = -48;
m.attr("ARK_DOMEIG_FAIL")            = -49;
m.attr("ARK_MAX_STAGE_LIMIT_FAIL")   = -50;
m.attr("ARK_SUNSTEPPER_ERR")         = -51;
m.attr("ARK_STEP_DIRECTION_ERR")     = -52;
m.attr("ARK_ADJ_CHECKPOINT_FAIL")    = -53;
m.attr("ARK_ADJ_RECOMPUTE_FAIL")     = -54;
m.attr("ARK_SUNADJSTEPPER_ERR")      = -55;
m.attr("ARK_DEE_FAIL")               = -56;
m.attr("ARK_UNRECOGNIZED_ERROR")     = -99;

auto pyEnumARKRelaxSolver = nb::enum_<ARKRelaxSolver>(m, "ARKRelaxSolver",
                                                      nb::is_arithmetic(), "")
                              .value("ARK_RELAX_BRENT", ARK_RELAX_BRENT, "")
                              .value("ARK_RELAX_NEWTON", ARK_RELAX_NEWTON, "")
                              .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyEnumARKAccumError =
  nb::enum_<ARKAccumError>(m, "ARKAccumError", nb::is_arithmetic(), "")
    .value("ARK_ACCUMERROR_NONE", ARK_ACCUMERROR_NONE, "")
    .value("ARK_ACCUMERROR_MAX", ARK_ACCUMERROR_MAX, "")
    .value("ARK_ACCUMERROR_SUM", ARK_ACCUMERROR_SUM, "")
    .value("ARK_ACCUMERROR_AVG", ARK_ACCUMERROR_AVG, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

m.def("ARKodeReset", ARKodeReset, nb::arg("arkode_mem"), nb::arg("tR"),
      nb::arg("yR"));

m.def(
  "ARKodeCreateMRIStepInnerStepper",
  [](void* arkode_mem)
    -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<MRIStepInnerStepper>>>
  {
    auto ARKodeCreateMRIStepInnerStepper_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, MRIStepInnerStepper>
    {
      MRIStepInnerStepper stepper_adapt_modifiable;

      int r = ARKodeCreateMRIStepInnerStepper(arkode_mem,
                                              &stepper_adapt_modifiable);
      return std::make_tuple(r, stepper_adapt_modifiable);
    };
    auto ARKodeCreateMRIStepInnerStepper_adapt_return_type_to_shared_ptr =
      [&ARKodeCreateMRIStepInnerStepper_adapt_modifiable_immutable_to_return](
        void* arkode_mem)
      -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<MRIStepInnerStepper>>>
    {
      auto lambda_result =
        ARKodeCreateMRIStepInnerStepper_adapt_modifiable_immutable_to_return(
          arkode_mem);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<std::remove_pointer_t<MRIStepInnerStepper>,
                                             MRIStepInnerStepperDeleter>(
                               std::get<1>(lambda_result)));
    };

    return ARKodeCreateMRIStepInnerStepper_adapt_return_type_to_shared_ptr(
      arkode_mem);
  },
  nb::arg("arkode_mem"), "Utility to wrap ARKODE as an MRIStepInnerStepper",
  nb::rv_policy::reference);

m.def("ARKodeSStolerances", ARKodeSStolerances, nb::arg("arkode_mem"),
      nb::arg("reltol"), nb::arg("abstol"));

m.def("ARKodeSVtolerances", ARKodeSVtolerances, nb::arg("arkode_mem"),
      nb::arg("reltol"), nb::arg("abstol"));

m.def("ARKodeResStolerance", ARKodeResStolerance, nb::arg("arkode_mem"),
      nb::arg("rabstol"));

m.def("ARKodeResVtolerance", ARKodeResVtolerance, nb::arg("arkode_mem"),
      nb::arg("rabstol"));

m.def(
  "ARKodeSetRootDirection",
  [](void* arkode_mem, std::vector<int> rootdir_1d) -> int
  {
    auto ARKodeSetRootDirection_adapt_arr_ptr_to_std_vector =
      [](void* arkode_mem, std::vector<int> rootdir_1d) -> int
    {
      int* rootdir_1d_ptr = rootdir_1d.empty() ? nullptr : rootdir_1d.data();

      auto lambda_result = ARKodeSetRootDirection(arkode_mem, rootdir_1d_ptr);
      return lambda_result;
    };

    return ARKodeSetRootDirection_adapt_arr_ptr_to_std_vector(arkode_mem,
                                                              rootdir_1d);
  },
  nb::arg("arkode_mem"), nb::arg("rootdir_1d"));

m.def("ARKodeSetNoInactiveRootWarn", ARKodeSetNoInactiveRootWarn,
      nb::arg("arkode_mem"));

m.def("ARKodeSetDefaults", ARKodeSetDefaults, nb::arg("arkode_mem"));

m.def("ARKodeSetOrder", ARKodeSetOrder, nb::arg("arkode_mem"), nb::arg("maxord"));

m.def("ARKodeSetInterpolantType", ARKodeSetInterpolantType,
      nb::arg("arkode_mem"), nb::arg("itype"));

m.def("ARKodeSetInterpolantDegree", ARKodeSetInterpolantDegree,
      nb::arg("arkode_mem"), nb::arg("degree"));

m.def("ARKodeSetMaxNumSteps", ARKodeSetMaxNumSteps, nb::arg("arkode_mem"),
      nb::arg("mxsteps"));

m.def("ARKodeSetInterpolateStopTime", ARKodeSetInterpolateStopTime,
      nb::arg("arkode_mem"), nb::arg("interp"));

m.def("ARKodeSetStopTime", ARKodeSetStopTime, nb::arg("arkode_mem"),
      nb::arg("tstop"));

m.def("ARKodeClearStopTime", ARKodeClearStopTime, nb::arg("arkode_mem"));

m.def("ARKodeSetFixedStep", ARKodeSetFixedStep, nb::arg("arkode_mem"),
      nb::arg("hfixed"));

m.def("ARKodeSetStepDirection", ARKodeSetStepDirection, nb::arg("arkode_mem"),
      nb::arg("stepdir"));

m.def("ARKodeSetNonlinearSolver", ARKodeSetNonlinearSolver,
      nb::arg("arkode_mem"), nb::arg("NLS"));

m.def("ARKodeSetLinear", ARKodeSetLinear, nb::arg("arkode_mem"),
      nb::arg("timedepend"));

m.def("ARKodeSetNonlinear", ARKodeSetNonlinear, nb::arg("arkode_mem"));

m.def("ARKodeSetAutonomous", ARKodeSetAutonomous, nb::arg("arkode_mem"),
      nb::arg("autonomous"));

m.def("ARKodeSetDeduceImplicitRhs", ARKodeSetDeduceImplicitRhs,
      nb::arg("arkode_mem"), nb::arg("deduce"));

m.def("ARKodeSetNonlinCRDown", ARKodeSetNonlinCRDown, nb::arg("arkode_mem"),
      nb::arg("crdown"));

m.def("ARKodeSetNonlinRDiv", ARKodeSetNonlinRDiv, nb::arg("arkode_mem"),
      nb::arg("rdiv"));

m.def("ARKodeSetDeltaGammaMax", ARKodeSetDeltaGammaMax, nb::arg("arkode_mem"),
      nb::arg("dgmax"));

m.def("ARKodeSetLSetupFrequency", ARKodeSetLSetupFrequency,
      nb::arg("arkode_mem"), nb::arg("msbp"));

m.def("ARKodeSetPredictorMethod", ARKodeSetPredictorMethod,
      nb::arg("arkode_mem"), nb::arg("method"));

m.def("ARKodeSetMaxNonlinIters", ARKodeSetMaxNonlinIters, nb::arg("arkode_mem"),
      nb::arg("maxcor"));

m.def("ARKodeSetMaxConvFails", ARKodeSetMaxConvFails, nb::arg("arkode_mem"),
      nb::arg("maxncf"));

m.def("ARKodeSetNonlinConvCoef", ARKodeSetNonlinConvCoef, nb::arg("arkode_mem"),
      nb::arg("nlscoef"));

m.def("ARKodeSetAdaptController", ARKodeSetAdaptController,
      nb::arg("arkode_mem"), nb::arg("C"));

m.def("ARKodeSetAdaptControllerByName", ARKodeSetAdaptControllerByName,
      nb::arg("arkode_mem"), nb::arg("cname"));

m.def("ARKodeSetAdaptivityAdjustment", ARKodeSetAdaptivityAdjustment,
      nb::arg("arkode_mem"), nb::arg("adjust"));

m.def("ARKodeSetCFLFraction", ARKodeSetCFLFraction, nb::arg("arkode_mem"),
      nb::arg("cfl_frac"));

m.def("ARKodeSetErrorBias", ARKodeSetErrorBias, nb::arg("arkode_mem"),
      nb::arg("bias"));

m.def("ARKodeSetSafetyFactor", ARKodeSetSafetyFactor, nb::arg("arkode_mem"),
      nb::arg("safety"));

m.def("ARKodeSetMaxGrowth", ARKodeSetMaxGrowth, nb::arg("arkode_mem"),
      nb::arg("mx_growth"));

m.def("ARKodeSetMinReduction", ARKodeSetMinReduction, nb::arg("arkode_mem"),
      nb::arg("eta_min"));

m.def("ARKodeSetFixedStepBounds", ARKodeSetFixedStepBounds,
      nb::arg("arkode_mem"), nb::arg("lb"), nb::arg("ub"));

m.def("ARKodeSetMaxFirstGrowth", ARKodeSetMaxFirstGrowth, nb::arg("arkode_mem"),
      nb::arg("etamx1"));

m.def("ARKodeSetMaxEFailGrowth", ARKodeSetMaxEFailGrowth, nb::arg("arkode_mem"),
      nb::arg("etamxf"));

m.def("ARKodeSetSmallNumEFails", ARKodeSetSmallNumEFails, nb::arg("arkode_mem"),
      nb::arg("small_nef"));

m.def("ARKodeSetMaxCFailGrowth", ARKodeSetMaxCFailGrowth, nb::arg("arkode_mem"),
      nb::arg("etacf"));

m.def("ARKodeSetMaxErrTestFails", ARKodeSetMaxErrTestFails,
      nb::arg("arkode_mem"), nb::arg("maxnef"));

m.def("ARKodeSetConstraints", ARKodeSetConstraints, nb::arg("arkode_mem"),
      nb::arg("constraints"));

m.def("ARKodeSetMaxHnilWarns", ARKodeSetMaxHnilWarns, nb::arg("arkode_mem"),
      nb::arg("mxhnil"));

m.def("ARKodeSetInitStep", ARKodeSetInitStep, nb::arg("arkode_mem"),
      nb::arg("hin"));

m.def("ARKodeSetMinStep", ARKodeSetMinStep, nb::arg("arkode_mem"),
      nb::arg("hmin"));

m.def("ARKodeSetMaxStep", ARKodeSetMaxStep, nb::arg("arkode_mem"),
      nb::arg("hmax"));

m.def("ARKodeSetMaxNumConstrFails", ARKodeSetMaxNumConstrFails,
      nb::arg("arkode_mem"), nb::arg("maxfails"));

m.def("ARKodeSetAdjointCheckpointScheme", ARKodeSetAdjointCheckpointScheme,
      nb::arg("arkode_mem"), nb::arg("checkpoint_scheme"));

m.def("ARKodeSetAdjointCheckpointIndex", ARKodeSetAdjointCheckpointIndex,
      nb::arg("arkode_mem"), nb::arg("step_index"));

m.def("ARKodeSetUseCompensatedSums", ARKodeSetUseCompensatedSums,
      nb::arg("arkode_mem"), nb::arg("onoff"));

m.def("ARKodeSetAccumulatedErrorType", ARKodeSetAccumulatedErrorType,
      nb::arg("arkode_mem"), nb::arg("accum_type"));

m.def("ARKodeResetAccumulatedError", ARKodeResetAccumulatedError,
      nb::arg("arkode_mem"));

m.def(
  "ARKodeEvolve",
  [](void* arkode_mem, sunrealtype tout, N_Vector yout,
     int itask) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeEvolve_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem, sunrealtype tout, N_Vector yout,
         int itask) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tret_adapt_modifiable;

      int r = ARKodeEvolve(arkode_mem, tout, yout, &tret_adapt_modifiable, itask);
      return std::make_tuple(r, tret_adapt_modifiable);
    };

    return ARKodeEvolve_adapt_modifiable_immutable_to_return(arkode_mem, tout,
                                                             yout, itask);
  },
  nb::arg("arkode_mem"), nb::arg("tout"), nb::arg("yout"), nb::arg("itask"),
  "Integrate the ODE over an interval in t");

m.def("ARKodeGetDky", ARKodeGetDky, nb::arg("arkode_mem"), nb::arg("t"),
      nb::arg("k"), nb::arg("dky"),
      "Computes the kth derivative of the y function at time t");

m.def("ARKodeComputeState", ARKodeComputeState, nb::arg("arkode_mem"),
      nb::arg("zcor"), nb::arg("z"),
      "Utility function to update/compute y based on zcor");

m.def(
  "ARKodeGetNumRhsEvals",
  [](void* arkode_mem, int partition_index) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem, int partition_index) -> std::tuple<int, long>
    {
      long num_rhs_evals_adapt_modifiable;

      int r = ARKodeGetNumRhsEvals(arkode_mem, partition_index,
                                   &num_rhs_evals_adapt_modifiable);
      return std::make_tuple(r, num_rhs_evals_adapt_modifiable);
    };

    return ARKodeGetNumRhsEvals_adapt_modifiable_immutable_to_return(arkode_mem,
                                                                     partition_index);
  },
  nb::arg("arkode_mem"), nb::arg("partition_index"));

m.def(
  "ARKodeGetNumStepAttempts",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumStepAttempts_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long step_attempts_adapt_modifiable;

      int r = ARKodeGetNumStepAttempts(arkode_mem,
                                       &step_attempts_adapt_modifiable);
      return std::make_tuple(r, step_attempts_adapt_modifiable);
    };

    return ARKodeGetNumStepAttempts_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumSteps",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumSteps_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nsteps_adapt_modifiable;

      int r = ARKodeGetNumSteps(arkode_mem, &nsteps_adapt_modifiable);
      return std::make_tuple(r, nsteps_adapt_modifiable);
    };

    return ARKodeGetNumSteps_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetLastStep",
  [](void* arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetLastStep_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hlast_adapt_modifiable;

      int r = ARKodeGetLastStep(arkode_mem, &hlast_adapt_modifiable);
      return std::make_tuple(r, hlast_adapt_modifiable);
    };

    return ARKodeGetLastStep_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetCurrentStep",
  [](void* arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetCurrentStep_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hcur_adapt_modifiable;

      int r = ARKodeGetCurrentStep(arkode_mem, &hcur_adapt_modifiable);
      return std::make_tuple(r, hcur_adapt_modifiable);
    };

    return ARKodeGetCurrentStep_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetStepDirection",
  [](void* arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetStepDirection_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype stepdir_adapt_modifiable;

      int r = ARKodeGetStepDirection(arkode_mem, &stepdir_adapt_modifiable);
      return std::make_tuple(r, stepdir_adapt_modifiable);
    };

    return ARKodeGetStepDirection_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def("ARKodeGetErrWeights", ARKodeGetErrWeights, nb::arg("arkode_mem"),
      nb::arg("eweight"));

m.def(
  "ARKodeGetNumGEvals",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumGEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long ngevals_adapt_modifiable;

      int r = ARKodeGetNumGEvals(arkode_mem, &ngevals_adapt_modifiable);
      return std::make_tuple(r, ngevals_adapt_modifiable);
    };

    return ARKodeGetNumGEvals_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetRootInfo",
  [](void* arkode_mem, std::vector<int> rootsfound_1d) -> int
  {
    auto ARKodeGetRootInfo_adapt_arr_ptr_to_std_vector =
      [](void* arkode_mem, std::vector<int> rootsfound_1d) -> int
    {
      int* rootsfound_1d_ptr = rootsfound_1d.empty() ? nullptr
                                                     : rootsfound_1d.data();

      auto lambda_result = ARKodeGetRootInfo(arkode_mem, rootsfound_1d_ptr);
      return lambda_result;
    };

    return ARKodeGetRootInfo_adapt_arr_ptr_to_std_vector(arkode_mem,
                                                         rootsfound_1d);
  },
  nb::arg("arkode_mem"), nb::arg("rootsfound_1d"));

m.def("ARKodePrintAllStats", ARKodePrintAllStats, nb::arg("arkode_mem"),
      nb::arg("outfile"), nb::arg("fmt"));

m.def("ARKodeGetReturnFlagName", ARKodeGetReturnFlagName, nb::arg("flag"));

m.def("ARKodeWriteParameters", ARKodeWriteParameters, nb::arg("arkode_mem"),
      nb::arg("fp"));

m.def(
  "ARKodeGetNumExpSteps",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumExpSteps_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long expsteps_adapt_modifiable;

      int r = ARKodeGetNumExpSteps(arkode_mem, &expsteps_adapt_modifiable);
      return std::make_tuple(r, expsteps_adapt_modifiable);
    };

    return ARKodeGetNumExpSteps_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumAccSteps",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumAccSteps_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long accsteps_adapt_modifiable;

      int r = ARKodeGetNumAccSteps(arkode_mem, &accsteps_adapt_modifiable);
      return std::make_tuple(r, accsteps_adapt_modifiable);
    };

    return ARKodeGetNumAccSteps_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumErrTestFails",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumErrTestFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long netfails_adapt_modifiable;

      int r = ARKodeGetNumErrTestFails(arkode_mem, &netfails_adapt_modifiable);
      return std::make_tuple(r, netfails_adapt_modifiable);
    };

    return ARKodeGetNumErrTestFails_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def("ARKodeGetEstLocalErrors", ARKodeGetEstLocalErrors, nb::arg("arkode_mem"),
      nb::arg("ele"));

m.def(
  "ARKodeGetActualInitStep",
  [](void* arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetActualInitStep_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype hinused_adapt_modifiable;

      int r = ARKodeGetActualInitStep(arkode_mem, &hinused_adapt_modifiable);
      return std::make_tuple(r, hinused_adapt_modifiable);
    };

    return ARKodeGetActualInitStep_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetTolScaleFactor",
  [](void* arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetTolScaleFactor_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tolsfac_adapt_modifiable;

      int r = ARKodeGetTolScaleFactor(arkode_mem, &tolsfac_adapt_modifiable);
      return std::make_tuple(r, tolsfac_adapt_modifiable);
    };

    return ARKodeGetTolScaleFactor_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumConstrFails",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumConstrFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nconstrfails_adapt_modifiable;

      int r = ARKodeGetNumConstrFails(arkode_mem, &nconstrfails_adapt_modifiable);
      return std::make_tuple(r, nconstrfails_adapt_modifiable);
    };

    return ARKodeGetNumConstrFails_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetStepStats",
  [](void* arkode_mem)
    -> std::tuple<int, long, sunrealtype, sunrealtype, sunrealtype, sunrealtype>
  {
    auto ARKodeGetStepStats_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem)
      -> std::tuple<int, long, sunrealtype, sunrealtype, sunrealtype, sunrealtype>
    {
      long nsteps_adapt_modifiable;
      sunrealtype hinused_adapt_modifiable;
      sunrealtype hlast_adapt_modifiable;
      sunrealtype hcur_adapt_modifiable;
      sunrealtype tcur_adapt_modifiable;

      int r = ARKodeGetStepStats(arkode_mem, &nsteps_adapt_modifiable,
                                 &hinused_adapt_modifiable,
                                 &hlast_adapt_modifiable,
                                 &hcur_adapt_modifiable, &tcur_adapt_modifiable);
      return std::make_tuple(r, nsteps_adapt_modifiable,
                             hinused_adapt_modifiable, hlast_adapt_modifiable,
                             hcur_adapt_modifiable, tcur_adapt_modifiable);
    };

    return ARKodeGetStepStats_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetAccumulatedError",
  [](void* arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetAccumulatedError_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype accum_error_adapt_modifiable;

      int r = ARKodeGetAccumulatedError(arkode_mem,
                                        &accum_error_adapt_modifiable);
      return std::make_tuple(r, accum_error_adapt_modifiable);
    };

    return ARKodeGetAccumulatedError_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumLinSolvSetups",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nlinsetups_adapt_modifiable;

      int r = ARKodeGetNumLinSolvSetups(arkode_mem, &nlinsetups_adapt_modifiable);
      return std::make_tuple(r, nlinsetups_adapt_modifiable);
    };

    return ARKodeGetNumLinSolvSetups_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetCurrentTime",
  [](void* arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetCurrentTime_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype tcur_adapt_modifiable;

      int r = ARKodeGetCurrentTime(arkode_mem, &tcur_adapt_modifiable);
      return std::make_tuple(r, tcur_adapt_modifiable);
    };

    return ARKodeGetCurrentTime_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetCurrentState",
  [](void* arkode_mem) -> std::tuple<int, N_Vector>
  {
    auto ARKodeGetCurrentState_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, N_Vector>
    {
      N_Vector state_adapt_modifiable;

      int r = ARKodeGetCurrentState(arkode_mem, &state_adapt_modifiable);
      return std::make_tuple(r, state_adapt_modifiable);
    };

    return ARKodeGetCurrentState_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

m.def(
  "ARKodeGetCurrentGamma",
  [](void* arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetCurrentGamma_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype gamma_adapt_modifiable;

      int r = ARKodeGetCurrentGamma(arkode_mem, &gamma_adapt_modifiable);
      return std::make_tuple(r, gamma_adapt_modifiable);
    };

    return ARKodeGetCurrentGamma_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumNonlinSolvIters",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nniters_adapt_modifiable;

      int r = ARKodeGetNumNonlinSolvIters(arkode_mem, &nniters_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable);
    };

    return ARKodeGetNumNonlinSolvIters_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumNonlinSolvConvFails",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nnfails_adapt_modifiable;

      int r = ARKodeGetNumNonlinSolvConvFails(arkode_mem,
                                              &nnfails_adapt_modifiable);
      return std::make_tuple(r, nnfails_adapt_modifiable);
    };

    return ARKodeGetNumNonlinSolvConvFails_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNonlinSolvStats",
  [](void* arkode_mem) -> std::tuple<int, long, long>
  {
    auto ARKodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long, long>
    {
      long nniters_adapt_modifiable;
      long nnfails_adapt_modifiable;

      int r = ARKodeGetNonlinSolvStats(arkode_mem, &nniters_adapt_modifiable,
                                       &nnfails_adapt_modifiable);
      return std::make_tuple(r, nniters_adapt_modifiable,
                             nnfails_adapt_modifiable);
    };

    return ARKodeGetNonlinSolvStats_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumStepSolveFails",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nncfails_adapt_modifiable;

      int r = ARKodeGetNumStepSolveFails(arkode_mem, &nncfails_adapt_modifiable);
      return std::make_tuple(r, nncfails_adapt_modifiable);
    };

    return ARKodeGetNumStepSolveFails_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetJac",
  [](void* arkode_mem) -> std::tuple<int, SUNMatrix>
  {
    auto ARKodeGetJac_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, SUNMatrix>
    {
      SUNMatrix J_adapt_modifiable;

      int r = ARKodeGetJac(arkode_mem, &J_adapt_modifiable);
      return std::make_tuple(r, J_adapt_modifiable);
    };

    return ARKodeGetJac_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"), "nb::rv_policy::reference", nb::rv_policy::reference);

m.def(
  "ARKodeGetJacTime",
  [](void* arkode_mem) -> std::tuple<int, sunrealtype>
  {
    auto ARKodeGetJacTime_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, sunrealtype>
    {
      sunrealtype t_J_adapt_modifiable;

      int r = ARKodeGetJacTime(arkode_mem, &t_J_adapt_modifiable);
      return std::make_tuple(r, t_J_adapt_modifiable);
    };

    return ARKodeGetJacTime_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetJacNumSteps",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetJacNumSteps_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nst_J_adapt_modifiable;

      int r = ARKodeGetJacNumSteps(arkode_mem, &nst_J_adapt_modifiable);
      return std::make_tuple(r, nst_J_adapt_modifiable);
    };

    return ARKodeGetJacNumSteps_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumJacEvals",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumJacEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long njevals_adapt_modifiable;

      int r = ARKodeGetNumJacEvals(arkode_mem, &njevals_adapt_modifiable);
      return std::make_tuple(r, njevals_adapt_modifiable);
    };

    return ARKodeGetNumJacEvals_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumPrecEvals",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumPrecEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long npevals_adapt_modifiable;

      int r = ARKodeGetNumPrecEvals(arkode_mem, &npevals_adapt_modifiable);
      return std::make_tuple(r, npevals_adapt_modifiable);
    };

    return ARKodeGetNumPrecEvals_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumPrecSolves",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumPrecSolves_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long npsolves_adapt_modifiable;

      int r = ARKodeGetNumPrecSolves(arkode_mem, &npsolves_adapt_modifiable);
      return std::make_tuple(r, npsolves_adapt_modifiable);
    };

    return ARKodeGetNumPrecSolves_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumLinIters",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumLinIters_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nliters_adapt_modifiable;

      int r = ARKodeGetNumLinIters(arkode_mem, &nliters_adapt_modifiable);
      return std::make_tuple(r, nliters_adapt_modifiable);
    };

    return ARKodeGetNumLinIters_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumLinConvFails",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumLinConvFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nlcfails_adapt_modifiable;

      int r = ARKodeGetNumLinConvFails(arkode_mem, &nlcfails_adapt_modifiable);
      return std::make_tuple(r, nlcfails_adapt_modifiable);
    };

    return ARKodeGetNumLinConvFails_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumJTSetupEvals",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long njtsetups_adapt_modifiable;

      int r = ARKodeGetNumJTSetupEvals(arkode_mem, &njtsetups_adapt_modifiable);
      return std::make_tuple(r, njtsetups_adapt_modifiable);
    };

    return ARKodeGetNumJTSetupEvals_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumJtimesEvals",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long njvevals_adapt_modifiable;

      int r = ARKodeGetNumJtimesEvals(arkode_mem, &njvevals_adapt_modifiable);
      return std::make_tuple(r, njvevals_adapt_modifiable);
    };

    return ARKodeGetNumJtimesEvals_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumLinRhsEvals",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nfevalsLS_adapt_modifiable;

      int r = ARKodeGetNumLinRhsEvals(arkode_mem, &nfevalsLS_adapt_modifiable);
      return std::make_tuple(r, nfevalsLS_adapt_modifiable);
    };

    return ARKodeGetNumLinRhsEvals_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetLastLinFlag",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetLastLinFlag_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long flag_adapt_modifiable;

      int r = ARKodeGetLastLinFlag(arkode_mem, &flag_adapt_modifiable);
      return std::make_tuple(r, flag_adapt_modifiable);
    };

    return ARKodeGetLastLinFlag_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def("ARKodeGetLinReturnFlagName", ARKodeGetLinReturnFlagName, nb::arg("flag"));

m.def(
  "ARKodeGetCurrentMassMatrix",
  [](void* arkode_mem) -> std::tuple<int, SUNMatrix>
  {
    auto ARKodeGetCurrentMassMatrix_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, SUNMatrix>
    {
      SUNMatrix M_adapt_modifiable;

      int r = ARKodeGetCurrentMassMatrix(arkode_mem, &M_adapt_modifiable);
      return std::make_tuple(r, M_adapt_modifiable);
    };

    return ARKodeGetCurrentMassMatrix_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"), " Optional output functions (non-identity mass matrices)\n\n nb::rv_policy::reference",
  nb::rv_policy::reference);

m.def("ARKodeGetResWeights", ARKodeGetResWeights, nb::arg("arkode_mem"),
      nb::arg("rweight"));

m.def(
  "ARKodeGetNumMassSetups",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassSetups_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmsetups_adapt_modifiable;

      int r = ARKodeGetNumMassSetups(arkode_mem, &nmsetups_adapt_modifiable);
      return std::make_tuple(r, nmsetups_adapt_modifiable);
    };

    return ARKodeGetNumMassSetups_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumMassMultSetups",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassMultSetups_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmvsetups_adapt_modifiable;

      int r = ARKodeGetNumMassMultSetups(arkode_mem, &nmvsetups_adapt_modifiable);
      return std::make_tuple(r, nmvsetups_adapt_modifiable);
    };

    return ARKodeGetNumMassMultSetups_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumMassMult",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassMult_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmvevals_adapt_modifiable;

      int r = ARKodeGetNumMassMult(arkode_mem, &nmvevals_adapt_modifiable);
      return std::make_tuple(r, nmvevals_adapt_modifiable);
    };

    return ARKodeGetNumMassMult_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumMassSolves",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassSolves_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmsolves_adapt_modifiable;

      int r = ARKodeGetNumMassSolves(arkode_mem, &nmsolves_adapt_modifiable);
      return std::make_tuple(r, nmsolves_adapt_modifiable);
    };

    return ARKodeGetNumMassSolves_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumMassPrecEvals",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassPrecEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmpevals_adapt_modifiable;

      int r = ARKodeGetNumMassPrecEvals(arkode_mem, &nmpevals_adapt_modifiable);
      return std::make_tuple(r, nmpevals_adapt_modifiable);
    };

    return ARKodeGetNumMassPrecEvals_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumMassPrecSolves",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassPrecSolves_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmpsolves_adapt_modifiable;

      int r = ARKodeGetNumMassPrecSolves(arkode_mem, &nmpsolves_adapt_modifiable);
      return std::make_tuple(r, nmpsolves_adapt_modifiable);
    };

    return ARKodeGetNumMassPrecSolves_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumMassIters",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassIters_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmiters_adapt_modifiable;

      int r = ARKodeGetNumMassIters(arkode_mem, &nmiters_adapt_modifiable);
      return std::make_tuple(r, nmiters_adapt_modifiable);
    };

    return ARKodeGetNumMassIters_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumMassConvFails",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMassConvFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmcfails_adapt_modifiable;

      int r = ARKodeGetNumMassConvFails(arkode_mem, &nmcfails_adapt_modifiable);
      return std::make_tuple(r, nmcfails_adapt_modifiable);
    };

    return ARKodeGetNumMassConvFails_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumMTSetups",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumMTSetups_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long nmtsetups_adapt_modifiable;

      int r = ARKodeGetNumMTSetups(arkode_mem, &nmtsetups_adapt_modifiable);
      return std::make_tuple(r, nmtsetups_adapt_modifiable);
    };

    return ARKodeGetNumMTSetups_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetLastMassFlag",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetLastMassFlag_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long flag_adapt_modifiable;

      int r = ARKodeGetLastMassFlag(arkode_mem, &flag_adapt_modifiable);
      return std::make_tuple(r, flag_adapt_modifiable);
    };

    return ARKodeGetLastMassFlag_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def("ARKodePrintMem", ARKodePrintMem, nb::arg("arkode_mem"), nb::arg("outfile"),
      "Output the ARKODE memory structure (useful when debugging)");

m.def("ARKodeSetRelaxEtaFail", ARKodeSetRelaxEtaFail, nb::arg("arkode_mem"),
      nb::arg("eta_rf"));

m.def("ARKodeSetRelaxLowerBound", ARKodeSetRelaxLowerBound,
      nb::arg("arkode_mem"), nb::arg("lower"));

m.def("ARKodeSetRelaxMaxFails", ARKodeSetRelaxMaxFails, nb::arg("arkode_mem"),
      nb::arg("max_fails"));

m.def("ARKodeSetRelaxMaxIters", ARKodeSetRelaxMaxIters, nb::arg("arkode_mem"),
      nb::arg("max_iters"));

m.def("ARKodeSetRelaxSolver", ARKodeSetRelaxSolver, nb::arg("arkode_mem"),
      nb::arg("solver"));

m.def("ARKodeSetRelaxResTol", ARKodeSetRelaxResTol, nb::arg("arkode_mem"),
      nb::arg("res_tol"));

m.def("ARKodeSetRelaxTol", ARKodeSetRelaxTol, nb::arg("arkode_mem"),
      nb::arg("rel_tol"), nb::arg("abs_tol"));

m.def("ARKodeSetRelaxUpperBound", ARKodeSetRelaxUpperBound,
      nb::arg("arkode_mem"), nb::arg("upper"));

m.def(
  "ARKodeGetNumRelaxFnEvals",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRelaxFnEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long r_evals_adapt_modifiable;

      int r = ARKodeGetNumRelaxFnEvals(arkode_mem, &r_evals_adapt_modifiable);
      return std::make_tuple(r, r_evals_adapt_modifiable);
    };

    return ARKodeGetNumRelaxFnEvals_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumRelaxJacEvals",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRelaxJacEvals_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long J_evals_adapt_modifiable;

      int r = ARKodeGetNumRelaxJacEvals(arkode_mem, &J_evals_adapt_modifiable);
      return std::make_tuple(r, J_evals_adapt_modifiable);
    };

    return ARKodeGetNumRelaxJacEvals_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumRelaxFails",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRelaxFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long relax_fails_adapt_modifiable;

      int r = ARKodeGetNumRelaxFails(arkode_mem, &relax_fails_adapt_modifiable);
      return std::make_tuple(r, relax_fails_adapt_modifiable);
    };

    return ARKodeGetNumRelaxFails_adapt_modifiable_immutable_to_return(arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumRelaxBoundFails",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRelaxBoundFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long fails_adapt_modifiable;

      int r = ARKodeGetNumRelaxBoundFails(arkode_mem, &fails_adapt_modifiable);
      return std::make_tuple(r, fails_adapt_modifiable);
    };

    return ARKodeGetNumRelaxBoundFails_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumRelaxSolveFails",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRelaxSolveFails_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long fails_adapt_modifiable;

      int r = ARKodeGetNumRelaxSolveFails(arkode_mem, &fails_adapt_modifiable);
      return std::make_tuple(r, fails_adapt_modifiable);
    };

    return ARKodeGetNumRelaxSolveFails_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeGetNumRelaxSolveIters",
  [](void* arkode_mem) -> std::tuple<int, long>
  {
    auto ARKodeGetNumRelaxSolveIters_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, long>
    {
      long iters_adapt_modifiable;

      int r = ARKodeGetNumRelaxSolveIters(arkode_mem, &iters_adapt_modifiable);
      return std::make_tuple(r, iters_adapt_modifiable);
    };

    return ARKodeGetNumRelaxSolveIters_adapt_modifiable_immutable_to_return(
      arkode_mem);
  },
  nb::arg("arkode_mem"));

m.def(
  "ARKodeCreateSUNStepper",
  [](void* arkode_mem)
    -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<SUNStepper>>>
  {
    auto ARKodeCreateSUNStepper_adapt_modifiable_immutable_to_return =
      [](void* arkode_mem) -> std::tuple<int, SUNStepper>
    {
      SUNStepper stepper_adapt_modifiable;

      int r = ARKodeCreateSUNStepper(arkode_mem, &stepper_adapt_modifiable);
      return std::make_tuple(r, stepper_adapt_modifiable);
    };
    auto ARKodeCreateSUNStepper_adapt_return_type_to_shared_ptr =
      [&ARKodeCreateSUNStepper_adapt_modifiable_immutable_to_return](
        void* arkode_mem)
      -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<SUNStepper>>>
    {
      auto lambda_result =
        ARKodeCreateSUNStepper_adapt_modifiable_immutable_to_return(arkode_mem);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<std::remove_pointer_t<SUNStepper>,
                                             SUNStepperDeleter>(
                               std::get<1>(lambda_result)));
    };

    return ARKodeCreateSUNStepper_adapt_return_type_to_shared_ptr(arkode_mem);
  },
  nb::arg("arkode_mem"), "SUNStepper functions", nb::rv_policy::reference);
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _ARKLS_H
//
// #ifdef __cplusplus
// #endif
//
m.attr("ARKLS_SUCCESS")          = 0;
m.attr("ARKLS_MEM_NULL")         = -1;
m.attr("ARKLS_LMEM_NULL")        = -2;
m.attr("ARKLS_ILL_INPUT")        = -3;
m.attr("ARKLS_MEM_FAIL")         = -4;
m.attr("ARKLS_PMEM_NULL")        = -5;
m.attr("ARKLS_MASSMEM_NULL")     = -6;
m.attr("ARKLS_JACFUNC_UNRECVR")  = -7;
m.attr("ARKLS_JACFUNC_RECVR")    = -8;
m.attr("ARKLS_MASSFUNC_UNRECVR") = -9;
m.attr("ARKLS_MASSFUNC_RECVR")   = -10;
m.attr("ARKLS_SUNMAT_FAIL")      = -11;
m.attr("ARKLS_SUNLS_FAIL")       = -12;

m.def(
  "ARKodeSetLinearSolver",
  [](void* arkode_mem, SUNLinearSolver LS,
     std::optional<SUNMatrix> A = std::nullopt) -> int
  {
    auto ARKodeSetLinearSolver_adapt_optional_arg_with_default_null =
      [](void* arkode_mem, SUNLinearSolver LS,
         std::optional<SUNMatrix> A = std::nullopt) -> int
    {
      SUNMatrix A_adapt_default_null = nullptr;
      if (A.has_value()) A_adapt_default_null = A.value();

      auto lambda_result = ARKodeSetLinearSolver(arkode_mem, LS,
                                                 A_adapt_default_null);
      return lambda_result;
    };

    return ARKodeSetLinearSolver_adapt_optional_arg_with_default_null(arkode_mem,
                                                                      LS, A);
  },
  nb::arg("arkode_mem"), nb::arg("LS"), nb::arg("A").none() = nb::none());

m.def("ARKodeSetJacEvalFrequency", ARKodeSetJacEvalFrequency,
      nb::arg("arkode_mem"), nb::arg("msbj"));

m.def("ARKodeSetLinearSolutionScaling", ARKodeSetLinearSolutionScaling,
      nb::arg("arkode_mem"), nb::arg("onoff"));

m.def("ARKodeSetEpsLin", ARKodeSetEpsLin, nb::arg("arkode_mem"),
      nb::arg("eplifac"));

m.def("ARKodeSetMassEpsLin", ARKodeSetMassEpsLin, nb::arg("arkode_mem"),
      nb::arg("eplifac"));

m.def("ARKodeSetLSNormFactor", ARKodeSetLSNormFactor, nb::arg("arkode_mem"),
      nb::arg("nrmfac"));

m.def("ARKodeSetMassLSNormFactor", ARKodeSetMassLSNormFactor,
      nb::arg("arkode_mem"), nb::arg("nrmfac"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _ARKODE_BUTCHER_H
//
// #ifdef __cplusplus
// #endif
//

auto pyClassARKodeButcherTableMem =
  nb::class_<ARKodeButcherTableMem>(m, "ARKodeButcherTableMem", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "ARKodeButcherTable_Create",
  [](int s, int q, int p, sundials4py::Array1d c_1d, sundials4py::Array1d A_1d,
     sundials4py::Array1d b_1d, sundials4py::Array1d d_1d)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
  {
    auto ARKodeButcherTable_Create_adapt_arr_ptr_to_std_vector =
      [](int s, int q, int p, sundials4py::Array1d c_1d,
         sundials4py::Array1d A_1d, sundials4py::Array1d b_1d,
         sundials4py::Array1d d_1d) -> ARKodeButcherTable
    {
      sunrealtype* c_1d_ptr = c_1d.size() == 0 ? nullptr : c_1d.data();
      sunrealtype* A_1d_ptr = A_1d.size() == 0 ? nullptr : A_1d.data();
      sunrealtype* b_1d_ptr = b_1d.size() == 0 ? nullptr : b_1d.data();
      sunrealtype* d_1d_ptr = d_1d.size() == 0 ? nullptr : d_1d.data();

      auto lambda_result = ARKodeButcherTable_Create(s, q, p, c_1d_ptr, A_1d_ptr,
                                                     b_1d_ptr, d_1d_ptr);
      return lambda_result;
    };
    auto ARKodeButcherTable_Create_adapt_return_type_to_shared_ptr =
      [&ARKodeButcherTable_Create_adapt_arr_ptr_to_std_vector](int s, int q,
                                                               int p,
                                                               sundials4py::Array1d c_1d,
                                                               sundials4py::Array1d A_1d,
                                                               sundials4py::Array1d b_1d,
                                                               sundials4py::Array1d d_1d)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
    {
      auto lambda_result =
        ARKodeButcherTable_Create_adapt_arr_ptr_to_std_vector(s, q, p, c_1d,
                                                              A_1d, b_1d, d_1d);

      return our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                             ARKodeButcherTableDeleter>(lambda_result);
    };

    return ARKodeButcherTable_Create_adapt_return_type_to_shared_ptr(s, q, p,
                                                                     c_1d, A_1d,
                                                                     b_1d, d_1d);
  },
  nb::arg("s"), nb::arg("q"), nb::arg("p"), nb::arg("c_1d"), nb::arg("A_1d"),
  nb::arg("b_1d"), nb::arg("d_1d"));

m.def(
  "ARKodeButcherTable_Copy",
  [](ARKodeButcherTable B)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
  {
    auto ARKodeButcherTable_Copy_adapt_return_type_to_shared_ptr =
      [](ARKodeButcherTable B)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
    {
      auto lambda_result = ARKodeButcherTable_Copy(B);

      return our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                             ARKodeButcherTableDeleter>(lambda_result);
    };

    return ARKodeButcherTable_Copy_adapt_return_type_to_shared_ptr(B);
  },
  nb::arg("B"));

m.def("ARKodeButcherTable_Write", ARKodeButcherTable_Write, nb::arg("B"),
      nb::arg("outfile"));

m.def("ARKodeButcherTable_IsStifflyAccurate",
      ARKodeButcherTable_IsStifflyAccurate, nb::arg("B"));

m.def(
  "ARKodeButcherTable_CheckOrder",
  [](ARKodeButcherTable B, FILE* outfile) -> std::tuple<int, int, int>
  {
    auto ARKodeButcherTable_CheckOrder_adapt_modifiable_immutable_to_return =
      [](ARKodeButcherTable B, FILE* outfile) -> std::tuple<int, int, int>
    {
      int q_adapt_modifiable;
      int p_adapt_modifiable;

      int r = ARKodeButcherTable_CheckOrder(B, &q_adapt_modifiable,
                                            &p_adapt_modifiable, outfile);
      return std::make_tuple(r, q_adapt_modifiable, p_adapt_modifiable);
    };

    return ARKodeButcherTable_CheckOrder_adapt_modifiable_immutable_to_return(B,
                                                                              outfile);
  },
  nb::arg("B"), nb::arg("outfile"));

m.def(
  "ARKodeButcherTable_CheckARKOrder",
  [](ARKodeButcherTable B1, ARKodeButcherTable B2,
     FILE* outfile) -> std::tuple<int, int, int>
  {
    auto ARKodeButcherTable_CheckARKOrder_adapt_modifiable_immutable_to_return =
      [](ARKodeButcherTable B1, ARKodeButcherTable B2,
         FILE* outfile) -> std::tuple<int, int, int>
    {
      int q_adapt_modifiable;
      int p_adapt_modifiable;

      int r = ARKodeButcherTable_CheckARKOrder(B1, B2, &q_adapt_modifiable,
                                               &p_adapt_modifiable, outfile);
      return std::make_tuple(r, q_adapt_modifiable, p_adapt_modifiable);
    };

    return ARKodeButcherTable_CheckARKOrder_adapt_modifiable_immutable_to_return(B1,
                                                                                 B2,
                                                                                 outfile);
  },
  nb::arg("B1"), nb::arg("B2"), nb::arg("outfile"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _ARKODE_ERK_TABLES_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumARKODE_ERKTableID =
  nb::enum_<ARKODE_ERKTableID>(m, "ARKODE_ERKTableID", nb::is_arithmetic(), "")
    .value("ARKODE_ERK_NONE", ARKODE_ERK_NONE, "ensure enum is signed int")
    .value("ARKODE_HEUN_EULER_2_1_2", ARKODE_HEUN_EULER_2_1_2, "")
    .value("ARKODE_MIN_ERK_NUM", ARKODE_MIN_ERK_NUM, "")
    .value("ARKODE_BOGACKI_SHAMPINE_4_2_3", ARKODE_BOGACKI_SHAMPINE_4_2_3, "")
    .value("ARKODE_ARK324L2SA_ERK_4_2_3", ARKODE_ARK324L2SA_ERK_4_2_3, "")
    .value("ARKODE_ZONNEVELD_5_3_4", ARKODE_ZONNEVELD_5_3_4, "")
    .value("ARKODE_ARK436L2SA_ERK_6_3_4", ARKODE_ARK436L2SA_ERK_6_3_4, "")
    .value("ARKODE_SAYFY_ABURUB_6_3_4", ARKODE_SAYFY_ABURUB_6_3_4, "")
    .value("ARKODE_CASH_KARP_6_4_5", ARKODE_CASH_KARP_6_4_5, "")
    .value("ARKODE_FEHLBERG_6_4_5", ARKODE_FEHLBERG_6_4_5, "")
    .value("ARKODE_DORMAND_PRINCE_7_4_5", ARKODE_DORMAND_PRINCE_7_4_5, "")
    .value("ARKODE_ARK548L2SA_ERK_8_4_5", ARKODE_ARK548L2SA_ERK_8_4_5, "")
    .value("ARKODE_VERNER_8_5_6", ARKODE_VERNER_8_5_6, "")
    .value("ARKODE_FEHLBERG_13_7_8", ARKODE_FEHLBERG_13_7_8, "")
    .value("ARKODE_KNOTH_WOLKE_3_3", ARKODE_KNOTH_WOLKE_3_3, "")
    .value("ARKODE_ARK437L2SA_ERK_7_3_4", ARKODE_ARK437L2SA_ERK_7_3_4, "")
    .value("ARKODE_ARK548L2SAb_ERK_8_4_5", ARKODE_ARK548L2SAb_ERK_8_4_5, "")
    .value("ARKODE_ARK2_ERK_3_1_2", ARKODE_ARK2_ERK_3_1_2, "")
    .value("ARKODE_SOFRONIOU_SPALETTA_5_3_4", ARKODE_SOFRONIOU_SPALETTA_5_3_4, "")
    .value("ARKODE_SHU_OSHER_3_2_3", ARKODE_SHU_OSHER_3_2_3, "")
    .value("ARKODE_VERNER_9_5_6", ARKODE_VERNER_9_5_6, "")
    .value("ARKODE_VERNER_10_6_7", ARKODE_VERNER_10_6_7, "")
    .value("ARKODE_VERNER_13_7_8", ARKODE_VERNER_13_7_8, "")
    .value("ARKODE_VERNER_16_8_9", ARKODE_VERNER_16_8_9, "")
    .value("ARKODE_FORWARD_EULER_1_1", ARKODE_FORWARD_EULER_1_1, "")
    .value("ARKODE_RALSTON_EULER_2_1_2", ARKODE_RALSTON_EULER_2_1_2, "")
    .value("ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2",
           ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2, "")
    .value("ARKODE_RALSTON_3_1_2", ARKODE_RALSTON_3_1_2, "")
    .value("ARKODE_TSITOURAS_7_4_5", ARKODE_TSITOURAS_7_4_5, "")
    .value("ARKODE_MAX_ERK_NUM", ARKODE_MAX_ERK_NUM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

m.def(
  "ARKodeButcherTable_LoadERK",
  [](ARKODE_ERKTableID emethod)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
  {
    auto ARKodeButcherTable_LoadERK_adapt_return_type_to_shared_ptr =
      [](ARKODE_ERKTableID emethod)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
    {
      auto lambda_result = ARKodeButcherTable_LoadERK(emethod);

      return our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                             ARKodeButcherTableDeleter>(lambda_result);
    };

    return ARKodeButcherTable_LoadERK_adapt_return_type_to_shared_ptr(emethod);
  },
  nb::arg("emethod"), "Accessor routine to load built-in ERK table");

m.def(
  "ARKodeButcherTable_LoadERKByName",
  [](const char* emethod)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
  {
    auto ARKodeButcherTable_LoadERKByName_adapt_return_type_to_shared_ptr =
      [](const char* emethod)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
    {
      auto lambda_result = ARKodeButcherTable_LoadERKByName(emethod);

      return our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                             ARKodeButcherTableDeleter>(lambda_result);
    };

    return ARKodeButcherTable_LoadERKByName_adapt_return_type_to_shared_ptr(
      emethod);
  },
  nb::arg("emethod"));

m.def("ARKodeButcherTable_ERKIDToName", ARKodeButcherTable_ERKIDToName,
      nb::arg("emethod"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _ARKODE_DIRK_TABLES_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumARKODE_DIRKTableID =
  nb::enum_<ARKODE_DIRKTableID>(m, "ARKODE_DIRKTableID", nb::is_arithmetic(), "")
    .value("ARKODE_DIRK_NONE", ARKODE_DIRK_NONE, "ensure enum is signed int")
    .value("ARKODE_SDIRK_2_1_2", ARKODE_SDIRK_2_1_2, "")
    .value("ARKODE_MIN_DIRK_NUM", ARKODE_MIN_DIRK_NUM, "")
    .value("ARKODE_BILLINGTON_3_3_2", ARKODE_BILLINGTON_3_3_2, "")
    .value("ARKODE_TRBDF2_3_3_2", ARKODE_TRBDF2_3_3_2, "")
    .value("ARKODE_KVAERNO_4_2_3", ARKODE_KVAERNO_4_2_3, "")
    .value("ARKODE_ARK324L2SA_DIRK_4_2_3", ARKODE_ARK324L2SA_DIRK_4_2_3, "")
    .value("ARKODE_CASH_5_2_4", ARKODE_CASH_5_2_4, "")
    .value("ARKODE_CASH_5_3_4", ARKODE_CASH_5_3_4, "")
    .value("ARKODE_SDIRK_5_3_4", ARKODE_SDIRK_5_3_4, "")
    .value("ARKODE_KVAERNO_5_3_4", ARKODE_KVAERNO_5_3_4, "")
    .value("ARKODE_ARK436L2SA_DIRK_6_3_4", ARKODE_ARK436L2SA_DIRK_6_3_4, "")
    .value("ARKODE_KVAERNO_7_4_5", ARKODE_KVAERNO_7_4_5, "")
    .value("ARKODE_ARK548L2SA_DIRK_8_4_5", ARKODE_ARK548L2SA_DIRK_8_4_5, "")
    .value("ARKODE_ARK437L2SA_DIRK_7_3_4", ARKODE_ARK437L2SA_DIRK_7_3_4, "")
    .value("ARKODE_ARK548L2SAb_DIRK_8_4_5", ARKODE_ARK548L2SAb_DIRK_8_4_5, "")
    .value("ARKODE_ESDIRK324L2SA_4_2_3", ARKODE_ESDIRK324L2SA_4_2_3, "")
    .value("ARKODE_ESDIRK325L2SA_5_2_3", ARKODE_ESDIRK325L2SA_5_2_3, "")
    .value("ARKODE_ESDIRK32I5L2SA_5_2_3", ARKODE_ESDIRK32I5L2SA_5_2_3, "")
    .value("ARKODE_ESDIRK436L2SA_6_3_4", ARKODE_ESDIRK436L2SA_6_3_4, "")
    .value("ARKODE_ESDIRK43I6L2SA_6_3_4", ARKODE_ESDIRK43I6L2SA_6_3_4, "")
    .value("ARKODE_QESDIRK436L2SA_6_3_4", ARKODE_QESDIRK436L2SA_6_3_4, "")
    .value("ARKODE_ESDIRK437L2SA_7_3_4", ARKODE_ESDIRK437L2SA_7_3_4, "")
    .value("ARKODE_ESDIRK547L2SA_7_4_5", ARKODE_ESDIRK547L2SA_7_4_5, "")
    .value("ARKODE_ESDIRK547L2SA2_7_4_5", ARKODE_ESDIRK547L2SA2_7_4_5, "")
    .value("ARKODE_ARK2_DIRK_3_1_2", ARKODE_ARK2_DIRK_3_1_2, "")
    .value("ARKODE_BACKWARD_EULER_1_1", ARKODE_BACKWARD_EULER_1_1, "")
    .value("ARKODE_IMPLICIT_MIDPOINT_1_2", ARKODE_IMPLICIT_MIDPOINT_1_2, "")
    .value("ARKODE_IMPLICIT_TRAPEZOIDAL_2_2", ARKODE_IMPLICIT_TRAPEZOIDAL_2_2, "")
    .value("ARKODE_MAX_DIRK_NUM", ARKODE_MAX_DIRK_NUM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

m.def(
  "ARKodeButcherTable_LoadDIRK",
  [](ARKODE_DIRKTableID imethod)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
  {
    auto ARKodeButcherTable_LoadDIRK_adapt_return_type_to_shared_ptr =
      [](ARKODE_DIRKTableID imethod)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
    {
      auto lambda_result = ARKodeButcherTable_LoadDIRK(imethod);

      return our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                             ARKodeButcherTableDeleter>(lambda_result);
    };

    return ARKodeButcherTable_LoadDIRK_adapt_return_type_to_shared_ptr(imethod);
  },
  nb::arg("imethod"), "Accessor routine to load built-in DIRK table");

m.def(
  "ARKodeButcherTable_LoadDIRKByName",
  [](const char* imethod)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
  {
    auto ARKodeButcherTable_LoadDIRKByName_adapt_return_type_to_shared_ptr =
      [](const char* imethod)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>
    {
      auto lambda_result = ARKodeButcherTable_LoadDIRKByName(imethod);

      return our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                             ARKodeButcherTableDeleter>(lambda_result);
    };

    return ARKodeButcherTable_LoadDIRKByName_adapt_return_type_to_shared_ptr(
      imethod);
  },
  nb::arg("imethod"), "Accessor routine to load built-in DIRK table");

m.def("ARKodeButcherTable_DIRKIDToName", ARKodeButcherTable_DIRKIDToName,
      nb::arg("imethod"));
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
// #ifndef _ARKODE_SPRKTABLE_H
//
// #ifdef __cplusplus
// #endif
//

auto pyEnumARKODE_SPRKMethodID =
  nb::enum_<ARKODE_SPRKMethodID>(m, "ARKODE_SPRKMethodID", nb::is_arithmetic(), "")
    .value("ARKODE_SPRK_NONE", ARKODE_SPRK_NONE, "ensure enum is signed int")
    .value("ARKODE_SPRK_EULER_1_1", ARKODE_SPRK_EULER_1_1, "")
    .value("ARKODE_MIN_SPRK_NUM", ARKODE_MIN_SPRK_NUM, "")
    .value("ARKODE_SPRK_LEAPFROG_2_2", ARKODE_SPRK_LEAPFROG_2_2, "")
    .value("ARKODE_SPRK_PSEUDO_LEAPFROG_2_2", ARKODE_SPRK_PSEUDO_LEAPFROG_2_2, "")
    .value("ARKODE_SPRK_RUTH_3_3", ARKODE_SPRK_RUTH_3_3, "")
    .value("ARKODE_SPRK_MCLACHLAN_2_2", ARKODE_SPRK_MCLACHLAN_2_2, "")
    .value("ARKODE_SPRK_MCLACHLAN_3_3", ARKODE_SPRK_MCLACHLAN_3_3, "")
    .value("ARKODE_SPRK_CANDY_ROZMUS_4_4", ARKODE_SPRK_CANDY_ROZMUS_4_4, "")
    .value("ARKODE_SPRK_MCLACHLAN_4_4", ARKODE_SPRK_MCLACHLAN_4_4, "")
    .value("ARKODE_SPRK_MCLACHLAN_5_6", ARKODE_SPRK_MCLACHLAN_5_6, "")
    .value("ARKODE_SPRK_YOSHIDA_6_8", ARKODE_SPRK_YOSHIDA_6_8, "")
    .value("ARKODE_SPRK_SUZUKI_UMENO_8_16", ARKODE_SPRK_SUZUKI_UMENO_8_16, "")
    .value("ARKODE_SPRK_SOFRONIOU_10_36", ARKODE_SPRK_SOFRONIOU_10_36, "")
    .value("ARKODE_MAX_SPRK_NUM", ARKODE_MAX_SPRK_NUM, "")
    .export_values();
// #ifndef SWIG
//
// #endif
//

auto pyClassARKodeSPRKTableMem =
  nb::class_<ARKodeSPRKTableMem>(m, "ARKodeSPRKTableMem", "")
    .def(nb::init<>()) // implicit default constructor
  ;

m.def(
  "ARKodeSPRKTable_Create",
  [](int s, int q, sundials4py::Array1d a_1d, sundials4py::Array1d ahat_1d)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
  {
    auto ARKodeSPRKTable_Create_adapt_arr_ptr_to_std_vector =
      [](int s, int q, sundials4py::Array1d a_1d,
         sundials4py::Array1d ahat_1d) -> ARKodeSPRKTable
    {
      sunrealtype* a_1d_ptr    = a_1d.size() == 0 ? nullptr : a_1d.data();
      sunrealtype* ahat_1d_ptr = ahat_1d.size() == 0 ? nullptr : ahat_1d.data();

      auto lambda_result = ARKodeSPRKTable_Create(s, q, a_1d_ptr, ahat_1d_ptr);
      return lambda_result;
    };
    auto ARKodeSPRKTable_Create_adapt_return_type_to_shared_ptr =
      [&ARKodeSPRKTable_Create_adapt_arr_ptr_to_std_vector](int s, int q,
                                                            sundials4py::Array1d a_1d,
                                                            sundials4py::Array1d ahat_1d)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
    {
      auto lambda_result =
        ARKodeSPRKTable_Create_adapt_arr_ptr_to_std_vector(s, q, a_1d, ahat_1d);

      return our_make_shared<std::remove_pointer_t<ARKodeSPRKTable>,
                             ARKodeSPRKTableDeleter>(lambda_result);
    };

    return ARKodeSPRKTable_Create_adapt_return_type_to_shared_ptr(s, q, a_1d,
                                                                  ahat_1d);
  },
  nb::arg("s"), nb::arg("q"), nb::arg("a_1d"), nb::arg("ahat_1d"));

m.def(
  "ARKodeSPRKTable_Load",
  [](ARKODE_SPRKMethodID id)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
  {
    auto ARKodeSPRKTable_Load_adapt_return_type_to_shared_ptr =
      [](ARKODE_SPRKMethodID id)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
    {
      auto lambda_result = ARKodeSPRKTable_Load(id);

      return our_make_shared<std::remove_pointer_t<ARKodeSPRKTable>,
                             ARKodeSPRKTableDeleter>(lambda_result);
    };

    return ARKodeSPRKTable_Load_adapt_return_type_to_shared_ptr(id);
  },
  nb::arg("id"));

m.def(
  "ARKodeSPRKTable_LoadByName",
  [](const char* method) -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
  {
    auto ARKodeSPRKTable_LoadByName_adapt_return_type_to_shared_ptr =
      [](const char* method)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
    {
      auto lambda_result = ARKodeSPRKTable_LoadByName(method);

      return our_make_shared<std::remove_pointer_t<ARKodeSPRKTable>,
                             ARKodeSPRKTableDeleter>(lambda_result);
    };

    return ARKodeSPRKTable_LoadByName_adapt_return_type_to_shared_ptr(method);
  },
  nb::arg("method"));

m.def(
  "ARKodeSPRKTable_Copy",
  [](ARKodeSPRKTable that_sprk_storage)
    -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
  {
    auto ARKodeSPRKTable_Copy_adapt_return_type_to_shared_ptr =
      [](ARKodeSPRKTable that_sprk_storage)
      -> std::shared_ptr<std::remove_pointer_t<ARKodeSPRKTable>>
    {
      auto lambda_result = ARKodeSPRKTable_Copy(that_sprk_storage);

      return our_make_shared<std::remove_pointer_t<ARKodeSPRKTable>,
                             ARKodeSPRKTableDeleter>(lambda_result);
    };

    return ARKodeSPRKTable_Copy_adapt_return_type_to_shared_ptr(that_sprk_storage);
  },
  nb::arg("that_sprk_storage"));

m.def("ARKodeSPRKTable_Write", ARKodeSPRKTable_Write, nb::arg("sprk_table"),
      nb::arg("outfile"));

m.def(
  "ARKodeSPRKTable_ToButcher",
  [](ARKodeSPRKTable sprk_storage)
    -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>,
                  std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>>
  {
    auto ARKodeSPRKTable_ToButcher_adapt_modifiable_immutable_to_return =
      [](ARKodeSPRKTable sprk_storage)
      -> std::tuple<int, ARKodeButcherTable, ARKodeButcherTable>
    {
      ARKodeButcherTable a_ptr_adapt_modifiable;
      ARKodeButcherTable b_ptr_adapt_modifiable;

      int r = ARKodeSPRKTable_ToButcher(sprk_storage, &a_ptr_adapt_modifiable,
                                        &b_ptr_adapt_modifiable);
      return std::make_tuple(r, a_ptr_adapt_modifiable, b_ptr_adapt_modifiable);
    };
    auto ARKodeSPRKTable_ToButcher_adapt_return_type_to_shared_ptr =
      [&ARKodeSPRKTable_ToButcher_adapt_modifiable_immutable_to_return](
        ARKodeSPRKTable sprk_storage)
      -> std::tuple<int, std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>,
                    std::shared_ptr<std::remove_pointer_t<ARKodeButcherTable>>>
    {
      auto lambda_result =
        ARKodeSPRKTable_ToButcher_adapt_modifiable_immutable_to_return(
          sprk_storage);

      return std::make_tuple(std::get<0>(lambda_result),
                             our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                                             ARKodeButcherTableDeleter>(
                               std::get<1>(lambda_result)),
                             our_make_shared<std::remove_pointer_t<ARKodeButcherTable>,
                                             ARKodeButcherTableDeleter>(
                               std::get<2>(lambda_result)));
    };

    return ARKodeSPRKTable_ToButcher_adapt_return_type_to_shared_ptr(sprk_storage);
  },
  nb::arg("sprk_storage"), nb::rv_policy::reference);
// #ifdef __cplusplus
//
// #endif
//
// #endif
//
