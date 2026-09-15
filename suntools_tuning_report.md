# SUNDIALS Serial Example Tuning Report

Run date: 2026-09-15

## Scope and method

This report summarizes fresh `suntools tune` searches for all 40 identified
serial examples from ARKODE, CVODE, CVODES, IDA, IDAS, and KINSOL. The
configurations use DeepHyper 0.13.2 with its default `ET` surrogate;
the previous `surrogate_model: DUMMY` override was removed.

Each search retained its existing solver stack, workload, parameter ranges,
objective, correctness checks, evaluation budget, and one-worker execution.
Deterministic work searches use one repetition; wall-time searches use five.

The campaign evaluated 814 sampled configurations: 774 were feasible
and 40 failed or violated a configured correctness constraint.
Every sampled winner subsequently passed three independent validation reruns.

The best-observed value below is the better of the unchanged baseline and the
best feasible sampled trial. The worst value is the worst feasible sampled
trial. Results are specific to this workload, build, machine, and search.

## Results

| Package | Example | Objective | Baseline | Best observed | Change | Sampled worst | Change | Feasible |
|---|---|---|---:|---:|---:|---:|---:|---:|
| ARKODE | `ark_analytic_lsrk_domeigest` | RHS and estimator work | 132464 | 132126 | -0.3% | 143786 | +8.5% | 32/32 |
| ARKODE | `ark_analytic_lsrk` | RHS evaluations | 142049 | 132404 | -6.8% | 156923 | +10.5% | 24/24 |
| ARKODE | `ark_analytic_lsrk_varjac` | Wall time | 0.012243 s | 0.010887 s | -11.1% | 0.012881 s | +5.2% | 16/16 |
| ARKODE | `ark_analytic_mels` | Wall time | 0.001475 s | 0.001368 s | -7.3% | 0.001522 s | +3.1% | 6/16 |
| ARKODE | `ark_analytic_ssprk` | RHS evaluations | 6489 | 2004 | -69.1% | 5222 | -19.5% | 12/12 |
| ARKODE | `ark_analytic` | Solver work | 355 | 310 | -12.7% | 3560 | +902.8% | 40/40 |
| ARKODE | `ark_brusselator_lsrk_domeigest` | Wall time | 0.003183 s | 0.003149 s | -1.1% | 0.003339 s | +4.9% | 16/16 |
| ARKODE | `ark_brusselator_lsrk_externaldomeigest` | Wall time | 0.003431 s | 0.003301 s | -3.8% | 0.003768 s | +9.8% | 16/16 |
| ARKODE | `ark_heat1D` | Wall time | 0.012633 s | 0.011680 s | -7.5% | 0.015891 s | +25.8% | 16/16 |
| ARKODE | `ark_kpr_mri` | RHS evaluations | 1570314 | 115581 | -92.6% | 3170581 | +101.9% | 18/18 |
| ARKODE | `ark_robertson_constraints` | Wall time | 0.003290 s | 0.002248 s | -31.7% | 0.003798 s | +15.5% | 12/16 |
| ARKODE | `ark_robertson` | Wall time | 0.002281 s | 0.002200 s | -3.6% | 0.004884 s | +114.1% | 12/16 |
| CVODE | `cv_kpr` | Wall time | 0.002714 s | 0.002619 s | -3.5% | 0.003374 s | +24.3% | 13/16 |
| CVODE | `cvAnalytic_mels` | RHS evaluations | 135 | 120 | -11.1% | 895 | +563.0% | 16/16 |
| CVODE | `cvParticle_dns` | Wall time | 0.002922 s | 0.002732 s | -6.5% | 0.009587 s | +228.1% | 16/16 |
| CVODE | `cvPendulum_dns` | Solver work | 3890 | 3498 | -10.1% | 23498 | +504.1% | 36/40 |
| CVODE | `cvRoberts_dns_constraints` | Wall time | 0.001512 s | 0.001411 s | -6.6% | 0.001572 s | +4.0% | 16/16 |
| CVODE | `cvRoberts_dns` | Solver work | 872 | 726 | -16.7% | 1168 | +33.9% | 32/32 |
| CVODE | `cvVdp_auto_nls` | Wall time | 0.002217 s | 0.001409 s | -36.4% | 0.003140 s | +41.7% | 16/16 |
| CVODES | `cvsAnalytic_mels` | RHS evaluations | 135 | 120 | -11.1% | 896 | +563.7% | 16/16 |
| CVODES | `cvsParticle_dns` | Wall time | 0.002998 s | 0.002760 s | -7.9% | 0.010206 s | +240.4% | 16/16 |
| CVODES | `cvsPendulum_dns` | Solver work | 3890 | 3498 | -10.1% | 20594 | +429.4% | 36/40 |
| CVODES | `cvsRoberts_dns_constraints` | Wall time | 0.001489 s | 0.001378 s | -7.4% | 0.001612 s | +8.2% | 16/16 |
| CVODES | `cvsRoberts_dns` | Solver work | 872 | 718 | -17.7% | 1335 | +53.1% | 32/32 |
| IDA | `idaAnalytic_mels` | Residual evaluations | 576 | 576 (baseline) | +0.0% | 576 | +0.0% | 16/16 |
| IDA | `idaFoodWeb_bnd` | Wall time | 0.071602 s | 0.054051 s | -24.5% | 0.075892 s | +6.0% | 16/16 |
| IDA | `idaHeat2D_bnd` | Wall time | 0.001757 s | 0.001571 s | -10.6% | 0.001882 s | +7.1% | 16/16 |
| IDA | `idaRoberts_dns` | Solver work | 657 | 558 | -15.1% | 1260 | +91.8% | 32/32 |
| IDA | `idaSlCrank_dns` | Wall time | 0.001660 s | 0.001539 s | -7.3% | 0.001762 s | +6.1% | 16/16 |
| IDAS | `idasAnalytic_mels` | Residual evaluations | 576 | 576 (baseline) | +0.0% | 576 | +0.0% | 16/16 |
| IDAS | `idasFoodWeb_bnd` | Wall time | 0.071076 s | 0.053687 s | -24.5% | 0.079715 s | +12.2% | 16/16 |
| IDAS | `idasHeat2D_bnd` | Wall time | 0.001811 s | 0.001636 s | -9.7% | 0.002003 s | +10.6% | 16/16 |
| IDAS | `idasRoberts_dns` | Solver work | 657 | 562 | -14.5% | 1155 | +75.8% | 32/32 |
| IDAS | `idasSlCrank_dns` | Wall time | 0.001789 s | 0.001747 s | -2.3% | 0.002342 s | +30.9% | 16/16 |
| KINSOL | `kinAnalytic_fp` | Function evaluations | 5 | 5 (baseline) | +0.0% | 6 | +20.0% | 23/24 |
| KINSOL | `kinFerTron_dns` | Wall time | 0.001329 s | 0.001195 s | -10.0% | 0.001354 s | +1.9% | 16/16 |
| KINSOL | `kinFoodWeb_kry` | Wall time | 0.005526 s | 0.003310 s | -40.1% | 0.005535 s | +0.2% | 16/16 |
| KINSOL | `kinLaplace_bnd` | Wall time | 0.002303 s | 0.002213 s | -3.9% | 0.002831 s | +22.9% | 16/16 |
| KINSOL | `kinRoberts_fp` | Function evaluations | 8 | 8 (baseline) | +0.0% | 16 | +100.0% | 19/24 |
| KINSOL | `kinRoboKin_dns` | Wall time | 0.002422 s | 0.001248 s | -48.5% | 0.002809 s | +16.0% | 11/16 |

## Best sampled settings

These are the best feasible sampled parameters even where the unchanged
baseline tied or beat the sampled result.

| Example | Parameters |
|---|---|
| `ark_analytic_lsrk_domeigest` | `arkode.dom_eig_frequency=11`, `arkode.dom_eig_safety_factor=1.0180220625015552`, `arkode.num_dom_eig_est_preprocess_iters=0`, `sundomeigestimator.max_iters=32`, `sundomeigestimator.rel_tol=0.0014520346596627438` |
| `ark_analytic_lsrk` | `arkode.dom_eig_safety_factor=1.0058701637045775`, `arkode.sts_method_name="ARKODE_LSRK_RKC_2"`, `arkode.use_analytic_stability_region=0` |
| `ark_analytic_lsrk_varjac` | `arkode.dom_eig_frequency=18`, `arkode.dom_eig_safety_factor=1.0244847543132465` |
| `ark_analytic_mels` | `arkode.order=3`, `arkode.safety_factor=0.8008465497334075` |
| `ark_analytic_ssprk` | `arkode.safety_factor=0.8315593368080384`, `arkode.ssp_method_name="ARKODE_LSRK_SSP_S_2"` |
| `ark_analytic` | `arkode.table_names="ARKODE_TRBDF2_3_3_2 ARKODE_ERK_NONE"` |
| `ark_brusselator_lsrk_domeigest` | `arkode.dom_eig_frequency=31`, `arkode.dom_eig_safety_factor=1.0454051525474497` |
| `ark_brusselator_lsrk_externaldomeigest` | `arkode.dom_eig_frequency=44`, `arkode.dom_eig_safety_factor=1.0469221623761091` |
| `ark_heat1D` | `arkode.lsetup_frequency=31`, `arkode.order=3` |
| `ark_kpr_mri` | `inner.fixed_step=0.00019964022382486455`, `inner.table_names="ARKODE_DIRK_NONE ARKODE_BOGACKI_SHAMPINE_4_2_3"`, `outer.coupling_table_name="ARKODE_MRI_GARK_IRK21a"` |
| `ark_robertson_constraints` | `arkode.jac_eval_frequency=30`, `arkode.lsetup_frequency=37`, `arkode.order=3` |
| `ark_robertson` | `arkode.jac_eval_frequency=28`, `arkode.lsetup_frequency=40`, `arkode.order=4` |
| `cv_kpr` | `cvode.jac_eval_frequency=54`, `cvode.lsetup_frequency=35`, `cvode.max_order=4` |
| `cvAnalytic_mels` | `cvode.max_order=5`, `cvode.nonlin_conv_coef=0.2608218393905527` |
| `cvParticle_dns` | `cvode.jac_eval_frequency=90`, `cvode.lsetup_frequency=35`, `cvode.max_order=5` |
| `cvPendulum_dns` | `cvode.jac_eval_frequency=85`, `cvode.lsetup_frequency=39`, `cvode.max_order=5` |
| `cvRoberts_dns_constraints` | `cvode.jac_eval_frequency=39`, `cvode.lsetup_frequency=30`, `cvode.max_order=5` |
| `cvRoberts_dns` | `cvode.jac_eval_frequency=20`, `cvode.lsetup_frequency=39`, `cvode.max_order=5` |
| `cvVdp_auto_nls` | `cvode.jac_eval_frequency=66`, `cvode.lsetup_frequency=22`, `cvode.max_order=3` |
| `cvsAnalytic_mels` | `cvodes.max_order=5`, `cvodes.nonlin_conv_coef=0.23295000527370746` |
| `cvsParticle_dns` | `cvodes.jac_eval_frequency=32`, `cvodes.lsetup_frequency=16`, `cvodes.max_order=5` |
| `cvsPendulum_dns` | `cvodes.jac_eval_frequency=85`, `cvodes.lsetup_frequency=39`, `cvodes.max_order=5` |
| `cvsRoberts_dns_constraints` | `cvodes.jac_eval_frequency=41`, `cvodes.lsetup_frequency=38`, `cvodes.max_order=3` |
| `cvsRoberts_dns` | `cvodes.jac_eval_frequency=26`, `cvodes.lsetup_frequency=8`, `cvodes.max_order=4` |
| `idaAnalytic_mels` | `ida.max_order=2`, `ida.nonlin_conv_coef=0.28366152731576105` |
| `idaFoodWeb_bnd` | `ida.delta_cj_lsetup=0.4314870625849645`, `ida.max_order=4`, `ida.nonlin_conv_coef=0.22003199437162826` |
| `idaHeat2D_bnd` | `ida.delta_cj_lsetup=0.4658315406062523`, `ida.max_order=2`, `ida.nonlin_conv_coef=0.2808957649163569` |
| `idaRoberts_dns` | `ida.delta_cj_lsetup=0.18019514670321335`, `ida.max_order=4`, `ida.nonlin_conv_coef=0.09289868471126751` |
| `idaSlCrank_dns` | `ida.delta_cj_lsetup=0.22316486880927022`, `ida.max_order=5`, `ida.nonlin_conv_coef=0.0775907765506927` |
| `idasAnalytic_mels` | `idas.max_order=5`, `idas.nonlin_conv_coef=0.07689651807894142` |
| `idasFoodWeb_bnd` | `idas.delta_cj_lsetup=0.4943400299375569`, `idas.max_order=4`, `idas.nonlin_conv_coef=0.04155455914962843` |
| `idasHeat2D_bnd` | `idas.delta_cj_lsetup=0.4161770320864516`, `idas.max_order=3`, `idas.nonlin_conv_coef=0.12058821202715875` |
| `idasRoberts_dns` | `idas.delta_cj_lsetup=0.19300326969484327`, `idas.max_order=4`, `idas.nonlin_conv_coef=0.17419669406363514` |
| `idasSlCrank_dns` | `idas.delta_cj_lsetup=0.3102160469474512`, `idas.max_order=5`, `idas.nonlin_conv_coef=0.20246575146515564` |
| `kinAnalytic_fp` | `kinsol.damping_aa=0.8950021784157111`, `kinsol.delay_aa=3`, `kinsol.m_aa=3` |
| `kinFerTron_dns` | `kinsol.eta_const_value=0.3135418953400806`, `kinsol.eta_form=3` |
| `kinFoodWeb_kry` | `kinsol.eta_const_value=0.3130322547517779`, `kinsol.eta_form=2`, `kinsol.max_setup_calls=1` |
| `kinLaplace_bnd` | `kinsol.max_setup_calls=72`, `kinsol.max_sub_setup_calls=2` |
| `kinRoberts_fp` | `kinsol.damping_aa=0.8457251546937141`, `kinsol.delay_aa=0`, `kinsol.m_aa=2` |
| `kinRoboKin_dns` | `kinsol.eta_const_value=0.3055570418734109`, `kinsol.eta_form=1`, `kinsol.max_setup_calls=15` |

## Worst feasible sampled settings

| Example | Parameters |
|---|---|
| `ark_analytic_lsrk_domeigest` | `arkode.dom_eig_frequency=13`, `arkode.dom_eig_safety_factor=1.1618865766041035`, `arkode.num_dom_eig_est_preprocess_iters=3`, `sundomeigestimator.max_iters=64`, `sundomeigestimator.rel_tol=0.008630791865915127` |
| `ark_analytic_lsrk` | `arkode.dom_eig_safety_factor=1.1666536593910672`, `arkode.sts_method_name="ARKODE_LSRK_RKL_2"`, `arkode.use_analytic_stability_region=0` |
| `ark_analytic_lsrk_varjac` | `arkode.dom_eig_frequency=38`, `arkode.dom_eig_safety_factor=1.1372674350632115` |
| `ark_analytic_mels` | `arkode.order=3`, `arkode.safety_factor=0.8078528184575422` |
| `ark_analytic_ssprk` | `arkode.safety_factor=0.8196655160712274`, `arkode.ssp_method_name="ARKODE_LSRK_SSP_S_3"` |
| `ark_analytic` | `arkode.table_names="ARKODE_CASH_5_2_4 ARKODE_ERK_NONE"` |
| `ark_brusselator_lsrk_domeigest` | `arkode.dom_eig_frequency=18`, `arkode.dom_eig_safety_factor=1.0761798987599809` |
| `ark_brusselator_lsrk_externaldomeigest` | `arkode.dom_eig_frequency=12`, `arkode.dom_eig_safety_factor=1.0632103976993121` |
| `ark_heat1D` | `arkode.lsetup_frequency=37`, `arkode.order=2` |
| `ark_kpr_mri` | `inner.fixed_step=9.544751300522987e-06`, `inner.table_names="ARKODE_DIRK_NONE ARKODE_ARK436L2SA_ERK_6_3_4"`, `outer.coupling_table_name="ARKODE_MRI_GARK_IRK21a"` |
| `ark_robertson_constraints` | `arkode.jac_eval_frequency=30`, `arkode.lsetup_frequency=40`, `arkode.order=3` |
| `ark_robertson` | `arkode.jac_eval_frequency=54`, `arkode.lsetup_frequency=14`, `arkode.order=2` |
| `cv_kpr` | `cvode.jac_eval_frequency=71`, `cvode.lsetup_frequency=40`, `cvode.max_order=3` |
| `cvAnalytic_mels` | `cvode.max_order=2`, `cvode.nonlin_conv_coef=0.046092784586293054` |
| `cvParticle_dns` | `cvode.jac_eval_frequency=83`, `cvode.lsetup_frequency=32`, `cvode.max_order=2` |
| `cvPendulum_dns` | `cvode.jac_eval_frequency=24`, `cvode.lsetup_frequency=10`, `cvode.max_order=2` |
| `cvRoberts_dns_constraints` | `cvode.jac_eval_frequency=35`, `cvode.lsetup_frequency=30`, `cvode.max_order=4` |
| `cvRoberts_dns` | `cvode.jac_eval_frequency=69`, `cvode.lsetup_frequency=12`, `cvode.max_order=2` |
| `cvVdp_auto_nls` | `cvode.jac_eval_frequency=62`, `cvode.lsetup_frequency=7`, `cvode.max_order=4` |
| `cvsAnalytic_mels` | `cvodes.max_order=2`, `cvodes.nonlin_conv_coef=0.04029919713927839` |
| `cvsParticle_dns` | `cvodes.jac_eval_frequency=21`, `cvodes.lsetup_frequency=6`, `cvodes.max_order=2` |
| `cvsPendulum_dns` | `cvodes.jac_eval_frequency=34`, `cvodes.lsetup_frequency=39`, `cvodes.max_order=2` |
| `cvsRoberts_dns_constraints` | `cvodes.jac_eval_frequency=60`, `cvodes.lsetup_frequency=35`, `cvodes.max_order=4` |
| `cvsRoberts_dns` | `cvodes.jac_eval_frequency=82`, `cvodes.lsetup_frequency=5`, `cvodes.max_order=2` |
| `idaAnalytic_mels` | `ida.max_order=2`, `ida.nonlin_conv_coef=0.28366152731576105` |
| `idaFoodWeb_bnd` | `ida.delta_cj_lsetup=0.2885293038235052`, `ida.max_order=2`, `ida.nonlin_conv_coef=0.11624666320077459` |
| `idaHeat2D_bnd` | `ida.delta_cj_lsetup=0.2840299010615597`, `ida.max_order=2`, `ida.nonlin_conv_coef=0.16538912175916204` |
| `idaRoberts_dns` | `ida.delta_cj_lsetup=0.31946938505313544`, `ida.max_order=2`, `ida.nonlin_conv_coef=0.030086744318393137` |
| `idaSlCrank_dns` | `ida.delta_cj_lsetup=0.20410024670670857`, `ida.max_order=2`, `ida.nonlin_conv_coef=0.0984308904586615` |
| `idasAnalytic_mels` | `idas.max_order=5`, `idas.nonlin_conv_coef=0.07689651807894142` |
| `idasFoodWeb_bnd` | `idas.delta_cj_lsetup=0.10564748304299437`, `idas.max_order=3`, `idas.nonlin_conv_coef=0.04439541988064023` |
| `idasHeat2D_bnd` | `idas.delta_cj_lsetup=0.27947314761312436`, `idas.max_order=5`, `idas.nonlin_conv_coef=0.23380604146559406` |
| `idasRoberts_dns` | `idas.delta_cj_lsetup=0.41827647638561427`, `idas.max_order=2`, `idas.nonlin_conv_coef=0.045042352297884794` |
| `idasSlCrank_dns` | `idas.delta_cj_lsetup=0.29102957599862356`, `idas.max_order=2`, `idas.nonlin_conv_coef=0.18116299812046258` |
| `kinAnalytic_fp` | `kinsol.damping_aa=0.8198036274578507`, `kinsol.delay_aa=0`, `kinsol.m_aa=2` |
| `kinFerTron_dns` | `kinsol.eta_const_value=0.10149145910622288`, `kinsol.eta_form=1` |
| `kinFoodWeb_kry` | `kinsol.eta_const_value=0.1947770701626817`, `kinsol.eta_form=1`, `kinsol.max_setup_calls=11` |
| `kinLaplace_bnd` | `kinsol.max_setup_calls=3`, `kinsol.max_sub_setup_calls=10` |
| `kinRoberts_fp` | `kinsol.damping_aa=0.554404157068713`, `kinsol.delay_aa=0`, `kinsol.m_aa=2` |
| `kinRoboKin_dns` | `kinsol.eta_const_value=0.01272469600224772`, `kinsol.eta_form=2`, `kinsol.max_setup_calls=17` |

## Validation and interpretation

- All 40 sampled winners completed successfully and satisfied their configured
  constraints in three independent validation reruns.
- The deterministic-work winners reproduced their objective exactly on every
  rerun. Wall-time results remain screening measurements because many examples
  finish in only a few milliseconds.
- The largest deterministic reduction was `ark_kpr_mri`: 1,570,314 to
  115,581 RHS evaluations (-92.6%). `ark_analytic_ssprk` fell from 6,489
  to 2,004 RHS evaluations (-69.1%).
- The default remained tied for best for `idaAnalytic_mels`,
  `idasAnalytic_mels`, `kinAnalytic_fp`, and `kinRoberts_fp`.
- Negative worst-case changes mean every feasible sampled setting measured
  better than that run's baseline; they do not eliminate timing noise.

## Comparison with the DUMMY-surrogate campaign

Across the 17 deterministic-objective searches, the default ET
surrogate produced a better best-observed result in 10,
the same result in 5, and a worse result in
2. The two regressions were
`ark_analytic_lsrk_domeigest` (132,126 versus 132,099, a 0.02% difference)
and `kinAnalytic_fp` (5 versus 4 function evaluations). The wall-time
campaigns are not compared as optimizer-quality evidence because their
millisecond-scale objectives vary with system noise.

## Reproduction and environment

Each search was run with:

```sh
export builddir="$PWD/builddir"
export LD_LIBRARY_PATH="<GCC 13 lib64>:$LD_LIBRARY_PATH"
.venv/bin/suntools tune --config <configuration.yaml>
```

Configurations are every file matching `examples/**/*_tune.yaml`. New raw
results are under `/tmp/suntools-tune-*-default-surrogate`; each directory
contains `baseline.json`, `best.json`, `worst.json`, `results.csv`,
`trials.jsonl`, `validation.json`, and DeepHyper logs.

- Git revision: `958e705f361d` plus working-tree tuning changes
- Host: `tux457.llnl.gov`, Linux 4.18.0-553.159.1.el8_10.x86_64
- CPU: Intel Xeon W-2255 at 3.70 GHz, 10 cores / 20 threads
- Build: `RelWithDebInfo`, double precision
- C compiler: Clang 17.0.4
- CMake: 3.30.2
- Tuning backend: DeepHyper 0.13.2, default ET surrogate
- Search execution: one worker; one repetition for deterministic objectives
  and five repetitions for wall-time objectives
