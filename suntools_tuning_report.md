# SUNDIALS Serial Example Tuning Report

Run date: 2026-09-14

## Scope and method

This report summarizes `suntools tune` searches for forty serial examples
from ARKODE, CVODE, CVODES, IDA, IDAS, and KINSOL. Each search used the YAML
configuration, the existing fixed solver stack and workload, the DeepHyper
backend, and one worker. The initial sixteen-example campaign used
deterministic solver-work metrics and one repetition. The subsequent
twenty-four-example candidate campaign used five-repetition wall time, except
for `cvsPendulum_dns`, which used the same constrained deterministic work
metric as `cvPendulum_dns`.

The reported **best observed** value is the better of the baseline and best
feasible sampled trial. The **worst** value is the worst feasible sampled
trial. A negative best change is an improvement. Objective values are
comparable only within one row because the examples use different work
metrics. Feasible trials completed successfully and satisfied any configured
accuracy constraint. Trials that returned a failure status were excluded.

These are best-observed settings for these workloads, search spaces, build,
and machine. They are not universal defaults.

Together the campaigns evaluated 814 sampled configurations. Of these, 763
were feasible. The initial campaign had sixteen infeasible or failed trials;
the candidate campaign had thirty-five, arising from existing answer checks,
solver failures, or the pendulum accuracy constraint.

## Results

| Package | Example | Objective | Baseline | Best observed | Best change | Sampled worst | Worst change | Feasible trials |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| ARKODE | `ark_analytic` | Solver work | 355 | 310 | -12.7% | 3560 | +902.8% | 40/40 |
| ARKODE | `ark_analytic_lsrk` | RHS evaluations | 142049 | 133332 | -6.1% | 159707 | +12.4% | 24/24 |
| ARKODE | `ark_analytic_lsrk_domeigest` | RHS and estimator work | 132464 | 132099 | -0.3% | 145618 | +9.9% | 32/32 |
| ARKODE | `ark_analytic_ssprk` | RHS evaluations | 6489 | 2004 | -69.1% | 5030 | -22.5% | 12/12 |
| ARKODE | `ark_kpr_mri` | RHS evaluations | 1570314 | 205581 | -86.9% | 4280581 | +172.6% | 18/18 |
| CVODE | `cvAnalytic_mels` | RHS evaluations | 135 | 135 (baseline) | 0.0% | 896 | +563.7% | 16/16 |
| CVODE | `cvPendulum_dns` | Solver work | 3890 | 3547 | -8.8% | 20917 | +437.7% | 34/40 |
| CVODE | `cvRoberts_dns` | Solver work | 872 | 752 | -13.8% | 1196 | +37.2% | 32/32 |
| CVODES | `cvsAnalytic_mels` | RHS evaluations | 135 | 122 | -9.6% | 897 | +564.4% | 16/16 |
| CVODES | `cvsRoberts_dns` | Solver work | 872 | 755 | -13.4% | 1169 | +34.1% | 32/32 |
| IDA | `idaAnalytic_mels` | Residual evaluations | 576 | 576 | 0.0% | 576 | 0.0% | 16/16 |
| IDA | `idaRoberts_dns` | Solver work | 657 | 575 | -12.5% | 1005 | +53.0% | 32/32 |
| IDAS | `idasAnalytic_mels` | Residual evaluations | 576 | 576 | 0.0% | 576 | 0.0% | 16/16 |
| IDAS | `idasRoberts_dns` | Solver work | 657 | 569 | -13.4% | 971 | +47.8% | 32/32 |
| KINSOL | `kinAnalytic_fp` | Function evaluations | 5 | 4 | -20.0% | 8 | +60.0% | 23/24 |
| KINSOL | `kinRoberts_fp` | Function evaluations | 8 | 8 | 0.0% | 14 | +75.0% | 15/24 |

The default CVODE analytic configuration remains better than every sampled
configuration. The IDA and IDAS analytic objectives were flat throughout their
feasible sampled spaces. The default KINSOL Robertson configuration tied the
best sampled result. Every other search found a lower-work setting.

## Best-observed sampled settings

| Example | Parameters |
|---|---|
| `ark_analytic` | `arkode.table_names="ARKODE_TRBDF2_3_3_2 ARKODE_ERK_NONE"` |
| `ark_analytic_lsrk` | `arkode.sts_method_name=ARKODE_LSRK_RKC_2`, `arkode.dom_eig_safety_factor=1.042904416960918`, `arkode.use_analytic_stability_region=1` |
| `ark_analytic_lsrk_domeigest` | `arkode.dom_eig_frequency=28`, `arkode.num_dom_eig_est_preprocess_iters=0`, `arkode.dom_eig_safety_factor=1.0196381357800368`, `sundomeigestimator.max_iters=59`, `sundomeigestimator.rel_tol=0.00465327266169915` |
| `ark_analytic_ssprk` | `arkode.ssp_method_name=ARKODE_LSRK_SSP_S_2`, `arkode.safety_factor=0.8563585214848191` |
| `ark_kpr_mri` | `outer.coupling_table_name=ARKODE_MRI_GARK_IRK21a`, `inner.table_names="ARKODE_DIRK_NONE ARKODE_DORMAND_PRINCE_7_4_5"`, `inner.fixed_step=0.00019257956450201874` |
| `cvAnalytic_mels` | `cvode.max_order=5`, `cvode.nonlin_conv_coef=0.07249570472138886` |
| `cvPendulum_dns` | `cvode.max_order=5`, `cvode.lsetup_frequency=39`, `cvode.jac_eval_frequency=74` |
| `cvRoberts_dns` | `cvode.max_order=5`, `cvode.lsetup_frequency=9`, `cvode.jac_eval_frequency=20` |
| `cvsAnalytic_mels` | `cvodes.max_order=5`, `cvodes.nonlin_conv_coef=0.16398055778872048` |
| `cvsRoberts_dns` | `cvodes.max_order=5`, `cvodes.lsetup_frequency=29`, `cvodes.jac_eval_frequency=28` |
| `idaAnalytic_mels` | `ida.max_order=3`, `ida.nonlin_conv_coef=0.08548916476509041` |
| `idaRoberts_dns` | `ida.max_order=4`, `ida.delta_cj_lsetup=0.25928084111462`, `ida.nonlin_conv_coef=0.25646078325217225` |
| `idasAnalytic_mels` | `idas.max_order=5`, `idas.nonlin_conv_coef=0.07926830518243651` |
| `idasRoberts_dns` | `idas.max_order=4`, `idas.delta_cj_lsetup=0.2212323365128911`, `idas.nonlin_conv_coef=0.15872689693217773` |
| `kinAnalytic_fp` | `kinsol.m_aa=4`, `kinsol.delay_aa=1`, `kinsol.damping_aa=0.9934936509521524` |
| `kinRoberts_fp` | `kinsol.m_aa=2`, `kinsol.delay_aa=0`, `kinsol.damping_aa=0.896804568504799` |

Every best sampled command was rerun three times after the search. Each rerun
completed successfully and reproduced the reported objective exactly.

## Worst feasible sampled settings

| Example | Parameters |
|---|---|
| `ark_analytic` | `arkode.table_names="ARKODE_CASH_5_2_4 ARKODE_ERK_NONE"` |
| `ark_analytic_lsrk` | `arkode.sts_method_name=ARKODE_LSRK_RKL_2`, `arkode.dom_eig_safety_factor=1.1931554123977857`, `arkode.use_analytic_stability_region=1` |
| `ark_analytic_lsrk_domeigest` | `arkode.dom_eig_frequency=50`, `arkode.num_dom_eig_est_preprocess_iters=4`, `arkode.dom_eig_safety_factor=1.1859629599558883`, `sundomeigestimator.max_iters=54`, `sundomeigestimator.rel_tol=0.00896592074893581` |
| `ark_analytic_ssprk` | `arkode.ssp_method_name=ARKODE_LSRK_SSP_S_3`, `arkode.safety_factor=0.8515804440455831` |
| `ark_kpr_mri` | `outer.coupling_table_name=ARKODE_MRI_GARK_IRK21a`, `inner.table_names="ARKODE_DIRK_NONE ARKODE_ARK436L2SA_ERK_6_3_4"`, `inner.fixed_step=7.051367978584598e-06` |
| `cvAnalytic_mels` | `cvode.max_order=2`, `cvode.nonlin_conv_coef=0.04108125367487805` |
| `cvPendulum_dns` | `cvode.max_order=2`, `cvode.lsetup_frequency=19`, `cvode.jac_eval_frequency=49` |
| `cvRoberts_dns` | `cvode.max_order=5`, `cvode.lsetup_frequency=22`, `cvode.jac_eval_frequency=94` |
| `cvsAnalytic_mels` | `cvodes.max_order=2`, `cvodes.nonlin_conv_coef=0.03844466595434592` |
| `cvsRoberts_dns` | `cvodes.max_order=2`, `cvodes.lsetup_frequency=9`, `cvodes.jac_eval_frequency=63` |
| `idaAnalytic_mels` | All feasible samples tied at 576 residual evaluations. |
| `idaRoberts_dns` | `ida.max_order=2`, `ida.delta_cj_lsetup=0.3781141230203494`, `ida.nonlin_conv_coef=0.09102443473874267` |
| `idasAnalytic_mels` | All feasible samples tied at 576 residual evaluations. |
| `idasRoberts_dns` | `idas.max_order=2`, `idas.delta_cj_lsetup=0.47346458666003965`, `idas.nonlin_conv_coef=0.10334356353655899` |
| `kinAnalytic_fp` | `kinsol.m_aa=2`, `kinsol.delay_aa=0`, `kinsol.damping_aa=0.7005573726807526` |
| `kinRoberts_fp` | `kinsol.m_aa=1`, `kinsol.delay_aa=3`, `kinsol.damping_aa=0.5960032443928074` |

## Interpretation

- The largest work reductions were 86.9% for the multirate KPR example and
  69.1% for the SSPRK example. Both winners remained within their explicit
  error constraints.
- The DIRK-table search favored TR-BDF2 and reduced its work proxy by 12.7%.
  The dominant-eigenvalue estimator search improved only 0.3%, which is too
  small to motivate changing defaults without broader validation.
- The Robertson integration examples consistently benefit from the sampled
  order/setup controls, reducing their work proxy by 12.5% to 13.8%.
- The pendulum search reduced aggregate projected and unprojected solver work
  by 8.8% at the default `rtol=atol=1e-5`. Its constraint sums the ten printed
  absolute solution and constraint errors and requires that sum to be at most
  8.0; six sampled configurations exceeded that bound.
- Capping BDF at order two was particularly costly in the analytic CVODE and
  CVODES searches and in both Robertson DAE searches for this workload.
- `idaAnalytic_mels` and `idasAnalytic_mels` remained at 343 steps and 576
  residual evaluations throughout the sampled space. Their current objective
  cannot distinguish the tested settings, so no tuned value should replace the
  defaults based on this experiment.
- Anderson acceleration reduced `kinAnalytic_fp` from five to four function
  evaluations. For `kinRoberts_fp`, the default remained tied for best and nine
  sampled configurations failed validation, indicating a narrower follow-up
  search would be prudent before changing defaults.

## Additional candidate results

These searches cover the examples identified in the follow-up serial-example
audit. Times are means of five executions. Most workloads complete in only a
few milliseconds, so these values are screening measurements rather than
recommendations to change example defaults. The selected commands were rerun
independently three times and all completed successfully.

| Package | Example | Baseline | Best observed | Best change | Sampled worst | Worst change | Feasible trials |
|---|---|---:|---:|---:|---:|---:|---:|
| ARKODE | `ark_analytic_mels` | 0.001369 s | 0.001300 s | -5.0% | 0.001433 s | +4.7% | 7/16 |
| ARKODE | `ark_analytic_lsrk_varjac` | 0.012206 s | 0.010680 s | -12.5% | 0.011690 s | -4.2% | 16/16 |
| ARKODE | `ark_robertson` | 0.003755 s | 0.002181 s | -41.9% | 0.003187 s | -15.1% | 7/16 |
| ARKODE | `ark_robertson_constraints` | 0.004363 s | 0.001822 s | -58.2% | 0.003012 s | -31.0% | 12/16 |
| CVODE | `cvRoberts_dns_constraints` | 0.001419 s | 0.001289 s | -9.1% | 0.001489 s | +5.0% | 16/16 |
| CVODE | `cvParticle_dns` | 0.002910 s | 0.002653 s | -8.9% | 0.009495 s | +226.2% | 16/16 |
| CVODES | `cvsRoberts_dns_constraints` | 0.001583 s | 0.001298 s | -18.0% | 0.001672 s | +5.6% | 16/16 |
| CVODES | `cvsParticle_dns` | 0.003052 s | 0.002757 s | -9.6% | 0.010019 s | +228.3% | 16/16 |
| CVODES | `cvsPendulum_dns` | 3890 work | 3727 work | -4.2% | 24194 work | +522.0% | 30/40 |
| IDA | `idaSlCrank_dns` | 0.001659 s | 0.001427 s | -14.0% | 0.001668 s | +0.5% | 16/16 |
| IDA | `idaHeat2D_bnd` | 0.001738 s | 0.001530 s | -12.0% | 0.001743 s | +0.2% | 16/16 |
| IDA | `idaFoodWeb_bnd` | 0.069875 s | 0.054631 s | -21.8% | 0.087385 s | +25.1% | 16/16 |
| IDAS | `idasSlCrank_dns` | 0.001687 s | 0.001608 s | -4.6% | 0.002436 s | +44.4% | 16/16 |
| IDAS | `idasHeat2D_bnd` | 0.001717 s | 0.001528 s | -11.0% | 0.001753 s | +2.1% | 16/16 |
| IDAS | `idasFoodWeb_bnd` | 0.070035 s | 0.051262 s | -26.8% | 0.087682 s | +25.2% | 16/16 |
| KINSOL | `kinFerTron_dns` | 0.001218 s | 0.001069 s | -12.2% | 0.001196 s | -1.8% | 16/16 |
| KINSOL | `kinRoboKin_dns` | 0.002780 s | 0.001123 s | -59.6% | 0.001301 s | -53.2% | 16/16 |
| KINSOL | `kinLaplace_bnd` | 0.002350 s | 0.002133 s | -9.2% | 0.002748 s | +17.0% | 16/16 |
| KINSOL | `kinFoodWeb_kry` | 0.005523 s | 0.003604 s | -34.7% | 0.005487 s | -0.6% | 16/16 |
| ARKODE | `ark_heat1D` | 0.012684 s | 0.011545 s | -9.0% | 0.016072 s | +26.7% | 16/16 |
| ARKODE | `ark_brusselator_lsrk_domeigest` | 0.003169 s | 0.003032 s | -4.3% | 0.003323 s | +4.9% | 16/16 |
| ARKODE | `ark_brusselator_lsrk_externaldomeigest` | 0.003303 s | 0.003157 s | -4.4% | 0.003577 s | +8.3% | 16/16 |
| CVODE | `cvVdp_auto_nls` | 0.002280 s | 0.001344 s | -41.1% | 0.011978 s | +425.4% | 16/16 |
| CVODE | `cv_kpr` | 0.002942 s | 0.002581 s | -12.2% | 0.003369 s | +14.5% | 13/16 |

The negative sampled-worst changes indicate a noisy baseline rather than a
guarantee that every sampled setting is faster. In particular, sub-10 ms
results should be confirmed with a longer representative workload before
using the sampled settings as defaults.

### Best observed candidate settings

| Example | Parameters |
|---|---|
| `ark_analytic_mels` | `arkode.order=3`, `arkode.safety_factor=0.8218310080658441` |
| `ark_analytic_lsrk_varjac` | `arkode.dom_eig_frequency=35`, `arkode.dom_eig_safety_factor=1.0144721343741945` |
| `ark_robertson` | `arkode.jac_eval_frequency=58`, `arkode.lsetup_frequency=36`, `arkode.order=4` |
| `ark_robertson_constraints` | `arkode.jac_eval_frequency=59`, `arkode.lsetup_frequency=30`, `arkode.order=2` |
| `cvRoberts_dns_constraints` | `cvode.jac_eval_frequency=52`, `cvode.lsetup_frequency=39`, `cvode.max_order=3` |
| `cvParticle_dns` | `cvode.jac_eval_frequency=34`, `cvode.lsetup_frequency=28`, `cvode.max_order=5` |
| `cvsRoberts_dns_constraints` | `cvodes.jac_eval_frequency=28`, `cvodes.lsetup_frequency=26`, `cvodes.max_order=3` |
| `cvsParticle_dns` | `cvodes.jac_eval_frequency=87`, `cvodes.lsetup_frequency=27`, `cvodes.max_order=5` |
| `cvsPendulum_dns` | `cvodes.jac_eval_frequency=43`, `cvodes.lsetup_frequency=32`, `cvodes.max_order=5` |
| `idaSlCrank_dns` | `ida.delta_cj_lsetup=0.30335962003578143`, `ida.max_order=4`, `ida.nonlin_conv_coef=0.15771736784605336` |
| `idaHeat2D_bnd` | `ida.delta_cj_lsetup=0.36383954458795087`, `ida.max_order=4`, `ida.nonlin_conv_coef=0.1908534702666119` |
| `idaFoodWeb_bnd` | `ida.delta_cj_lsetup=0.4306394881927077`, `ida.max_order=3`, `ida.nonlin_conv_coef=0.19121134025364464` |
| `idasSlCrank_dns` | `idas.delta_cj_lsetup=0.33271468535105575`, `idas.max_order=4`, `idas.nonlin_conv_coef=0.15770878414355732` |
| `idasHeat2D_bnd` | `idas.delta_cj_lsetup=0.46795086318053614`, `idas.max_order=4`, `idas.nonlin_conv_coef=0.2957092663888814` |
| `idasFoodWeb_bnd` | `idas.delta_cj_lsetup=0.46019042629604523`, `idas.max_order=5`, `idas.nonlin_conv_coef=0.1848333979202713` |
| `kinFerTron_dns` | `kinsol.eta_const_value=0.10180309394016727`, `kinsol.eta_form=2` |
| `kinRoboKin_dns` | `kinsol.eta_const_value=0.2578731866872149`, `kinsol.eta_form=3`, `kinsol.max_setup_calls=17` |
| `kinLaplace_bnd` | `kinsol.max_setup_calls=22`, `kinsol.max_sub_setup_calls=10` |
| `kinFoodWeb_kry` | `kinsol.eta_const_value=0.30349888989492`, `kinsol.eta_form=3`, `kinsol.max_setup_calls=3` |
| `ark_heat1D` | `arkode.lsetup_frequency=32`, `arkode.order=3` |
| `ark_brusselator_lsrk_domeigest` | `arkode.dom_eig_frequency=30`, `arkode.dom_eig_safety_factor=1.0451862616076188` |
| `ark_brusselator_lsrk_externaldomeigest` | `arkode.dom_eig_frequency=50`, `arkode.dom_eig_safety_factor=1.0603950089979146` |
| `cvVdp_auto_nls` | `cvode.jac_eval_frequency=51`, `cvode.lsetup_frequency=27`, `cvode.max_order=3` |
| `cv_kpr` | `cvode.jac_eval_frequency=93`, `cvode.lsetup_frequency=35`, `cvode.max_order=5` |

## Reproduction and environment

Each search was run with:

```sh
export builddir="$PWD/builddir"
.venv/bin/suntools tune --config <configuration.yaml>
```

Configurations:

- `examples/arkode/C_serial/ark_analytic_tune.yaml`
- `examples/arkode/C_serial/ark_analytic_lsrk_tune.yaml`
- `examples/arkode/C_serial/ark_analytic_lsrk_domeigest_tune.yaml`
- `examples/arkode/C_serial/ark_analytic_ssprk_tune.yaml`
- `examples/arkode/C_serial/ark_kpr_mri_tune.yaml`
- `examples/cvode/serial/cvAnalytic_mels_tune.yaml`
- `examples/cvode/serial/cvPendulum_dns_tune.yaml`
- `examples/cvode/serial/cvRoberts_dns_tune.yaml`
- `examples/cvodes/serial/cvsAnalytic_mels_tune.yaml`
- `examples/cvodes/serial/cvsRoberts_dns_tune.yaml`
- `examples/ida/serial/idaAnalytic_mels_tune.yaml`
- `examples/ida/serial/idaRoberts_dns_tune.yaml`
- `examples/idas/serial/idasAnalytic_mels_tune.yaml`
- `examples/idas/serial/idasRoberts_dns_tune.yaml`
- `examples/kinsol/serial/kinAnalytic_fp_tune.yaml`
- `examples/kinsol/serial/kinRoberts_fp_tune.yaml`

Additional candidate configurations:

- `examples/arkode/C_serial/ark_analytic_mels_tune.yaml`
- `examples/arkode/C_serial/ark_analytic_lsrk_varjac_tune.yaml`
- `examples/arkode/C_serial/ark_robertson_tune.yaml`
- `examples/arkode/C_serial/ark_robertson_constraints_tune.yaml`
- `examples/arkode/C_serial/ark_heat1D_tune.yaml`
- `examples/arkode/C_serial/ark_brusselator_lsrk_domeigest_tune.yaml`
- `examples/arkode/C_serial/ark_brusselator_lsrk_externaldomeigest_tune.yaml`
- `examples/cvode/serial/cvRoberts_dns_constraints_tune.yaml`
- `examples/cvode/serial/cvParticle_dns_tune.yaml`
- `examples/cvode/serial/cvVdp_auto_nls_tune.yaml`
- `examples/cvode/CXX_serial/cv_kpr_tune.yaml`
- `examples/cvodes/serial/cvsRoberts_dns_constraints_tune.yaml`
- `examples/cvodes/serial/cvsParticle_dns_tune.yaml`
- `examples/cvodes/serial/cvsPendulum_dns_tune.yaml`
- `examples/ida/serial/idaSlCrank_dns_tune.yaml`
- `examples/ida/serial/idaHeat2D_bnd_tune.yaml`
- `examples/ida/serial/idaFoodWeb_bnd_tune.yaml`
- `examples/idas/serial/idasSlCrank_dns_tune.yaml`
- `examples/idas/serial/idasHeat2D_bnd_tune.yaml`
- `examples/idas/serial/idasFoodWeb_bnd_tune.yaml`
- `examples/kinsol/serial/kinFerTron_dns_tune.yaml`
- `examples/kinsol/serial/kinRoboKin_dns_tune.yaml`
- `examples/kinsol/serial/kinLaplace_bnd_tune.yaml`
- `examples/kinsol/serial/kinFoodWeb_kry_tune.yaml`

Environment:

- Git revision: `c6d5d9fc4f35` plus the working-tree tuning changes
- Host: Linux 4.18.0-553.159.1.el8_10.x86_64, x86_64
- CPU: Intel Xeon W-2255 at 3.70 GHz
- Build: `RelWithDebInfo`, double precision
- C compiler: Clang 17.0.4
- CMake: 3.30.2
- Tuning backend: DeepHyper 0.13.2
- Search execution: one worker; one repetition for deterministic-work searches
  and five repetitions for candidate wall-time searches
- The `cv_kpr` run used the GCC 13 runtime library matching the configured
  Clang toolchain through `LD_LIBRARY_PATH`.

The raw result directories were written under `/tmp/suntools-tune-*-yaml` and
contain `baseline.json`, `best.json`, `worst.json`, `results.csv`, and
`trials.jsonl` for each search.
