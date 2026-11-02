# ExtSTS: incorporating the STS temporal error estimate into MRI step-size control

Planning document for SUNDIALS issue
[#1075](https://github.com/llnl/sundials/issues/1075), against `develop` at
`8e87ef1e4` (branch `feature/extsts-stserror`).  Paths are relative to the SUNDIALS root.

**Revision (2026-11-02).**  The STS error estimate is now opt-in: it is disabled by
default and enabled with a new public routine `MRIStepSetUseExtSTSErrorEstimate` and the
matching command-line key `arkid.use_extsts_error_estimate` (§4.6, §5).  The original
design, in which the estimate was always on, is implemented in `93d05f4f4`; §10 lists what
must change in that commit and in the uncommitted documentation and answer-file edits.

## 1. Problem

An ExtSTS method (`MRIStepCreateExtSTS`) controls its step size `H` using only the
MRI embedding, `dsm_S = ||y_n - ỹ_n||_WRMS`, which measures the error of the *slow*
operators `f^E`/`f^I`.  The error committed by the inner LSRKStep (RKC/RKL) solves of
the fast IVPs for `f^D` never reaches the controller, because:

1. `extSTSInnerStepper_Evolve` runs LSRKStep with `ark_mem->fixedstep = SUNTRUE`
   ([arkode_mristep_extsts.c:258-260](src/arkode/arkode_mristep_extsts.c#L258-L260)),
   and the RKC/RKL step routines only form the embedding when `!ark_mem->fixedstep`
   ([arkode_lsrkstep.c:969](src/arkode/arkode_lsrkstep.c#L969),
   [:1402](src/arkode/arkode_lsrkstep.c#L1402)).
2. MRIStep only retrieves an inner error estimate for `SUN_ADAPTCONTROLLER_MRI_H_TOL`
   controllers ([arkode_mristep.c:3746-3764](src/arkode/arkode_mristep.c#L3746-L3764)),
   which are not applicable to ExtSTS (§3).

When `f^E`, `f^I` are small, `dsm_S ≈ 0`, `H` grows by `etamax` every step, and the
STS error is uncontrolled (issue reproducer: 2 steps, error 100× the solution).

## 2. Key observations

- ExtSTS does not "evolve" LSRKStep in the usual sense: `extSTSInnerStepper_Evolve`
  resets the inner object and calls `ark_mem->step` once per fast interval.  There is no
  inner adaptivity, no inner tolerance, and the inner `ewt` is unused.
- MRIStep and LSRKStep are both ARKODE modules, and `arkode_mristep_extsts.c` already
  includes `arkode_lsrkstep_impl.h`, so the MRI side can read LSRKStep's memory directly.
- After an RKC/RKL step with the embedding enabled, the local error *vector*
  `0.8(y_n - y_{n+1}) + 0.4 h (f_n + f_{n+1})` is left in the inner `ark_mem->tempv1`
  (this is what `lsrkStep_GetEstLocalErrors` returns).  The end-of-step RHS `f_{n+1}`
  (`tempv2`) is already evaluated unconditionally, so forming the vector costs one
  `N_VLinearCombination` of 4 vectors.  `arkCompleteStep` does not touch `tempv1`:
  ExtSTS sets `ARK_INTERP_NONE`, so the interpolation update (the only vector work in
  `arkCompleteStep`) is skipped.
- The STS error estimate over a fast interval of length `c_i H` is controlled only by
  `H`; it is a local error estimate of the composite step, i.e. the same kind of
  quantity as `dsm_S`.  Like the MRI embedding, the RKC/RKL embedding is one order
  lower than the method itself, so the STS estimate scales as `H^2`.

Added cost per fast solve: one 4-vector linear combination (inner) and one WRMS norm
(outer), versus the `s`-stage STS solve itself (`s` RHS evaluations, `s` ≈ 10–100).

## 3. Why not the H-Tol path

The `MRI_H_TOL` controllers model the fast error as `c(t) H RTOL^F` and adapt `RTOL^F`
for an *adaptive* inner integrator.  ExtSTS has no `RTOL^F`, so the STS estimate must be
combined with `dsm_S` for the ordinary `SUN_ADAPTCONTROLLER_H` controllers.
`MRI_H_TOL` is currently rejected for ExtSTS (the inner stepper has no `setrtol`; see
`mriStepInnerStepper_SupportsRTolAdaptivity`) and must stay rejected; the new scalar
below is deliberately kept separate from `step_mem->inner_dsm`, which feeds the
H-Tol controller wrapper in `arkode_mristep_controller.c`.

## 4. Design

### 4.1 LSRKStep: an "estimate error only" flag

Add `sunbooleantype compute_err_est` to `ARKodeLSRKStepMem` (default `SUNFALSE`;
set in `lsrkStep_Create*`/defaults).  Semantics: compute the embedded error vector
during the step, but leave every other aspect of fixed-step operation unchanged.

Changes in `src/arkode/arkode_lsrkstep.c`:

1. `lsrkStep_TakeStepRKC`, `lsrkStep_TakeStepRKL`: form the embedding when
   `!ark_mem->fixedstep || step_mem->compute_err_est`.  Leave the inner `*dsmPtr` at
   `ZERO` in fixed-step mode: the inner `ewt` is meaningless there (`arkEwtSetSmallReal`),
   the inner step is always accepted, and nothing downstream of the inner step
   (`arkCompleteStep`, inner controller bookkeeping) should change.
2. `lsrkStep_DomEigUpdateLogic` keeps receiving `*dsmPtr` (`ZERO` in fixed-step mode), so
   the dominant-eigenvalue update cadence is unchanged.  The outer method owns step
   rejection.
3. The stage-limit retry logic (`if (!ark_mem->fixedstep)` at lines 639, 707, 1090, 1158)
   is unaffected.
4. The SSP routines need no change (ExtSTS uses RKC/RKL only); the flag is documented as
   RKC/RKL-specific in the impl header.

No public LSRKStep API is added; the flag is set by ExtSTS through the impl header.

### 4.2 MRIStep: a new ExtSTS error accumulator

Add three scalars to `ARKodeMRIStepMem` (initialized to `ZERO` where `inner_dsm` is):

```c
sunrealtype extsts_dsm;     /* STS error estimate for the current step attempt */
sunrealtype extsts_acc;     /* running sum of |dt_k| * ||e_k||^2 over fast solves */
sunrealtype extsts_dt_acc;  /* running sum of |dt_k| over fast solves */
```

and three helpers in `src/arkode/arkode_mristep_extsts.c`, declared in
`arkode_mristep_impl.h`, so `arkode_mristep.c` needs no LSRKStep internals:

```c
/* Enable/disable the STS error estimate and zero the accumulators */
int extSTS_BeginStepError(ARKodeMRIStepMem step_mem, sunbooleantype on);

/* After one fast solve over [t0, tf]: with e_k = inner tempv1 and the OUTER ewt,
   extsts_acc += |tf - t0| * ||e_k||_WRMS^2;  extsts_dt_acc += |tf - t0| */
int extSTS_AccumulateStageError(ARKodeMem ark_mem, ARKodeMRIStepMem step_mem,
                                sunrealtype t0, sunrealtype tf);

/* At the end of the step attempt: extsts_dsm = sqrt(extsts_acc / extsts_dt_acc)
   (ZERO if no fast solve contributed) */
int extSTS_FinalizeStepError(ARKodeMRIStepMem step_mem);
```

`extSTS_AccumulateStageError` reads `((ARKodeMem)step_mem->stepper->content)->tempv1`
and uses `N_VWrmsNorm` with the outer `ark_mem->ewt`, so the STS and MRI estimates are
measured in the same norm with no tolerance forwarding.  The absolute values keep the
weights positive for negative-direction stepping.

Changes in `src/arkode/arkode_mristep.c`:

Below, `use_sts_err` is shorthand for
`step_mem->extsts_method && step_mem->extsts_use_err_est`, where `extsts_use_err_est` is
the user option of §4.6.

1. `mriStep_TakeStepMRIGARK`, `mriStep_TakeStepMRISR`, `mriStep_TakeStepMERK`, at the
   top of each step attempt, when `step_mem->extsts_method`:
   - `extSTS_BeginStepError(step_mem, use_sts_err && do_embedding)` where
     `do_embedding = !ark_mem->fixedstep || ark_mem->AccumErrorType != ARK_ACCUMERROR_NONE`
     (MRIGARK already names it; MRISR/MERK compute this condition inline today, so
     name it there too).  Either way the accumulators are zeroed, and the inner LSRKStep
     `compute_err_est` flag is set to exactly that condition, so an ExtSTS run with the
     option disabled (the default), or a fixed-step ExtSTS run without error
     accumulation, pays nothing and is bit-identical to the pre-#1075 code;
   - `need_inner_dsm = use_sts_err && do_embedding` (keep the
     `ResetAccumulatedError`/`SetRTol` calls inside the `MRI_H_TOL` branch only).
2. `mriStep_StageERKFast`: after a successful `mriStepInnerStepper_Evolve`, when
   `get_inner_dsm && use_sts_err`, call
   `extSTS_AccumulateStageError(ark_mem, step_mem, t0, tf)`.  The existing `MRI_H_TOL`
   retrieval branch is untouched.
3. Embedding stages: the MRI embedding is still computed, but the STS error from the fast
   solve *inside* the embedding stage must not be added to `extsts_dsm`, since it would
   double-count the final fast interval.  No change is needed: MRIGARK computes its
   embedding stage before the final stage and already passes `SUNFALSE` to
   `mriStep_StageERKFast` there
   ([arkode_mristep.c:2186-2187](src/arkode/arkode_mristep.c#L2186-L2187)), while the
   final solution stage passes `need_inner_dsm`
   ([:2266](src/arkode/arkode_mristep.c#L2266)); MRISR and MERK pass
   `need_inner_dsm && !embedding`.
4. After `dsm_S` is formed in each `TakeStep*`, when `use_sts_err`:
   `extSTS_FinalizeStepError(step_mem)`, then
   `*dsmPtr = SUNRsqrt((*dsmPtr) * (*dsmPtr) + extsts_dsm * extsts_dsm)`, logging
   `dsm_S`, `extsts_dsm`, and the combined value.

### 4.3 Accumulation and combination

Write `e_k` for the STS error vector of fast solve `k` of the current step attempt,
`Δt_k = |t_f − t_0|` for the length of its subinterval, and `||·||` for the WRMS norm in
the outer weights.

- **Across the fast solves of one slow step (time-weighted RMS):**

  ```text
  extsts_dsm = sqrt( Σ_k w_k ||e_k||² ),   w_k = Δt_k / Σ_j Δt_j .
  ```

  The `K` stage estimates sample the same quantity (the STS error density along the
  step), so a weighted mean of squares is the appropriate reduction, and the result is a
  weighted WRMS norm of the stacked vector `(e_1, …, e_K)`, bounded by `max_k ||e_k||`.
  Weighting by subinterval length lets longer fast solves contribute more.  Normalizing
  by `Σ_j Δt_j` rather than by `H` keeps the weights summing to 1 for every method
  family: for MRI-GARK the fast solves tile `[t_n, t_n + H]`, so the two agree, but for
  MRISR and MERK the solves all start at `t_n` with lengths `c_i H`, so `Σ_j Δt_j ≠ H`.
  A plain sum of norms was rejected: it is a triangle-inequality bound, not a norm, and
  grows with `K`.

- **Between slow and STS estimates (root-sum-square):**

  ```text
  dsm = sqrt( dsm_S² + extsts_dsm² ).
  ```

  Here `e_S` and `e_D` are two distinct additive components of one local error vector,
  `e = e_S + e_D`, and `||e||² = ||e_S||² + ||e_D||² + 2⟨e_S, e_D⟩_w`; dropping the
  cross term (different operators, different embeddings, no reason to expect correlation)
  gives the RSS.  It has the right limits: when `f^D` is negligible ExtSTS controls
  exactly like MRIStep, and when `f^E`, `f^I` are negligible (the issue's case) it
  controls approximately like LSRKStep alone (exactly, for a single fast solve), which
  is what the reporter compared against.
  A `1/2` inside the root would make both limits `1/√2` too small (steps `2^{1/6} ≈ 12%`
  too long relative to the single-rate methods); `max` has the right limits but is not
  the norm of anything and ignores the smaller component when both are comparable.

- **Caveat (constant factor only).**  For smooth solutions the leading term of the STS
  error estimate `c(t) Δt²` varies slowly, so adjacent-stage errors tend to share a sign
  and `||Σ_k e_k||` is nearer `Σ_k ||e_k||` than the RMS.  That is a factor of at most
  `K ≈ 2–3` in the error constant, absorbed by the controller as a `K^{1/(p+1)}` shift in
  `H`; it does not affect the proportionality between tolerance and achieved error.

- **Controller order.**  The `H` controller is driven with the MRI embedding order
  `hadapt_mem->p = MRIC->p`.  Both embeddings are one order below their methods, so for
  the second-order ExtSTS tables the MRI and STS estimates scale alike (`H^2`).  A single
  exponent for a combined estimate is standard for partitioned methods.  No change.

### 4.4 Rejected steps and retries

Nothing new.  A rejected slow step reduces `H`; the next attempt resets the inner
stepper to `(t_n, y_n)` and `extSTS_BeginStepError` zeroes the accumulators.  The `ARK_RETRY_STEP` path from the
inner stage-limit check is unaffected.

### 4.5 `extSTSInnerStepper_Evolve`

Only cosmetic: the `end-step-attempt` log line hard-codes `dsm = 0`; keep it, since the
inner `dsm` genuinely stays zero, but add a debug log of the stage error norm in
`extSTS_AccumulateStageError`.

### 4.6 User control: `MRIStepSetUseExtSTSErrorEstimate`

New public routine, declared in `include/arkode/arkode_mristep.h` next to
`MRIStepSetCoupling` and defined in `src/arkode/arkode_mristep_io.c`:

```c
SUNDIALS_EXPORT int MRIStepSetUseExtSTSErrorEstimate(void* arkode_mem,
                                                     sunbooleantype onoff);
```

- Stores `onoff` in a new `sunbooleantype extsts_use_err_est` field of
  `ARKodeMRIStepMem`.
- Default `SUNFALSE`, set in `mriStep_SetDefaults`, so `ARKodeSetDefaults` also clears
  it.  The option persists across `MRIStepReInitExtSTS`/`ARKodeReset`, as other
  optional inputs do.
- Returns `ARK_MEM_NULL` for a `NULL` or non-MRIStep memory block (via
  `mriStep_AccessStepMem`), and `ARK_ILL_INPUT` with an error message if the object was
  not created by `MRIStepCreateExtSTS` (`!step_mem->extsts_method`).  Otherwise it
  returns `ARK_SUCCESS`.
- It may be called at any time; it takes effect at the next step attempt, since
  `extSTS_BeginStepError` reads the flag at the start of every attempt.
- `mriStep_WriteParameters` reports the setting for ExtSTS objects only, e.g.
  `"  ExtSTS STS error estimate = enabled|disabled"`.

**Command-line control.**  `mriStep_SetOptions` gains a `sunKeyIntPair` table, following
`lsrkStep_SetOptions`:

```c
static const struct sunKeyIntPair int_pairs[] = {
  {"use_extsts_error_estimate", MRIStepSetUseExtSTSErrorEstimate}};
```

It is processed with `sunCheckAndSetIntArgs` before the existing
`coupling_table_name` branch.  This is the same `key 0|1` form used by booleans such as
`use_compensated_sums`.  On a non-ExtSTS object the key fails with `ARK_ILL_INPUT`,
consistent with the setter.

## 5. API impact

One new public function, `MRIStepSetUseExtSTSErrorEstimate`, plus one new
command-line key, `arkid.use_extsts_error_estimate`.  The default is **disabled**, so
existing ExtSTS users see no change in results: with the option off, every ExtSTS run,
adaptive or fixed-step and with or without error accumulation, is bit-identical to the
pre-#1075 code.  When the option is enabled:
- adaptive ExtSTS runs respond to the STS error;
- fixed-step runs with outer error accumulation include the STS error in
  `ARKodeGetAccumulatedError`.  The fixed-step solution itself is unchanged.

This is a new feature rather than a bug fix.  It is announced under "New Features and
Enhancements" in `CHANGELOG.md` and `doc/shared/RecentChanges.rst`.

Bindings: Fortran (SWIG) and sundials4py are generated from the public headers.
- Regenerate the Fortran module sources (`src/arkode/fmod_int{32,64}/farkode_mristep_mod.*`).
- Regenerate `bindings/sundials4py/arkode/arkode_mristep_generated.hpp`.

Check that the new function is not excluded by `bindings/sundials4py/arkode/generate.yaml`.
It takes no callbacks, so no hand-written wrapper is needed.

`ARKodeGetEstLocalErrors` on the ExtSTS object continues to return the slow embedding
vector only (there is no single combined error vector); document this.

## 6. Documentation

- `doc/arkode/guide/source/Mathematics.rst`, ExtSTS section: a paragraph on temporal
  error estimation, containing:
  - by default only the MRI embedding is used;
  - when enabled with `MRIStepSetUseExtSTSErrorEstimate`, a per-stage STS estimate from
    the RKC/RKL embedding, measured in the MRIStep weights;
  - these estimates are combined over the fast stage solves as a subinterval-weighted
    RMS, then combined with the MRI embedding as a root-sum-square (both formulas from
    §4.3);
  - when to enable it: the STS error may dominate, e.g. when `f^E`, `f^I` are small
    relative to `f^D`;
  - H-Tol controllers are not applicable to ExtSTS.
- `doc/arkode/guide/source/Usage/MRIStep/User_callable.rst`:
  - Note under `MRIStepCreateExtSTS`:
    - tolerances are set on the ExtSTS object only;
    - the STS error estimate is off by default (link to the new function);
    - `ARKodeGetEstLocalErrors` returns the slow component;
    - `MRI_H_TOL` controllers return `ARK_ILL_INPUT`.
  - A row in the "Optional inputs for MRIStep" section
    (`ARKODE.Usage.MRIStep.MRIStepInput`).
  - A full `.. c:function::` entry for `MRIStepSetUseExtSTSErrorEstimate`, containing:
    - its parameters and return values;
    - a reference to the Mathematics section;
    - the note "This routine will be called by :c:func:`ARKodeSetOptions` when using
      the key "arkid.use_extsts_error_estimate"" (the wording used for other options);
    - `.. versionadded:: x.y.z`.
  - Drop the `versionchanged:: x.y.z` block under `MRIStepCreateExtSTS` that the
    uncommitted edit adds: default behavior no longer changes.

## 7. Testing

All tests in `ark_test_extsts_stserror.c` that rely on the STS estimate (1, 2, 4, and 3
for the ExtSTS run) enable it explicitly with `MRIStepSetUseExtSTSErrorEstimate`.

1. **Regression for #1075.**  Turn the issue's reproducer (`extsts_zero_slow.c`: 1D
   diffusion with `f^E = c u_x`, RKL-2, analytic dominant eigenvalue,
   `rtol = atol = 1e-6`, `T = 5`) into
   `test/unit_tests/arkode/C_serial/ark_test_extsts_stserror.c`.  For `c = 0` compare to
   the exact semi-discrete solution and assert the error is a small multiple of the
   tolerance (e.g. `max error < 1e-4`, versus `7.2e-1` with the option off) and the step
   count is within a factor of ~2 of LSRKStep alone (111 steps).  Run `c = 1e-2` to
   confirm the intended regime still behaves.  Run with RKC as well as RKL.
2. **Example: default and opt-in runs.**
   - In `examples/arkode/CXX_serial/ark_adr1d_extsts.cpp`, add
     `ARKodeSetOptions(arkode_mem, NULL, NULL, argc, argv)` (default `arkode.` prefix)
     after all other
     programmatic `ARKodeSet*`/`MRIStepSet*` calls on the ExtSTS object, so command-line
     keys override them.  The example's own `--key value` parser ignores the
     `arkode.`-prefixed arguments, and `ARKodeSetOptions` skips the example's own
     arguments, since they lack the prefix.
   - Document the new key in `InputHelp()`, e.g. a line pointing to
     `arkode.<ARKODE option>` (in particular `arkode.use_extsts_error_estimate 1`).
   - Do not print the setting in `PrintSetup`, so that the default output is unchanged.
   - In `examples/arkode/CXX_serial/CMakeLists.txt`, keep the existing default entry and
     add a second one:

     ```cmake
     "ark_adr1d_extsts.cpp\;\;develop"
     "ark_adr1d_extsts.cpp\;arkode.use_extsts_error_estimate 1\;develop"
     ```

   - Answer files (the second name follows the existing
     `ark_analytic_arkode.scalar_tolerances_…out` precedent):
     - `ark_adr1d_extsts.out`: revert to the committed version (`git checkout`), since
       the default run is again bit-identical to the pre-#1075 code.  The `.hpp` changes
       in `93d05f4f4` (default MRI table with `--no-reaction`, `--build_mri_table` check)
       do not affect the default run.  Verify by rerunning.
     - `ark_adr1d_extsts_arkode.use_extsts_error_estimate_1.out` (new): the content of
       the current modified `ark_adr1d_extsts.out`.  Regenerate it from the new build
       rather than copying, and confirm it matches.
   - The Python example `examples/python/arkode/ark_adr1d_extsts.py` is unchanged.  It
     has no answer file, and its argparse parser does not pass through ARKODE keys.
3. **Fixed-step invariance.**  With the option enabled, `ark_adr1d_extsts` in
   fixed-step mode (`--fixed_h`) without error accumulation must produce bit-identical
   output to the default run.  With error accumulation enabled, only the accumulated
   error estimate may change.
4. **Embedding-stage accounting.**  With the option enabled and logging on, verify that
   the accumulators receive one contribution per non-embedding fast solve (2 for the
   GKC21 tables) and none from the embedding stage, for MRIGARK, MRISR, and MERK tables;
   for MRIGARK check `extsts_dt_acc == H` at finalize.
5. **Limit checks.**  (a) `f^D ≡ 0` (zero diffusion RHS): the ExtSTS step sequence must
   match the same MRI method run as plain MRIStep with an exact/trivial inner stepper,
   confirming `extsts_dsm = 0` leaves `dsm_S` unchanged.  (b) The issue's `c = 0` case
   with a 2-stage table (e.g. GKC21): `dsm` should approximately equal the embedding norm
   from LSRKStep run directly on the same step, confirming the RSS reduces to
   (approximately) the LSRKStep estimate.
6. **H-Tol rejection.**  Attaching an `MRI_H_TOL` controller to an ExtSTS object still
   fails at init with `ARK_ILL_INPUT`.
7. **Option handling** (new unit test in `ark_test_extsts_stserror.c`):
   - (a) *Default off.*  For the `c = 0` problem, an adaptive ExtSTS run that never
     calls the setter and one that calls
     `MRIStepSetUseExtSTSErrorEstimate(mem, SUNFALSE)` must have identical step counts
     and bitwise-identical solutions.  Together with (c), whose enabled run must differ
     from the default, this checks both directions without relying on the #1075
     symptoms.
   - (b) *Setter.*  `MRIStepSetUseExtSTSErrorEstimate` returns:
     - `ARK_ILL_INPUT` on a plain `MRIStepCreate` object;
     - `ARK_MEM_NULL` on `NULL`;
     - `ARK_SUCCESS` on an ExtSTS object.
   - (c) *Command line.*  `ARKodeSetOptions(mem, NULL, NULL, argc, argv)` with
     `arkode.use_extsts_error_estimate 1` gives the same step count and solution as
     calling the setter with `SUNTRUE`, and a different step count from the default
     run in (a).
   - (d) *Reset.*  `ARKodeSetDefaults` turns the option off again.
8. **Full CTest** with `-DSUNDIALS_TEST_ENABLE_UNIT_TESTS=ON` and
   `-DSUNDIALS_TEST_ENABLE_DEV_TESTS=ON`.  Non-ExtSTS tests, including the MRIStep H-Tol
   tests, must be unaffected.  Build the Fortran and Python interfaces with the
   regenerated bindings, and check that `sundials4py` exposes the new function.

## 8. Risks and open points

- **Numerics unchanged by default.**  Only users who opt in see new results.  The risk
  shifts to discoverability: users who hit #1075 must find the option.  So the
  Mathematics and `MRIStepCreateExtSTS` docs say when to enable it.
- **Accumulation across stages** (subinterval-weighted RMS) and **slow/STS
  combination** (root-sum-square) are settled per §4.3; the constant-factor caveat there
  is accepted.
- **Inner `PostStepFn`.**  The outer reads inner `tempv1` after `arkCompleteStep` has
  called the inner post-step hook.  User code has no legitimate access to that
  workspace, so this is a documentation note, not a design constraint.

## 9. Implementation order

1. LSRKStep flag and embedding guard (§4.1).  No output changes.  *(done in
   `93d05f4f4`)*
2. MRIStep accumulator fields, ExtSTS helpers, `TakeStep*`/`StageERKFast` wiring
   (§4.2).  *(done in `93d05f4f4`; needs the `use_sts_err` gating)*
3. Public setter, defaults, `WriteParameters`, and command-line key (§4.6); gate the
   wiring of §4.2 on the option.  Afterwards `ark_adr1d_extsts` must reproduce the
   committed `.out` exactly.
4. Unit-test updates: explicit enabling in the existing tests 1–4 (test 5, H-Tol
   rejection, is unaffected) and the new option tests (§7.1, §7.4–7.7).
5. Example, CMake entry, and the two answer files (§7.2–7.3).
6. Regenerate the Fortran and sundials4py bindings (§5).
7. Docs, CHANGELOG/RecentChanges (§5, §6, §10), then full CTest (§7.8).

## 10. Changes relative to the current branch state

Committed in `93d05f4f4` (always-on design):
- `src/arkode/arkode_mristep_impl.h`: add the `extsts_use_err_est` field.
- `src/arkode/arkode_mristep.c`: gate `extSTS_BeginStepError`, `need_inner_dsm`,
  `extSTS_AccumulateStageError`, and `extSTS_FinalizeStepError`/RSS on `use_sts_err`
  (§4.2).
- `src/arkode/arkode_mristep_io.c`: setter, default, `WriteParameters` line, CLI key.
- `include/arkode/arkode_mristep.h`: prototype.
- `test/unit_tests/arkode/C_serial/ark_test_extsts_stserror.c`: enable the option where
  needed; add §7.7; update the header comment's test list.
- LSRKStep changes (§4.1): unchanged.

Uncommitted edits:
- `CHANGELOG.md`, `doc/shared/RecentChanges.rst`:
  - Move the ExtSTS entry from "Bug Fixes" to "New Features and Enhancements" and
    reword it: "Added `MRIStepSetUseExtSTSErrorEstimate` (command-line key
    `arkid.use_extsts_error_estimate`) to include the STS methods' embedded error
    estimates in the temporal error estimate of ExtSTS methods, for use when the STS
    error dominates.  The option is disabled by default."
  - Keep the `ark_adr1d_extsts` example fix under "Bug Fixes".
- `doc/arkode/guide/source/Mathematics.rst`: make the paragraph conditional on the
  option (§6).
- `doc/arkode/guide/source/Usage/MRIStep/User_callable.rst`: per §6.
- `examples/arkode/CXX_serial/ark_adr1d_extsts.out`: revert; its current content becomes
  the new opt-in answer file (§7.2).
