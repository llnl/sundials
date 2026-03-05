#ifndef SUNDIALS_NONLINSOLAUTO_H_
#define SUNDIALS_NONLINSOLAUTO_H_

#include <sundials/sundials_context.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#ifdef __cplusplus
extern "C" {
#endif

enum SUNNonlinSolAutoType
{
  SUNNONLINSOL_AUTO_FIXEDPOINT = 0,
  SUNNONLINSOL_AUTO_NEWTON     = 1
};

#ifndef SWIG
typedef enum SUNNonlinSolAutoType SUNNonlinSolAutoType;
#endif

struct _SUNNonlinearSolverContent_Auto
{
  SUNNonlinSolAutoType type;
  SUNNonlinearSolver fp_solver;
  SUNNonlinearSolver newton_solver;
  SUNNonlinSolConvTestFn user_ctest_fn;
  void* user_ctest_data;
  int maxiters;
  int curiter;
  long int fp_to_newt_delay;
  long int newt_to_fp_delay;
  long int nsolves_since_switch;
  sunrealtype newt_to_fp_threshold;
  sunrealtype fp_to_newt_threshold;
  long int niters;
  long int nconvfails;
  long int switch_count;
  void* auto_ctest_data;
};

typedef struct _SUNNonlinearSolverContent_Auto* SUNNonlinearSolverContent_Auto;

SUNDIALS_EXPORT
SUNNonlinearSolver SUNNonlinSol_Auto(N_Vector y, int m, SUNNonlinSolAutoType type,
                                     SUNContext sunctx);

SUNDIALS_EXPORT
SUNNonlinearSolver_Type SUNNonlinSolGetType_Auto(SUNNonlinearSolver NLS);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolInitialize_Auto(SUNNonlinearSolver NLS);

SUNDIALS_EXPORT
int SUNNonlinSolSolve_Auto(SUNNonlinearSolver NLS, N_Vector y0, N_Vector ycor,
                           N_Vector w, sunrealtype tol,
                           sunbooleantype callSetup, void* mem);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolFree_Auto(SUNNonlinearSolver NLS);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetSysFn_Auto(SUNNonlinearSolver NLS,
                                     SUNNonlinSolSysFn SysFn);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetSysFns_Auto(SUNNonlinearSolver NLS,
                                      SUNNonlinSolSysFn root_fn,
                                      SUNNonlinSolSysFn fixed_point_fn);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetConvTestFn_Auto(SUNNonlinearSolver NLS,
                                          SUNNonlinSolConvTestFn CTestFn,
                                          void* ctest_data);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetLSetupFn_Auto(SUNNonlinearSolver NLS,
                                        SUNNonlinSolLSetupFn LSetupFn);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetLSolveFn_Auto(SUNNonlinearSolver NLS,
                                        SUNNonlinSolLSolveFn LSolveFn);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetMaxIters_Auto(SUNNonlinearSolver NLS, int maxiters);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolSetSwitchingParameters_Auto(
  SUNNonlinearSolver NLS, sunrealtype newt_to_fp_threshold,
  long int newt_to_fp_delay, sunrealtype fp_to_newt_threshold,
  long int fp_to_newt_delay);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetNumIters_Auto(SUNNonlinearSolver NLS, long int* niters);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetCurIter_Auto(SUNNonlinearSolver NLS, int* iter);

SUNDIALS_EXPORT
SUNErrCode SUNNonlinSolGetNumConvFails_Auto(SUNNonlinearSolver NLS,
                                            long int* nconvfails);

#ifdef __cplusplus
}
#endif

#endif /* SUNDIALS_NONLINSOLAUTO_H_ */
