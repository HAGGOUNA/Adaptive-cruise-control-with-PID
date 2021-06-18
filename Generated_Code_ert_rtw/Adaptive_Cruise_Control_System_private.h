/*
 * File: Adaptive_Cruise_Control_System_private.h
 *
 * Code generated for Simulink model 'Adaptive_Cruise_Control_System'.
 *
 * Model version                  : 1.28
 * Simulink Coder version         : 8.11 (R2016b) 25-Aug-2016
 * C/C++ source code generated on : Fri Jun 18 13:49:13 2021
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_Adaptive_Cruise_Control_System_private_h_
#define RTW_HEADER_Adaptive_Cruise_Control_System_private_h_
#include "rtwtypes.h"
#include "zero_crossing_types.h"
#include "Adaptive_Cruise_Control_System.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmSetFirstInitCond
# define rtmSetFirstInitCond(rtm, val) ((rtm)->Timing.firstInitCondFlag = (val))
#endif

#ifndef rtmIsFirstInitCond
# define rtmIsFirstInitCond(rtm)       ((rtm)->Timing.firstInitCondFlag)
#endif

#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

extern real_T rt_powd_snf(real_T u0, real_T u1);
extern real_T look1_binlcpw(real_T u0, const real_T bp0[], const real_T table[], uint32_T maxIndex);
extern void Adaptive_Cruise_SimpleMagicTire(real_T rtu_Re, real_T rtu_Fz, real_T rtu_omega, real_T rtu_Vx, real_T rtu_lam_mux, real_T rtu_D, real_T rtu_C, real_T rtu_B, real_T rtu_E, const real_T rtu_kappaFx[3], const real_T rtu_FzFx[3], const real_T rtu_FxMap[9], real_T rtu_press, real_T rtu_FNOMIN, real_T rtu_NOMPRES, real_T rtu_QSY1, real_T rtu_QSY2, real_T rtu_QSY3, real_T rtu_QSY4, real_T rtu_QSY5, real_T rtu_QSY6, real_T rtu_QSY7, real_T rtu_QSY8, real_T rtu_gamma, real_T rtu_lam_My, real_T rtu_UNLOADED_RADIUS, real_T rtu_PRESMIN, real_T rtu_PRESMAX, const real_T rtu_VxMy[3], const real_T rtu_FzMy[3], const real_T rtu_MyMap[9], real_T rtu_FxType, real_T rtu_rollingType, B_SimpleMagicTire_Adaptive_Cr_T *localB, real_T rtp_FZMAX, real_T rtp_VXLOW, real_T rtp_kappamax);
extern void Adaptive_Cruise_Control__Locked(RT_MODEL_Adaptive_Cruise_Cont_T * const Adaptive_Cruise_Control_Syst_M, real_T *rty_locked_wout);

/* private model entry point functions */
extern void Adaptive_Cruise_Control_System_derivatives(void);

#endif                                 /* RTW_HEADER_Adaptive_Cruise_Control_System_private_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
