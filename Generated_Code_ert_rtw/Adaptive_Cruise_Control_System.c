/*
 * File: Adaptive_Cruise_Control_System.c
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

#include "Adaptive_Cruise_Control_System.h"
#include "Adaptive_Cruise_Control_System_private.h"

/* Block signals (auto storage) */
B_Adaptive_Cruise_Control_Sys_T Adaptive_Cruise_Control_Syste_B;

/* Continuous states */
X_Adaptive_Cruise_Control_Sys_T Adaptive_Cruise_Control_Syste_X;

/* Block states (auto storage) */
DW_Adaptive_Cruise_Control_Sy_T Adaptive_Cruise_Control_Syst_DW;

/* Previous zero-crossings (trigger) states */
PrevZCX_Adaptive_Cruise_Contr_T Adaptive_Cruise_Control_PrevZCX;

/* Real-time model */
RT_MODEL_Adaptive_Cruise_Cont_T Adaptive_Cruise_Control_Syst_M_;
RT_MODEL_Adaptive_Cruise_Cont_T *const Adaptive_Cruise_Control_Syst_M = &Adaptive_Cruise_Control_Syst_M_;

/* Forward declaration for local functions */
static void Adaptive_Cruise_Control_S_power(const real_T a_data[], const int32_T a_sizes[2], real_T y_data[], int32_T y_sizes[2]);
static void Adaptive_Cruis_automltirelongFx(real_T Re, real_T Fz, real_T omega, real_T Vx, real_T lam_mux, real_T D, real_T C, real_T B, real_T E, real_T b_FZMAX, real_T b_VXLOW, real_T b_kappamax, real_T *Fx, real_T *kappa);
static real_T Adaptive_Cruise_Control_interp2(const real_T varargin_1[3], const real_T varargin_2[3], const real_T varargin_3[9], real_T varargin_4, real_T varargin_5);
static void Adaptive_automltirelongFxMapped(real_T Re, real_T Fz, real_T omega, real_T Vx, const real_T kappaFx[3], const real_T FzFx[3], const real_T FxMap[9], real_T b_FZMAX, real_T b_VXLOW, real_T b_kappamax, real_T *Fx, real_T *kappa);
static real_T Adaptive_Cruis_automltirelongMy(real_T Fx, real_T Fz, real_T omega, real_T Vx, real_T press, real_T FNOMIN, real_T NOMPRES, real_T QSY1, real_T QSY2, real_T QSY3, real_T QSY4, real_T QSY5, real_T QSY6, real_T QSY7, real_T QSY8, real_T b_gamma, real_T lam_My, real_T UNLOADED_RADIUS, real_T b_FZMAX, real_T PRESMIN, real_T PRESMAX);
real_T look1_binlcpw(real_T u0, const real_T bp0[], const real_T table[], uint32_T maxIndex)
{
  real_T frac;
  uint32_T iRght;
  uint32_T iLeft;
  uint32_T bpIdx;

  /* Lookup 1-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear'
     Extrapolation method: 'Clip'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Clip'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = 0.0;
  } else if (u0 < bp0[maxIndex]) {
    /* Binary Search */
    bpIdx = maxIndex >> 1U;
    iLeft = 0U;
    iRght = maxIndex;
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex - 1U;
    frac = 1.0;
  }

  /* Interpolation 1-D
     Interpolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  return (table[iLeft + 1U] - table[iLeft]) * frac + table[iLeft];
}

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 15;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  Adaptive_Cruise_Control_System_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  Adaptive_Cruise_Control_System_step();
  Adaptive_Cruise_Control_System_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  Adaptive_Cruise_Control_System_step();
  Adaptive_Cruise_Control_System_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Function for MATLAB Function: '<S12>/Simple Magic Tire' */
static void Adaptive_Cruise_Control_S_power(const real_T a_data[], const int32_T a_sizes[2], real_T y_data[], int32_T y_sizes[2])
{
  int32_T k;
  y_sizes[0] = 1;
  y_sizes[1] = (int8_T)a_sizes[1];
  k = 1;
  while (k <= a_sizes[1]) {
    y_data[0] = a_data[0] * a_data[0];
    k = 2;
  }
}

/* Function for MATLAB Function: '<S12>/Simple Magic Tire' */
static void Adaptive_Cruis_automltirelongFx(real_T Re, real_T Fz, real_T omega, real_T Vx, real_T lam_mux, real_T D, real_T C, real_T B, real_T E, real_T b_FZMAX, real_T b_VXLOW, real_T b_kappamax, real_T *Fx, real_T *kappa)
{
  real_T Vxpabs;
  real_T d;
  int32_T d_trueCount;
  int32_T i;
  real_T y_data;
  real_T Vxpabs_data;
  int32_T Vxpabs_sizes[2];
  real_T tmp_data;
  int32_T tmp_sizes[2];
  real_T b_idx_0;
  Vxpabs = fabs(Vx);
  d_trueCount = 0;
  if (Vxpabs < b_VXLOW) {
    d_trueCount = 1;
  }

  Vxpabs_sizes[0] = 1;
  Vxpabs_sizes[1] = d_trueCount;
  for (i = 0; i < d_trueCount; i++) {
    Vxpabs_data = Vxpabs / b_VXLOW;
  }

  Adaptive_Cruise_Control_S_power(&Vxpabs_data, Vxpabs_sizes, &tmp_data, tmp_sizes);
  for (i = 0; i < d_trueCount; i++) {
    y_data = 2.0 * b_VXLOW / (3.0 - tmp_data);
  }

  if (Vxpabs < b_VXLOW) {
    Vxpabs = y_data;
  }

  b_idx_0 = Fz;
  if (Fz < 0.0) {
    b_idx_0 = 0.0;
  }

  if (b_idx_0 > b_FZMAX) {
    b_idx_0 = b_FZMAX;
  }

  *kappa = (Re * omega - Vx) / Vxpabs;
  d = *kappa;
  d_trueCount = 0;
  if (*kappa < -b_kappamax) {
    d_trueCount = 1;
  }

  for (i = 0; i < d_trueCount; i++) {
    d = -b_kappamax;
  }

  *kappa = d;
  if (d > b_kappamax) {
    *kappa = b_kappamax;
  }

  *Fx = sin(atan(B * *kappa - (B * *kappa - atan(B * *kappa)) * E) * C) * D * (b_idx_0 * lam_mux);
}

/* Function for MATLAB Function: '<S12>/Simple Magic Tire' */
static real_T Adaptive_Cruise_Control_interp2(const real_T varargin_1[3], const real_T varargin_2[3], const real_T varargin_3[9], real_T varargin_4, real_T varargin_5)
{
  real_T Vq;
  int32_T low_i;
  int32_T low_ip1;
  int32_T high_i;
  int32_T b_high_i;
  real_T qx1;
  real_T rx;
  if ((varargin_4 >= varargin_1[0]) && (varargin_4 <= varargin_1[2]) && (varargin_5 >= varargin_2[0]) && (varargin_5 <= varargin_2[2])) {
    low_i = 0;
    low_ip1 = 2;
    high_i = 3;
    while (high_i > low_ip1) {
      if (varargin_4 >= varargin_1[1]) {
        low_i = 1;
        low_ip1 = 3;
      } else {
        high_i = 2;
      }
    }

    low_ip1 = 0;
    high_i = 2;
    b_high_i = 3;
    while (b_high_i > high_i) {
      if (varargin_5 >= varargin_2[1]) {
        low_ip1 = 1;
        high_i = 3;
      } else {
        b_high_i = 2;
      }
    }

    if (varargin_4 == varargin_1[low_i]) {
      qx1 = varargin_3[3 * low_i + low_ip1];
      Vq = varargin_3[(3 * low_i + low_ip1) + 1];
    } else if (varargin_1[low_i + 1] == varargin_4) {
      qx1 = varargin_3[(low_i + 1) * 3 + low_ip1];
      Vq = varargin_3[((low_i + 1) * 3 + low_ip1) + 1];
    } else {
      rx = (varargin_4 - varargin_1[low_i]) / (varargin_1[low_i + 1] - varargin_1[low_i]);
      if (varargin_3[(low_i + 1) * 3 + low_ip1] == varargin_3[3 * low_i + low_ip1]) {
        qx1 = varargin_3[3 * low_i + low_ip1];
      } else {
        qx1 = varargin_3[(low_i + 1) * 3 + low_ip1] * rx + varargin_3[3 * low_i + low_ip1] * (1.0 - rx);
      }

      if (varargin_3[((low_i + 1) * 3 + low_ip1) + 1] == varargin_3[(3 * low_i + low_ip1) + 1]) {
        Vq = varargin_3[(3 * low_i + low_ip1) + 1];
      } else {
        Vq = varargin_3[((low_i + 1) * 3 + low_ip1) + 1] * rx + varargin_3[(3 * low_i + low_ip1) + 1] * (1.0 - rx);
      }
    }

    if ((varargin_5 == varargin_2[low_ip1]) || (qx1 == Vq)) {
      Vq = qx1;
    } else {
      if (!(varargin_2[low_ip1 + 1] == varargin_5)) {
        rx = (varargin_5 - varargin_2[low_ip1]) / (varargin_2[low_ip1 + 1] - varargin_2[low_ip1]);
        Vq = (1.0 - rx) * qx1 + rx * Vq;
      }
    }
  } else {
    Vq = 0.0;
  }

  return Vq;
}

/* Function for MATLAB Function: '<S12>/Simple Magic Tire' */
static void Adaptive_automltirelongFxMapped(real_T Re, real_T Fz, real_T omega, real_T Vx, const real_T kappaFx[3], const real_T FzFx[3], const real_T FxMap[9], real_T b_FZMAX, real_T b_VXLOW, real_T b_kappamax, real_T *Fx, real_T *kappa)
{
  real_T Vxpabs;
  real_T d;
  int32_T d_trueCount;
  real_T FxMap_0[9];
  int32_T i;
  real_T y_data;
  real_T Vxpabs_data;
  int32_T Vxpabs_sizes[2];
  real_T tmp_data;
  int32_T tmp_sizes[2];
  real_T b_idx_0;
  Vxpabs = fabs(Vx);
  d_trueCount = 0;
  if (Vxpabs < b_VXLOW) {
    d_trueCount = 1;
  }

  Vxpabs_sizes[0] = 1;
  Vxpabs_sizes[1] = d_trueCount;
  for (i = 0; i < d_trueCount; i++) {
    Vxpabs_data = Vxpabs / b_VXLOW;
  }

  Adaptive_Cruise_Control_S_power(&Vxpabs_data, Vxpabs_sizes, &tmp_data, tmp_sizes);
  for (i = 0; i < d_trueCount; i++) {
    y_data = 2.0 * b_VXLOW / (3.0 - tmp_data);
  }

  if (Vxpabs < b_VXLOW) {
    Vxpabs = y_data;
  }

  b_idx_0 = Fz;
  if (Fz < 0.0) {
    b_idx_0 = 0.0;
  }

  if (b_idx_0 > b_FZMAX) {
    b_idx_0 = b_FZMAX;
  }

  *kappa = (Re * omega - Vx) / Vxpabs;
  d = *kappa;
  d_trueCount = 0;
  if (*kappa < -b_kappamax) {
    d_trueCount = 1;
  }

  for (i = 0; i < d_trueCount; i++) {
    d = -b_kappamax;
  }

  *kappa = d;
  if (d > b_kappamax) {
    *kappa = b_kappamax;
  }

  for (i = 0; i < 3; i++) {
    FxMap_0[3 * i] = FxMap[i];
    FxMap_0[1 + 3 * i] = FxMap[i + 3];
    FxMap_0[2 + 3 * i] = FxMap[i + 6];
  }

  *Fx = Adaptive_Cruise_Control_interp2(kappaFx, FzFx, FxMap_0, *kappa, b_idx_0);
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  real_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = (rtNaN);
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S12>/Simple Magic Tire' */
static real_T Adaptive_Cruis_automltirelongMy(real_T Fx, real_T Fz, real_T omega, real_T Vx, real_T press, real_T FNOMIN, real_T NOMPRES, real_T QSY1, real_T QSY2, real_T QSY3, real_T QSY4, real_T QSY5, real_T QSY6, real_T QSY7, real_T QSY8, real_T b_gamma, real_T lam_My, real_T UNLOADED_RADIUS, real_T b_FZMAX, real_T PRESMIN, real_T PRESMAX)
{
  real_T My;
  real_T b_idx_0;
  real_T d_idx_0;
  b_idx_0 = press;
  if (press < PRESMIN) {
    b_idx_0 = PRESMIN;
  }

  if (b_idx_0 > PRESMAX) {
    b_idx_0 = PRESMAX;
  }

  d_idx_0 = Fz;
  if (Fz < 0.0) {
    d_idx_0 = 0.0;
  }

  if (d_idx_0 > b_FZMAX) {
    d_idx_0 = b_FZMAX;
  }

  My = ((((QSY2 * Fx / FNOMIN + QSY1) + fabs(Vx / 16.7) * QSY3) + rt_powd_snf(Vx / 16.7, 4.0) * QSY4) + (QSY6 * d_idx_0 / FNOMIN + QSY5) * (b_gamma * b_gamma)) * (tanh(omega) * d_idx_0 * UNLOADED_RADIUS) * (rt_powd_snf(d_idx_0 / FNOMIN, QSY7) * rt_powd_snf(b_idx_0 / NOMPRES, QSY8)) * lam_My;
  return My;
}

/*
 * Output and update for atomic system:
 *    '<S12>/Simple Magic Tire'
 *    '<S48>/Simple Magic Tire'
 */
void Adaptive_Cruise_SimpleMagicTire(real_T rtu_Re, real_T rtu_Fz, real_T rtu_omega, real_T rtu_Vx, real_T rtu_lam_mux, real_T rtu_D, real_T rtu_C, real_T rtu_B, real_T rtu_E, const real_T rtu_kappaFx[3], const real_T rtu_FzFx[3], const real_T rtu_FxMap[9], real_T rtu_press, real_T rtu_FNOMIN, real_T rtu_NOMPRES, real_T rtu_QSY1, real_T rtu_QSY2, real_T rtu_QSY3, real_T rtu_QSY4, real_T rtu_QSY5, real_T rtu_QSY6, real_T rtu_QSY7, real_T rtu_QSY8, real_T rtu_gamma, real_T rtu_lam_My, real_T rtu_UNLOADED_RADIUS, real_T rtu_PRESMIN, real_T rtu_PRESMAX, const real_T rtu_VxMy[3], const real_T rtu_FzMy[3], const real_T rtu_MyMap[9], real_T rtu_FxType, real_T rtu_rollingType, B_SimpleMagicTire_Adaptive_Cr_T *localB, real_T rtp_FZMAX, real_T rtp_VXLOW, real_T rtp_kappamax)
{
  real_T b_kappa;
  real_T c_Fx;
  real_T rtu_MyMap_0[9];
  int32_T i;
  real_T d_idx_0;

  /* MATLAB Function 'Longitudinal Wheel/Longitudinal Basic Magic Tire/Simple Magic Tire': '<S16>:1' */
  /* '<S16>:1:4' coder.allowpcode('plain') */
  /* '<S16>:1:6' switch FxType */
  switch ((int32_T)rtu_FxType) {
   case 0:
    /* '<S16>:1:7' case 0 */
    /* '<S16>:1:8' [Fx,kappa] = automltirelongFx(Re,Fz,omega,Vx,lam_mux,D,C,B,E,FZMAX,VXLOW,kappamax); */
    Adaptive_Cruis_automltirelongFx(rtu_Re, rtu_Fz, rtu_omega, rtu_Vx, rtu_lam_mux, rtu_D, rtu_C, rtu_B, rtu_E, rtp_FZMAX, rtp_VXLOW, rtp_kappamax, &c_Fx, &b_kappa);
    break;

   case 3:
    /* '<S16>:1:9' case 3 */
    /* '<S16>:1:10' [Fx,kappa] = automltirelongFxMapped(Re,Fz,omega,Vx,kappaFx,FzFx,FxMap,FZMAX,VXLOW,kappamax); */
    Adaptive_automltirelongFxMapped(rtu_Re, rtu_Fz, rtu_omega, rtu_Vx, rtu_kappaFx, rtu_FzFx, rtu_FxMap, rtp_FZMAX, rtp_VXLOW, rtp_kappamax, &c_Fx, &b_kappa);
    break;

   default:
    /* '<S16>:1:11' otherwise */
    /* '<S16>:1:12' Fx = zeros(size(Fz)); */
    c_Fx = 0.0;

    /* '<S16>:1:13' kappa = zeros(size(Fz)); */
    break;
  }

  /* '<S16>:1:17' switch rollingType */
  switch ((int32_T)rtu_rollingType) {
   case 0:
    /* '<S16>:1:18' case 0 %'None' */
    /* 'None' */
    /* '<S16>:1:19' My = zeros(size(Fz)); */
    b_kappa = 0.0;
    break;

   case 1:
    /* '<S16>:1:20' case 1 %'Simple' */
    /* 'Simple' */
    /* '<S16>:1:21' My = automltirelongMySAE(Fz,omega,Vx,press,QSY1,QSY2,... */
    /* '<S16>:1:22'                                  QSY3, QSY7,QSY8, UNLOADED_RADIUS,FZMAX,PRESMIN,PRESMAX); */
    b_kappa = rtu_press;
    if (rtu_press < rtu_PRESMIN) {
      b_kappa = rtu_PRESMIN;
    }

    if (b_kappa > rtu_PRESMAX) {
      b_kappa = rtu_PRESMAX;
    }

    d_idx_0 = rtu_Fz;
    if (rtu_Fz > rtp_FZMAX) {
      d_idx_0 = rtp_FZMAX;
    }

    b_kappa = ((rtu_QSY2 * fabs(rtu_Vx) + rtu_QSY1) + rtu_Vx * rtu_Vx * rtu_QSY3) * (tanh(rtu_omega) * d_idx_0 * rtu_UNLOADED_RADIUS) * (rt_powd_snf(d_idx_0, rtu_QSY7) * rt_powd_snf(b_kappa, rtu_QSY8));
    break;

   case 2:
    /* '<S16>:1:23' case 2 %'Magic Formula' */
    /* 'Magic Formula' */
    /* '<S16>:1:24' My = automltirelongMy(Fx,Fz,omega,Vx,press,FNOMIN, NOMPRES, QSY1, QSY2,... */
    /* '<S16>:1:25'             QSY3,QSY4,QSY5,QSY6,QSY7,QSY8,gamma,lam_My,UNLOADED_RADIUS,FZMAX,PRESMIN, PRESMAX); */
    b_kappa = Adaptive_Cruis_automltirelongMy(c_Fx, rtu_Fz, rtu_omega, rtu_Vx, rtu_press, rtu_FNOMIN, rtu_NOMPRES, rtu_QSY1, rtu_QSY2, rtu_QSY3, rtu_QSY4, rtu_QSY5, rtu_QSY6, rtu_QSY7, rtu_QSY8, rtu_gamma, rtu_lam_My, rtu_UNLOADED_RADIUS, rtp_FZMAX, rtu_PRESMIN, rtu_PRESMAX);
    break;

   case 3:
    /* '<S16>:1:26' case 3 %'Mapped Torque' */
    /* 'Mapped Torque' */
    /* '<S16>:1:27' My = automltirelongMyMapped(omega,Fz,Vx,VxMy,FzMy,MyMap,FZMAX); */
    b_kappa = rtu_Fz;
    if (rtu_Fz < 0.0) {
      b_kappa = 0.0;
    }

    if (b_kappa > rtp_FZMAX) {
      b_kappa = rtp_FZMAX;
    }

    for (i = 0; i < 3; i++) {
      rtu_MyMap_0[3 * i] = rtu_MyMap[i];
      rtu_MyMap_0[1 + 3 * i] = rtu_MyMap[i + 3];
      rtu_MyMap_0[2 + 3 * i] = rtu_MyMap[i + 6];
    }

    b_kappa = tanh(rtu_omega) * Adaptive_Cruise_Control_interp2(rtu_VxMy, rtu_FzMy, rtu_MyMap_0, rtu_Vx, b_kappa);
    break;

   default:
    /* '<S16>:1:28' otherwise */
    /* '<S16>:1:29' My = zeros(size(Fz)); */
    b_kappa = 0.0;
    break;
  }

  localB->Fx = c_Fx;
  localB->My = b_kappa;
}

/*
 * Output and update for action system:
 *    '<S15>/Locked'
 *    '<S51>/Locked'
 */
void Adaptive_Cruise_Control__Locked(RT_MODEL_Adaptive_Cruise_Cont_T * const Adaptive_Cruise_Control_Syst_M, real_T *rty_locked_wout)
{
  if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
    /* SignalConversion: '<S27>/OutportBufferForlocked_wout' incorporates:
     *  Constant: '<S27>/locked'
     */
    *rty_locked_wout = 0.0;
  }
}

/* Model step function */
void Adaptive_Cruise_Control_System_step(void)
{
  /* local block i/o variables */
  uint16_T rtb_FixPtSwitch;
  if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
    /* set solver stop time */
    rtsiSetSolverStopTime(&Adaptive_Cruise_Control_Syst_M->solverInfo,((Adaptive_Cruise_Control_Syst_M->Timing.clockTick0+1)*Adaptive_Cruise_Control_Syst_M->Timing.stepSize0));
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
    Adaptive_Cruise_Control_Syst_M->Timing.t[0] = rtsiGetT(&Adaptive_Cruise_Control_Syst_M->solverInfo);
  }

  {
    real_T *lastU;
    ZCEventType zcEvent;
    uint16_T rtb_FixPtSum1;
    real_T rtb_Sum3;
    real_T rtb_Sum2_j;
    real_T rtb_Sum_j;
    real_T rtb_Switch_a;
    real_T rtb_product;
    real_T rtb_Divide1;
    real_T rtb_Switch;
    real_T rtb_Signconvention;
    real_T rtb_Sum1;
    int8_T rtPrevAction;
    int8_T rtAction;
    real_T rtb_Signconvention_n;
    real_T rtb_Product1;
    if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
      /* UnitDelay: '<S83>/Output' */
      rtb_FixPtSum1 = Adaptive_Cruise_Control_Syst_DW.Output_DSTATE;

      /* Gain: '<Root>/Gain' incorporates:
       *  DataTypeConversion: '<S4>/Data Type Conversion'
       *  Lookup_n-D: '<S4>/Lookup'
       *  SampleTimeMath: '<S4>/Sample Time Math'
       *  UnitDelay: '<S83>/Output'
       *
       * About '<S4>/Sample Time Math':
       *  y = u * K where K = ( w * Ts )
       */
      Adaptive_Cruise_Control_Syste_B.velocity_lead_carms = look1_binlcpw((real_T)Adaptive_Cruise_Control_Syst_DW.Output_DSTATE * 0.01, Adaptive_Cruise_Control__ConstP.Lookup_bp01Data, Adaptive_Cruise_Control__ConstP.Lookup_tableData, 8U) * 0.27777777777777779;
    }

    /* Derivative: '<Root>/Derivative' */
    if ((Adaptive_Cruise_Control_Syst_DW.TimeStampA >= Adaptive_Cruise_Control_Syst_M->Timing.t[0]) && (Adaptive_Cruise_Control_Syst_DW.TimeStampB >= Adaptive_Cruise_Control_Syst_M->Timing.t[0])) {
      rtb_Sum3 = 0.0;
    } else {
      rtb_Sum3 = Adaptive_Cruise_Control_Syst_DW.TimeStampA;
      lastU = &Adaptive_Cruise_Control_Syst_DW.LastUAtTimeA;
      if (Adaptive_Cruise_Control_Syst_DW.TimeStampA < Adaptive_Cruise_Control_Syst_DW.TimeStampB) {
        if (Adaptive_Cruise_Control_Syst_DW.TimeStampB < Adaptive_Cruise_Control_Syst_M->Timing.t[0]) {
          rtb_Sum3 = Adaptive_Cruise_Control_Syst_DW.TimeStampB;
          lastU = &Adaptive_Cruise_Control_Syst_DW.LastUAtTimeB;
        }
      } else {
        if (Adaptive_Cruise_Control_Syst_DW.TimeStampA >= Adaptive_Cruise_Control_Syst_M->Timing.t[0]) {
          rtb_Sum3 = Adaptive_Cruise_Control_Syst_DW.TimeStampB;
          lastU = &Adaptive_Cruise_Control_Syst_DW.LastUAtTimeB;
        }
      }

      rtb_Sum3 = (Adaptive_Cruise_Control_Syste_B.velocity_lead_carms - *lastU) / (Adaptive_Cruise_Control_Syst_M->Timing.t[0] - rtb_Sum3);
    }

    /* End of Derivative: '<Root>/Derivative' */

    /* Product: '<Root>/Product1' incorporates:
     *  Constant: '<Root>/r'
     *  Gain: '<Root>/Gain4'
     */
    rtb_Sum2_j = 750.0 * rtb_Sum3 * 0.3;
    if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
      /* Gain: '<Root>/Gain3' incorporates:
       *  UnitDelay: '<Root>/Unit Delay'
       */
      Adaptive_Cruise_Control_Syste_B.Gain3 = 0.3 * Adaptive_Cruise_Control_Syst_DW.UnitDelay_DSTATE;

      /* UnitDelay: '<Root>/Unit Delay1' */
      Adaptive_Cruise_Control_Syste_B.UnitDelay1 = Adaptive_Cruise_Control_Syst_DW.UnitDelay1_DSTATE;
    }

    /* Sum: '<Root>/Sum' */
    rtb_Sum_j = rtb_Sum2_j - Adaptive_Cruise_Control_Syste_B.Gain3;

    /* Gain: '<S6>/Filter Coefficient' incorporates:
     *  Gain: '<S6>/Derivative Gain'
     *  Integrator: '<S6>/Filter'
     *  Sum: '<S6>/SumD'
     */
    Adaptive_Cruise_Control_Syste_B.FilterCoefficient = (0.0 * rtb_Sum_j - Adaptive_Cruise_Control_Syste_X.Filter_CSTATE) * 100.0;

    /* Sum: '<S6>/Sum' incorporates:
     *  Gain: '<S6>/Proportional Gain'
     *  Integrator: '<S6>/Integrator'
     */
    rtb_product = (-2.0 * rtb_Sum_j + Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE) + Adaptive_Cruise_Control_Syste_B.FilterCoefficient;

    /* Product: '<S1>/Divide1' incorporates:
     *  Gain: '<S1>/Gain2'
     */
    rtb_Divide1 = 0.2 * rtb_product / Adaptive_Cruise_Control__ConstB.Product;

    /* Sum: '<Root>/Sum3' */
    rtb_Sum3 -= Adaptive_Cruise_Control_Syste_B.UnitDelay1;

    /* Gain: '<S7>/Filter Coefficient' incorporates:
     *  Gain: '<S7>/Derivative Gain'
     *  Integrator: '<S7>/Filter'
     *  Sum: '<S7>/SumD'
     */
    Adaptive_Cruise_Control_Syste_B.FilterCoefficient_j = (0.0 * rtb_Sum3 - Adaptive_Cruise_Control_Syste_X.Filter_CSTATE_p) * 100.0;

    /* Product: '<S3>/Divide' incorporates:
     *  Constant: '<S3>/H'
     *  Constant: '<S3>/L'
     *  Constant: '<S3>/M'
     *  Gain: '<S7>/Proportional Gain'
     *  Integrator: '<S7>/Integrator'
     *  Sum: '<S7>/Sum'
     */
    rtb_Switch_a = ((rtb_Sum3 + Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE_c) + Adaptive_Cruise_Control_Syste_B.FilterCoefficient_j) * 0.1 * 0.55 * 750.0 / 2.6;

    /* Product: '<S3>/Product1' incorporates:
     *  Constant: '<S3>/g'
     *  Sum: '<S3>/Subtract'
     */
    rtb_Product1 = (Adaptive_Cruise_Control__ConstB.Product_p - rtb_Switch_a) * 9.81;
    if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
      /* HitCross: '<S38>/Velocities Match' */
      zcEvent = rt_ZCFcn(ANY_ZERO_CROSSING,&Adaptive_Cruise_Control_PrevZCX.VelocitiesMatch_Input_ZCE,
                         (Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE_j));
      if (Adaptive_Cruise_Control_Syst_DW.VelocitiesMatch_MODE == 0) {
        if (zcEvent != NO_ZCEVENT) {
          Adaptive_Cruise_Control_Syste_B.VelocitiesMatch = !Adaptive_Cruise_Control_Syste_B.VelocitiesMatch;
          Adaptive_Cruise_Control_Syst_DW.VelocitiesMatch_MODE = 1;
        } else if (Adaptive_Cruise_Control_Syste_B.VelocitiesMatch) {
          if (Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE_j != 0.0) {
            Adaptive_Cruise_Control_Syste_B.VelocitiesMatch = false;
          }
        } else {
          if (Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE_j == 0.0) {
            Adaptive_Cruise_Control_Syste_B.VelocitiesMatch = true;
          }
        }
      } else {
        if (Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE_j != 0.0) {
          Adaptive_Cruise_Control_Syste_B.VelocitiesMatch = false;
        }

        Adaptive_Cruise_Control_Syst_DW.VelocitiesMatch_MODE = 0;
      }

      /* End of HitCross: '<S38>/Velocities Match' */
    }

    /* Sum: '<Root>/Sum2' */
    rtb_Sum2_j -= Adaptive_Cruise_Control_Syste_B.Gain3;

    /* Gain: '<S5>/Filter Coefficient' incorporates:
     *  Gain: '<S5>/Derivative Gain'
     *  Integrator: '<S5>/Filter'
     *  Sum: '<S5>/SumD'
     */
    Adaptive_Cruise_Control_Syste_B.FilterCoefficient_n = (0.0 * rtb_Sum2_j - Adaptive_Cruise_Control_Syste_X.Filter_CSTATE_a) * 100.0;

    /* Integrator: '<S43>/Integrator' */
    rtb_Switch = Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE_o;

    /* Gain: '<S15>/Sign convention' incorporates:
     *  Constant: '<S29>/Constant12'
     *  Gain: '<S5>/Proportional Gain'
     *  Integrator: '<S43>/Integrator'
     *  Integrator: '<S5>/Integrator'
     *  Product: '<S29>/Product1'
     *  Sum: '<S15>/Add1'
     *  Sum: '<S5>/Sum'
     */
    rtb_Signconvention = -(((2.0 * rtb_Sum2_j + Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE_j) + Adaptive_Cruise_Control_Syste_B.FilterCoefficient_n) - Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE_o * 0.3);

    /* Product: '<S34>/product' incorporates:
     *  Constant: '<S34>/Disk brake actuator bore'
     *  Constant: '<S34>/Number of brake pads'
     *  Gain: '<S1>/Gain1'
     *  Product: '<S1>/Divide'
     */
    rtb_Sum1 = 0.8 * rtb_product / Adaptive_Cruise_Control__ConstB.Product * Adaptive_Cruise_Control__ConstB.TorqueConversion1 * 0.05 * 2.0;

    /* Saturate: '<S34>/Disallow Negative Brake Torque' */
    if (rtb_Sum1 <= 2.2204460492503131E-16) {
      rtb_Sum1 = 2.2204460492503131E-16;
    }

    /* End of Saturate: '<S34>/Disallow Negative Brake Torque' */

    /* Gain: '<S34>/Torque Conversion' */
    rtb_Sum1 *= 0.03556;
    if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
      /* Memory: '<S39>/Memory' */
      Adaptive_Cruise_Control_Syste_B.Memory = Adaptive_Cruise_Control_Syst_DW.Memory_PreviousInput;

      /* InitialCondition: '<S26>/IC' */
      if (Adaptive_Cruise_Control_Syst_DW.IC_FirstOutputTime) {
        Adaptive_Cruise_Control_Syst_DW.IC_FirstOutputTime = false;
      }

      /* End of InitialCondition: '<S26>/IC' */

      /* Switch: '<S26>/Switch' incorporates:
       *  Constant: '<S15>/locked wo'
       */
      Adaptive_Cruise_Control_Syste_B.Switch = 0.0;
    }

    /* CombinatorialLogic: '<S39>/Combinatorial  Logic' incorporates:
     *  Abs: '<S37>/Abs'
     *  Abs: '<S42>/Abs'
     *  Logic: '<S38>/Logic'
     *  RelationalOperator: '<S37>/Relational Operator'
     *  RelationalOperator: '<S42>/Relational Operator'
     *  Sum: '<S40>/Sum1'
     *  Sum: '<S40>/Sum2'
     *  UnaryMinus: '<S41>/Unary Minus'
     */
    Adaptive_Cruise_Control_Syste_B.CombinatorialLogic = Adaptive_Cruise_Control__ConstP.pooled27[(((fabs(((0.0 - rtb_Signconvention) - Adaptive_Cruise_Control__ConstB.OutputDamping) + Adaptive_Cruise_Control__ConstB.OutputDamping) >= rtb_Sum1) + ((uint32_T)(Adaptive_Cruise_Control_Syste_B.VelocitiesMatch && (fabs(-rtb_Signconvention) <= rtb_Sum1)) << 1)) << 1) + Adaptive_Cruise_Control_Syste_B.Memory];

    /* If: '<S15>/If' */
    rtPrevAction = Adaptive_Cruise_Control_Syst_DW.If_ActiveSubsystem;
    if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
      rtAction = (int8_T)!Adaptive_Cruise_Control_Syste_B.CombinatorialLogic;
      Adaptive_Cruise_Control_Syst_DW.If_ActiveSubsystem = rtAction;
    } else {
      rtAction = Adaptive_Cruise_Control_Syst_DW.If_ActiveSubsystem;
    }

    switch (rtAction) {
     case 0:
      /* Outputs for IfAction SubSystem: '<S15>/Locked' incorporates:
       *  ActionPort: '<S27>/Action'
       */
      Adaptive_Cruise_Control__Locked(Adaptive_Cruise_Control_Syst_M, &Adaptive_Cruise_Control_Syste_B.omega);

      /* End of Outputs for SubSystem: '<S15>/Locked' */
      break;

     case 1:
      if (rtAction != rtPrevAction) {
        /* InitializeConditions for IfAction SubSystem: '<S15>/Unlocked' incorporates:
         *  InitializeConditions for ActionPort: '<S28>/Action'
         */
        /* InitializeConditions for If: '<S15>/If' incorporates:
         *  InitializeConditions for Integrator: '<S28>/Output Integrator'
         */
        if (rtmIsFirstInitCond(Adaptive_Cruise_Control_Syst_M)) {
          Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE_j = 0.0;
        }

        Adaptive_Cruise_Control_Syst_DW.OutputIntegrator_IWORK_b = 1;

        /* End of InitializeConditions for If: '<S15>/If' */
        /* End of InitializeConditions for SubSystem: '<S15>/Unlocked' */
      }

      /* Outputs for IfAction SubSystem: '<S15>/Unlocked' incorporates:
       *  ActionPort: '<S28>/Action'
       */
      /* Integrator: '<S28>/Output Integrator' */
      if (Adaptive_Cruise_Control_Syst_DW.OutputIntegrator_IWORK_b != 0) {
        Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE_j = Adaptive_Cruise_Control_Syste_B.Switch;
        rtsiSetBlkStateChange(&Adaptive_Cruise_Control_Syst_M->solverInfo, true);
      }

      /* Gain: '<S28>/Output Inertia' incorporates:
       *  Gain: '<S28>/-4'
       *  Gain: '<S28>/Output Damping'
       *  Integrator: '<S28>/Output Integrator'
       *  Product: '<S28>/Max Dynamic Friction Torque'
       *  Sum: '<S28>/Output Sum'
       *  Trigonometry: '<S28>/Trigonometric Function'
       */
      Adaptive_Cruise_Control_Syste_B.OutputInertia_p = ((tanh(-4.0 * Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE_j) * rtb_Sum1 - rtb_Signconvention) - 0.001 * Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE_j) * 1.25;

      /* SignalConversion: '<S28>/Signal Conversion1' incorporates:
       *  Integrator: '<S28>/Output Integrator'
       */
      Adaptive_Cruise_Control_Syste_B.omega = Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE_j;

      /* End of Outputs for SubSystem: '<S15>/Unlocked' */
      break;
    }

    /* End of If: '<S15>/If' */

    /* Integrator: '<Root>/Integrator1' */
    Adaptive_Cruise_Control_Syste_B.velocity_ego_carms = Adaptive_Cruise_Control_Syste_X.Integrator1_CSTATE;

    /* Sum: '<Root>/Sum1' */
    rtb_Sum1 = Adaptive_Cruise_Control_Syste_B.velocity_lead_carms - Adaptive_Cruise_Control_Syste_B.velocity_ego_carms;

    /* Gain: '<S8>/Filter Coefficient' incorporates:
     *  Gain: '<S8>/Derivative Gain'
     *  Integrator: '<S8>/Filter'
     *  Sum: '<S8>/SumD'
     */
    Adaptive_Cruise_Control_Syste_B.FilterCoefficient_c = (0.0 * rtb_Sum1 - Adaptive_Cruise_Control_Syste_X.Filter_CSTATE_ao) * 100.0;

    /* Gain: '<S8>/Proportional Gain' incorporates:
     *  Integrator: '<S8>/Integrator'
     *  Sum: '<S8>/Sum'
     */
    rtb_Signconvention = ((rtb_Sum1 + Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE_p) + Adaptive_Cruise_Control_Syste_B.FilterCoefficient_c) * -15.0;

    /* MATLAB Function: '<S12>/Simple Magic Tire' incorporates:
     *  Constant: '<S11>/Constant9'
     *  Constant: '<S11>/FxType'
     *  Constant: '<S11>/rollType'
     *  Constant: '<S17>/Constant'
     *  Constant: '<S17>/Constant1'
     *  Constant: '<S17>/Constant12'
     *  Constant: '<S17>/Constant14'
     *  Constant: '<S17>/Constant19'
     *  Constant: '<S17>/Constant3'
     *  Constant: '<S17>/Constant6'
     *  Constant: '<S17>/Constant7'
     *  Constant: '<S21>/Constant1'
     *  Constant: '<S21>/Constant10'
     *  Constant: '<S21>/Constant11'
     *  Constant: '<S21>/Constant12'
     *  Constant: '<S21>/Constant13'
     *  Constant: '<S21>/Constant14'
     *  Constant: '<S21>/Constant15'
     *  Constant: '<S21>/Constant16'
     *  Constant: '<S21>/Constant17'
     *  Constant: '<S21>/Constant18'
     *  Constant: '<S21>/Constant19'
     *  Constant: '<S21>/Constant2'
     *  Constant: '<S21>/Constant3'
     *  Constant: '<S21>/Constant4'
     *  Constant: '<S21>/Constant5'
     *  Constant: '<S21>/Constant6'
     *  Constant: '<S21>/Constant7'
     *  Constant: '<S21>/Constant8'
     *  Constant: '<S21>/Constant9'
     */
    Adaptive_Cruise_SimpleMagicTire(0.3, rtb_Product1, Adaptive_Cruise_Control_Syste_B.omega, rtb_Signconvention, 1.0, 1.0, 1.65, 10.0, 0.01, Adaptive_Cruise_Control__ConstP.pooled22, Adaptive_Cruise_Control__ConstP.pooled22, Adaptive_Cruise_Control__ConstP.pooled23, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, Adaptive_Cruise_Control__ConstP.pooled22, Adaptive_Cruise_Control__ConstP.pooled22, Adaptive_Cruise_Control__ConstP.pooled23, 0.0, 0.0, &Adaptive_Cruise_Control_Syste_B.sf_SimpleMagicTire, 10000.0, 1.0, 1.5);

    /* Product: '<S3>/Product3' incorporates:
     *  Constant: '<S3>/g'
     *  Sum: '<S3>/Subtract1'
     */
    rtb_Switch_a = (Adaptive_Cruise_Control__ConstB.Product2 - rtb_Switch_a) * 9.81;
    if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
      /* HitCross: '<S74>/Velocities Match' */
      zcEvent = rt_ZCFcn(ANY_ZERO_CROSSING,&Adaptive_Cruise_Control_PrevZCX.VelocitiesMatch_Input_ZCE_a,
                         (Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE));
      if (Adaptive_Cruise_Control_Syst_DW.VelocitiesMatch_MODE_j == 0) {
        if (zcEvent != NO_ZCEVENT) {
          Adaptive_Cruise_Control_Syste_B.VelocitiesMatch_a = !Adaptive_Cruise_Control_Syste_B.VelocitiesMatch_a;
          Adaptive_Cruise_Control_Syst_DW.VelocitiesMatch_MODE_j = 1;
        } else if (Adaptive_Cruise_Control_Syste_B.VelocitiesMatch_a) {
          if (Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE != 0.0) {
            Adaptive_Cruise_Control_Syste_B.VelocitiesMatch_a = false;
          }
        } else {
          if (Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE == 0.0) {
            Adaptive_Cruise_Control_Syste_B.VelocitiesMatch_a = true;
          }
        }
      } else {
        if (Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE != 0.0) {
          Adaptive_Cruise_Control_Syste_B.VelocitiesMatch_a = false;
        }

        Adaptive_Cruise_Control_Syst_DW.VelocitiesMatch_MODE_j = 0;
      }

      /* End of HitCross: '<S74>/Velocities Match' */
    }

    /* Integrator: '<S79>/Integrator' */
    rtb_Product1 = Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE_g;

    /* Gain: '<S51>/Sign convention' incorporates:
     *  Constant: '<S2>/Constant'
     *  Constant: '<S65>/Constant12'
     *  Integrator: '<S79>/Integrator'
     *  Product: '<S65>/Product1'
     *  Sum: '<S51>/Add1'
     */
    rtb_Signconvention_n = -(0.0 - Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE_g * 0.3);

    /* Product: '<S70>/product' incorporates:
     *  Constant: '<S70>/Disk brake actuator bore'
     *  Constant: '<S70>/Number of brake pads'
     */
    rtb_product = rtb_Divide1 * Adaptive_Cruise_Control__ConstB.TorqueConversion1_c * 0.05 * 2.0;

    /* Saturate: '<S70>/Disallow Negative Brake Torque' */
    if (rtb_product <= 2.2204460492503131E-16) {
      rtb_product = 2.2204460492503131E-16;
    }

    /* End of Saturate: '<S70>/Disallow Negative Brake Torque' */

    /* Gain: '<S70>/Torque Conversion' */
    rtb_Divide1 = 0.03556 * rtb_product;
    if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
      /* Memory: '<S75>/Memory' */
      Adaptive_Cruise_Control_Syste_B.Memory_b = Adaptive_Cruise_Control_Syst_DW.Memory_PreviousInput_a;

      /* InitialCondition: '<S62>/IC' */
      if (Adaptive_Cruise_Control_Syst_DW.IC_FirstOutputTime_g) {
        Adaptive_Cruise_Control_Syst_DW.IC_FirstOutputTime_g = false;
      }

      /* End of InitialCondition: '<S62>/IC' */

      /* Switch: '<S62>/Switch' incorporates:
       *  Constant: '<S51>/locked wo'
       */
      Adaptive_Cruise_Control_Syste_B.Switch_n = 0.0;
    }

    /* CombinatorialLogic: '<S75>/Combinatorial  Logic' incorporates:
     *  Abs: '<S73>/Abs'
     *  Abs: '<S78>/Abs'
     *  Logic: '<S74>/Logic'
     *  RelationalOperator: '<S73>/Relational Operator'
     *  RelationalOperator: '<S78>/Relational Operator'
     *  Sum: '<S76>/Sum1'
     *  Sum: '<S76>/Sum2'
     *  UnaryMinus: '<S77>/Unary Minus'
     */
    Adaptive_Cruise_Control_Syste_B.CombinatorialLogic_c = Adaptive_Cruise_Control__ConstP.pooled27[(((fabs(((0.0 - rtb_Signconvention_n) - Adaptive_Cruise_Control__ConstB.OutputDamping_a) + Adaptive_Cruise_Control__ConstB.OutputDamping_a) >= rtb_Divide1) + ((uint32_T)(Adaptive_Cruise_Control_Syste_B.VelocitiesMatch_a && (fabs(-rtb_Signconvention_n) <= rtb_Divide1)) << 1)) << 1) + Adaptive_Cruise_Control_Syste_B.Memory_b];

    /* If: '<S51>/If' */
    rtPrevAction = Adaptive_Cruise_Control_Syst_DW.If_ActiveSubsystem_i;
    if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
      rtAction = (int8_T)!Adaptive_Cruise_Control_Syste_B.CombinatorialLogic_c;
      Adaptive_Cruise_Control_Syst_DW.If_ActiveSubsystem_i = rtAction;
    } else {
      rtAction = Adaptive_Cruise_Control_Syst_DW.If_ActiveSubsystem_i;
    }

    switch (rtAction) {
     case 0:
      /* Outputs for IfAction SubSystem: '<S51>/Locked' incorporates:
       *  ActionPort: '<S63>/Action'
       */
      Adaptive_Cruise_Control__Locked(Adaptive_Cruise_Control_Syst_M, &Adaptive_Cruise_Control_Syste_B.omega_c);

      /* End of Outputs for SubSystem: '<S51>/Locked' */
      break;

     case 1:
      if (rtAction != rtPrevAction) {
        /* InitializeConditions for IfAction SubSystem: '<S51>/Unlocked' incorporates:
         *  InitializeConditions for ActionPort: '<S64>/Action'
         */
        /* InitializeConditions for If: '<S51>/If' incorporates:
         *  InitializeConditions for Integrator: '<S64>/Output Integrator'
         */
        if (rtmIsFirstInitCond(Adaptive_Cruise_Control_Syst_M)) {
          Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE = 0.0;
        }

        Adaptive_Cruise_Control_Syst_DW.OutputIntegrator_IWORK = 1;

        /* End of InitializeConditions for If: '<S51>/If' */
        /* End of InitializeConditions for SubSystem: '<S51>/Unlocked' */
      }

      /* Outputs for IfAction SubSystem: '<S51>/Unlocked' incorporates:
       *  ActionPort: '<S64>/Action'
       */
      /* Integrator: '<S64>/Output Integrator' */
      if (Adaptive_Cruise_Control_Syst_DW.OutputIntegrator_IWORK != 0) {
        Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE = Adaptive_Cruise_Control_Syste_B.Switch_n;
        rtsiSetBlkStateChange(&Adaptive_Cruise_Control_Syst_M->solverInfo, true);
      }

      /* Gain: '<S64>/Output Inertia' incorporates:
       *  Gain: '<S64>/-4'
       *  Gain: '<S64>/Output Damping'
       *  Integrator: '<S64>/Output Integrator'
       *  Product: '<S64>/Max Dynamic Friction Torque'
       *  Sum: '<S64>/Output Sum'
       *  Trigonometry: '<S64>/Trigonometric Function'
       */
      Adaptive_Cruise_Control_Syste_B.OutputInertia = ((tanh(-4.0 * Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE) * rtb_Divide1 - rtb_Signconvention_n) - 0.001 * Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE) * 1.25;

      /* SignalConversion: '<S64>/Signal Conversion1' incorporates:
       *  Integrator: '<S64>/Output Integrator'
       */
      Adaptive_Cruise_Control_Syste_B.omega_c = Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE;

      /* End of Outputs for SubSystem: '<S51>/Unlocked' */
      break;
    }

    /* End of If: '<S51>/If' */

    /* MATLAB Function: '<S48>/Simple Magic Tire' incorporates:
     *  Constant: '<S47>/Constant9'
     *  Constant: '<S47>/FxType'
     *  Constant: '<S47>/rollType'
     *  Constant: '<S53>/Constant'
     *  Constant: '<S53>/Constant1'
     *  Constant: '<S53>/Constant12'
     *  Constant: '<S53>/Constant14'
     *  Constant: '<S53>/Constant19'
     *  Constant: '<S53>/Constant3'
     *  Constant: '<S53>/Constant6'
     *  Constant: '<S53>/Constant7'
     *  Constant: '<S57>/Constant1'
     *  Constant: '<S57>/Constant10'
     *  Constant: '<S57>/Constant11'
     *  Constant: '<S57>/Constant12'
     *  Constant: '<S57>/Constant13'
     *  Constant: '<S57>/Constant14'
     *  Constant: '<S57>/Constant15'
     *  Constant: '<S57>/Constant16'
     *  Constant: '<S57>/Constant17'
     *  Constant: '<S57>/Constant18'
     *  Constant: '<S57>/Constant19'
     *  Constant: '<S57>/Constant2'
     *  Constant: '<S57>/Constant3'
     *  Constant: '<S57>/Constant4'
     *  Constant: '<S57>/Constant5'
     *  Constant: '<S57>/Constant6'
     *  Constant: '<S57>/Constant7'
     *  Constant: '<S57>/Constant8'
     *  Constant: '<S57>/Constant9'
     */
    Adaptive_Cruise_SimpleMagicTire(0.3, rtb_Switch_a, Adaptive_Cruise_Control_Syste_B.omega_c, rtb_Signconvention, 1.0, 1.0, 1.65, 10.0, 0.01, Adaptive_Cruise_Control__ConstP.pooled22, Adaptive_Cruise_Control__ConstP.pooled22, Adaptive_Cruise_Control__ConstP.pooled23, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, Adaptive_Cruise_Control__ConstP.pooled22, Adaptive_Cruise_Control__ConstP.pooled22, Adaptive_Cruise_Control__ConstP.pooled23, 0.0, 0.0, &Adaptive_Cruise_Control_Syste_B.sf_SimpleMagicTire_g, 10000.0, 1.0, 1.5);

    /* Sum: '<S2>/Sum' */
    Adaptive_Cruise_Control_Syste_B.Sum = Adaptive_Cruise_Control_Syste_B.sf_SimpleMagicTire.Fx + Adaptive_Cruise_Control_Syste_B.sf_SimpleMagicTire_g.Fx;

    /* Product: '<Root>/Divide' incorporates:
     *  Constant: '<Root>/M'
     */
    Adaptive_Cruise_Control_Syste_B.Acceleration_ego_car = Adaptive_Cruise_Control_Syste_B.Sum / 750.0;

    /* Product: '<S29>/Product3' incorporates:
     *  Constant: '<S29>/Constant12'
     */
    rtb_Divide1 = Adaptive_Cruise_Control_Syste_B.omega * 0.3;

    /* Switch: '<S44>/Switch' incorporates:
     *  Abs: '<S44>/Abs'
     *  Constant: '<S45>/Constant'
     *  Constant: '<S46>/Constant'
     *  Fcn: '<S44>/Fcn'
     *  Logic: '<S44>/Logical Operator'
     *  RelationalOperator: '<S45>/Compare'
     *  RelationalOperator: '<S46>/Compare'
     */
    if ((rtb_Divide1 >= -1.0) && (rtb_Divide1 <= 1.0)) {
      rtb_Switch_a = 2.0 / (3.0 - rt_powd_snf(rtb_Divide1, 2.0));
    } else {
      rtb_Switch_a = fabs(rtb_Divide1);
    }

    /* End of Switch: '<S44>/Switch' */

    /* Product: '<S43>/Product' incorporates:
     *  Constant: '<S29>/Constant12'
     *  Constant: '<S29>/Constant2'
     *  Product: '<S29>/Product'
     *  Product: '<S29>/Product2'
     *  Sum: '<S29>/Add'
     *  Sum: '<S43>/Sum'
     */
    Adaptive_Cruise_Control_Syste_B.Product = ((Adaptive_Cruise_Control_Syste_B.sf_SimpleMagicTire.My / 0.3 + Adaptive_Cruise_Control_Syste_B.sf_SimpleMagicTire.Fx) - rtb_Switch) * (rtb_Switch_a / 0.5);

    /* Product: '<S65>/Product3' incorporates:
     *  Constant: '<S65>/Constant12'
     */
    rtb_Switch = Adaptive_Cruise_Control_Syste_B.omega_c * 0.3;

    /* Switch: '<S80>/Switch' incorporates:
     *  Abs: '<S80>/Abs'
     *  Constant: '<S81>/Constant'
     *  Constant: '<S82>/Constant'
     *  Fcn: '<S80>/Fcn'
     *  Logic: '<S80>/Logical Operator'
     *  RelationalOperator: '<S81>/Compare'
     *  RelationalOperator: '<S82>/Compare'
     */
    if ((rtb_Switch >= -1.0) && (rtb_Switch <= 1.0)) {
      rtb_Switch = 2.0 / (3.0 - rt_powd_snf(rtb_Switch, 2.0));
    } else {
      rtb_Switch = fabs(rtb_Switch);
    }

    /* End of Switch: '<S80>/Switch' */

    /* Product: '<S79>/Product' incorporates:
     *  Constant: '<S65>/Constant12'
     *  Constant: '<S65>/Constant2'
     *  Product: '<S65>/Product'
     *  Product: '<S65>/Product2'
     *  Sum: '<S65>/Add'
     *  Sum: '<S79>/Sum'
     */
    Adaptive_Cruise_Control_Syste_B.Product_d = ((Adaptive_Cruise_Control_Syste_B.sf_SimpleMagicTire_g.My / 0.3 + Adaptive_Cruise_Control_Syste_B.sf_SimpleMagicTire_g.Fx) - rtb_Product1) * (rtb_Switch / 0.5);
    if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
      /* Sum: '<S84>/FixPt Sum1' incorporates:
       *  Constant: '<S84>/FixPt Constant'
       */
      rtb_FixPtSum1++;

      /* Switch: '<S85>/FixPt Switch' incorporates:
       *  Constant: '<S85>/Constant'
       */
      if (rtb_FixPtSum1 > 12000) {
        rtb_FixPtSwitch = 0U;
      } else {
        rtb_FixPtSwitch = rtb_FixPtSum1;
      }

      /* End of Switch: '<S85>/FixPt Switch' */
    }

    /* Gain: '<S5>/Integral Gain' */
    Adaptive_Cruise_Control_Syste_B.IntegralGain = 0.0 * rtb_Sum2_j;

    /* Gain: '<S6>/Integral Gain' */
    Adaptive_Cruise_Control_Syste_B.IntegralGain_o = 0.0 * rtb_Sum_j;

    /* Gain: '<S7>/Integral Gain' */
    Adaptive_Cruise_Control_Syste_B.IntegralGain_e = 0.0 * rtb_Sum3;

    /* Gain: '<S8>/Integral Gain' */
    Adaptive_Cruise_Control_Syste_B.IntegralGain_b = 0.0 * rtb_Sum1;
  }

  if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
    real_T *lastU;
    if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
      /* Update for UnitDelay: '<S83>/Output' */
      Adaptive_Cruise_Control_Syst_DW.Output_DSTATE = rtb_FixPtSwitch;
    }

    /* Update for Derivative: '<Root>/Derivative' */
    if (Adaptive_Cruise_Control_Syst_DW.TimeStampA == (rtInf)) {
      Adaptive_Cruise_Control_Syst_DW.TimeStampA = Adaptive_Cruise_Control_Syst_M->Timing.t[0];
      lastU = &Adaptive_Cruise_Control_Syst_DW.LastUAtTimeA;
    } else if (Adaptive_Cruise_Control_Syst_DW.TimeStampB == (rtInf)) {
      Adaptive_Cruise_Control_Syst_DW.TimeStampB = Adaptive_Cruise_Control_Syst_M->Timing.t[0];
      lastU = &Adaptive_Cruise_Control_Syst_DW.LastUAtTimeB;
    } else if (Adaptive_Cruise_Control_Syst_DW.TimeStampA < Adaptive_Cruise_Control_Syst_DW.TimeStampB) {
      Adaptive_Cruise_Control_Syst_DW.TimeStampA = Adaptive_Cruise_Control_Syst_M->Timing.t[0];
      lastU = &Adaptive_Cruise_Control_Syst_DW.LastUAtTimeA;
    } else {
      Adaptive_Cruise_Control_Syst_DW.TimeStampB = Adaptive_Cruise_Control_Syst_M->Timing.t[0];
      lastU = &Adaptive_Cruise_Control_Syst_DW.LastUAtTimeB;
    }

    *lastU = Adaptive_Cruise_Control_Syste_B.velocity_lead_carms;

    /* End of Update for Derivative: '<Root>/Derivative' */
    if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
      /* Update for UnitDelay: '<Root>/Unit Delay' */
      Adaptive_Cruise_Control_Syst_DW.UnitDelay_DSTATE = Adaptive_Cruise_Control_Syste_B.Sum;

      /* Update for UnitDelay: '<Root>/Unit Delay1' */
      Adaptive_Cruise_Control_Syst_DW.UnitDelay1_DSTATE = Adaptive_Cruise_Control_Syste_B.Acceleration_ego_car;

      /* Update for Memory: '<S39>/Memory' */
      Adaptive_Cruise_Control_Syst_DW.Memory_PreviousInput = Adaptive_Cruise_Control_Syste_B.CombinatorialLogic;

      /* Update for Memory: '<S75>/Memory' */
      Adaptive_Cruise_Control_Syst_DW.Memory_PreviousInput_a = Adaptive_Cruise_Control_Syste_B.CombinatorialLogic_c;
    }

    /* Update for If: '<S15>/If' */
    if (Adaptive_Cruise_Control_Syst_DW.If_ActiveSubsystem == 1) {
      /* Update for IfAction SubSystem: '<S15>/Unlocked' incorporates:
       *  Update for ActionPort: '<S28>/Action'
       */
      /* Update for Integrator: '<S28>/Output Integrator' */
      Adaptive_Cruise_Control_Syst_DW.OutputIntegrator_IWORK_b = 0;

      /* End of Update for SubSystem: '<S15>/Unlocked' */
    }

    /* End of Update for If: '<S15>/If' */

    /* Update for If: '<S51>/If' */
    if (Adaptive_Cruise_Control_Syst_DW.If_ActiveSubsystem_i == 1) {
      /* Update for IfAction SubSystem: '<S51>/Unlocked' incorporates:
       *  Update for ActionPort: '<S64>/Action'
       */
      /* Update for Integrator: '<S64>/Output Integrator' */
      Adaptive_Cruise_Control_Syst_DW.OutputIntegrator_IWORK = 0;

      /* End of Update for SubSystem: '<S51>/Unlocked' */
    }

    /* End of Update for If: '<S51>/If' */

    /* BlkStateChangeFlag is set, need to run a minor output */
    if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
      if (rtsiGetBlkStateChange(&Adaptive_Cruise_Control_Syst_M->solverInfo)) {
        rtsiSetSimTimeStep(&Adaptive_Cruise_Control_Syst_M->solverInfo,MINOR_TIME_STEP);
        rtsiSetBlkStateChange(&Adaptive_Cruise_Control_Syst_M->solverInfo, false);
        Adaptive_Cruise_Control_System_step();
        rtsiSetSimTimeStep(&Adaptive_Cruise_Control_Syst_M->solverInfo, MAJOR_TIME_STEP);
      }
    }
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(Adaptive_Cruise_Control_Syst_M)) {
    rt_ertODEUpdateContinuousStates(&Adaptive_Cruise_Control_Syst_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     */
    ++Adaptive_Cruise_Control_Syst_M->Timing.clockTick0;
    Adaptive_Cruise_Control_Syst_M->Timing.t[0] = rtsiGetSolverStopTime(&Adaptive_Cruise_Control_Syst_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.01s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.01, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       */
      Adaptive_Cruise_Control_Syst_M->Timing.clockTick1++;
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void Adaptive_Cruise_Control_System_derivatives(void)
{
  XDot_Adaptive_Cruise_Control__T *_rtXdot;
  _rtXdot = ((XDot_Adaptive_Cruise_Control__T *) Adaptive_Cruise_Control_Syst_M->derivs);

  /* Derivatives for Integrator: '<S6>/Integrator' */
  _rtXdot->Integrator_CSTATE = Adaptive_Cruise_Control_Syste_B.IntegralGain_o;

  /* Derivatives for Integrator: '<S6>/Filter' */
  _rtXdot->Filter_CSTATE = Adaptive_Cruise_Control_Syste_B.FilterCoefficient;

  /* Derivatives for Integrator: '<S7>/Integrator' */
  _rtXdot->Integrator_CSTATE_c = Adaptive_Cruise_Control_Syste_B.IntegralGain_e;

  /* Derivatives for Integrator: '<S7>/Filter' */
  _rtXdot->Filter_CSTATE_p = Adaptive_Cruise_Control_Syste_B.FilterCoefficient_j;

  /* Derivatives for Integrator: '<S5>/Integrator' */
  _rtXdot->Integrator_CSTATE_j = Adaptive_Cruise_Control_Syste_B.IntegralGain;

  /* Derivatives for Integrator: '<S5>/Filter' */
  _rtXdot->Filter_CSTATE_a = Adaptive_Cruise_Control_Syste_B.FilterCoefficient_n;

  /* Derivatives for Integrator: '<S43>/Integrator' */
  _rtXdot->Integrator_CSTATE_o = Adaptive_Cruise_Control_Syste_B.Product;

  /* Derivatives for If: '<S15>/If' */
  ((XDot_Adaptive_Cruise_Control__T *) Adaptive_Cruise_Control_Syst_M->derivs)->OutputIntegrator_CSTATE_j = 0.0;
  if (Adaptive_Cruise_Control_Syst_DW.If_ActiveSubsystem == 1) {
    /* Derivatives for IfAction SubSystem: '<S15>/Unlocked' incorporates:
     *  Derivatives for ActionPort: '<S28>/Action'
     */
    /* Derivatives for Integrator: '<S28>/Output Integrator' */
    _rtXdot->OutputIntegrator_CSTATE_j = Adaptive_Cruise_Control_Syste_B.OutputInertia_p;

    /* End of Derivatives for SubSystem: '<S15>/Unlocked' */
  }

  /* End of Derivatives for If: '<S15>/If' */

  /* Derivatives for Integrator: '<Root>/Integrator1' */
  _rtXdot->Integrator1_CSTATE = Adaptive_Cruise_Control_Syste_B.Acceleration_ego_car;

  /* Derivatives for Integrator: '<S8>/Integrator' */
  _rtXdot->Integrator_CSTATE_p = Adaptive_Cruise_Control_Syste_B.IntegralGain_b;

  /* Derivatives for Integrator: '<S8>/Filter' */
  _rtXdot->Filter_CSTATE_ao = Adaptive_Cruise_Control_Syste_B.FilterCoefficient_c;

  /* Derivatives for Integrator: '<S79>/Integrator' */
  _rtXdot->Integrator_CSTATE_g = Adaptive_Cruise_Control_Syste_B.Product_d;

  /* Derivatives for If: '<S51>/If' */
  ((XDot_Adaptive_Cruise_Control__T *) Adaptive_Cruise_Control_Syst_M->derivs)->OutputIntegrator_CSTATE = 0.0;
  if (Adaptive_Cruise_Control_Syst_DW.If_ActiveSubsystem_i == 1) {
    /* Derivatives for IfAction SubSystem: '<S51>/Unlocked' incorporates:
     *  Derivatives for ActionPort: '<S64>/Action'
     */
    /* Derivatives for Integrator: '<S64>/Output Integrator' */
    _rtXdot->OutputIntegrator_CSTATE = Adaptive_Cruise_Control_Syste_B.OutputInertia;

    /* End of Derivatives for SubSystem: '<S51>/Unlocked' */
  }

  /* End of Derivatives for If: '<S51>/If' */

  /* Derivatives for Integrator: '<Root>/Integrator2' */
  _rtXdot->Integrator2_CSTATE = Adaptive_Cruise_Control_Syste_B.velocity_ego_carms;

  /* Derivatives for Integrator: '<Root>/Integrator3' */
  _rtXdot->Integrator3_CSTATE = Adaptive_Cruise_Control_Syste_B.velocity_lead_carms;
}

/* Model initialize function */
void Adaptive_Cruise_Control_System_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Adaptive_Cruise_Control_Syst_M->solverInfo, &Adaptive_Cruise_Control_Syst_M->Timing.simTimeStep);
    rtsiSetTPtr(&Adaptive_Cruise_Control_Syst_M->solverInfo, &rtmGetTPtr(Adaptive_Cruise_Control_Syst_M));
    rtsiSetStepSizePtr(&Adaptive_Cruise_Control_Syst_M->solverInfo, &Adaptive_Cruise_Control_Syst_M->Timing.stepSize0);
    rtsiSetdXPtr(&Adaptive_Cruise_Control_Syst_M->solverInfo, &Adaptive_Cruise_Control_Syst_M->derivs);
    rtsiSetContStatesPtr(&Adaptive_Cruise_Control_Syst_M->solverInfo, (real_T **) &Adaptive_Cruise_Control_Syst_M->contStates);
    rtsiSetNumContStatesPtr(&Adaptive_Cruise_Control_Syst_M->solverInfo, &Adaptive_Cruise_Control_Syst_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&Adaptive_Cruise_Control_Syst_M->solverInfo, &Adaptive_Cruise_Control_Syst_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&Adaptive_Cruise_Control_Syst_M->solverInfo, &Adaptive_Cruise_Control_Syst_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&Adaptive_Cruise_Control_Syst_M->solverInfo, &Adaptive_Cruise_Control_Syst_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&Adaptive_Cruise_Control_Syst_M->solverInfo, (&rtmGetErrorStatus(Adaptive_Cruise_Control_Syst_M)));
    rtsiSetRTModelPtr(&Adaptive_Cruise_Control_Syst_M->solverInfo, Adaptive_Cruise_Control_Syst_M);
  }

  rtsiSetSimTimeStep(&Adaptive_Cruise_Control_Syst_M->solverInfo, MAJOR_TIME_STEP);
  Adaptive_Cruise_Control_Syst_M->intgData.y = Adaptive_Cruise_Control_Syst_M->odeY;
  Adaptive_Cruise_Control_Syst_M->intgData.f[0] = Adaptive_Cruise_Control_Syst_M->odeF[0];
  Adaptive_Cruise_Control_Syst_M->intgData.f[1] = Adaptive_Cruise_Control_Syst_M->odeF[1];
  Adaptive_Cruise_Control_Syst_M->intgData.f[2] = Adaptive_Cruise_Control_Syst_M->odeF[2];
  Adaptive_Cruise_Control_Syst_M->contStates = ((X_Adaptive_Cruise_Control_Sys_T *) &Adaptive_Cruise_Control_Syste_X);
  rtsiSetSolverData(&Adaptive_Cruise_Control_Syst_M->solverInfo, (void *)&Adaptive_Cruise_Control_Syst_M->intgData);
  rtsiSetSolverName(&Adaptive_Cruise_Control_Syst_M->solverInfo,"ode3");
  rtmSetTPtr(Adaptive_Cruise_Control_Syst_M, &Adaptive_Cruise_Control_Syst_M->Timing.tArray[0]);
  Adaptive_Cruise_Control_Syst_M->Timing.stepSize0 = 0.01;
  rtmSetFirstInitCond(Adaptive_Cruise_Control_Syst_M, 1);

  /* Start for InitialCondition: '<S26>/IC' */
  Adaptive_Cruise_Control_Syst_DW.IC_FirstOutputTime = true;

  /* Start for If: '<S15>/If' */
  Adaptive_Cruise_Control_Syst_DW.If_ActiveSubsystem = -1;

  /* Start for InitialCondition: '<S62>/IC' */
  Adaptive_Cruise_Control_Syst_DW.IC_FirstOutputTime_g = true;

  /* Start for If: '<S51>/If' */
  Adaptive_Cruise_Control_Syst_DW.If_ActiveSubsystem_i = -1;
  Adaptive_Cruise_Control_PrevZCX.VelocitiesMatch_Input_ZCE = UNINITIALIZED_ZCSIG;
  Adaptive_Cruise_Control_PrevZCX.VelocitiesMatch_Input_ZCE_a = UNINITIALIZED_ZCSIG;

  /* InitializeConditions for Derivative: '<Root>/Derivative' */
  Adaptive_Cruise_Control_Syst_DW.TimeStampA = (rtInf);
  Adaptive_Cruise_Control_Syst_DW.TimeStampB = (rtInf);

  /* InitializeConditions for Integrator: '<S6>/Integrator' */
  Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S6>/Filter' */
  Adaptive_Cruise_Control_Syste_X.Filter_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S7>/Integrator' */
  Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE_c = 0.0;

  /* InitializeConditions for Integrator: '<S7>/Filter' */
  Adaptive_Cruise_Control_Syste_X.Filter_CSTATE_p = 0.0;

  /* InitializeConditions for Integrator: '<S5>/Integrator' */
  Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE_j = 0.0;

  /* InitializeConditions for Integrator: '<S5>/Filter' */
  Adaptive_Cruise_Control_Syste_X.Filter_CSTATE_a = 0.0;

  /* InitializeConditions for Integrator: '<S43>/Integrator' */
  Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE_o = 0.0;

  /* InitializeConditions for Integrator: '<Root>/Integrator1' */
  Adaptive_Cruise_Control_Syste_X.Integrator1_CSTATE = 8.33;

  /* InitializeConditions for Integrator: '<S8>/Integrator' */
  Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE_p = 0.0;

  /* InitializeConditions for Integrator: '<S8>/Filter' */
  Adaptive_Cruise_Control_Syste_X.Filter_CSTATE_ao = 0.0;

  /* InitializeConditions for Integrator: '<S79>/Integrator' */
  Adaptive_Cruise_Control_Syste_X.Integrator_CSTATE_g = 0.0;

  /* InitializeConditions for Integrator: '<Root>/Integrator2' */
  Adaptive_Cruise_Control_Syste_X.Integrator2_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<Root>/Integrator3' */
  Adaptive_Cruise_Control_Syste_X.Integrator3_CSTATE = 25.0;

  /* SystemInitialize for IfAction SubSystem: '<S51>/Unlocked' */
  /* SystemInitialize for IfAction SubSystem: '<S15>/Unlocked' */
  /* InitializeConditions for Integrator: '<S28>/Output Integrator' incorporates:
   *  InitializeConditions for Integrator: '<S64>/Output Integrator'
   */
  if (rtmIsFirstInitCond(Adaptive_Cruise_Control_Syst_M)) {
    Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE_j = 0.0;
    Adaptive_Cruise_Control_Syste_X.OutputIntegrator_CSTATE = 0.0;
  }

  /* End of SystemInitialize for SubSystem: '<S51>/Unlocked' */
  Adaptive_Cruise_Control_Syst_DW.OutputIntegrator_IWORK_b = 1;

  /* End of InitializeConditions for Integrator: '<S28>/Output Integrator' */
  /* End of SystemInitialize for SubSystem: '<S15>/Unlocked' */

  /* SystemInitialize for IfAction SubSystem: '<S51>/Unlocked' */
  /* InitializeConditions for Integrator: '<S64>/Output Integrator' */
  Adaptive_Cruise_Control_Syst_DW.OutputIntegrator_IWORK = 1;

  /* End of SystemInitialize for SubSystem: '<S51>/Unlocked' */

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond(Adaptive_Cruise_Control_Syst_M)) {
    rtmSetFirstInitCond(Adaptive_Cruise_Control_Syst_M, 0);
  }
}

/* Model terminate function */
void Adaptive_Cruise_Control_System_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
