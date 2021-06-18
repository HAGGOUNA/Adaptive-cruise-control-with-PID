/*
 * File: Adaptive_Cruise_Control_System_data.c
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

/* Invariant block signals (auto storage) */
const ConstB_Adaptive_Cruise_Contro_T Adaptive_Cruise_Control__ConstB = {
  0.00013939520000000002,              /* '<S1>/Product' */
  450.0,                               /* '<S3>/Product' */
  0.039269908169872414,                /* '<S34>/Torque Conversion1' */
  0.0,                                 /* '<S40>/Output Damping' */
  300.0,                               /* '<S3>/Product2' */
  0.039269908169872414,                /* '<S70>/Torque Conversion1' */
  0.0                                  /* '<S76>/Output Damping' */
};

/* Constant parameters (auto storage) */
const ConstP_Adaptive_Cruise_Contro_T Adaptive_Cruise_Control__ConstP = {
  /* Expression: OutValues
   * Referenced by: '<S4>/Lookup'
   */
  { 30.0, 30.0, 100.0, 100.0, 60.0, 60.0, 100.0, 100.0, 180.0 },

  /* Expression: TimeValues
   * Referenced by: '<S4>/Lookup'
   */
  { 0.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 110.0, 120.0 },

  /* Pooled Parameter (Expression: zeros(1,3))
   * Referenced by:
   *   '<S17>/Constant12'
   *   '<S17>/Constant19'
   *   '<S21>/Constant12'
   *   '<S21>/Constant19'
   *   '<S53>/Constant12'
   *   '<S53>/Constant19'
   *   '<S57>/Constant12'
   *   '<S57>/Constant19'
   */
  { 0.0, 0.0, 0.0 },

  /* Pooled Parameter (Expression: zeros(3,3))
   * Referenced by:
   *   '<S17>/Constant14'
   *   '<S21>/Constant14'
   *   '<S53>/Constant14'
   *   '<S57>/Constant14'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /* Pooled Parameter (Expression: [0;1;0;0;1;1;1;0])
   * Referenced by:
   *   '<S39>/Combinatorial  Logic'
   *   '<S75>/Combinatorial  Logic'
   */
  { 0, 1, 0, 0, 1, 1, 1, 0 }
};

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
