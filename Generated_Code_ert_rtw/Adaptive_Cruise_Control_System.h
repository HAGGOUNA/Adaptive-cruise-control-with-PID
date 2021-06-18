/*
 * File: Adaptive_Cruise_Control_System.h
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

#ifndef RTW_HEADER_Adaptive_Cruise_Control_System_h_
#define RTW_HEADER_Adaptive_Cruise_Control_System_h_
#include <math.h>
#include <string.h>
#ifndef Adaptive_Cruise_Control_System_COMMON_INCLUDES_
# define Adaptive_Cruise_Control_System_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "zero_crossing_types.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* Adaptive_Cruise_Control_System_COMMON_INCLUDES_ */

#include "Adaptive_Cruise_Control_System_types.h"
#include "rtGetInf.h"
#include "rt_nonfinite.h"
#include "rt_zcfcn.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

/* Block signals for system '<S12>/Simple Magic Tire' */
typedef struct {
  real_T Fx;                           /* '<S12>/Simple Magic Tire' */
  real_T My;                           /* '<S12>/Simple Magic Tire' */
} B_SimpleMagicTire_Adaptive_Cr_T;

/* Block signals (auto storage) */
typedef struct {
  real_T velocity_lead_carms;          /* '<Root>/Gain' */
  real_T Gain3;                        /* '<Root>/Gain3' */
  real_T FilterCoefficient;            /* '<S6>/Filter Coefficient' */
  real_T UnitDelay1;                   /* '<Root>/Unit Delay1' */
  real_T FilterCoefficient_j;          /* '<S7>/Filter Coefficient' */
  real_T FilterCoefficient_n;          /* '<S5>/Filter Coefficient' */
  real_T Switch;                       /* '<S26>/Switch' */
  real_T omega;                        /* '<S15>/Merge' */
  real_T velocity_ego_carms;           /* '<Root>/Integrator1' */
  real_T FilterCoefficient_c;          /* '<S8>/Filter Coefficient' */
  real_T Switch_n;                     /* '<S62>/Switch' */
  real_T omega_c;                      /* '<S51>/Merge' */
  real_T Sum;                          /* '<S2>/Sum' */
  real_T Acceleration_ego_car;         /* '<Root>/Divide' */
  real_T Product;                      /* '<S43>/Product' */
  real_T Product_d;                    /* '<S79>/Product' */
  real_T IntegralGain;                 /* '<S5>/Integral Gain' */
  real_T IntegralGain_o;               /* '<S6>/Integral Gain' */
  real_T IntegralGain_e;               /* '<S7>/Integral Gain' */
  real_T IntegralGain_b;               /* '<S8>/Integral Gain' */
  real_T OutputInertia;                /* '<S64>/Output Inertia' */
  real_T OutputInertia_p;              /* '<S28>/Output Inertia' */
  boolean_T VelocitiesMatch;           /* '<S38>/Velocities Match' */
  boolean_T Memory;                    /* '<S39>/Memory' */
  boolean_T CombinatorialLogic;        /* '<S39>/Combinatorial  Logic' */
  boolean_T VelocitiesMatch_a;         /* '<S74>/Velocities Match' */
  boolean_T Memory_b;                  /* '<S75>/Memory' */
  boolean_T CombinatorialLogic_c;      /* '<S75>/Combinatorial  Logic' */
  B_SimpleMagicTire_Adaptive_Cr_T sf_SimpleMagicTire_g;/* '<S48>/Simple Magic Tire' */
  B_SimpleMagicTire_Adaptive_Cr_T sf_SimpleMagicTire;/* '<S12>/Simple Magic Tire' */
} B_Adaptive_Cruise_Control_Sys_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T UnitDelay_DSTATE;             /* '<Root>/Unit Delay' */
  real_T UnitDelay1_DSTATE;            /* '<Root>/Unit Delay1' */
  real_T TimeStampA;                   /* '<Root>/Derivative' */
  real_T LastUAtTimeA;                 /* '<Root>/Derivative' */
  real_T TimeStampB;                   /* '<Root>/Derivative' */
  real_T LastUAtTimeB;                 /* '<Root>/Derivative' */
  int_T OutputIntegrator_IWORK;        /* '<S64>/Output Integrator' */
  int_T OutputIntegrator_IWORK_b;      /* '<S28>/Output Integrator' */
  int_T VelocitiesMatch_MODE;          /* '<S38>/Velocities Match' */
  int_T VelocitiesMatch_MODE_j;        /* '<S74>/Velocities Match' */
  uint16_T Output_DSTATE;              /* '<S83>/Output' */
  int8_T If_ActiveSubsystem;           /* '<S15>/If' */
  int8_T If_ActiveSubsystem_i;         /* '<S51>/If' */
  boolean_T Memory_PreviousInput;      /* '<S39>/Memory' */
  boolean_T IC_FirstOutputTime;        /* '<S26>/IC' */
  boolean_T Memory_PreviousInput_a;    /* '<S75>/Memory' */
  boolean_T IC_FirstOutputTime_g;      /* '<S62>/IC' */
} DW_Adaptive_Cruise_Control_Sy_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T Integrator_CSTATE;            /* '<S6>/Integrator' */
  real_T Filter_CSTATE;                /* '<S6>/Filter' */
  real_T Integrator_CSTATE_c;          /* '<S7>/Integrator' */
  real_T Filter_CSTATE_p;              /* '<S7>/Filter' */
  real_T Integrator_CSTATE_j;          /* '<S5>/Integrator' */
  real_T Filter_CSTATE_a;              /* '<S5>/Filter' */
  real_T Integrator_CSTATE_o;          /* '<S43>/Integrator' */
  real_T Integrator1_CSTATE;           /* '<Root>/Integrator1' */
  real_T Integrator_CSTATE_p;          /* '<S8>/Integrator' */
  real_T Filter_CSTATE_ao;             /* '<S8>/Filter' */
  real_T Integrator_CSTATE_g;          /* '<S79>/Integrator' */
  real_T Integrator2_CSTATE;           /* '<Root>/Integrator2' */
  real_T Integrator3_CSTATE;           /* '<Root>/Integrator3' */
  real_T OutputIntegrator_CSTATE;      /* '<S64>/Output Integrator' */
  real_T OutputIntegrator_CSTATE_j;    /* '<S28>/Output Integrator' */
} X_Adaptive_Cruise_Control_Sys_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T Integrator_CSTATE;            /* '<S6>/Integrator' */
  real_T Filter_CSTATE;                /* '<S6>/Filter' */
  real_T Integrator_CSTATE_c;          /* '<S7>/Integrator' */
  real_T Filter_CSTATE_p;              /* '<S7>/Filter' */
  real_T Integrator_CSTATE_j;          /* '<S5>/Integrator' */
  real_T Filter_CSTATE_a;              /* '<S5>/Filter' */
  real_T Integrator_CSTATE_o;          /* '<S43>/Integrator' */
  real_T Integrator1_CSTATE;           /* '<Root>/Integrator1' */
  real_T Integrator_CSTATE_p;          /* '<S8>/Integrator' */
  real_T Filter_CSTATE_ao;             /* '<S8>/Filter' */
  real_T Integrator_CSTATE_g;          /* '<S79>/Integrator' */
  real_T Integrator2_CSTATE;           /* '<Root>/Integrator2' */
  real_T Integrator3_CSTATE;           /* '<Root>/Integrator3' */
  real_T OutputIntegrator_CSTATE;      /* '<S64>/Output Integrator' */
  real_T OutputIntegrator_CSTATE_j;    /* '<S28>/Output Integrator' */
} XDot_Adaptive_Cruise_Control__T;

/* State disabled  */
typedef struct {
  boolean_T Integrator_CSTATE;         /* '<S6>/Integrator' */
  boolean_T Filter_CSTATE;             /* '<S6>/Filter' */
  boolean_T Integrator_CSTATE_c;       /* '<S7>/Integrator' */
  boolean_T Filter_CSTATE_p;           /* '<S7>/Filter' */
  boolean_T Integrator_CSTATE_j;       /* '<S5>/Integrator' */
  boolean_T Filter_CSTATE_a;           /* '<S5>/Filter' */
  boolean_T Integrator_CSTATE_o;       /* '<S43>/Integrator' */
  boolean_T Integrator1_CSTATE;        /* '<Root>/Integrator1' */
  boolean_T Integrator_CSTATE_p;       /* '<S8>/Integrator' */
  boolean_T Filter_CSTATE_ao;          /* '<S8>/Filter' */
  boolean_T Integrator_CSTATE_g;       /* '<S79>/Integrator' */
  boolean_T Integrator2_CSTATE;        /* '<Root>/Integrator2' */
  boolean_T Integrator3_CSTATE;        /* '<Root>/Integrator3' */
  boolean_T OutputIntegrator_CSTATE;   /* '<S64>/Output Integrator' */
  boolean_T OutputIntegrator_CSTATE_j; /* '<S28>/Output Integrator' */
} XDis_Adaptive_Cruise_Control__T;

/* Zero-crossing (trigger) state */
typedef struct {
  ZCSigState VelocitiesMatch_Input_ZCE;/* '<S38>/Velocities Match' */
  ZCSigState VelocitiesMatch_Input_ZCE_a;/* '<S74>/Velocities Match' */
} PrevZCX_Adaptive_Cruise_Contr_T;

/* Invariant block signals (auto storage) */
typedef struct {
  const real_T Product;                /* '<S1>/Product' */
  const real_T Product_p;              /* '<S3>/Product' */
  const real_T TorqueConversion1;      /* '<S34>/Torque Conversion1' */
  const real_T OutputDamping;          /* '<S40>/Output Damping' */
  const real_T Product2;               /* '<S3>/Product2' */
  const real_T TorqueConversion1_c;    /* '<S70>/Torque Conversion1' */
  const real_T OutputDamping_a;        /* '<S76>/Output Damping' */
} ConstB_Adaptive_Cruise_Contro_T;

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
} ODE3_IntgData;

#endif

/* Constant parameters (auto storage) */
typedef struct {
  /* Expression: OutValues
   * Referenced by: '<S4>/Lookup'
   */
  real_T Lookup_tableData[9];

  /* Expression: TimeValues
   * Referenced by: '<S4>/Lookup'
   */
  real_T Lookup_bp01Data[9];

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
  real_T pooled22[3];

  /* Pooled Parameter (Expression: zeros(3,3))
   * Referenced by:
   *   '<S17>/Constant14'
   *   '<S21>/Constant14'
   *   '<S53>/Constant14'
   *   '<S57>/Constant14'
   */
  real_T pooled23[9];

  /* Pooled Parameter (Expression: [0;1;0;0;1;1;1;0])
   * Referenced by:
   *   '<S39>/Combinatorial  Logic'
   *   '<S75>/Combinatorial  Logic'
   */
  boolean_T pooled27[8];
} ConstP_Adaptive_Cruise_Contro_T;

/* Real-time Model Data Structure */
struct tag_RTM_Adaptive_Cruise_Contr_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_Adaptive_Cruise_Control_Sys_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T blkStateChange;
  real_T odeY[15];
  real_T odeF[3][15];
  ODE3_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    time_T stepSize0;
    uint32_T clockTick1;
    boolean_T firstInitCondFlag;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

/* Block signals (auto storage) */
extern B_Adaptive_Cruise_Control_Sys_T Adaptive_Cruise_Control_Syste_B;

/* Continuous states (auto storage) */
extern X_Adaptive_Cruise_Control_Sys_T Adaptive_Cruise_Control_Syste_X;

/* Block states (auto storage) */
extern DW_Adaptive_Cruise_Control_Sy_T Adaptive_Cruise_Control_Syst_DW;
extern const ConstB_Adaptive_Cruise_Contro_T Adaptive_Cruise_Control__ConstB;/* constant block i/o */

/* Constant parameters (auto storage) */
extern const ConstP_Adaptive_Cruise_Contro_T Adaptive_Cruise_Control__ConstP;

/* Model entry point functions */
extern void Adaptive_Cruise_Control_System_initialize(void);
extern void Adaptive_Cruise_Control_System_step(void);
extern void Adaptive_Cruise_Control_System_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Adaptive_Cruise_Cont_T *const Adaptive_Cruise_Control_Syst_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S11>/Constant2' : Unused code path elimination
 * Block '<S27>/Unary Minus' : Unused code path elimination
 * Block '<S28>/Signal Conversion' : Unused code path elimination
 * Block '<S47>/Constant2' : Unused code path elimination
 * Block '<S63>/Unary Minus' : Unused code path elimination
 * Block '<S64>/Signal Conversion' : Unused code path elimination
 * Block '<S4>/Data Type Propagation' : Unused code path elimination
 * Block '<S83>/Data Type Propagation' : Unused code path elimination
 * Block '<S84>/FixPt Data Type Duplicate' : Unused code path elimination
 * Block '<S85>/FixPt Data Type Duplicate1' : Unused code path elimination
 * Block '<Root>/Scope1' : Unused code path elimination
 * Block '<Root>/Subtract' : Unused code path elimination
 * Block '<Root>/Scope' : Unused code path elimination
 * Block '<S24>/Ratio of static to kinetic' : Eliminated nontunable gain of 1
 * Block '<S60>/Ratio of static to kinetic' : Eliminated nontunable gain of 1
 * Block '<S4>/Output' : Eliminate redundant signal conversion block
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Adaptive_Cruise_Control_System'
 * '<S1>'   : 'Adaptive_Cruise_Control_System/Brk_pres_calc'
 * '<S2>'   : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics'
 * '<S3>'   : 'Adaptive_Cruise_Control_System/FZ_calc'
 * '<S4>'   : 'Adaptive_Cruise_Control_System/Leading vehicle  velocity [km//h]'
 * '<S5>'   : 'Adaptive_Cruise_Control_System/PID_Axial_trq'
 * '<S6>'   : 'Adaptive_Cruise_Control_System/PID_BRK'
 * '<S7>'   : 'Adaptive_Cruise_Control_System/PID_acceleration'
 * '<S8>'   : 'Adaptive_Cruise_Control_System/PID_long_velocity'
 * '<S9>'   : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake'
 * '<S10>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1'
 * '<S11>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake'
 * '<S12>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Longitudinal Basic Magic Tire'
 * '<S13>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Longitudinal Parameters'
 * '<S14>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Rolling Parameters'
 * '<S15>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module'
 * '<S16>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Longitudinal Basic Magic Tire/Simple Magic Tire'
 * '<S17>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Longitudinal Parameters/Magic Formula Peak Value'
 * '<S18>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Longitudinal Parameters/Mapped Force'
 * '<S19>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Rolling Parameters/Magic'
 * '<S20>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Rolling Parameters/Mapped'
 * '<S21>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Rolling Parameters/None'
 * '<S22>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Rolling Parameters/Simple'
 * '<S23>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Brakes'
 * '<S24>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Friction Model'
 * '<S25>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic'
 * '<S26>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/IC tunable'
 * '<S27>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Locked'
 * '<S28>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Unlocked'
 * '<S29>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/relaxation'
 * '<S30>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Brakes/Disk Brake'
 * '<S31>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Brakes/Drum Brake'
 * '<S32>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Brakes/Mapped Brake'
 * '<S33>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Brakes/No Brake'
 * '<S34>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Brakes/Disk Brake/Disk Brake'
 * '<S35>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Brakes/Drum Brake/Drum Brake'
 * '<S36>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Brakes/Drum Brake/Drum Brake/Drum Brake Torque Calculation'
 * '<S37>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic/Break Apart Detection'
 * '<S38>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic/Lockup Detection'
 * '<S39>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic/Lockup FSM'
 * '<S40>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic/Requisite Friction'
 * '<S41>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic/Lockup Detection/Friction Calc'
 * '<S42>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic/Lockup Detection/Required Friction for Lockup'
 * '<S43>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/relaxation/Cont LPF Dyn'
 * '<S44>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/relaxation/div0protect - abs poly'
 * '<S45>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/relaxation/div0protect - abs poly/Compare To Constant'
 * '<S46>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake/Longitudinal Wheel with Brake/Wheel Module/relaxation/div0protect - abs poly/Compare To Constant1'
 * '<S47>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake'
 * '<S48>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Longitudinal Basic Magic Tire'
 * '<S49>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Longitudinal Parameters'
 * '<S50>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Rolling Parameters'
 * '<S51>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module'
 * '<S52>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Longitudinal Basic Magic Tire/Simple Magic Tire'
 * '<S53>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Longitudinal Parameters/Magic Formula Peak Value'
 * '<S54>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Longitudinal Parameters/Mapped Force'
 * '<S55>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Rolling Parameters/Magic'
 * '<S56>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Rolling Parameters/Mapped'
 * '<S57>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Rolling Parameters/None'
 * '<S58>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Rolling Parameters/Simple'
 * '<S59>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Brakes'
 * '<S60>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Friction Model'
 * '<S61>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic'
 * '<S62>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/IC tunable'
 * '<S63>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Locked'
 * '<S64>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Unlocked'
 * '<S65>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/relaxation'
 * '<S66>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Brakes/Disk Brake'
 * '<S67>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Brakes/Drum Brake'
 * '<S68>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Brakes/Mapped Brake'
 * '<S69>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Brakes/No Brake'
 * '<S70>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Brakes/Disk Brake/Disk Brake'
 * '<S71>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Brakes/Drum Brake/Drum Brake'
 * '<S72>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Brakes/Drum Brake/Drum Brake/Drum Brake Torque Calculation'
 * '<S73>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic/Break Apart Detection'
 * '<S74>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic/Lockup Detection'
 * '<S75>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic/Lockup FSM'
 * '<S76>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic/Requisite Friction'
 * '<S77>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic/Lockup Detection/Friction Calc'
 * '<S78>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/Friction Mode Logic/Lockup Detection/Required Friction for Lockup'
 * '<S79>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/relaxation/Cont LPF Dyn'
 * '<S80>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/relaxation/div0protect - abs poly'
 * '<S81>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/relaxation/div0protect - abs poly/Compare To Constant'
 * '<S82>'  : 'Adaptive_Cruise_Control_System/Ego vehicle wheel dynamics/Longitudinal Wheel - Disk Brake1/Longitudinal Wheel with Brake/Wheel Module/relaxation/div0protect - abs poly/Compare To Constant1'
 * '<S83>'  : 'Adaptive_Cruise_Control_System/Leading vehicle  velocity [km//h]/LimitedCounter'
 * '<S84>'  : 'Adaptive_Cruise_Control_System/Leading vehicle  velocity [km//h]/LimitedCounter/Increment Real World'
 * '<S85>'  : 'Adaptive_Cruise_Control_System/Leading vehicle  velocity [km//h]/LimitedCounter/Wrap To Zero'
 */
#endif                                 /* RTW_HEADER_Adaptive_Cruise_Control_System_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
