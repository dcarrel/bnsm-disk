#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     POTENTIAL
#define  COOLING                        TABULATED
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            7

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      EXPLICIT
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  MBH                            0
#define  NGAM                           1
#define  TSTART                         2
#define  ALPHA                          3
#define  ESRC                           4
#define  RHOCRIT                        5
#define  EXP_LIM                        6

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   1.e8
#define  UNIT_LENGTH                    1.e7
#define  UNIT_VELOCITY                  1.e8
#define  SHOCK_FLATTENING               MULTID
#define  FAILSAFE                       YES
#define  EPS_SHOCK_FLATTENING           1.0

/* [End] user-defined constants (do not change this line) */
