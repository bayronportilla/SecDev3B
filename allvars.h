#include <stdio.h>
#include <gsl/gsl_spline.h>

#define AU     149.6e9
#define MS     1.989e30
#define RS     6.957e8
#define MJ     1.898e27
#define RJ     6.9911e7
#define ME     5.972e24
#define RE     6.371e6
#define YEARS  (365.25*86400)
#define Myr    (1.0e6*YEARS)
#define DAYS   86400.0 
#define G      1.0
#define PI     3.14159265358979323
#define InDeg  (180.0/PI)
#define cspeed 300.0e6
#define Gyr    (1.0e9*YEARS)

struct Inpar_params{
  
  // Inner A properties
  double m_A;
  double R_A;
  double k_A;
  double tv_A;
  double gyr_rad_A;
  double a_in;
  double e_in;
  double I_in;
  double W_in;
  double w_in;
  double P_rot_A;
  double alpha_A;
  double beta_A;

  // Inner B properties
  double m_B ;
  double R_B ;
  double k_B ;
  double tv_B;
  double gyr_rad_B;
  double P_rot_B;
  double alpha_B;
  double beta_B;

  // Inner C properties
  double m_C;
  double R_C;
  double k_C;
  double tv_C;
  double gyr_rad_C;
  double a_out;
  double e_out;
  double I_out;
  double W_out;
  double w_out;

  // General properties
  const char *sim_name;
  double t_ini;
  double t_end;
  double h_output;
  double I_tot;
  int q_orb;
  int q_tid;
  int q_GR;
  int rcu;
  int rse;

  
  // Rotational velocities
  double Om_Ax;
  double Om_Ay;
  double Om_Az;
  double Om_Bx;
  double Om_By;
  double Om_Bz;
  
  
  // Canonical units info.
  double uM;
  double uL;
  double uT;
 
};

typedef struct Inpar_params Inpar;
extern Inpar st;
extern int Nlines;
extern double *tim_A;
extern double *tim_B;
extern double *rad_A;
extern double *rad_B;

extern gsl_interp_accel *acc;
extern gsl_spline       *spline;





/*  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spline.h>
#include <libconfig.h>
#include "Units.h"
#include "matrix.h"

#include "proto.h"
#include "Params.h"
#include "interpol.h"
#include "bulk.h"
#include "generic.h"

#include "FetchInfo.h"
#include "ModOct.h"
//#include "ModQuad.h"
*/




