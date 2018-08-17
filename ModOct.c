#include "allvars.h"
#include "proto.h"
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spline.h>


////////////////////////////////////////////////////////////
//
// Differential equations
//
////////////////////////////////////////////////////////////


double da_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params){
  ////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }

  
  double da_in_dt_orb;
  double da_in_dt_tid;
 
  da_in_dt_orb = 0.0;
  
  da_in_dt_tid = -2.0*a_in*( W_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)
			     + W_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) ) - 
    2.0*a_in*pow(e_in,2)/(1.0-e_in*e_in) * ( V_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + 
					     V_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) );

  return (params.q_orb*da_in_dt_orb) + (params.q_tid*da_in_dt_tid);
  
}


double da_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  
  double da_out_dt_orb;
  double da_out_dt_tid;
 
  da_out_dt_orb = 0.0;
  da_out_dt_tid = 0.0;
  
  return (params.q_orb*da_out_dt_orb) + (params.q_tid*da_out_dt_tid);
}



double de_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params){
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  
  double de_in_dt_orb;
  double de_in_dt_tid;
  double I_tot;
  double B;
  double A;
  double cphi;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
 
  de_in_dt_orb = C2(m_A,m_B,m_C,a_in,a_out,e_out)*(1-e_in*e_in)/G_1(m_A,m_B,a_in,e_in) * 30.0*e_in*pow(sin(I_tot),2) * sin(2.0*w_in) +
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_out*(1-e_in*e_in)/G_1(m_A,m_B,a_in,e_in)*(35.0*cphi*pow(sin(I_tot),2)*e_in*e_in*sin(2.0*w_in) -
										 10.0*cos(I_tot)*pow(sin(I_tot),2)*cos(w_in)*sin(w_out)*(1.0-e_in*e_in) -
										 A*(sin(w_in)*cos(w_out)-cos(I_tot)*cos(w_in)*sin(w_out)));
  
  de_in_dt_tid = -(V_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + V_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz))*e_in;
  
  return (params.q_orb*de_in_dt_orb) + (params.q_tid*de_in_dt_tid);
}



double de_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  
  double de_out_dt_orb;
  double de_out_dt_tid;
  double I_tot;
  double B;
  double A;
  double cphi;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
 
  de_out_dt_orb = -C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*(1.0-e_out*e_out)/G_2(m_A,m_B,m_C,a_out,e_out) *
    (10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-e_in*e_in)*sin(w_in)*cos(w_out) +
     A*(cos(w_in)*sin(w_out)-cos(I_tot)*sin(w_in)*cos(w_out)));
  de_out_dt_tid = 0.0;
  
  return (params.q_orb*de_out_dt_orb) + (params.q_tid*de_out_dt_tid);
}




double dI_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params){
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  
  double dI_in_dt_orb;
  double dI_in_dt_tid;
  double I_tot;
  double B;
  double A;
  double cphi;
  double dGin_dt;
  double dGout_dt;
  double dHin_dt;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  

  dGin_dt  = -C2(m_A,m_B,m_C,a_in,a_out,e_out)*30.0*e_in*e_in*sin(2.0*w_in)*pow(sin(I_tot),2) +
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*e_out*(-35.0*e_in*e_in*pow(sin(I_tot),2)*sin(2.0*w_in)*cphi +
						 A*(sin(w_in)*cos(w_out)-cos(I_tot)*cos(w_in)*sin(w_out)) +
						 10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-e_in*e_in)*cos(w_in)*sin(w_out));
  
  dGout_dt = C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*e_out*(A*(cos(w_in)*sin(w_out)-cos(I_tot)*sin(w_in)*cos(w_out)) +
							  10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-e_in*e_in)*sin(w_in)*cos(w_out));
  dHin_dt  = (sin(I_out)*dGin_dt - sin(I_in)*dGout_dt)/sin(I_tot);

  dI_in_dt_orb = -1.0 * ( dHin_dt - dGin_dt * cos(I_in) ) / (sin(I_in) * G_1(m_A,m_B,a_in,e_in) );
  dI_in_dt_tid = ( X_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + X_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) ) * cos(w_in) - 
                 ( Y_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + Y_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) ) * sin(w_in);

  
  return (params.q_orb*dI_in_dt_orb) + (params.q_tid*dI_in_dt_tid);
}



double dI_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  
  double dI_out_dt_orb;
  double dI_out_dt_tid;
  double I_tot;
  double B;
  double A;
  double cphi;
  double dGin_dt;
  double dGout_dt;
  double dHin_dt;
  double dHout_dt;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
  
  dGin_dt  = -C2(m_A,m_B,m_C,a_in,a_out,e_out)*30.0*e_in*e_in*sin(2.0*w_in)*pow(sin(I_tot),2) +
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*e_out*(-35.0*e_in*e_in*pow(sin(I_tot),2)*sin(2.0*w_in)*cphi +
						 A*(sin(w_in)*cos(w_out)-cos(I_tot)*cos(w_in)*sin(w_out)) +
						 10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-e_in*e_in)*cos(w_in)*sin(w_out));
  dGout_dt = C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*e_out*(A*(cos(w_in)*sin(w_out)-cos(I_tot)*sin(w_in)*cos(w_out)) +
							  10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-e_in*e_in)*sin(w_in)*cos(w_out));
  dHin_dt  = (sin(I_out)*dGin_dt - sin(I_in)*dGout_dt)/sin(I_tot);
  dHout_dt = -1.0*dHin_dt;
    
  dI_out_dt_orb = -1.0 * (dHout_dt - dGout_dt * cos(I_out))/(sin(I_out)*G_2(m_A,m_B,m_C,a_out,e_out));
  dI_out_dt_tid = 0.0;
  
  return (params.q_orb*dI_out_dt_orb) + (params.q_tid*dI_out_dt_tid);
}




double dW_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params){
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  
  double dW_in_dt_orb;
  double dW_in_dt_tid;
  double I_tot;
  double B;
  double A;
  double cphi;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
 
  dW_in_dt_orb = -3.0*C2(m_A,m_B,m_C,a_in,a_out,e_out)/(G_1(m_A,m_B,a_in,e_in)*sin(I_in)) * (2.0+3.0*e_in*e_in-5.0*e_in*e_in*cos(2.0*w_in))*sin(2.0*I_tot) -
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*e_out*(5.0*B*cos(I_tot)*cphi - A*sin(w_in)*sin(w_out) +
						 10.0*(1.0-3.0*pow(cos(I_tot),2))*(1.0-e_in*e_in)*sin(w_in)*sin(w_out))*sin(I_tot)/(G_1(m_A,m_B,a_in,e_in)*sin(I_in));

  dW_in_dt_tid = ( X_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + X_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) )
    * sin(w_in)/sin(I_tot) +				
    ( Y_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + Y_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) )
    * cos(w_in)/sin(I_tot);

  return (params.q_orb * dW_in_dt_orb) + (params.q_tid * dW_in_dt_tid);

}


// Ojo aqui
double dW_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  double dW_out_dt_orb;
  double I_tot;
  double B;
  double A;
  double cphi;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
  
  dW_out_dt_orb = -3.0*C2(m_A,m_B,m_C,a_in,a_out,e_out)/(G_1(m_A,m_B,a_in,e_in)*sin(I_in)) * (2.0+3.0*e_in*e_in-5.0*e_in*e_in*cos(2.0*w_in))*sin(2.0*I_tot) -
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*e_out*(5.0*B*cos(I_tot)*cphi - A*sin(w_in)*sin(w_out) +
						 10.0*(1.0-3.0*pow(cos(I_tot),2))*(1.0-e_in*e_in)*sin(w_in)*sin(w_out))*sin(I_tot)/(G_1(m_A,m_B,a_in,e_in)*sin(I_in));

  return (params.q_orb * dW_out_dt_orb);
  
}


double dw_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params){
  
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }

  double cvel = 300.0e6 * params.uT/params.uL;
  
  
  double dw_in_dt_orb;
  double dw_in_dt_tid;
  double dw_in_dt_GR;
  double I_tot;
  double B;
  double A;
  double cphi;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
  
  dw_in_dt_orb = 6.0*C2(m_A,m_B,m_C,a_in,a_out,e_out)*(1.0/G_1(m_A,m_B,a_in,e_in)*(4.0*pow(cos(I_tot),2)
										   + (5.0*cos(2.0*w_in)-1.0)*(1.0-e_in*e_in-pow(cos(I_tot),2))) + 
						       cos(I_tot)/G_2(m_A,m_B,m_C,a_out,e_out) * (2.0+e_in*e_in*(3.0-5.0*cos(2.0*w_in)))) -
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_out*(e_in*(1.0/G_2(m_A,m_B,m_C,a_out,e_out) + cos(I_tot)/G_1(m_A,m_B,a_in,e_in)) *
					    (sin(w_in)*sin(w_out)*(10.0*(3.0*pow(cos(I_tot),2)-1.0)*(1.0-e_in*e_in)+A) - 5.0*B*cos(I_tot)*cphi) -
					    (1.0-e_in*e_in)/(e_in*G_1(m_A,m_B,a_in,e_in))*(sin(w_in)*sin(w_out)*10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-3.0*e_in*e_in) +
											   cphi*(3.0*A-10.0*pow(cos(I_tot),2)+2.0))); 

  
  dw_in_dt_tid = Z_A(R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + Z_B(R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)  - 
    ( X_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + X_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) ) * sin(w_in)*cos(I_tot)/sin(I_tot) - 
    ( Y_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az) + Y_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz) ) * cos(w_in)*cos(I_tot)/sin(I_tot);
		

  dw_in_dt_GR = 3.0*pow(G*(m_A+m_B),1.5)/(pow(a_in,2.5)*cvel*cvel*(1.0-e_in*e_in));

  return (params.q_orb * dw_in_dt_orb) + (params.q_tid * dw_in_dt_tid) + (params.q_GR * dw_in_dt_GR);


}




double dw_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  
  double dw_out_dt_orb;
  double dw_out_dt_tid;
  double I_tot;
  double B;
  double A;
  double cphi;

  I_tot = I_in + I_out;
  B     = 2.0 + 5.0*e_in*e_in - 7.0*e_in*e_in*cos(2.0*w_in);
  A     = 4.0 + 3.0*e_in*e_in - 2.5*B*pow(sin(I_tot),2);
  cphi  = -cos(w_in)*cos(w_out) - cos(I_tot)*sin(w_in)*sin(w_out);  
  
  dw_out_dt_orb = 3.0*C2(m_A,m_B,m_C,a_in,a_out,e_out)*(2.0*cos(I_tot)/G_1(m_A,m_B,a_in,e_in) * (2.0+e_in*e_in*(3.0-5.0*cos(2.0*w_in))) + 
							1.0/G_2(m_A,m_B,m_C,a_out,e_out) * (4.0+6.0*e_in*e_in+(5.0*pow(cos(I_tot),2)-3.0)*
											    (2.0+e_in*e_in*(3.0-5.0*cos(2.0*w_in))))) +
    C3(m_A,m_B,m_C,a_in,a_out,e_out)*e_in*(sin(w_in)*sin(w_out)*((4.0*e_out*e_out + 1.0)/(e_out*G_2(m_A,m_B,m_C,a_out,e_out))*
								 10.0*cos(I_tot)*pow(sin(I_tot),2)*(1.0-e_in*e_in) -
								 e_out*(1.0/G_1(m_A,m_B,a_in,e_in) + cos(I_tot)/G_2(m_A,m_B,m_C,a_out,e_out)) *
								 (A+10.0*(3.0*pow(cos(I_tot),2)-1.0)*(1.0-e_in*e_in))) +
					   cphi*(5.0*B*cos(I_tot)*e_out*(1.0/G_1(m_A,m_B,a_in,e_in) + cos(I_tot)/G_2(m_A,m_B,m_C,a_out,e_out)) +
						 (4.0*e_out*e_out + 1.0)/(e_out*G_2(m_A,m_B,m_C,a_out,e_out))*A));
  
  dw_out_dt_tid = 0.0;
  
  return (params.q_orb * dw_out_dt_orb) + (params.q_tid * dw_out_dt_tid);

}





double dOm_Ax_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  
  double mu;

  mu = m_A*m_B/(m_A+m_B);

  
  return params.q_tid*( mu*pow(G*(m_A+m_B)*a_in*(1.0-e_in*e_in),0.5)/(gyr_rad_A*m_A*pow(R_A,2)) * 
    ( X_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(-cos(W_in)*sin(w_in) - sin(W_in)*cos(w_in)*cos(I_in)) + 
      W_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(sin(I_in)*sin(W_in)) - 
      Y_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(cos(W_in)*cos(w_in)-sin(W_in)*sin(w_in)*cos(I_in))));
  
}




double dOm_Ay_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
  
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  double mu;

  mu = m_A*m_B/(m_A+m_B);

  return params.q_tid * ( mu*pow(G*(m_A+m_B)*a_in*(1.0-e_in*e_in),0.5)/(gyr_rad_A*m_A*pow(R_A,2)) * 
    ( X_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(-sin(W_in)*sin(w_in) + cos(W_in)*cos(w_in)*cos(I_in)) + 
      W_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(-sin(I_in)*cos(W_in)) - 
      Y_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(sin(W_in)*cos(w_in)+cos(W_in)*sin(w_in)*cos(I_in))) ) ;

}





double dOm_Az_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  
  double mu;

  mu = m_A*m_B/(m_A+m_B);
  
  return params.q_tid * ( mu*pow(G*(m_A+m_B)*a_in*(1.0-e_in*e_in),0.5)/(gyr_rad_A*m_A*pow(R_A,2)) * 
			  ( X_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(sin(I_in)*cos(w_in)) + 
			    W_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*cos(I_in) - 
			    Y_A(tv_A,R_A,k_A,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Ax,Om_Ay,Om_Az)*(sin(I_in)*sin(w_in))) );

}





double dOm_Bx_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  
  double mu;

  mu = m_A*m_B/(m_A+m_B);

  return params.q_tid*( mu*pow(G*(m_A+m_B)*a_in*(1.0-e_in*e_in),0.5)/(gyr_rad_B*m_B*pow(R_B,2)) * 
    ( X_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(-cos(W_in)*sin(w_in) - sin(W_in)*cos(w_in)*cos(I_in)) + 
      W_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(sin(I_in)*sin(W_in)) - 
      Y_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(cos(W_in)*cos(w_in)-sin(W_in)*sin(w_in)*cos(I_in))));

}




double dOm_By_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  
  double mu;

  mu = m_A*m_B/(m_A+m_B);

  return params.q_tid * ( mu*pow(G*(m_A+m_B)*a_in*(1.0-e_in*e_in),0.5)/(gyr_rad_B*m_B*pow(R_B,2)) * 
    ( X_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(-sin(W_in)*sin(w_in) + cos(W_in)*cos(w_in)*cos(I_in)) + 
      W_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(-sin(I_in)*cos(W_in)) - 
      Y_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(sin(W_in)*cos(w_in)+cos(W_in)*sin(w_in)*cos(I_in))) ) ;

}





double dOm_Bz_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params){
////////////////////////////////////////////////////////////
  // bulk properties
  double tv_A;
  double tv_B;
  double tv_C;
  double R_A;  
  double R_B;
  double R_C;  
  double k_A;  
  double k_B;
  double k_C;  
  double m_A;  
  double m_B;
  double m_C;
  double gyr_rad_A;
  double gyr_rad_B;
  double gyr_rad_C;
  
  if (params.rse==0){
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = params.R_A;
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  else{
    tv_A      = params.tv_A;
    tv_B      = params.tv_B;
    tv_C      = params.tv_C;
    R_A       = fn_R_A(params,t,tim_A,rad_A);
    R_B       = params.R_B;
    R_C       = params.R_C;
    k_A       = params.k_A;
    k_B       = params.k_B;
    k_C       = params.k_C;
    m_A       = params.m_A;
    m_B       = params.m_B;
    m_C       = params.m_C;
    gyr_rad_A = params.gyr_rad_A;
    gyr_rad_B = params.gyr_rad_B;
    gyr_rad_C = params.gyr_rad_C;
  }
  
  
  double mu;

  mu = m_A*m_B/(m_A+m_B);
  
  return params.q_tid * ( mu*pow(G*(m_A+m_B)*a_in*(1.0-e_in*e_in),0.5)/(gyr_rad_B*m_B*pow(R_B,2)) * 
			  ( X_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(sin(I_in)*cos(w_in)) + 
			    W_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*cos(I_in) - 
			    Y_B(tv_B,R_B,k_B,m_A,m_B,a_in,e_in,W_in,I_in,w_in,Om_Bx,Om_By,Om_Bz)*(sin(I_in)*sin(w_in))) );

}




////////////////////////////////////////////////////////////
//
// State vector function
//
////////////////////////////////////////////////////////////


int func (double t, const double y[], double f[], void *params){
  (void)(t); // avoid unused parameter warning

  Inpar parameters = *(Inpar *)params;
    
  double a_in  = y[0];
  double a_out = y[1];
  double e_in  = y[2];
  double e_out = y[3];
  double I_in  = y[4];
  double I_out = y[5];
  double W_in  = y[6];
  double W_out = y[7];
  double w_in  = y[8];
  double w_out = y[9];
  double Om_Ax = y[10];
  double Om_Ay = y[11];
  double Om_Az = y[12];
  double Om_Bx = y[13];
  double Om_By = y[14];
  double Om_Bz = y[15];

  f[0]  = da_in_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[1]  = da_out_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[2]  = de_in_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[3]  = de_out_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[4]  = dI_in_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[5]  = dI_out_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[6]  = dW_in_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[7]  = dW_out_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[8]  = dw_in_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[9]  = dw_out_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[10] = dOm_Ax_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[11] = dOm_Ay_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[12] = dOm_Az_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[13] = dOm_Bx_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[14] = dOm_By_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  f[15] = dOm_Bz_dt(a_in,a_out,e_in,e_out,I_in,I_out,W_in,W_out,w_in,w_out,Om_Ax,Om_Ay,Om_Az,Om_Bx,Om_By,Om_Bz,t,parameters);
  
  return GSL_SUCCESS;
  
}


