#include <math.h>
#include "allvars.h"
#include "proto.h"
//#include "matrix.c"


////////////////////////////////////////////////////////////
//
// Definition of generic functions
//
////////////////////////////////////////////////////////////


int sign(double x){
  if(x > 0.0){
    return 1;
  }
  else if(x<0.0){
    return -1;
  }
  else{
    return 0;
  }
}


////////////////////////////////////////////////////////////
// Following equations are from Naoz 2016 and Fabricky and
// Tremaine 2007.

double f_1(double e){
  return ( 1.0 + 3.0*pow(e,2) + (3.0/8.0)*pow(e,4) ) / pow(1.0-e*e,5);
}


double f_2(double e){
  return ( 1.0 + (15.0/2.0)*pow(e,2) + (45.0/8.0)*pow(e,4) + (5.0/16.0)*pow(e,6) ) / pow(1.0-e*e,13.0/2.0); 
}

double f_3(double e){
  return ( 1.0 + (9.0/2.0)*pow(e,2) + (5.0/8.0)*pow(e,4) ) / pow(1.0-e*e,5);
}

double f_4(double e){
  return ( 1.0 + (3.0/2.0)*pow(e,2) + (1.0/8.0)*pow(e,4) ) / pow(1.0-e*e,5);
}

double f_5(double e){
  return ( 1.0 + (15.0/4.0)*pow(e,2) + (15.0/8.0)*pow(e,4) + (5.0/64.0)*pow(e,6) ) / pow(1.0-e*e,13.0/2.0);
}

double n(double m_1, double m_2, double a){
  return pow(G*(m_1+m_2)/pow(a,3),0.5);
}

double L_1(double m_1, double m_2, double a_1){
  return m_1*m_2/(m_1+m_2) * pow(G*(m_1+m_2)*a_1,0.5);
}

double L_2(double m_1,double m_2,double m_3, double a_2){
  return m_3*(m_1+m_2)/(m_1+m_2+m_3) * pow(G*(m_1+m_2+m_3)*a_2,0.5);
}

double G_1(double m_1,double m_2,double a_1, double e_1){
  return L_1(m_1,m_2,a_1)*pow(1.0-e_1*e_1,0.5);
}

double G_2(double m_1,double m_2,double m_3, double a_2, double e_2){
  return L_2(m_1,m_2,m_3,a_2)*pow(1.0-e_2*e_2,0.5);
}

double H_1(double m_1,double m_2,double a_1,double e_1,double I_1){
  return G_1(m_1,m_2,a_1,e_1)*cos(I_1);
}

double H_2(double m_1,double m_2,double m_3,double a_2,double e_2,double I_2 ){
  return G_2(m_1,m_2,m_3,a_2,e_2)*cos(I_2);
}

double C2(double m_1,double m_2,double m_3,double a_1,double a_2,double e_2){
return (G*G/16.0) * pow(m_1+m_2,7)/pow(m_1+m_2+m_3,3) * pow(m_3,7)/pow(m_1*m_2,3) * pow(L_1(m_1,m_2,a_1),4)/(pow(L_2(m_1,m_2,m_3,a_2),3) * pow(G_2(m_1,m_2,m_3,a_2,e_2),3));
}

double epsilon_m(double m_1,double m_2,double a_1,double a_2,double e_2){
return (m_1-m_2)/(m_1+m_2) * (a_1/a_2) * e_2/(1.0-e_2*e_2);
}

// Ojo, si es negativo? ver eq 24 de Naoz 2013
double C3(double m_1,double m_2,double m_3,double a_1,double a_2,double e_2){
return -C2(m_1,m_2,m_3,a_1,a_2,e_2)*(15.0/4.0) * epsilon_m(m_1,m_2,a_1,a_2,e_2)/e_2;
}

double tf_A(double tv_A, double R_A, double k_A, double a, double m_A, double m_B){
  return tv_A/9.0 * pow(a/R_A,8) * m_A*m_A/((m_A+m_B)*m_B) * 1.0/pow(1.0+2.0*k_A,2);
}

double tf_B(double tv_B, double R_B, double k_B, double a, double m_A, double m_B){
  return tv_B/9.0 * pow(a/R_B,8) * m_B*m_B/((m_A+m_B)*m_A) * 1.0/pow(1.0+2.0*k_B,2);
}

double V_A(double tv_A, double R_A, double k_A, double m_A, double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay, double Om_Az){
  
  double Om_A_in[3] = {Om_Ax,Om_Ay,Om_Az};
  double Om_Az_orb = RotMat(W,I,w,Om_A_in,"InOrb")[2];
  
return 9.0/tf_A(tv_A,R_A,k_A,a,m_A,m_B) * (f_5(e) - (11.0/18.0)*Om_Az_orb/n(m_A,m_B,a) * f_4(e));
}

double V_B(double tv_B, double R_B, double k_B, double m_A, double m_B, 
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By, double Om_Bz){
  
  double Om_B_in[3] = {Om_Bx,Om_By,Om_Bz};
  double Om_Bz_orb = RotMat(W,I,w,Om_B_in,"InOrb")[2];
  
return 9.0/tf_B(tv_B,R_B,k_B,a,m_A,m_B) * (f_5(e) - (11.0/18.0)*Om_Bz_orb/n(m_A,m_B,a) * f_4(e));
}


double W_A(double tv_A, double R_A, double k_A, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay, double Om_Az){
  
  double Om_A_in[3] = {Om_Ax,Om_Ay,Om_Az};
  double Om_Az_orb = RotMat(W,I,w,Om_A_in,"InOrb")[2];

  return 1.0/tf_A(tv_A,R_A,k_A,a,m_A,m_B) * ( f_2(e) - Om_Az_orb/n(m_A,m_B,a) * f_1(e) );
}


double W_B(double tv_B, double R_B, double k_B, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By, double Om_Bz){
  
  double Om_B_in[3] = {Om_Bx,Om_By,Om_Bz};
  double Om_Bz_orb = RotMat(W,I,w,Om_B_in,"InOrb")[2];

  return 1.0/tf_B(tv_B,R_B,k_B,a,m_A,m_B) * ( f_2(e) - Om_Bz_orb/n(m_A,m_B,a) * f_1(e) );
}


double X_A(double tv_A, double R_A, double k_A, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay, double Om_Az){


  double Om_A_in[3] = {Om_Ax,Om_Ay,Om_Az};
  double Om_Ax_orb = RotMat(W,I,w,Om_A_in,"InOrb")[0];
  double Om_Ay_orb = RotMat(W,I,w,Om_A_in,"InOrb")[1];
  double Om_Az_orb = RotMat(W,I,w,Om_A_in,"InOrb")[2];
  double mu;
  mu = m_A*m_B/(m_A+m_B);
  

return -m_B*k_A*pow(R_A,5)/(mu*n(m_A,m_B,a)*pow(a,5)) * Om_Az_orb*Om_Ax_orb/pow(1.0-e*e,2)- Om_Ay_orb/(2.0*n(m_A,m_B,a)*tf_A(tv_A,R_A,k_A,a,m_A,m_B)) * f_3(e);
}

double X_B(double tv_B, double R_B, double k_B, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By,double Om_Bz){

  double Om_B_in[3] = {Om_Bx,Om_By,Om_Bz};
  double Om_Bx_orb = RotMat(W,I,w,Om_B_in,"InOrb")[0];
  double Om_By_orb = RotMat(W,I,w,Om_B_in,"InOrb")[1];
  double Om_Bz_orb = RotMat(W,I,w,Om_B_in,"InOrb")[2];
  double mu;
  mu = m_A*m_B/(m_A+m_B);
  
return -m_A*k_B*pow(R_B,5)/(mu*n(m_A,m_B,a)*pow(a,5)) * Om_Bz_orb*Om_Bx_orb/pow(1.0-e*e,2) - Om_By_orb/(2.0*n(m_A,m_B,a)*tf_B(tv_B,R_B,k_B,a,m_A,m_B)) * f_3(e);
}


double Y_A(double tv_A, double R_A, double k_A, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay,double Om_Az){

  double Om_A_in[3] = {Om_Ax,Om_Ay,Om_Az};
  double Om_Ax_orb = RotMat(W,I,w,Om_A_in,"InOrb")[0];
  double Om_Ay_orb = RotMat(W,I,w,Om_A_in,"InOrb")[1];
  double Om_Az_orb = RotMat(W,I,w,Om_A_in,"InOrb")[2];
  double mu;
  mu = m_A*m_B/(m_A+m_B);
 
return -m_B*k_A*pow(R_A,5)/(mu*n(m_A,m_B,a)*pow(a,5)) * Om_Az_orb*Om_Ay_orb/pow(1.0-e*e,2) + Om_Ax_orb/(2.0*n(m_A,m_B,a)*tf_A(tv_A,R_A,k_A,a,m_A,m_B)) * f_4(e);
  
}

double Y_B(double tv_B, double R_B, double k_B, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By,double Om_Bz){

  double Om_B_in[3] = {Om_Bx,Om_By,Om_Bz};
  double Om_Bx_orb = RotMat(W,I,w,Om_B_in,"InOrb")[0];
  double Om_By_orb = RotMat(W,I,w,Om_B_in,"InOrb")[1];
  double Om_Bz_orb = RotMat(W,I,w,Om_B_in,"InOrb")[2];
  double mu;
  mu = m_A*m_B/(m_A+m_B);
 
  return -m_A*k_B*pow(R_B,5)/(mu*n(m_A,m_B,a)*pow(a,5)) * Om_Bz_orb*Om_By_orb/pow(1.0-e*e,2) + Om_Bx_orb/(2.0*n(m_A,m_B,a)*tf_B(tv_B,R_B,k_B,a,m_A,m_B)) * f_4(e);
  
}

double Z_A(double R_A, double k_A, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay,double Om_Az){

  double Om_A_in[3] = {Om_Ax,Om_Ay,Om_Az};
  double Om_Ax_orb = RotMat(W,I,w,Om_A_in,"InOrb")[0];
  double Om_Ay_orb = RotMat(W,I,w,Om_A_in,"InOrb")[1];
  double Om_Az_orb = RotMat(W,I,w,Om_A_in,"InOrb")[2];
  double mu;
  mu = m_A*m_B/(m_A+m_B);
  
return m_B*k_A*pow(R_A,5)/(mu*n(m_A,m_B,a)*pow(a,5)) * ( (2.0*pow(Om_Az_orb,2) - pow(Om_Ay_orb,2) - pow(Om_Ax_orb,2))/(2.0*pow(1.0-e*e,2)) + 15.0*G*m_B/pow(a,3) * f_4(e) );
  
}


double Z_B(double R_B, double k_B, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By, double Om_Bz){

  double Om_B_in[3] = {Om_Bx,Om_By,Om_Bz};
  double Om_Bx_orb = RotMat(W,I,w,Om_B_in,"InOrb")[0];
  double Om_By_orb = RotMat(W,I,w,Om_B_in,"InOrb")[1];
  double Om_Bz_orb = RotMat(W,I,w,Om_B_in,"InOrb")[2];
  double mu;
  mu = m_A*m_B/(m_A+m_B);
  
  return m_A*k_B*pow(R_B,5)/(mu*n(m_A,m_B,a)*pow(a,5)) * ( (2.0*pow(Om_Bz_orb,2) - pow(Om_By_orb,2) - pow(Om_Bx_orb,2))/(2.0*pow(1.0-e*e,2)) + 15.0*G*m_A/pow(a,3) * f_4(e) );
  
}



