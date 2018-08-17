//#include "allvars.h"

int sign(double x);

int FetchInfo(Inpar st);

double da_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params);


double da_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params);


double de_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params);



double de_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params);


double dI_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params);


double dI_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params);


double dW_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params);


// Ojo aqui
double dW_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params);


double dw_in_dt(double a_in, double a_out, double e_in,
		double e_out, double I_in, double I_out,
		double W_in, double W_out, double w_in,
		double w_out, double Om_Ax, double Om_Ay,
		double Om_Az, double Om_Bx, double Om_By,
		double Om_Bz, double t, Inpar params);


double dw_out_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params);


double dOm_Ax_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params);


double dOm_Ay_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params);


double dOm_Az_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params);


double dOm_Bx_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params);


double dOm_By_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params);


double dOm_Bz_dt(double a_in, double a_out, double e_in,
		 double e_out, double I_in, double I_out,
		 double W_in, double W_out, double w_in,
		 double w_out, double Om_Ax, double Om_Ay,
		 double Om_Az, double Om_Bx, double Om_By,
		 double Om_Bz, double t, Inpar params);

int func (double t, const double y[], double f[], void *params);

double *units(double a, double b, char var[]);
double fn_R_A(Inpar st, double t, double x[], double y[]);
double fn_R_B(Inpar st, double t, double x[], double y[]);
double interpol(Inpar st, double t, double x[], double y[]);
Inpar params();

double f_1(double e);

double f_2(double e);

double f_3(double e);

double f_4(double e);

double f_5(double e);

double n(double m_1, double m_2, double a);

double L_1(double m_1, double m_2, double a_1);

double L_2(double m_1,double m_2,double m_3, double a_2);

double G_1(double m_1,double m_2,double a_1, double e_1);

double G_2(double m_1,double m_2,double m_3, double a_2, double e_2);

double H_1(double m_1,double m_2,double a_1,double e_1,double I_1);

double H_2(double m_1,double m_2,double m_3,double a_2,double e_2,double I_2 );

double C2(double m_1,double m_2,double m_3,double a_1,double a_2,double e_2);

double epsilon_m(double m_1,double m_2,double a_1,double a_2,double e_2);

// Ojo, si es negativo? ver eq 24 de Naoz 2013
double C3(double m_1,double m_2,double m_3,double a_1,double a_2,double e_2);

double tf_A(double tv_A, double R_A, double k_A, double a, double m_A, double m_B);

double tf_B(double tv_B, double R_B, double k_B, double a, double m_A, double m_B);

double V_A(double tv_A, double R_A, double k_A, double m_A, double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay, double Om_Az);

double V_B(double tv_B, double R_B, double k_B, double m_A, double m_B, 
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By, double Om_Bz);

double W_A(double tv_A, double R_A, double k_A, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay, double Om_Az);

double W_B(double tv_B, double R_B, double k_B, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By, double Om_Bz);

double X_A(double tv_A, double R_A, double k_A, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay, double Om_Az);

double X_B(double tv_B, double R_B, double k_B, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By,double Om_Bz);

double Y_A(double tv_A, double R_A, double k_A, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay,double Om_Az);

double Y_B(double tv_B, double R_B, double k_B, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By,double Om_Bz);

double Z_A(double R_A, double k_A, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Ax, double Om_Ay,double Om_Az);

double Z_B(double R_B, double k_B, double m_A,double m_B,
	   double a, double e, double W, double I, double w,
	   double Om_Bx, double Om_By, double Om_Bz);

double *RotMat(double W, double I, double w, double *in_array, char dir[]);

int counterLines(char *infile);
