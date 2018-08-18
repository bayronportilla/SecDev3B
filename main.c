#include "allvars.h"
#include "proto.h"
#include <gsl/gsl_odeiv2.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <math.h>


int main (void){

  Inpar st;
  st = params();
  
  ////////////////////////////////////////////////////////////
  //
  // Charging data files for interpolation
  //
  ////////////////////////////////////////////////////////////

  /*
  FILE *file;
  file   = fopen("red_dwarf_radius.dat","r");
  Nlines = counterLines("red_dwarf_radius.dat");
  tim_A  = (double *)malloc(Nlines*sizeof(double));
  rad_A  = (double *)malloc(Nlines*sizeof(double));
  
  int j=0;
  double col1,col2,col3;
  while(fscanf(file,"%lf %lf %lf",&col1,&col2,&col3)!=EOF){
    tim_A[j] = col1*Gyr/st.uT;
    rad_A[j] = col3*RS/st.uL;
    j++;
  }

  acc    = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_linear, Nlines);
  gsl_spline_init(spline,tim_A,rad_A,Nlines);
  */
  /*
  FILE *files;
  files = fopen("interpolated_radius.dat","w");
  double time;
  for(time=1.0e6*YEARS/st.uT;time<5.0e9*YEARS/st.uT;time+=5.0e5*YEARS/st.uT){
    fprintf(files,"%1.4e %1.4e \n",time*st.uT/Gyr,fn_R_A(st,time,tim_A,rad_A)*st.uL/RS);
  }
  */
  
  //gsl_odeiv2_driver_alloc_y_new(const gsl_odeiv2_system * sys,
  //                              const gsl_odeiv2_step_type * T,
  //                              const double hstart,
  //                              const double epsabs,
  //                              const double epsrel)

  //printf("%f\n",tim_B[0]*st.uT/YEARS);

  ////////////////////////////////////////////////////////////
  // This is the only line must be modified by the user
  //
  const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
  //
  ////////////////////////////////////////////////////////////

  double hstart = 1.0e0*YEARS/st.uT;
  gsl_odeiv2_system sys = { func, NULL, 16, &st}; // Define sistema de ecuaciones
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T, hstart, 1e-10, 1e-10);

  double y[16] = { st.a_in, st.a_out, st.e_in, st.e_out, 
		   st.I_in, st.I_out, st.W_in, st.W_out,
		   st.w_in, st.w_out, st.Om_Ax, st.Om_Ay,
		   st.Om_Az, st.Om_Bx, st.Om_By, st.Om_Bz}; // y[number of entries of the array] = {}
  
  int s;
  double t = st.t_ini;
  double progress;
  
  FetchInfo(st);
  
   FILE *fp;
   char src[100];
   char dest[100];
   char name_files[100];
   strcpy(src,st.sim_name);
   strcpy(dest,".dat");
   strcpy(name_files,strcat(src,dest));
   fp  = fopen(name_files,"w");
   
   FILE *fp2;
   fp2  = fopen("final_values.dat","w");

   FILE *fp3;
   fp3  = fopen("derivada.dat","w");

   
   double ti = t;
   double a_stop,e_stop,d_roche,t_stop,a_min;
   a_stop = 0.005*st.a_in;
   a_min  = 0.001*AU/st.uL;
   e_stop = 0.001;
   t_stop = 5.0e6*YEARS/st.uT;
   d_roche = 1.66*st.R_A*pow((st.m_A+st.m_B)/st.m_B , 1/3.);

   char stop_reason[100];
   int twrite=0;




   while(t<=st.t_end){
     s = gsl_odeiv2_driver_apply(d, &t, ti, y); // from t to ti
     if (s != GSL_SUCCESS){
       printf ("error: driver returned %d\n", s);
       break;
     }
     
     ti += st.h_output;     
     progress = (t/st.t_end)*100.0;

     /*
     if ( (int)progress%10 == 0){
       printf("Progress: %d per cent \n",(int)progress);
     }
     */
     
     printf("t=%.4f Myr a_in=%.16e e_in=%.16e \n",t*st.uT/Myr,y[0]*st.uL/AU,y[2]);
     
     
     fprintf(fp,"%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e \
%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e \n",
	     t,y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10],y[11],y[12],y[13],y[14],y[15],
	     st.m_A,st.m_B,fn_R_A(st,t,tim_A,rad_A),st.R_B,st.k_A,st.k_B,st.tv_A,st.tv_B,st.gyr_rad_A,st.gyr_rad_B);
     
     
     fprintf(fp3,"%.5e %.5e \n",
	     t,de_in_dt(y[0],y[1],y[2],y[3],
			y[4],y[5],y[6],y[7],
			y[8],y[9],y[10],y[11],
			y[12],y[13],y[14],y[15],
			t,st));
     
     
     if(y[0]<a_min){
       strcpy(stop_reason,"a_min_or_a_stop_reached");
       break;
     }
     
   } // end while
      


   ////////////////////////////////////////////////////////////
   // write decay time vs a_ini,e_ini,I_ini.
   
   double Itot_ini=(st.I_in + st.I_out)*InDeg;
   
   if(strcmp(stop_reason,"")==0){
     fprintf(fp2,"%s %.5e %.5e %.5e %.5e %s \n",st.sim_name,t*st.uT/Myr,st.a_in*st.uL/AU,st.e_in,Itot_ini,"time_exhausted");
   }
   else{
     fprintf(fp2,"%s %.5e %.5e %.5e %.5e %s \n",st.sim_name,t*st.uT/Myr,st.a_in*st.uL/AU,st.e_in,Itot_ini,stop_reason);
   }  

   

   ////////////////////////////////////////////////////////////
   // Releasing memory
   
   gsl_odeiv2_driver_free (d);
   fclose(fp);
   fclose(fp2);
   //free(e_extreme);
   //free(t_extreme);
   gsl_spline_free (spline);
   gsl_interp_accel_free (acc);
  
   
   
  return 0;
   
}

