#include "allvars.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spline.h>
#include <string.h>
#include <stdio.h>
#include "proto.h"

double interpol(Inpar st, double t, double x[], double y[]){

  ////////////////////////////////////////////////////////////
  // Interpolation for Body A
  
  //acc    = gsl_interp_accel_alloc ();
  //spline = gsl_spline_alloc (gsl_interp_linear, Nlines);
  
  //gsl_spline_init(spline,x,y,Nlines);

  
  return gsl_spline_eval (spline, t, acc);
    


  /*  
  ////////////////////////////////////////////////////////////
  //
  // Finding max and min of time entries
  //
  ////////////////////////////////////////////////////////////

  int i;
  gsl_vector *time_A_gsl = gsl_vector_alloc (Nlines);
  gsl_vector *time_B_gsl = gsl_vector_alloc (Nlines);
  gsl_vector *time_compa = gsl_vector_alloc (2);
  
  for (i = 0; i < Nlines; i++){
      gsl_vector_set (time_A_gsl, i, tim_A[i]);
      gsl_vector_set (time_B_gsl, i, tim_B[i]);
  }

  double min_time_A, min_time_B;
  double max_time_A, max_time_B;

  min_time_A = gsl_vector_min(time_A_gsl);
  min_time_B = gsl_vector_min(time_B_gsl);
  max_time_A = gsl_vector_max(time_A_gsl);
  max_time_B = gsl_vector_max(time_B_gsl);
  
  if (min_time_A<min_time_B){
    gsl_vector_set (time_compa,0,min_time_A);
  }
  else{
    gsl_vector_set (time_compa,0,min_time_B);
  }

  if (max_time_A>max_time_B){
    gsl_vector_set (time_compa,1,max_time_A);
  }
  elzse{
    gsl_vector_set (time_compa,1,max_time_B);
  }


  
  ////////////////////////////////////////////////////////////
  //
  // Printing info
  //
  ////////////////////////////////////////////////////////////

  //printf("t_ini can't be lower than:   %1.3e yr\n",gsl_vector_get(time_compa,0)*1e9);
  //printf("t_end can't be greater than: %1.3e yr\n",gsl_vector_get(time_compa,1)*1e9);
  

  


  ////////////////////////////////////////////////////////////
  //
  // Releasing memory
  //
  ////////////////////////////////////////////////////////////
  
  gsl_spline_free (spline_rad_A);
  gsl_spline_free (spline_rad_B);
  gsl_interp_accel_free (acc_rad_A);
  gsl_interp_accel_free (acc_rad_B);
  free(tim_A);
  free(tim_B);
  free(rad_A);
  free(rad_B);  
  */
  
  //return 0;
  
}
