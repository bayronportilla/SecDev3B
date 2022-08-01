#include "allvars.h"
#include "proto.h"

////////////////////////////////////////////////////////////
//
// Definition of functions for bulk changes
//
////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
//
// Body A
//
////////////////////////////////////////////////////////////

/*
// fictitious
double fn_R_A(Inpar st, double t, double x[], double y[]){

  double m,R_ini,R_end,t_ini,t_end,R;

  R_end = 1.1*RS/st.uL;
  R_ini = 2.62*RS/st.uL;
  t_end = 19.0e6*YEARS/st.uT;
  t_ini = st.t_ini*YEARS/st.uT;
  m = (R_end - R_ini)/(t_end - t_ini);
  R = m*t+R_ini;

  if(t<t_end){
    return R;
  }
  else{
    return R_end;
  }
}
*/


// realistic
double fn_R_A(Inpar st, double t, double x[], double y[]){
  //return interpol(st,t,x,y);
  return st.R_A;
}


double fn_R_B(Inpar st, double t, double x[], double y[]){
  //return interpol(st,t,x,y);
  return st.R_B;
}
