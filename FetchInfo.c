#include <stdio.h>
#include <string.h>
#include "allvars.h"


int FetchInfo(Inpar st){
  FILE *fp;
  char src[100];
  char dest[100];
  char name_files[100];
  strcpy(src,st.sim_name);
  strcpy(dest,".log");
  strcpy(name_files,strcat(src,dest));

  fp = fopen(name_files,"w");

  fprintf(fp,"0   m_A          %1.17e \n\
1   m_B          %1.17e \n\
2   m_C          %1.17e \n\
3   R_A          %1.17e \n\
4   R_B          %1.17e \n\
5   k_A          %1.17e \n\
6   k_B          %1.17e \n\
7   tv_A         %1.17e \n\
8   tv_A         %1.17e \n\
9   gyr_rad_A    %1.17e \n\
10  gyr_rad_B    %1.17e \n\
11  a_in         %1.17e \n\
12  a_out        %1.17e \n\
13  e_in         %1.17e \n\
14  e_out        %1.17e \n\
15  I_in         %1.17e \n\
16  I_out        %1.17e \n\
17  W_in         %1.17e \n\
18  W_out        %1.17e \n\
19  w_in         %1.17e \n\
20  w_out        %1.17e \n\
21  P_rot_A      %1.17e \n\
22  P_rot_B      %1.17e \n\
23  alpha_A      %1.17e \n\
24  alpha_B      %1.17e \n\
25  beta_A       %1.17e \n\
26  beta_B       %1.17e \n\
27  t_ini        %1.17e \n\
28  t_end        %1.17e \n\
29  q_orb        %d     \n\
30  q_tid        %d     \n\
31  q_GR         %d     \n\
32  uM           %1.17e \n\
33  uL           %1.17e \n\
34  uT           %1.17e \n\
35  h_output     %1.17e \n\
36  sim_name     %s     \n",
	  st.m_A,st.m_B,st.m_C,st.R_A,st.R_B,st.k_A,st.k_B,st.tv_A,st.tv_B,
	  st.gyr_rad_A,st.gyr_rad_B,st.a_in,st.a_out,st.e_in,st.e_out,
	  st.I_in,st.I_out,st.W_in,st.W_out,st.w_in,st.w_out,st.P_rot_A,
	  st.P_rot_B,st.alpha_A,st.alpha_B,st.beta_A,st.beta_B,st.t_ini,
	  st.t_end,st.q_orb,st.q_tid,st.q_GR,st.uM,st.uL,st.uT,st.h_output,st.sim_name);


  return 0;
}
