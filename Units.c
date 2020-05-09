#include "allvars.h"
#include "proto.h"
#include <math.h>
#include <string.h>

/*
Este modulo calcula el sistema de unidades canonicas
apropiado para un problema determinado. Para usarlo,
identifique primero dos unidades fundamentales que
serviran para deducir la tercera.  Supongamos que 
fijamos las unidades a y b. Para determinar la unidad c
es necesario invocar a la funcion *units* siguiendo la 
sigueinte convencion units(a,b,"c") donde c es una cadena
que puede adoptar los siguientes valores "uM", "uL" y "uT". 
a y b son variables tipo double que se deben introducir 
de manera ciclica asi "uM" ---> "uL" ---> "uT". 

  Ej 1. Si uM = 1.980e30 y uL = 149.6e9. La unidad a derivar
  es uT. La funcion se invoca units como units(1.989e30,149.6e9,"uT")
  Note que los parametros introducidos en la funcion respetan 
  el orden ciclico.  

  Ej 2. Si uM = 5.94e24 y uT = 86400. La unidad a derivar es 
  uL. La funcion se invoca como units(5.94e24,86400,"uL"). Una 
  vez mas, el orden en que se han introducido los parametros 
  NO es arbitrario, debe respetar el orden ciclico. 

El retorno de la funcion es un puntero a un arreglo de tres componentes
(uM,uL,uT), en ese orden estricto.
*/

 
double *units(double a, double b, char var[]){

  char   var1[] = "uM";
  char   var2[] = "uL";
  char   var3[] = "uT";
  double G_univ      = 6.6740831e-11;
  double uM,uL,uT;
   
  static double array[3];
   
  if ( strcmp(var,var1) == 0 ) {  
    // I will find the unit of mass
    uL = a;
    uT = b;
    uM = uL*uL*uL/(uT*uT*G_univ);
  }
  
  if ( strcmp(var,var2) == 0 ) {  
    // I will find the unit of longitude
    uT = a;
    uM = b;
    uL = pow(G_univ*uM*uT*uT,1.0/3.0);
  }
  
  if ( strcmp(var,var3) == 0) {  
    // I will find the unit of time
    uM = a;
    uL = b;
    uT = pow(uL*uL*uL/(G_univ*uM),0.5);
  }

  array[0] = uM;
  array[1] = uL;
  array[2] = uT;
    
  return array;
}

 
/*
double *units(double a, double b, char var[]){

  char   var1[] = "uM";
  char   var2[] = "uL";
  char   var3[] = "uT";
  double B3_univ      = 351.6e-120;
  double uM,uL,uT;
   
  static double array[3];
   
  // I will find the unit of time
  uM = a;
  uL = b;
  uT = pow(uL*uL*uL*uM*uM*uM*uM/(B3_univ),0.5);


  array[0] = uM;
  array[1] = uL;
  array[2] = uT;
    
  return array;
}
*/
