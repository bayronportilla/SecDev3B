#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_vector.h>
double *RotMat(double W, double I, double w, double *in_array, char dir[]){

  /*
    Esta funcion permite hacer la conversion de un arreglo unidimensional de 
    tres entradas medidas respecto a un sistema de referencia inicial a otro 
    arreglo de la misma dimension con entradas medidas en otro sistema que se
    ha construido mediante tres rotaciones sucesivas parametrizadas por los 
    angulos W,I,w. Tambien permite hacer la conversion inversa. En el contexto 
    de este programa, las posibles direcciones son: de un sistema inercial a uno 
    orbital "InOrb" y de uno orbital a uno inercial "OrbIn". Para usarlo, los
    argumentos a introducir (respetando el orden son) W, I, w, un puntero a un vector
    que contiene las entradas en el plano de partida (inercial u orbital) y una cadena
    de caracteres que indica la direccion de la conversion "InOrb" o "OrbIn". El return 
    de esta funcion es un vector unidimensional de tres entradas, medidas en el plano de 
    llegada. 

    Ej 1. 

    #include <ModQuad.h>
    ... y demas includes de GSL y math.h

    int main(void){
       double vector[3]={0.1,0.2,0.3};
       printf("%f n",RotMat(2,3,4,vector,"InOrb")[0]);
    }
    
   */

  ////////////////////////////////////////////////////////////
  //
  // vectors
  //
  ////////////////////////////////////////////////////////////

  static double out_array[3]; // the output array with the elements in the orbital plane
  
  gsl_vector *v = gsl_vector_alloc (3);

  gsl_vector_set (v, 0, in_array[0]);
  gsl_vector_set (v, 1, in_array[1]);
  gsl_vector_set (v, 2, in_array[2]);
  
  ////////////////////////////////////////////////////////////
  //
  // matrices
  //
  ////////////////////////////////////////////////////////////
  
  gsl_matrix *m = gsl_matrix_alloc (3, 3);
  
  gsl_matrix_set (m, 0, 0,  cos(W)*cos(w) - sin(W)*sin(w)*cos(I));
  gsl_matrix_set (m, 0, 1,  sin(W)*cos(w) + cos(W)*sin(w)*cos(I));
  gsl_matrix_set (m, 0, 2,  sin(I)*sin(w));
  gsl_matrix_set (m, 1, 0, -cos(W)*sin(w) - sin(W)*cos(w)*cos(I));
  gsl_matrix_set (m, 1, 1, -sin(W)*sin(w) + cos(W)*cos(w)*cos(I));
  gsl_matrix_set (m, 1, 2,  sin(I)*cos(w));
  gsl_matrix_set (m, 2, 0,  sin(I)*sin(W));
  gsl_matrix_set (m, 2, 1, -sin(I)*cos(W));
  gsl_matrix_set (m, 2, 2,  cos(I));
 
  ////////////////////////////////////////////////////////////
  //
  // multiplication
  //
  ////////////////////////////////////////////////////////////
  
  gsl_vector *y = gsl_vector_alloc (3); // Result will be stored in y

  char str1[] = "InOrb";
  char str2[] = "OrbIn";
  
  if( strcmp(dir,str1) == 0){
      gsl_blas_dgemv(CblasNoTrans, 1.0, m, v, 0.0, y);
    }

  if( strcmp(dir,str2) == 0){
      gsl_blas_dgemv(CblasTrans, 1.0, m, v, 0.0, y);
    }
  
  out_array[0] = gsl_vector_get (y, 0);
  out_array[1] = gsl_vector_get (y, 1);
  out_array[2] = gsl_vector_get (y, 2);

  gsl_vector_free (v);
  gsl_vector_free (y);
  gsl_matrix_free (m);
 
  return out_array;
  
}



//double *RotMatt(double W, double I, double w, double *in_array, char dir[]){

  /*
    Esta funcion permite hacer la conversion de un arreglo unidimensional de 
    tres entradas medidas respecto a un sistema de referencia inicial a otro 
    arreglo de la misma dimension con entradas medidas en otro sistema que se
    ha construido mediante tres rotaciones sucesivas parametrizadas por los 
    angulos W,I,w. Tambien permite hacer la conversion inversa. En el contexto 
    de este programa, las posibles direcciones son: de un sistema inercial a uno 
    orbital "InOrb" y de uno orbital a uno inercial "OrbIn". Para usarlo, los
    argumentos a introducir (respetando el orden son) W, I, w, un puntero a un vector
    que contiene las entradas en el plano de partida (inercial u orbital) y una cadena
    de caracteres que indica la direccion de la conversion "InOrb" o "OrbIn". El return 
    de esta funcion es un vector unidimensional de tres entradas, medidas en el plano de 
    llegada. 

    Ej 1. 

    #include <ModQuad.h>
    ... y demas includes de GSL y math.h

    int main(void){
       double vector[3]={0.1,0.2,0.3};
       printf("%f n",RotMat(2,3,4,vector,"InOrb")[0]);
    }
    
   */
  /*
  ////////////////////////////////////////////////////////////
  //
  // vectors
  //
  ////////////////////////////////////////////////////////////

  static double out_array[3]; // the output array with the elements in the orbital plane
  //out_array=(double*)malloc(3*sizeof(double));
  
  gsl_vector *v = gsl_vector_alloc (3);

  
  gsl_vector_set (v, 0, in_array[0]);
  gsl_vector_set (v, 1, in_array[1]);
  gsl_vector_set (v, 2, in_array[2]);
  


  ////////////////////////////////////////////////////////////
  //
  // matrices
  //
  ////////////////////////////////////////////////////////////
  
  gsl_matrix *m = gsl_matrix_alloc (3, 3);
  
  gsl_matrix_set (m, 0, 0,  cos(W)*cos(w) - sin(W)*sin(w)*cos(I));
  gsl_matrix_set (m, 0, 1,  sin(W)*cos(w) + cos(W)*sin(w)*cos(I));
  gsl_matrix_set (m, 0, 2,  sin(I)*sin(w));
  gsl_matrix_set (m, 1, 0, -cos(W)*sin(w) - sin(W)*cos(w)*cos(I));
  gsl_matrix_set (m, 1, 1, -sin(W)*sin(w) + cos(W)*cos(w)*cos(I));
  gsl_matrix_set (m, 1, 2,  sin(I)*cos(w));
  gsl_matrix_set (m, 2, 0,  sin(I)*sin(W));
  gsl_matrix_set (m, 2, 1, -sin(I)*cos(W));
  gsl_matrix_set (m, 2, 2,  cos(I));
 
  ////////////////////////////////////////////////////////////
  //
  // multiplication
  //
  ////////////////////////////////////////////////////////////
  
  gsl_vector *y = gsl_vector_alloc (3); // Result will be stored in y

  char str1[] = "InOrb";
  char str2[] = "OrbIn";
  
  if( strcmp(dir,str1) == 0){
      gsl_blas_dgemv(CblasNoTrans, 1.0, m, v, 0.0, y);
    }

  if( strcmp(dir,str2) == 0){
      gsl_blas_dgemv(CblasTrans, 1.0, m, v, 0.0, y);
    }
  
  out_array[0] = gsl_vector_get (y, 0);
  out_array[1] = gsl_vector_get (y, 1);
  out_array[2] = gsl_vector_get (y, 2);

  gsl_vector_free (v);
  gsl_vector_free (y);
  gsl_matrix_free (m);
 
  return out_array;
  
}
*/
