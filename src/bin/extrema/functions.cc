/*##################################################################
#
#  small_functions.cc
#
#  various small functions for extrema
#                                                  J.Kenny 7-22-00
###################################################################*/						  

#define EXTERN
#include "extrema.h"

/*this needs to be in C*/
extern "C" {
char *gprgid()
{
   char *prgid = "EXTREMA";
 
   return(prgid);
   }
}                      


void punt(char *mess)
{
  fprintf(outfile, "  error: %s\n", mess);
  fprintf(stderr, "  EXTREMA error: %s\n", mess);
  // stop_io();
  exit(1);
}



double dot_pdt(double *vec1, double *vec2, int num) {

  int i;
  double result=0;

  for(i=0;i<num;++i) {
      result += vec1[i] * vec2[i];
    }

  return result;
}



/*----------------------------------------------------------------------------

       **SYM_MATRIX_INVERT

       inverts a matrix by diagonalization

       parameters:
             **A = matrix to be inverted
             dim = dimension of A
             print_det = print determinant if 1, nothing if 0
             redundant = zero eigenvalues allowed if 1

       returns:
             **inv_A = inverse of A

                                                      written by Rollin King
----------------------------------------------------------------------------*/

double **symm_matrix_invert(double **A, int dim, int print_det, int redundant) {
  int i;
  double **A_inv, **A_vects, *A_vals, **A_temp, det=1.0;

  A_inv   = init_matrix(dim,dim);
  A_temp  = init_matrix(dim,dim);
  A_vects = init_matrix(dim,dim);
  A_vals  = init_array(dim);

  sq_rsp(dim,dim,A,A_vals,1,A_vects,1E-10);

  if (redundant == 0) {
     for (i=0;i<dim;++i) {
        det *= A_vals[i];
        A_inv[i][i] = 1.0/A_vals[i];
     }
     if (print_det)
        fprintf(outfile,"Determinant: %10.6e\n",det);
     if (fabs(det) < 1E-10) {
        fprintf(outfile,"Determinant: %10.6e\n",det);
        fprintf(outfile,"Determinant is too small...aborting.\n");
        fclose(outfile);
        exit(2);
     }
  }
  else {
     for (i=0;i<dim;++i) {
        det *= A_vals[i];
        if (fabs(A_vals[i]) > 1E-10)
           A_inv[i][i] = 1.0/A_vals[i];
     }
     if (print_det)
        fprintf(outfile,"Determinant: %10.6e\n",det);
  }

  mmult(A_inv,0,A_vects,1,A_temp,0,dim,dim,dim,0);
  mmult(A_vects,0,A_temp,0,A_inv,0,dim,dim,dim,0);

  free(A_vals);
  free_matrix(A_vects,dim);
  free_matrix(A_temp,dim);
  return A_inv;
}

