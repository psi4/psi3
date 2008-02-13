/*! \file
    \ingroup OPTKING
    \brief miscellaneous little matrix and print functions
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <psifiles.h>

#define EXTERN
#include "opt.h"
#undef EXTERN

namespace psi { namespace optking {

/*** PRINT_MAT2   prints a matrix to output file ***/
void print_mat2(double **matrix, int rows, int cols, FILE *of) {
  int i,j,col;
  for (i=0;i<rows;++i) {
    col = 0;
    for (j=0;j<cols;++j) {
      if (col == 9) {
        fprintf(outfile,"\n");
        col = 0;
      }
      fprintf(of,"%15.10lf",matrix[i][j]);
      ++col;
    }
    fprintf(outfile,"\n");
  } 
  return;
} 

/*** PRINT_MAT5   prints a matrix to output file ***/
void print_mat5(double **matrix, int rows, int cols, FILE *of) {
  int i,j,col;
  for (i=0;i<rows;++i) {
    col = 0;
    for (j=0;j<cols;++j) {
      if (col == 14) {
        fprintf(outfile,"\n");
        col = 0;
      }
      fprintf(of,"%9.5lf",matrix[i][j]);
      ++col;
    }
    fprintf(outfile,"\n");
  } 
  return;
} 

/*** CROSS_PRODUCT   computes cross product of two vectors ***/
void cross_product(double *u,double *v,double *out)
{
  out[0] = u[1]*v[2]-u[2]*v[1];
  out[1] = -1.0*(u[0]*v[2]-u[2]*v[0]);
  out[2] = u[0]*v[1]-u[1]*v[0];
  return;
} 

/*** SCALAR_MULT   performs scalar multiplication of a vector***/ 
void scalar_mult(double a, double *vect, int dim) {
  int i;
  for (i=0;i<dim;++i)
    vect[i] *= a;
  return;
}

/*** SCALAR_DIV   performs scalar division of a vector ***/
void scalar_div(double a, double *vect) {
  int i;
  for (i=0;i<3;++i)
    vect[i] /= a;
  return;
}


void punt(char *message) {
  fprintf(outfile,"\nerror: %s\n", message);
  fprintf(outfile,"         *** stopping execution ***\n");
  fprintf(stderr,"\n OPTKING error: %s\n", message);
  fprintf(stderr,"                 *** stopping execution ***\n");
  /* fclose(outfile); */
  /* exit(1); */
  exit_io();
}

void open_PSIF(void) {
  psio_open(PSIF_OPTKING, PSIO_OPEN_OLD);
  return;
}

void close_PSIF(void) {
  psio_close(PSIF_OPTKING, 1);
  return;
}

void exit_io(void) {
  fprintf(outfile,"\n******** OPTKING execution completed ********\n\n");
  psio_done();
  psi_stop(infile,outfile,psi_file_prefix);
}

/*** SWAP_TORS -- canonical torsion order is atom a < atom d ***/
void swap_tors(int *a, int *b, int *c, int *d) {
  int p;
  if (*a > *d) {
    p = *a;
    *a = *d;
    *d = p;
    p = *b;
    *b = *c;
    *c = p;
  }
}

/*** SWAP -- canonical order of bonds is atom a < atom b ***/
void swap(int *a, int *b) {
  int c;
  if (*a > *b) {
    c = *a;
    *a = *b;
    *b = c;
  }
  return;
}

/*** DIV_INT Does little go into big? ***/
int div_int(int big, int little) {
  if (little > big) return 0;
  for(;;) {
    big -= little;
    if (big == 0) return 1;
    if (big < 0) return 0;
  }
}

/*** SYM_MATRIX_INVERT inverts a matrix by diagonalization
 *  **A = matrix to be inverted
 *  dim = dimension of A
 *  print_det = print determinant if 1, nothing if 0
 *  redundant = zero eigenvalues allowed if 1
 *  returns: inverse of A ***/
double **symm_matrix_invert(double **A, int dim, int print_det, int redundant) {
  int i;
  double **A_inv, **A_vects, *A_vals, **A_temp, det=1.0;

  A_inv   = block_matrix(dim,dim);
  A_temp  = block_matrix(dim,dim);
  A_vects = block_matrix(dim,dim);
  A_vals  = init_array(dim);

  sq_rsp(dim,dim,A,A_vals,1,A_vects,EVAL_TOL);

  if (redundant == 0) {
    for (i=0;i<dim;++i) {
      det *= A_vals[i];
      A_inv[i][i] = 1.0/A_vals[i];
    }
    if (print_det)
      fprintf(outfile,"Determinant: %10.6e\n",det);
    if (fabs(det) < 1E-10)
      fprintf(outfile,"Warning, determinant is small: %10.6e\n",det);
  }
  else {
    for (i=0;i<dim;++i) {
      det *= A_vals[i];
      if (fabs(A_vals[i]) > REDUNDANT_EVAL_TOL)
        A_inv[i][i] = 1.0/A_vals[i];
      else
        A_inv[i][i] = 0.0;
    }
    if (print_det)
      fprintf(outfile,"Determinant: %10.6e\n",det);
  }

  mmult(A_inv,0,A_vects,1,A_temp,0,dim,dim,dim,0);
  mmult(A_vects,0,A_temp,0,A_inv,0,dim,dim,dim,0);

  free(A_vals);
  free_block(A_vects);
  free_block(A_temp);
  return A_inv;
} 

double energy_chkpt(void) {
  double energy;

  chkpt_init(PSIO_OPEN_OLD);
  /* energy = chkpt_rd_escf(); */
  energy = chkpt_rd_etot();
  chkpt_close();
  return energy;
}
     

/* This routine transposes matrices and calls lapack dgeev() in libqt *
 * to diagonalize a square nonsymmetric matrix.  The eigenvalues      *
 * are returned in random order.                                      */
void dgeev_optking(int L, double **G, double *lambda, double **alpha) {

  int i, j, lwork, info;
  double *evals_i, *work, **left_evects, tval, temp;

  evals_i = init_array(L); 
  left_evects = block_matrix(L,L);

  work = init_array(20*L);
  lwork = 20*L;          

  for (i=0; i<L; ++i)
    for (j=0; j<i; ++j) {
      temp = G[i][j];
      G[i][j] = G[j][i];
      G[j][i] = temp;
    }

  i = C_DGEEV(L, G, L, lambda, evals_i, left_evects,
    L, alpha, L, work, lwork, info);

  for (i=0; i<L; ++i)
    for (j=0; j<i; ++j) {
      temp = alpha[i][j];
      alpha[i][j] = alpha[j][i];
      alpha[j][i] = temp;
    }

  free(work);
  free(evals_i);
  free_block(left_evects);
  return;
}

double **mass_mat(double *masses) {
    int i, dim;
    double **u;

    dim = 3*optinfo.natom;
    u = block_matrix(dim,dim);

    for (i=0; i<dim; ++i) {
      if (masses[i] == 0.0)
        u[i][i] = 0.0;
      else
        u[i][i] = 1.0/masses[i];
    }
    return u;
}

double **unit_mat(int dim) {
    int i;
    double **u;
    u = block_matrix(dim,dim);

    for (i=0; i<dim; ++i)
        u[i][i] = 1.0;

    return u;
}

void print_evects(double **evects, double *evals, int nrow, int ncol, FILE *out)
{
      int ii,jj,kk,nn;
      int i,j;

      ii=0;jj=0;
L200:
      ii++;
      jj++;
      kk=10*jj;
      nn=ncol;
      if (nn > kk) nn=kk;
      fprintf (out,"\n");
      for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
      fprintf (out,"\n");
      for (i=0; i < nrow; i++) {
         fprintf (out,"\n%5d",i+1);
         for (j=ii-1; j < nn; j++) {
            fprintf (out,"%12.7f",evects[i][j]);
            }
         }
      fprintf (out,"\n");
      fprintf (out,"\n     ");
      for (j=ii-1; j < nn; j++) {
         fprintf(out,"%12.7f",evals[j]);
         }
      fprintf (out,"\n");
      if (ncol <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
}

}} /* namespace psi::optking */

