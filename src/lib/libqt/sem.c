/*
** Simultaneous Expansion Method for the Iterative Solution of
** Several of the Lowest Eigenvalues and Corresponding Eivenvectors of
** Large Real-Symmetric Matrices
**
** Algorithm due to Bowen Liu
** IBM Research Laboratory
**
** Implemented for Schaefer Group by David Sherrill
** Center for Computational Quantum Chemistry, UGA
**
** In-core version for now!
** February 1994
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr.h>
#include "qt.h"

#define MAX_B_ROWS 200



/*** do a little test routine
main()
{
double **A ;
double *evals, **evecs ;
int i, j, used;
void sem() ;
FILE *outfile ;

   ffile(&outfile, "output.dat", 0) ;
   tstart(outfile) ;

   A = init_matrix(50,50) ;
   evals = init_array(4) ;
   evecs = init_matrix(4, 50) ;
   for (i=0; i<50; i++) {
      for (j=0; j<=i; j++) {
         if (i!=j) {A[i][j] = 1.0; A[j][i] = 1.0; }
         else if (i<5) A[i][j] = 1.0 + 0.1 * (double) i ;
         else A[i][j] = 2.0 * (double) i + 1.0 ;
         }
      }
   sem(A, 50, 4, evecs, evals, 1.0E-10, 6, &used) ;


   fprintf(outfile, "Ok, the eigenvectors are sideways!\n") ;
   eivout(evecs, evals, 4, 50, outfile) ;
   fprintf(outfile, "\nused %d expansion vectors\n", used) ;

   tstop(outfile) ;
   fclose(outfile) ;
}
***/


void sem(A, N, M, evecs, evals, conv_tol, maxiter, vu, offset, outfile)
      double **A ;
      int N, M ;
      double **evecs, *evals, conv_tol;
      int maxiter, *vu ;
      double offset ;
      FILE *outfile ;
{
double *tmp_vec, **tmp_mat ;
int sm_tridim ;
double *sm_mat ;
double *sm_evals ;
double **sm_evecs ;
int i, j, ij, k, I, L; 
double **G, **d ;
double *lambda, **alpha, **b, **f ;
double tval ;
int converged = 0, iter = 1 ;
double lastroot = 0.0 ;

   printf("offset = %f\n", offset) ;
   
/* check parameters */
   if (evecs == NULL || evals == NULL) {
      printf("(sem): passed uncallocated pointers for evecs or evals\n") ;
      return ;
      }

/* make space for temp vector */
   tmp_vec = init_array(N) ;

/* obtain a set of L orthonormal trial vectors, L > M */
   L = 5 + M ;   /* just a guess */
   sm_tridim = L * (L + 1) / 2 ;
   sm_mat = init_array(sm_tridim) ;
   sm_evals = init_array(L) ;
   sm_evecs = init_matrix(L, L) ;
   for (i=0, ij=0; i<L; i++) 
      for (j=0; j<=i; j++, ij++)
         sm_mat[ij] = A[i][j] ;
   rsp(L, L, sm_tridim, sm_mat, sm_evals, 1, sm_evecs, 1E-14) ;
   
/* need to fill out sm_evecs into b (pad w/ 0's) */
   b = (double **) malloc (MAX_B_ROWS * sizeof(double *)) ;
   for (i=0; i<MAX_B_ROWS; i++) {
      if (i<L) b[i] = init_array(N) ;
      else b[i] = NULL ;
      }

   for (i=0; i<L; i++) 
      for (j=0; j<L; j++)
         b[i][j] = sm_evecs[i][j] ;

/* allocate other arrays with ~fixed dimensions during iteration */
   d = init_matrix(M, N) ;    /* M and N are both fixed */
   f = init_matrix(M, N) ;
 

/* ITERATE */
   while (!converged && iter <= maxiter) {

   /* form G matrix */
      G = init_matrix(L, L) ;
      tmp_mat = init_matrix(L, N) ;
      mmult(b, 0, A, 0, tmp_mat, 0, L, N, N, 0) ; /* tmp = B * A    */
      mmult(tmp_mat, 0, b, 1, G, 0, L, N, L, 0) ; /* G = tmp * B(T) */
      free_matrix(tmp_mat, L) ;

   /* solve the L x L eigenvalue problem G a = lambda a for M roots */
      lambda = init_array(L) ;
      alpha = init_matrix(L, L) ;
      sq_rsp(L, L, G, lambda, 1, alpha, 1E-14) ;
      free_matrix(G, L) ;


   /* form the d part of the correction vector */
      zero_mat(d, M, N) ;
      for (k=0; k<M; k++) {
         for (i=0; i<L; i++) {
            mmult(A,0,&(b[i]),1,&(tmp_vec),1,N,N,1,0) ; /* tmp=A*b[i] */
            for (I=0; I<N; I++) {
               d[k][I] += alpha[i][k] * (tmp_vec[I] - lambda[k] * b[i][I]) ;
               }
            }
         }

   /* check for convergence */
      dot_arr(d[M-1], d[M-1], N, &tval) ;
      tval = sqrt(tval) ;
      fprintf(outfile, "Iter %3d  Root 1 = %12.7f  RMS Corr = %12.7f\n", 
         iter, (lambda[0]+offset), tval) ;
      if (fabs(lambda[0]-lastroot) <= conv_tol) {
         converged = 1 ;
         for (i=0; i<M; i++) {
            evals[i] = lambda[i] ;
            for (j=0; j<L; j++) {
               tval = alpha[j][i] ;
               for (I=0; I<N; I++) 
                  evecs[i][I] += tval * b[j][I] ;
               }
            } 
         free_matrix(alpha, L) ;
         free(lambda) ;
         break ;
         }
      else lastroot = lambda[0] ;

   /* form the correction vector and normalize */
      for (k=0; k<M; k++) {
         for (I=0; I<N; I++) {
            f[k][I] = d[k][I] / (lambda[k] - A[I][I]) ;
            }
         }
      normalize(f, M, N) ;
      free_matrix(alpha, L) ;
      free(lambda) ;

   /* Schmidt orthog and append f's to b */
      for (i=0; i<M; i++) 
         if (schmidt_add(b, L, N, f[i])) L++;

   /* Again Schmidt orthog b's (minimize numerical error) */
      schmidt(b, L, N, outfile);
      iter++ ;
      }

   *vu = L ;
 
   /* printf("(sem): Used %d basis vectors\n", L) ; */
   free(tmp_vec) ;
   free_matrix(d, M) ;
   free_matrix(f, M) ;
   free_matrix(b, L) ;
}
         
