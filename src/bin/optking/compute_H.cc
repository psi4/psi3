// This function reads in hessian (H) from fconst.dat, inverts
// it to H_inv, and then BFGS updates it.  It writes the new
// Hessian to fconst.dat and returns H_inv

#include <math.h>

extern "C" {
  #include <stdio.h>
  #include <file30.h>
  #include <stdlib.h>
  #include <string.h>
  #include <libciomr.h>
  #include <physconst.h>
}

#define EXTERN
#include "opt.h"
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

double **compute_H(salc_set &symm, double *q, double *f_intcos, double **P) {
  double **F, **H, **H_inv, **H_inv_new, **H_new, **temp_mat, **temp_mat2;
  double a, b;
  int i,j,col,count,dim,update;
  char buffer[MAX_LINELENGTH];

  dim = symm.get_num();
  F = init_matrix(dim,dim);

  ffile(&fp_fconst,"fconst.dat",2);
  for (i=0;i<dim;++i) {
    fgets(buffer,MAX_LINELENGTH,fp_fconst);
    count = 0;
    for (j=0;j<=i;++j) {
      if ( div_int(j,8) ) {
        fgets(buffer,MAX_LINELENGTH,fp_fconst);
        count = 0;
      }
      if (sscanf(buffer+count,"%lf",&(F[i][j])) != 1) {
        fprintf(outfile,"\nProblem reading force constants.\n");
        exit(2);
      }
      count += 10;
    }
  }
  fclose(fp_fconst);
  for (i=0;i<dim;++i)
    for (j=0;j<i;++j)
      F[j][i] = F[i][j];

 // Form H_inv = P((PFP)^-1)P
 // H_inv is same as regular inverse if nonredundant coordinates are used
  temp_mat  = init_matrix(dim,dim);
  temp_mat2 = init_matrix(dim,dim);
  H_inv     = init_matrix(dim,dim);
   
  mmult(F,0,P,0,temp_mat,0,dim,dim,dim,0);
  mmult(P,0,temp_mat,0,temp_mat2,0,dim,dim,dim,0);
  temp_mat = symm_matrix_invert(temp_mat2,dim,0,optinfo.redundant);
  mmult(temp_mat,0,P,0,temp_mat2,0,dim,dim,dim,0);
  mmult(P,0,temp_mat2,0,H_inv,0,dim,dim,dim,0);
  free_matrix(F,dim);

  update = optinfo.bfgs;
  if (update) {
    fp_opt_aux = fopen("opt.aux","r");
    if (fp_opt_aux == NULL) update = 0;
  }
  if (update) {
     if (fgets(buffer,MAX_LINELENGTH,fp_opt_aux) == NULL) {
       update = 0;
       fclose(fp_opt_aux);
     }
  }
  if (!update) {
    fprintf(outfile,"\nNo BFGS update performed.\n");
    free_matrix(temp_mat,dim);
    free_matrix(temp_mat2,dim);
    return H_inv;
  }

 // The BFGS update is performed on H_inv
  fprintf(outfile,"\nPerforming BFGS Hessian update.\n");
  double *q_old, *f_old, *dq, *df;
  q_old = init_array(dim);
  f_old = init_array(dim);
  dq    = init_array(dim);
  df    = init_array(dim);

  for (i=0;i<dim;++i)
    if (fgets(buffer,MAX_LINELENGTH,fp_opt_aux) != NULL)
      if (sscanf(buffer, "%lf %lf", &(q_old[i]), &(f_old[i])) != 2) {
        fprintf(outfile,"compute_H: Trouble reading opt.aux\n") ;
        exit(0);
      }

  for (i=0;i<dim;++i) {
    dq[i] = q[i] - q_old[i];
    df[i] = (-1.0) * (f_intcos[i] - f_old[i]); // gradients -- not forces!
  }
  free(q_old);
  free(f_old);

  H_inv_new = init_matrix(dim,dim);
  mmult(&df,0,H_inv,0,temp_mat,0,1,dim,dim,0);
  dot_arr(temp_mat[0],df,dim,&b);
  dot_arr(dq,df,dim,&a);
  mmult(&dq,1,&dq,0,temp_mat,0,dim,1,dim,0);
  for (i=0;i<dim;++i)
    for (j=0;j<dim;++j)
      H_inv_new[i][j] = H_inv[i][j] + ((a + b)/(a*a)) * temp_mat[i][j]; 

  mmult(&df,0,H_inv,0,temp_mat,0,1,dim,dim,0);
  mmult(&dq,1,temp_mat,0,temp_mat2,0,dim,1,dim,0);
  for (i=0;i<dim;++i)
    for (j=0;j<dim;++j)
      H_inv_new[i][j] -= temp_mat2[i][j]/a; 

  mmult(&df,1,&dq,0,temp_mat,0,dim,1,dim,0);
  mmult(H_inv,0,temp_mat,0,temp_mat2,0,dim,dim,dim,0);
  for (i=0;i<dim;++i)
    for (j=0;j<dim;++j)
      H_inv_new[i][j] -= temp_mat2[i][j]/a; 

// Is is necessary to project Hessian after update?
// mmult(H_inv_new,0,P,0,temp_mat,0,dim,dim,dim,0);
// mmult(P,0,temp_mat,0,H_inv_new,0,dim,dim,dim,0);

  free(dq);
  free(df);

 // I think you should do this inversion the same way as above
  H_new = init_matrix(dim,dim);
  mmult(H_inv_new,0,P,0,temp_mat,0,dim,dim,dim,0);
  mmult(P,0,temp_mat,0,temp_mat2,0,dim,dim,dim,0);
  temp_mat = symm_matrix_invert(temp_mat2,dim,0,optinfo.redundant);
  mmult(temp_mat,0,P,0,temp_mat2,0,dim,dim,dim,0);
  mmult(P,0,temp_mat2,0,H_new,0,dim,dim,dim,0);

 //fprintf(outfile,"The Updated Hession Matrix:\n");
 //print_mat2(H_new,dim,dim,outfile);

 // Writing new force constants to fconst.dat
  fp_fconst = fopen("fconst.dat","w");
  for (i=0;i<symm.get_num();++i) {
    col = 0;
    for (j=0; j<=i ; ++j) {
      if (col == 8) {
        fprintf(fp_fconst,"\n");
        col = 0;
      }
      fprintf(fp_fconst,"%10.6f",H_new[i][j]);
      ++col;
    }
    fprintf(fp_fconst,"\n");
  }
  fclose(fp_fconst);

  free_matrix(H_inv,dim);
  free_matrix(H_new,dim);
  free_matrix(temp_mat,dim);
  free_matrix(temp_mat2,dim);
  return H_inv_new;
}

