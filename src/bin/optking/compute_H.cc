// This function reads in 
// Force Constants H from PSIF_OPTKING (in redundant internal coodinates)
// does a BFGS update on H
// inverts H to form H_inv and returns H_inv

#include <cmath>
extern "C" {
#include <stdio.h>
#include <libchkpt/chkpt.h>
#include <stdlib.h>
#include <string.h>
#include <libciomr/libciomr.h>
#include <physconst.h>
#include <libipv1/ip_lib.h>
#include <psifiles.h>
#include <libpsio/psio.h>
}

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

void bfgs(double **H, internals &simples, salc_set &symm, cartesians &carts);
extern double *compute_q(internals &simples, salc_set &symm);

double **compute_H(internals &simples, salc_set &symm, double **P, cartesians &carts) {
  double **H, **H_inv, **H_inv_new, **H_new, **temp_mat;
  int i,j,dim, nbfgs;
  char buffer[MAX_LINELENGTH];

  dim = symm.get_num();
  H = block_matrix(dim,dim);

  /*** Read in force constants from PSIF_OPTKING ***/
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "Force Constants",
      (char *) &(H[0][0]),dim*dim*sizeof(double));
  close_PSIF();
  fprintf(outfile,"\nForce Constants read from PSIF_OPTKING\n");

  // Do BFGS update on H if desired and possible
  nbfgs = 0;
  open_PSIF();
  if ( psio_tocscan(PSIF_OPTKING, "Num. of BFGS Entries") != NULL)
    psio_read_entry(PSIF_OPTKING, "Num. of BFGS Entries", (char *) &nbfgs, sizeof(int));
  close_PSIF();

  if (optinfo.bfgs && (nbfgs > 1) )
    bfgs(H, simples, symm, carts);
  else
    fprintf(outfile,"\nNo BFGS update performed.\n");

  if (optinfo.print_hessian) {
    fprintf(outfile,"The Hessian (Second Derivative) Matrix\n");
    print_mat5(H,dim,dim,outfile);
  }

  // Project redundancies out of H according to Peng JCC 1996
  // and produce invertible Hessian
  // Form H <= PHP + 1000 (1-P)
  if (optinfo.redundant) {
    temp_mat = block_matrix(dim,dim);
    mmult(P,0,H,0,temp_mat,0,dim,dim,dim,0);
    mmult(temp_mat,0,P,0,H,0,dim,dim,dim,0);
    free(temp_mat);

    temp_mat = unit_mat(dim);
    for (i=0;i<dim;++i)
      for (j=0;j<dim;++j) {
        H[i][j] += 1000 * (temp_mat[i][j] - P[i][j]);
      }
    free(temp_mat);
  }

  H_inv = symm_matrix_invert(H,dim,0,0);

  /*** write new force constants H to PSIF_OPTKING ***/
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Force Constants",
      (char *) &(H[0][0]),dim*dim*sizeof(double));
  close_PSIF();
  free_block(H);
  return H_inv;
}

/* This functions performs a BFGS update on H_inv */
void bfgs(double **H, internals &simples, salc_set &symm, cartesians &carts) {
  int i,j,dim,nbfgs,i_bfgs, natom, bfgs_start;
  double a, *q, *f, *q_old, *f_old;
  double *dq, *df, **X, **temp_mat, tval, *x, *x_old;
  char force_string[30], x_string[30];

  fprintf(outfile,"\nPerforming BFGS Hessian update");
  natom = carts.get_natom();

  dim = symm.get_num();
  f = init_array(dim);
  f_old = init_array(dim);
  dq    = init_array(dim);
  df    = init_array(dim);
  x = init_array(3*natom);
  x_old = init_array(3*natom);

  /*** read old internals and forces from PSIF_OPTKING ***/
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "Num. of BFGS Entries", (char *) &nbfgs, sizeof(int));
  close_PSIF();

  sprintf(force_string,"BFGS Internal Forces %d", nbfgs-1);
  sprintf(x_string,"BFGS Cartesian Values %d", nbfgs-1);
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, force_string, (char *) &(f[0]), dim * sizeof(double));
  psio_read_entry(PSIF_OPTKING, x_string, (char *) &(x[0]), 3*natom*sizeof(double));
  close_PSIF();

  /*
  fprintf(outfile,"bfgs internal forces\n");
  for (i=0;i<dim;++i)
    fprintf(outfile,"%15.10lf\n",f[i]);
  fprintf(outfile,"bfgs cartesian values\n");
  for (i=0;i<3*natom;++i)
    fprintf(outfile,"%15.10lf\n",x[i]);
    */

  simples.compute_internals(natom,x);
  simples.fix_near_lin(); // fix configuration for torsions
  q = compute_q(simples, symm);

  if (optinfo.bfgs_use_last == 0) { /* include all available gradients */
    bfgs_start = 0;
  }
  else {
    bfgs_start = nbfgs - optinfo.bfgs_use_last - 1;
    if (bfgs_start < 0)
      bfgs_start = 0;
  }

  fprintf(outfile," with previous %d gradient(s).\n", nbfgs-1-bfgs_start);

  for (i_bfgs=bfgs_start; i_bfgs<(nbfgs-1); ++i_bfgs) {
    sprintf(force_string,"BFGS Internal Forces %d", i_bfgs);
    sprintf(x_string,"BFGS Cartesian Values %d", i_bfgs);
    open_PSIF();
    psio_read_entry(PSIF_OPTKING, force_string, (char *) &(f_old[0]), dim * sizeof(double));
    psio_read_entry(PSIF_OPTKING, x_string, (char *) &(x_old[0]), 3*natom*sizeof(double));
    close_PSIF();

    simples.compute_internals(natom, x_old);
    q_old = compute_q(simples, symm);

    for (i=0;i<dim;++i) {
      dq[i] = q[i] - q_old[i];
      df[i] = (-1.0) * (f[i] - f_old[i]); // gradients -- not forces!
    }

    // Let a = DQT.DG
    // Let X = (I - DQ*DGT/a)
    // Then Hk = X * H_(k-1) * XT + DQ*DQT/a
    // Schlegel 1987 Ab Initio Methods in Quantum Chemistry 
    // as noted in article, to make formula work for Hessian (Schlegel's B)
    // you have to switch DQ and DG in the equation

    dot_arr(dq,df,dim,&a);
    X = unit_mat(dim);
    for (i=0;i<dim;++i)
      for (j=0;j<dim;++j)
        X[i][j] -= (df[i] * dq[j]) / a ; 

    temp_mat = block_matrix(dim,dim);
    mmult(X,0,H,0,temp_mat,0,dim,dim,dim,0);
    mmult(temp_mat,0,X,1,H,0,dim,dim,dim,0);
    free_block(temp_mat);
    free_block(X);
  
    for (i=0;i<dim;++i)
      for (j=0;j<dim;++j)
        H[i][j] += (df[i] * df[j]) / a ; 
  }
  free(q);
  free(f);
  free(q_old);
  free(f_old);
  free(dq);
  free(df);
  free(x);
  free(x_old);
  return;
}
