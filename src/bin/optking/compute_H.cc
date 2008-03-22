/*! \file
    \ingroup OPTKING
    \brief This function reads in Force Constants H from
    PSIF_OPTKING (in redundant internal coodinates) does a
    BFGS update on H inverts H to form H_inv and returns H_inv.
*/

#include <cmath>
#include <cstdio>
#include <libchkpt/chkpt.h>
#include <cstdlib>
#include <cstring>
#include <libciomr/libciomr.h>
#include <physconst.h>
#include <libipv1/ip_lib.h>
#include <psifiles.h>
#include <libpsio/psio.h>

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

namespace psi { namespace optking {

void bfgs(double **H, internals &simples, salc_set &symm, cartesians &carts);
extern double *compute_q(internals &simples, salc_set &symm);
extern void empirical_H(internals &simples, salc_set &symm, cartesians &carts);

double **compute_H(internals &simples, salc_set &symm, double **P, cartesians &carts) {
  double **H, **H_inv, **H_inv_new, **H_new, **temp_mat;
  int i,j,dim, nbfgs;
  char buffer[MAX_LINELENGTH];

  dim = symm.get_num();
  H = block_matrix(dim,dim);

  /*** Read in force constants from PSIF_OPTKING ***/
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "Symmetric Force Constants",
      (char *) &(H[0][0]),dim*dim*sizeof(double));
  close_PSIF();
  fprintf(outfile,"\nForce Constants read from PSIF_OPTKING\n");

  // Do BFGS update on H if desired and possible
  nbfgs = 0;
  open_PSIF();
  if ( psio_tocscan(PSIF_OPTKING, "Num. of Previous Entries") != NULL)
    psio_read_entry(PSIF_OPTKING, "Num. of Previous Entries", (char *) &nbfgs, sizeof(int));
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
  psio_write_entry(PSIF_OPTKING, "Symmetric Force Constants",
      (char *) &(H[0][0]),dim*dim*sizeof(double));
  close_PSIF();
  free_block(H);
  return H_inv;
}

/* This functions performs a BFGS update on H_inv */
void bfgs(double **H, internals &simples, salc_set &symm, cartesians &carts) {
  int i,j,dim,nbfgs,i_bfgs, natom, bfgs_start, skip;
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
  psio_read_entry(PSIF_OPTKING, "Num. of Previous Entries", (char *) &nbfgs, sizeof(int));
  close_PSIF();

  sprintf(force_string,"Previous Internal Forces %d", nbfgs-1);
  sprintf(x_string,"Previous Cartesian Values %d", nbfgs-1);
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
    sprintf(force_string,"Previous Internal Forces %d", i_bfgs);
    sprintf(x_string,"Previous Cartesian Values %d", i_bfgs);
    open_PSIF();
    psio_read_entry(PSIF_OPTKING, force_string, (char *) &(f_old[0]), dim * sizeof(double));
    psio_read_entry(PSIF_OPTKING, x_string, (char *) &(x_old[0]), 3*natom*sizeof(double));
    close_PSIF();

    simples.compute_internals(natom, x_old);
    q_old = compute_q(simples, symm);

    skip=0;
    for (i=0;i<dim;++i) {
      if (q[i] * q_old[i] < 0.0)
        skip = 1;
    }
    if (skip) {
      fprintf(outfile,"Warning a coordinate has passed through 0. Skipping a BFGS update.\n");
      skip = 0;
      continue;
    }

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

/*! \file fconst_init()
    \ingroup OPTKING
    \brief 
  FCONST_INIT -- make sure there are _some_ force constants in PSIF_OPTKING
  1) Confirm PSIF_OPTKING has them
  2) read them from FCONST: section of input
  3) read them from fconst.dat
  4) generate empirical force constants
*/

void fconst_init(cartesians &carts, internals &simples, salc_set &symm) {
  int i, j, dim, count, constants_in_PSIF, cnt;
  char *buffer;
  double **F, **temp_mat;
  buffer = new char[MAX_LINELENGTH];

  open_PSIF();
  if (psio_tocscan(PSIF_OPTKING, "Symmetric Force Constants") != NULL) {
    close_PSIF();
    return;
  }
  close_PSIF();

   /* read force constants from fconst section of input */
  if (ip_exist(":FCONST",0) ) {
    ip_cwk_add(":FCONST");
    fprintf(outfile,"Reading force constants from FCONST: \n");
    dim = symm.get_num();
    temp_mat = block_matrix(dim,dim);
    ip_count("SYMM_FC",&i,0);
    if (i != (dim*(dim+1))/2) {
      fprintf(outfile,"fconst has wrong number of entries\n");
      exit(2);
    }
    cnt = -1;
    for (i=0;i<dim;++i) {
      for (j=0;j<=i;++j) {
        ++cnt;
        ip_data("SYMM_FC","%lf",&(temp_mat[i][j]), 1, cnt);
      }
    }
    for (i=0;i<dim;++i)
      for (j=0;j<=i;++j)
        temp_mat[j][i] = temp_mat[i][j];

    /*** write to PSIF_OPTKING ***/
    open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Symmetric Force Constants",
        (char *) &(temp_mat[0][0]),dim*dim*sizeof(double));
    close_PSIF();

    free_block(temp_mat);
    return;
  }

  ffile_noexit(&fp_fconst, "fconst.dat",2);
  if (fp_fconst == NULL) { // generate empirical Hessian
    fprintf(outfile, "\nGenerating empirical Hessian.\n");
    empirical_H(simples,symm,carts);
    return;
  }
  else { // read force constants from fconst.dat
    ip_append(fp_fconst, outfile);
    if (ip_exist(":FCONST",0) ) { // file has libipv1 format
      ip_cwk_add(":FCONST");
      fprintf(outfile,"Reading force constants from FCONST: \n");
      dim = symm.get_num();
      temp_mat = block_matrix(dim,dim);
      ip_count("SYMM_FC",&i,0);
      if (i != (dim*(dim+1))/2) {
        fprintf(outfile,"fconst has wrong number of entries\n");
        exit(2);
      }
      cnt = -1;
      for (i=0;i<dim;++i) {
        for (j=0;j<=i;++j) {
          ++cnt;
          ip_data("SYMM_FC","%lf",&(temp_mat[i][j]), 1, cnt);
        }
      }
    }
  
/*    else {
      fprintf(outfile,"Reading force constants from fconst.dat\n");
      dim = symm.get_num();
      temp_mat = block_matrix(dim,dim);
      for (i=0;i<dim;i++) {
        for (j=0;j<=i;j++) {
          fscanf(fp_fconst,"%lf",&(temp_mat[i][j]));
          temp_mat[j][i] = temp_mat[i][j];
        }
      }
    } */
    fclose(fp_fconst);

    for (i=0;i<dim;++i)
      for (j=0;j<i;++j)
        temp_mat[j][i] = temp_mat[i][j];

    open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Symmetric Force Constants",
        (char *) &(temp_mat[0][0]),dim*dim*sizeof(double));
    close_PSIF();
    free_block(temp_mat);
  }
  delete [] buffer;
}

}} /* namespace psi::optking */

