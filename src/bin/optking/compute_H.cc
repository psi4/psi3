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

void H_update(double **H, internals &simples, salc_set &symm, cartesians &carts);
extern double *compute_q(internals &simples, salc_set &symm);
extern void empirical_H(internals &simples, salc_set &symm, cartesians &carts);

double **compute_H(internals &simples, salc_set &symm, double **P, cartesians &carts) {
  double **H, **H_inv, **H_inv_new, **H_new, **temp_mat;
  int i,j,dim, n_previous;
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
  // Current point has already been put in PSIF_OPTKING by opt_step()
  n_previous = 0;
  open_PSIF();
  if ( psio_tocscan(PSIF_OPTKING, "Num. of Previous Entries") != NULL)
    psio_read_entry(PSIF_OPTKING, "Num. of Previous Entries", (char *) &n_previous, sizeof(int));
  close_PSIF();

  if (optinfo.H_update == OPTInfo::NONE || n_previous < 2 )
    fprintf(outfile,"\nNo Hessian update performed.\n");
  else
    H_update(H, simples, symm, carts);

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

/* This functions performs update of Hessian */
void H_update(double **H, internals &simples, salc_set &symm, cartesians &carts) {
  int i,j,dim,n_previous,i_step, natom, step_start, skip;
  double qq, qg, qz, zz, *q, *f, *q_old, *f_old, *Z;
  double *dq, *dg, **X, **temp_mat, *x, *x_old, phi;
  char force_string[30], x_string[30];

  if (optinfo.H_update == OPTInfo::BFGS)
    fprintf(outfile,"\nPerforming BFGS update");
  else if (optinfo.H_update == OPTInfo::MS)
    fprintf(outfile,"\nPerforming Murtagh/Sargent update");
  else if (optinfo.H_update == OPTInfo::POWELL)
    fprintf(outfile,"\nPerforming Powell update");
  else if (optinfo.H_update == OPTInfo::BOFILL)
    fprintf(outfile,"\nPerforming Bofill update");

  natom = carts.get_natom();
  dim = symm.get_num();

  f = init_array(dim);
  f_old = init_array(dim);
  dq    = init_array(dim);
  dg    = init_array(dim);
  x = init_array(3*natom);
  x_old = init_array(3*natom);

  /*** read/compute current internals and forces from PSIF_OPTKING ***/
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "Num. of Previous Entries", (char *) &n_previous, sizeof(int));
  close_PSIF();

  sprintf(force_string,"Previous Internal Forces %d", n_previous-1);
  sprintf(x_string,"Previous Cartesian Values %d", n_previous-1);
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, force_string, (char *) &(f[0]), dim * sizeof(double));
  psio_read_entry(PSIF_OPTKING, x_string, (char *) &(x[0]), 3*natom*sizeof(double));
  close_PSIF();

  simples.compute_internals(natom,x);
  simples.fix_near_lin(); // fix configuration for torsions
  q = compute_q(simples, symm);

  if (optinfo.H_update_use_last == 0) { /* include all available gradients */
    step_start = 0;
  }
  else {
    step_start = n_previous - optinfo.H_update_use_last - 1;
    if (step_start < 0)
      step_start = 0;
  }

  fprintf(outfile," with previous %d gradient(s).\n", n_previous-1-step_start);

  for (i_step=step_start; i_step<(n_previous-1); ++i_step) {
    /* read/compute old internals and forces from PSIF_OPTKING ***/
    sprintf(force_string,"Previous Internal Forces %d", i_step);
    sprintf(x_string,"Previous Cartesian Values %d", i_step);
    open_PSIF();
    psio_read_entry(PSIF_OPTKING, force_string, (char *) &(f_old[0]), dim * sizeof(double));
    psio_read_entry(PSIF_OPTKING, x_string, (char *) &(x_old[0]), 3*natom*sizeof(double));
    close_PSIF();

    simples.compute_internals(natom, x_old);
    q_old = compute_q(simples, symm);

    /* check for scary passages through 0 */
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

    // compute delta(coordinate) and delta(gradient)
    for (i=0;i<dim;++i) {
      dq[i] = q[i] - q_old[i];
      dg[i] = (-1.0) * (f[i] - f_old[i]); // gradients -- not forces!
    }

    if (optinfo.H_update == OPTInfo::BFGS) {
      // Do BFGS update: Schlegel 1987 Ab Initio Methods in Quantum Chemistry 
      // Let a = DQT.DG and X = (I - DQ*DGT/a)
      // Then Hk = X * H_(k-1) * XT + DQ*DQT/a
      // To make formula work for Hessian (some call it "B")
      // you have to switch DQ and DG in the equation
  
      dot_arr(dq,dg,dim,&qg);
      X = unit_mat(dim);
      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          X[i][j] -= (dg[i] * dq[j]) / qg ; 
  
      temp_mat = block_matrix(dim,dim);
      mmult(X,0,H,0,temp_mat,0,dim,dim,dim,0);
      mmult(temp_mat,0,X,1,H,0,dim,dim,dim,0);
      free_block(temp_mat);
      free_block(X);
    
      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          H[i][j] += dg[i] * dg[j] / qg ; 
    }
    else if (optinfo.H_update == OPTInfo::MS) {
      // Equations taken from Bofill article below
      Z = init_array(dim);
      mmult(H,0,&dq,1,&Z,1,dim,dim,1,0);
      for (i=0;i<dim;++i)
        Z[i] = dg[i] - Z[i];

      dot_arr(dq,Z,dim,&qz);

      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          H[i][j] += Z[i] * Z[j] / qz ;

      free(Z);
    }
    else if (optinfo.H_update == OPTInfo::POWELL) {
      // Equations taken from Bofill article below
      Z = init_array(dim);
      mmult(H,0,&dq,1,&Z,1,dim,dim,1,0);
      for (i=0;i<dim;++i)
        Z[i] = dg[i] - Z[i];

      dot_arr(dq,Z,dim,&qz);
      dot_arr(dq,dq,dim,&qq);

      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          H[i][j] += -1.0*qz/(qq*qq)*dq[i]*dq[j] + (Z[i]*dq[j] + dq[i]*Z[j])/qq;

      free(Z);
    }
    else if (optinfo.H_update == OPTInfo::BOFILL) {
      /* This functions performs a Bofill update on the Hessian according to
      J. M. Bofill, J. Comp. Chem., Vol. 15, pages 1-11 (1994). */
      // Bofill = (1-phi) * MS + phi * Powell
      Z = init_array(dim);
      mmult(H,0,&dq,1,&Z,1,dim,dim,1,0);
      for (i=0;i<dim;++i)
        Z[i] = dg[i] - Z[i];

      dot_arr(dq,Z,dim,&qz);
      dot_arr(dq,dq,dim,&qq);
      dot_arr(Z,Z,dim,&zz);

      phi = 1.0 - qz*qz/(qq*zz);
      if (phi < 0.0) phi = 0.0;
      if (phi > 1.0) phi = 1.0;

      for (i=0;i<dim;++i) 
        for (j=0;j<dim;++j)
          H[i][j] += (1.0-phi) * Z[i] * Z[j] / qz ;

      for (i=0;i<dim;++i)
        for (j=0;j<dim;++j)
          H[i][j] += phi*(-1.0*qz/(qq*qq)*dq[i]*dq[j] + (Z[i]*dq[j] + dq[i]*Z[j])/qq);
      free(Z);
    }
  } //end over old steps

  free(q);
  free(f);
  free(q_old);
  free(f_old);
  free(dq);
  free(dg);
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

