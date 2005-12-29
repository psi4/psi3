#include <stdio.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"
#include <math.h>

/** AEL() computes the approximate excitation level according to
 ** Stanton and Bartlett, JCP, 98, 1993, 7034.
 ** Trace [rho(excited) - rho(ground)] = AEL
 ** where both densities are expressed in the bases that diagonalizes
 ** the ground-state CCSD density. **/

void ael(struct RHO_Params *rho_params)
{
  int dim,i,j;
  double **rho_g, *evals, **evects, **tmat, **rho_x, ael;

  dim = moinfo.nmo - moinfo.nfzv;
  rho_g = block_matrix(dim,dim);

  /* read in and diagonalize the ground-state rho */
  psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
  psio_read_entry(PSIF_MO_OPDM, rho_params[0].opdm_lbl, (char *) &(rho_g[0][0]),
		   sizeof(double)*dim*dim);
  psio_close(PSIF_MO_OPDM, 1);

  evals = init_array(dim);
  evects = block_matrix(dim,dim);
  tmat = block_matrix(dim,dim);

  sq_rsp(dim, dim, rho_g, evals, 3, evects, 1.0E-14);
  C_DGEMM('t', 'n', dim, dim, dim, 1.0, &(evects[0][0]), dim,
    &(rho_g[0][0]), dim, 0.0, &(tmat[0][0]), dim);
  C_DGEMM('n', 'n', dim, dim, dim, 1.0,   &(tmat[0][0]), dim,
    &(evects[0][0]), dim, 0.0, &(rho_g[0][0]), dim);

  fprintf(outfile,"\n");
  rho_x = block_matrix(dim,dim);
  for (i=1; i<params.nstates; ++i) {
    /* read in and transform the ith density */
    psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
    psio_read_entry(PSIF_MO_OPDM, rho_params[i].opdm_lbl, (char *) &(rho_x[0][0]),
		     sizeof(double)*dim*dim);
    psio_close(PSIF_MO_OPDM, 1);

    C_DGEMM('t', 'n', dim, dim, dim, 1.0, &(evects[0][0]), dim,
      &(rho_x[0][0]), dim, 0.0, &(tmat[0][0]), dim);
    C_DGEMM('n', 'n', dim, dim, dim, 1.0,   &(tmat[0][0]), dim,
      &(evects[0][0]), dim, 0.0, &(rho_x[0][0]), dim);
   /* mat_print(rho_x, dim, dim, outfile); */

    /* compute the ith AEL */
    ael = 0.0;
    for (j=0; j<dim; ++j) {
      ael += 0.5 * fabs(rho_x[j][j] - rho_g[j][j]);
    }
    fprintf(outfile, "\tAEL (approximate excitation level) for excited state %d: %10.7lf\n",
      i, ael);
  }

  free(evals);
  free_block(evects);
  free_block(tmat);
  free_block(rho_g);
  free_block(rho_x);
  return;
}

