#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

/* build_hessian(): Compute the vibrational (Cartesian) hessian using
** the skeleton second derivatives, the derivative Fock integrals, and
** the derivative overlap integrals computed by "cints --deriv2", as
** well as the nuclear-perturbation CPHF coefficients computed in
** cphf_X().  The equation for the electronic contribution to the
** Hessian is (MO basis):
**
** d E/da db = E^ab - 2 (S^ab)_ii eps_i - 2 (eta^ab)_ii eps_i
**     + 4 [ (U^b)_pj (F^a)_pj + (U^a)_pj (F^b)_pj ]
**     + 4 (U^a)_pj (U^b)_pj eps_i
**     + 4 (U^a)_pj (U^b)_ql [ 4(pj|ql) - (pq|jl) - (pl|jq) ]
**
** where E^ab is the energy expression evaluated using
** second-derivative integrals, S^ab is the overlap second derivative
** integrals, eps_p is the energy of orbital p, U^a and U^b are CPHF
** coefficients, F^a and F^b are derivative Fock integrals, and eta^ab
** is defined as
**
** (eta^ab)_pq = [(U^a)_pr (U^b)_qr + (U^a)_pr (U^b)_qr 
**                    - (S^a)_pr (S^b)_qr - (S^a)_pr (S^b)_qr ]
**
** For more details (and clearer notation) see:
**
** Y. Yamaguchi et al., "A New Dimension to Quantum Chemistry:
** Analytic Derivative Methods in Ab Initio Molecular Electronic
** Structure Theory", Oxford Press, New York, 1994.  Ch.4,
** esp. pp. 60-63.
**
** TDC, December 2001 (revised October 2002)
*/

void build_hessian(double ***U, double **hessian)
{
  double ***F, ***S;
  int coord, coord_a, coord_b;
  int i, isym, ifirst, ilast;
  int j, jsym, jfirst, jlast;
  int l, lsym, lfirst, llast;
  int p, psym, pfirst, plast;
  int q, qsym, qfirst, qlast;
  int m, ij, pj, pq, jq, ql, jl, pl, pjql, pqjl, pljq;
  char *label;
  double *inbuf, *eta;
  FILE *file15;

  /* Allocate space for the F, S, and B vectors */
  F = (double ***) malloc(natom*3 * sizeof(double **));
  S = (double ***) malloc(natom*3 * sizeof(double **));
  for(coord=0; coord < natom*3; coord++) {
    F[coord] = block_matrix(nmo, nmo);
    S[coord] = block_matrix(nmo, nmo);
  }

  /* Grab the MO-basis overlap and Fock derivative integrals from disk */
  label = (char *) malloc(PSIO_KEYLEN * sizeof(char));
  inbuf = init_array(ntri);
  for(coord=0; coord < natom*3; coord++) {

    sprintf(label, "MO-basis Fock Derivs (%d)", coord);
    iwl_rdone(PSIF_OEI, label, inbuf, ntri, 0, 0, NULL);
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';

    for(i=0, ij=0; i < nmo; i++)
      for(j=0; j <= i; j++, ij++)
	F[coord][i][j] = F[coord][j][i] = inbuf[ij];
  }

  for(coord=0; coord < natom*3; coord++) {
    sprintf(label, "MO-basis Overlap Derivs (%d)", coord);
    iwl_rdone(PSIF_OEI, label, inbuf, ntri, 0, 0, NULL);
    for(i=0; i < PSIO_KEYLEN; i++) label[i] = '\0';

    for(i=0, ij=0; i < nmo; i++)
      for(j=0; j <= i; j++, ij++)
	S[coord][i][j] = S[coord][j][i] = inbuf[ij];
  }
  free(inbuf);
  free(label);

  /* grab skeleton derivatives from cints --deriv2 **/
  psio_open(PSIF_DERINFO, PSIO_OPEN_OLD);
  psio_read_entry(PSIF_DERINFO, "Skeleton Hessian", (char *) hessian[0], natom*3*natom*3*sizeof(double));
  psio_close(PSIF_DERINFO, 1);

  eta = init_array(nmo);
  for(coord_a=0; coord_a < natom*3; coord_a++) {
    for(coord_b=0; coord_b < natom*3; coord_b++) {

      /* build eta intermediate */
      for(isym=0; isym < nirreps; isym++) {
	ifirst = ofirst[isym];
	ilast = olast[isym];
	for(i=ifirst; i <= ilast; i++) {
	  eta[i] = 0.0;
	  for(m=0; m < nmo; m++) {
	    eta[i] += (2.0 * U[coord_a][i][m] * U[coord_b][i][m] - 
		       2.0 * S[coord_a][i][m] * S[coord_b][i][m]);
	  }
	}
      }

      /* eta contribution to vibrational hessian */
      for(isym=0; isym < nirreps; isym++) {
	ifirst = ofirst[isym];
	ilast = olast[isym];
	for(i=ifirst; i <= ilast; i++) 
	  hessian[coord_a][coord_b] -= 2.0 * eta[i] * evals[i];
      }

      for(psym=0; psym < nirreps; psym++) {

	pfirst = first[psym];
	plast = last[psym];

	for(p=pfirst; p <= plast; p++) {

	  for(jsym=0; jsym < nirreps; jsym++) {
	    jfirst = ofirst[jsym];
	    jlast = olast[jsym];

	    for(j=jfirst; j <= jlast; j++) {

	      pj = INDEX(p,j);

	      hessian[coord_a][coord_b] += 4.0 * (U[coord_b][p][j] * F[coord_a][p][j] + 
						  U[coord_a][p][j] * F[coord_b][p][j]);

	      hessian[coord_a][coord_b] += 4.0 * U[coord_a][p][j] * U[coord_b][p][j] * evals[p];

	      for(qsym=0; qsym < nirreps; qsym++) {
		qfirst = first[qsym];
		qlast = last[qsym];

		for(q=qfirst; q <= qlast; q++) {

		  pq = INDEX(p,q);
		  jq = INDEX(j,q);

		  for(lsym=0; lsym < nirreps; lsym++) {

		    lfirst = ofirst[lsym];
		    llast = olast[lsym];

		    for(l=lfirst; l <= llast; l++) {

		      ql = INDEX(q,l);
		      jl = INDEX(j,l);
		      pl = INDEX(p,l);

		      pjql = INDEX(pj,ql);
		      pqjl = INDEX(pq,jl);
		      pljq = INDEX(pl,jq);

		      hessian[coord_a][coord_b] += 4.0 * U[coord_a][p][j] * U[coord_b][q][l] *
			(4.0 * ints[pjql] - ints[pqjl] - ints[pljq]);
		    }
		  }
		}
	      }

	    }
	  }
	}
      }

    }
  }

  free(eta);
  for(coord=0; coord < natom*3; coord++) {
    free_block(F[coord]);
    free_block(S[coord]);
  }
  free(F); free(S);

  /*
  fprintf(outfile, "SCF Molecular Hessian:\n");
  print_mat(hessian, natom*3, natom*3, outfile);
  */

  /* write the hessian to file15 in the PSI2 standard format */
  ffile(&file15, "file15.dat", 0);
  fprintf(file15, "%5d%5d\n", natom, natom*6);
  for(i=0; i < natom*3; i++) {
    for(j=0; j < natom; j++) {
      fprintf(file15, "%20.10f%20.10f%20.10f\n", 
	      hessian[i][j*3], hessian[i][j*3+1], hessian[i][j*3+2]);
    }
  }
  fclose(file15);
}
