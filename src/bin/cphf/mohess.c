#include <stdio.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

/* mohess(): Builds the molecular orbital hessian matrix, which is
** needed for solving the CPHF equations and for Hartree-Fock wave
** function stability analysis.
**
** The MO hessian may be expressed as:
**
** A_ai,bj = delta_ab delta_ij (eps_i - eps_a) 
**            - [4 (ai|bj) - (ab|ij) - (aj|ib) ]
**
** For more details (and clearer notation) see:
**
** Y. Yamaguchi et al., "A New Dimension to Quantum Chemistry:
** Analytic Derivative Methods in Ab Initio Molecular Electronic
** Structure Theory", Oxford Press, New York, 1994.  Ch.10,
** esp. pp. 128-132.
**
** Note that, although the loop structure here involves irreps,
** symmetry is only used because of the Pitzer ordering of the MO's
** (all orbitals together in each irrep).  Symmetry is not actually
** used to streamline the calculation or storage.  This could be done
** (and *should* be) by computing offsets for symmetry blocks similar
** to the direct-product decomposition approach used in the PSI3
** coupled cluster codes.  However, when the CPHF equations are
** solved, the perturbations themselves must be symmetrized or the
** irrep structure of the CPHF coefficients will be lost.
**
** TDC, December 2001 (revised October 2002)
*/

void mohess(double **A)
{
  int asym, isym, bsym, jsym;
  int a, i, b, j;
  int ai, bj, ab, ij, ib, aj, aibj, abij, ajib;
  int AI, BJ;
  int afirst, alast, ifirst, ilast, bfirst, blast, jfirst, jlast;

  for(asym=0,AI=0; asym < nirreps; asym++) {

    afirst = vfirst[asym];
    alast = vlast[asym];

    for(a=afirst; a <= alast; a++) {

      for(isym=0; isym < nirreps; isym++) {
	ifirst = ofirst[isym];
	ilast = olast[isym];

	for(i=ifirst; i <= ilast; i++,AI++) {
	  ai = INDEX(a,i);

	  for(bsym=0,BJ=0; bsym < nirreps; bsym++) {

	    bfirst = vfirst[bsym];
	    blast = vlast[bsym];

	    for(b=bfirst; b <= blast; b++) {
	      ab = INDEX(a,b);
	      ib = INDEX(i,b);

	      for(jsym=0; jsym < nirreps; jsym++) {

		jfirst = ofirst[jsym];
		jlast = olast[jsym];

		for(j=jfirst; j <= jlast; j++,BJ++) {
		  bj = INDEX(b,j);
		  ij = INDEX(i,j);
		  aj = INDEX(a,j);

		  aibj = INDEX(ai,bj);
		  abij = INDEX(ab,ij);
		  ajib = INDEX(aj,ib);

		  A[AI][BJ] = (a==b) * (i==j) * (evals[i] - evals[a]);
		  A[AI][BJ] -= 4.0 * ints[aibj] - ints[abij] - ints[ajib];

		}
	      }
	    }
	  }
	}
      }
    }
  }

  /* dump the hessian to disk */
  psio_open(PSIF_CPHF, 1);
  psio_write_entry(PSIF_CPHF, "RHF MO Hessian", (char *) A[0], num_ai*num_ai*sizeof(double));
  psio_close(PSIF_CPHF, 1);

  if(print_lvl & 16) {
    fprintf(outfile, "\nMO Hessian:\n");
    print_mat(A, num_ai, num_ai, outfile);
  }

  return;
}
