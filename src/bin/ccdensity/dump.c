#include <stdio.h>
#include <libciomr.h>
#include <iwl.h>
#include <dpd.h>
#include <psifiles.h>
#define EXTERN
#include "globals.h"

/* DUMP(): Mulliken-order the two-electron density and dump it to a
** file for subsequent backtransformation.   Basically all we have to
** do is swap indices two and three, e.g.
**
** G'(pr,qs) = G(pq,rs)
**
** In order for the Mulliken-ordered density to be valid for the
** backtransformation algorithm used in TRANSQT, the final density
** must have eight-fold permutational symmetry like the original
** integrals.  Unfortunately, there are a couple of complications
** introduced by the redundant storage I use for open-shell orbitals
** (useful for spin-restricted references --- see the CCSORT code). In
** particular, if the Mulliken-ordered density is not bra-ket
** symmetric, specific elements of the final density must be
** multiplied by two or they will not appear with the correct
** prefactor in the backtransformation.  This only affects the IJKA,
** IAJB, and ABCI Dirac-ordered densities, since the remaining three
** components are bra-ket symmetric in Mulliken order.
**
** I really need to give an example of this problem using specific
** elements of GIJKA so that the code below will be clearer.*/

void dump(struct iwlbuf *OutBuf)
{
  int nirreps, nmo, nfzv;
  int *qt_occ, *qt_vir;
  int h, row, col, p, q, r, s;
  PSI_FPTR next;
  struct dpdbuf G;

  qt_occ = moinfo.qt_occ;  qt_vir = moinfo.qt_vir;
  nirreps = moinfo.nirreps;
  nmo = moinfo.nmo;
  nfzv = moinfo.nfzv;

  rfile(PSIF_MO_OPDM);
  next = 0;
  for(p=0; p < (nmo-nfzv); p++)
      wwritw(PSIF_MO_OPDM, (char *) (moinfo.opdm[p]), 
              sizeof(double)*(nmo-nfzv), next, &next);
  rclose(PSIF_MO_OPDM, 3);

  if(!moinfo.nfzc && !moinfo.nfzv) {  /* Can't deal with frozen orbitals yet */
      rfile(PSIF_MO_LAG);
      next = 0;
      for(p=0; p < (nmo-nfzv); p++)
	  wwritw(PSIF_MO_LAG, (char *) (moinfo.I[p]), 
		 sizeof(double)*(nmo), next, &next);
      rclose(PSIF_MO_LAG, 3);
    }


  dpd_buf_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, "GIjKl", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, prqs, 0, 0, "G(IK,JL)", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 0, 0, 0, 0, 0, "G(IK,JL)", 0, outfile);
  dpd_buf_dump(&G, OutBuf, qt_occ, qt_occ, qt_occ, qt_occ, 1, 0, 0, outfile);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GIjKa", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, prqs, 0, 10, "G(IK,JA)", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 0, 10, 0, 10, 0, "G(IK,JA)", 0, outfile);
  
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);

      for(row=0; row < G.params->rowtot[h]; row++) {
	  p = G.params->roworb[h][row][0];
	  q = G.params->roworb[h][row][1];
	  for(col=0; col < G.params->coltot[h]; col++) {
	      r = G.params->colorb[h][col][0];
	      s = G.params->colorb[h][col][1];

	      if((qt_occ[q] == qt_vir[s]) && (p == r))
		  G.matrix[h][row][col] *= 2;
	    }
	}

      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }

  dpd_buf_dump(&G, OutBuf, qt_occ, qt_occ, qt_occ, qt_vir, 0, 0, 0, outfile);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  dpd_buf_sort(&G, CC_TMP9, prqs, 10, 10, "G(IA,JB)", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP9, 10, 10, 10, 10, 0, "G(IA,JB)", 0, outfile);
  dpd_buf_symm(&G);
  dpd_buf_dump(&G, OutBuf, qt_occ, qt_vir, qt_occ, qt_vir, 1, 0, 0, outfile);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIBJA", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, prqs, 0, 5, "G(IJ,AB)", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 0, 5, 0, 5, 0, "G(IJ,AB)", 0, outfile);
  dpd_scm(&G, 0.5, 0, outfile);

  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);

      for(row=0; row < G.params->rowtot[h]; row++) {
	  p = G.params->roworb[h][row][0];
	  q = G.params->roworb[h][row][1];
	  for(col=0; col < G.params->coltot[h]; col++) {
	      r = G.params->colorb[h][col][0];
	      s = G.params->colorb[h][col][1];

	      if((qt_occ[p] == qt_vir[r]) && (qt_occ[q] == qt_vir[s]))
		  G.matrix[h][row][col] *= 2;
	    }
	}

      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  
  dpd_buf_dump(&G, OutBuf, qt_occ, qt_occ, qt_vir, qt_vir, 0, 0, 0, outfile);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 11, 5, 11, 5, 0, "GCiAb", 0, outfile); 
  dpd_buf_sort(&G, CC_TMP0, prqs, 5, 10, "G(ca,IB)", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 5, 10, 5, 10, 0, "G(ca,IB)", 0, outfile);

  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);

      for(row=0; row < G.params->rowtot[h]; row++) {
	  p = G.params->roworb[h][row][0];
	  q = G.params->roworb[h][row][1];
	  for(col=0; col < G.params->coltot[h]; col++) {
	      r = G.params->colorb[h][col][0];
	      s = G.params->colorb[h][col][1];

	      if((qt_vir[p] == qt_occ[r]) && (q == s))
		  G.matrix[h][row][col] *= 2;
	    }
	}

      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }

  dpd_buf_dump(&G, OutBuf, qt_vir, qt_vir, qt_occ, qt_vir, 0, 0, 0, outfile);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 5, 5, 5, 5, 0, "GAbCd", 0, outfile);
  dpd_buf_sort(&G, CC_TMP0, prqs, 5, 5, "G(AC,BD)", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_TMP0, 5, 5, 5, 5, 0, "G(AC,BD)", 0, outfile);
  dpd_buf_dump(&G, OutBuf, qt_vir, qt_vir, qt_vir, qt_vir, 1, 0, 0, outfile);
  dpd_buf_close(&G);
}
