#include <stdio.h>
#include <dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

void denom(void)
{
  int nirreps;
  int h, i, j, a, b, ij, ab;
  int I, J, A, B;
  int isym, jsym, asym, bsym;
  int *occpi, *virtpi;
  int *occ_off, *vir_off;
  int *openpi;
  double fii, fjj, faa, fbb;
  struct oe_dpdfile fIJ, fij, fAB, fab;
  struct oe_dpdfile dIA, dia;
  struct dpdfile dIJAB, dijab, dIjAb;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;

  /* Grab Fock matrices from disk */
  dpd_oe_file_init(&fIJ, CC_OEI, 0, 0, "fIJ", 0, outfile);
  dpd_oe_file_mat_init(&fIJ);
  dpd_oe_file_mat_rd(&fIJ, 0, outfile);
  
  dpd_oe_file_init(&fij, CC_OEI, 0, 0, "fij", 0, outfile);
  dpd_oe_file_mat_init(&fij);
  dpd_oe_file_mat_rd(&fij, 0, outfile);
  
  dpd_oe_file_init(&fAB, CC_OEI, 1, 1, "fAB", 0, outfile);
  dpd_oe_file_mat_init(&fAB);
  dpd_oe_file_mat_rd(&fAB, 0, outfile);
  
  dpd_oe_file_init(&fab, CC_OEI, 1, 1, "fab", 0, outfile);
  dpd_oe_file_mat_init(&fab);
  dpd_oe_file_mat_rd(&fab, 0, outfile);

  /* Alpha one-electron denominator */
  dpd_oe_file_init(&dIA, CC_OEI, 0, 1, "dIA", 0, outfile);
  dpd_oe_file_mat_init(&dIA);

  for(h=0; h < nirreps; h++) {
      
      for(i=0; i < occpi[h]; i++) {
	  fii = fIJ.matrix[h][i][i];

	  for(a=0; a < (virtpi[h] - openpi[h]); a++) {
	      faa = fAB.matrix[h][a][a];

	      dIA.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
      
    }

  dpd_oe_file_mat_wrt(&dIA, 0, outfile);
  dpd_oe_file_mat_close(&dIA);
  dpd_oe_file_close(&dIA);

  /* Beta one-electron denominator */
  dpd_oe_file_init(&dia, CC_OEI, 0, 1, "dia", 0, outfile);
  dpd_oe_file_mat_init(&dia);

  for(h=0; h < nirreps; h++) {
      
      for(i=0; i < (occpi[h] - openpi[h]); i++) {
	  fii = fij.matrix[h][i][i];
 
	  for(a=0; a < virtpi[h]; a++) {
	      faa = fab.matrix[h][a][a];

	      dia.matrix[h][i][a] = 1.0/(fii - faa);
	    }
	}
      
    }
  
  dpd_oe_file_mat_wrt(&dia, 0, outfile);
  dpd_oe_file_mat_close(&dia);
  dpd_oe_file_close(&dia);

  /* Alpha-alpha two-electron denominator */
  dpd_file_init(&dIJAB, CC_DENOM, 1, 6, "dIJAB", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_file_mat_irrep_init(&dIJAB, h);

      /* Loop over the rows */
      for(ij=0; ij < dIJAB.params->rowtot[h]; ij++) {
	  i = dIJAB.params->roworb[h][ij][0];
	  j = dIJAB.params->roworb[h][ij][1];
	  isym = dIJAB.params->psym[i];
	  jsym = dIJAB.params->qsym[j];

	  /* Convert to relative orbital index */
	  I = i - occ_off[isym];
	  J = j - occ_off[jsym];

	  fii = fIJ.matrix[isym][I][I];
	  fjj = fIJ.matrix[jsym][J][J];
	  
	  /* Loop over the columns */
	  for(ab=0; ab < dIJAB.params->coltot[h]; ab++) {
	      a = dIJAB.params->colorb[h][ab][0];
	      b = dIJAB.params->colorb[h][ab][1];
	      asym = dIJAB.params->rsym[a];
	      bsym = dIJAB.params->ssym[b];

	      /* Convert to relative orbital index */
	      A = a - vir_off[asym];
	      B = b - vir_off[bsym];

	      faa = fAB.matrix[asym][A][A];
	      fbb = fAB.matrix[bsym][B][B];

	      dIJAB.matrix[h][ij][ab] =
                ((A >= (virtpi[asym] - openpi[asym])) ||
		 (B >= (virtpi[bsym] - openpi[bsym])) ?
		 0.0 : 1.0/(fii + fjj - faa - fbb));
	    }
	}

      dpd_file_mat_irrep_wrt(&dIJAB, h, 0, outfile);
      dpd_file_mat_irrep_close(&dIJAB, h);

    }

  dpd_file_close(&dIJAB);

  /* Beta-beta two-electron denominator */
  dpd_file_init(&dijab, CC_DENOM, 1, 6, "dijab", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_file_mat_irrep_init(&dijab, h);

      /* Loop over the rows */
      for(ij=0; ij < dijab.params->rowtot[h]; ij++) {
	  i = dijab.params->roworb[h][ij][0];
	  j = dijab.params->roworb[h][ij][1];
	  isym = dijab.params->psym[i];
	  jsym = dijab.params->qsym[j];

	  /* Convert to relative orbital index */
	  I = i - occ_off[isym];
	  J = j - occ_off[jsym];

	  fii = fij.matrix[isym][I][I];
	  fjj = fij.matrix[jsym][J][J];
	  
	  /* Loop over the columns */
	  for(ab=0; ab < dijab.params->coltot[h]; ab++) {
	      a = dijab.params->colorb[h][ab][0];
	      b = dijab.params->colorb[h][ab][1];
	      asym = dijab.params->rsym[a];
	      bsym = dijab.params->ssym[b];

	      /* Convert to relative orbital index */
	      A = a - vir_off[asym];
	      B = b - vir_off[bsym];

	      faa = fab.matrix[asym][A][A];
	      fbb = fab.matrix[bsym][B][B];

	      dijab.matrix[h][ij][ab] =
                ((I >= (occpi[isym] - openpi[isym])) ||
		 (J >= (occpi[jsym] - openpi[jsym])) ?
		 0.0 : 1.0/(fii + fjj - faa - fbb));
	    }
	}

      dpd_file_mat_irrep_wrt(&dijab, h, 0, outfile);
      dpd_file_mat_irrep_close(&dijab, h);

    }

  dpd_file_close(&dijab);

  /* Alpha-beta two-electron denominator */
  dpd_file_init(&dIjAb, CC_DENOM, 0, 5, "dIjAb", 0, outfile);

  for(h=0; h < nirreps; h++) {

      dpd_file_mat_irrep_init(&dIjAb, h);

      /* Loop over the rows */
      for(ij=0; ij < dIjAb.params->rowtot[h]; ij++) {
	  i = dIjAb.params->roworb[h][ij][0];
	  j = dIjAb.params->roworb[h][ij][1];
	  isym = dIjAb.params->psym[i];
	  jsym = dIjAb.params->qsym[j];

	  /* Convert to relative orbital index */
	  I = i - occ_off[isym];
	  J = j - occ_off[jsym];

	  fii = fIJ.matrix[isym][I][I];
	  fjj = fij.matrix[jsym][J][J];
	  
	  /* Loop over the columns */
	  for(ab=0; ab < dIjAb.params->coltot[h]; ab++) {
	      a = dIjAb.params->colorb[h][ab][0];
	      b = dIjAb.params->colorb[h][ab][1];
	      asym = dIjAb.params->rsym[a];
	      bsym = dIjAb.params->ssym[b];

	      /* Convert to relative orbital index */
	      A = a - vir_off[asym];
	      B = b - vir_off[bsym];

	      faa = fAB.matrix[asym][A][A];
	      fbb = fab.matrix[bsym][B][B];

	      dIjAb.matrix[h][ij][ab] =
                ((A >= (virtpi[asym] - openpi[asym])) ||
		 (J >= (occpi[jsym] - openpi[jsym])) ?
		 0.0 : 1.0/(fii + fjj - faa - fbb));
	    }
	}

      dpd_file_mat_irrep_wrt(&dIjAb, h, 0, outfile);
      dpd_file_mat_irrep_close(&dIjAb, h);

    }

  dpd_file_close(&dIjAb);

  dpd_oe_file_mat_close(&fIJ);
  dpd_oe_file_mat_close(&fij);
  dpd_oe_file_mat_close(&fAB);
  dpd_oe_file_mat_close(&fab);
  dpd_oe_file_close(&fIJ);
  dpd_oe_file_close(&fij);
  dpd_oe_file_close(&fAB);
  dpd_oe_file_close(&fab);

}
