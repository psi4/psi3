#include <stdio.h>
#include <math.h>
#define EXTERN
#include "globals.h"

void form_diagonal(int irrep) {
  dpdfile2 DIA, Dia, FAE, FMI, Fae, Fmi;
  dpdbuf4 DIJAB, Dijab, DIjAb, Cmnef, CMnEf;
  int *occpi, *virtpi, *occ_off, *vir_off, *openpi;
  int nirreps, *occ_sym, *vir_sym;
  int ij, ab, i, j, a, b, h, I, J, A, B;
  int isym, jsym, asym, bsym;
  char lbl[32];
  double tval;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;
  occ_sym = moinfo.occ_sym; vir_sym = moinfo.vir_sym;
  openpi = moinfo.openpi;

  dpd_file2_init(&FAE, CC_OEI, irrep, 1, 1, "FAE");
  dpd_file2_init(&FMI, CC_OEI, irrep, 0, 0, "FMI");
  dpd_file2_init(&Fae, CC_OEI, irrep, 1, 1, "Fae");
  dpd_file2_init(&Fmi, CC_OEI, irrep, 0, 0, "Fmi");
  dpd_file2_mat_init(&FAE);
  dpd_file2_mat_init(&FMI);
  dpd_file2_mat_init(&Fae);
  dpd_file2_mat_init(&Fmi);
  dpd_file2_mat_rd(&FAE);
  dpd_file2_mat_rd(&FMI);
  dpd_file2_mat_rd(&Fae);
  dpd_file2_mat_rd(&Fmi);

  dpd_file2_init(&DIA, EOM_D, irrep, 0, 1, "DIA");
  dpd_file2_mat_init(&DIA);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < occpi[h]; i++)
      for(a=0; a < (virtpi[h]-openpi[h]); a++)
	DIA.matrix[h][i][a] = FAE.matrix[h][a][a] - FMI.matrix[h][i][i];
  }
  dpd_file2_mat_wrt(&DIA);
  dpd_file2_close(&DIA);

  dpd_file2_init(&Dia, EOM_D, irrep, 0, 1, "Dia");
  dpd_file2_mat_init(&Dia);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < (occpi[h]-openpi[h]); i++)
      for(a=0; a < virtpi[h]; a++) 
	Dia.matrix[h][i][a] = Fae.matrix[h][a][a] - Fmi.matrix[h][i][i];
  }
  dpd_file2_mat_wrt(&Dia);
  dpd_file2_close(&Dia);


  dpd_buf4_init(&DIJAB, EOM_D, irrep, 2, 7, 2, 7, 0, "DIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&DIJAB, h);
    for(ij=0; ij < DIJAB.params->rowtot[h]; ij++) {
      i = DIJAB.params->roworb[h][ij][0];
      j = DIJAB.params->roworb[h][ij][1];
      isym = DIJAB.params->psym[i];
      jsym = DIJAB.params->qsym[j];
      I = i - occ_off[isym];
      J = j - occ_off[jsym];
      for(ab=0; ab < DIJAB.params->coltot[h]; ab++) {
	a = DIJAB.params->colorb[h][ab][0];
	b = DIJAB.params->colorb[h][ab][1];
	asym = DIJAB.params->rsym[a];
	bsym = DIJAB.params->ssym[b];
	A = a - vir_off[asym];
	B = b - vir_off[bsym];
	tval = FAE.matrix[asym][A][A] + FAE.matrix[bsym][B][B]
	  - FMI.matrix[isym][I][I] - FMI.matrix[jsym][J][J];
	DIJAB.matrix[h][ij][ab] =
	  ((A >= (virtpi[asym] - openpi[asym])) ||
	   (B >= (virtpi[bsym] - openpi[bsym])) ?
	   0.0 : tval);
      }
    }
    dpd_buf4_mat_irrep_wrt(&DIJAB, h);
    dpd_buf4_mat_irrep_close(&DIJAB, h);
  }
  dpd_buf4_close(&DIJAB);

  dpd_buf4_init(&Dijab, EOM_D, irrep, 2, 7, 2, 7, 0, "Dijab");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&Dijab, h);
    for(ij=0; ij < Dijab.params->rowtot[h]; ij++) {
      i = Dijab.params->roworb[h][ij][0];
      j = Dijab.params->roworb[h][ij][1];
      isym = Dijab.params->psym[i];
      jsym = Dijab.params->qsym[j];
      I = i - occ_off[isym];
      J = j - occ_off[jsym];
      for(ab=0; ab < Dijab.params->coltot[h]; ab++) {
	a = Dijab.params->colorb[h][ab][0];
	b = Dijab.params->colorb[h][ab][1];
	asym = Dijab.params->rsym[a];
	bsym = Dijab.params->ssym[b];
	A = a - vir_off[asym];
	B = b - vir_off[bsym];
	tval = Fae.matrix[asym][A][A] + Fae.matrix[bsym][B][B]
	  - Fmi.matrix[isym][I][I] - Fmi.matrix[jsym][J][J];
	Dijab.matrix[h][ij][ab] =
	  ((I >= (occpi[isym] - openpi[isym])) ||
	   (J >= (occpi[jsym] - openpi[jsym])) ?
	   0.0 : tval);
      }
    }
    dpd_buf4_mat_irrep_wrt(&Dijab, h);
    dpd_buf4_mat_irrep_close(&Dijab, h);
  }
  dpd_buf4_close(&Dijab);

  dpd_buf4_init(&DIjAb, EOM_D, irrep, 0, 5, 0, 5, 0, "DIjAb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&DIjAb, h);
    for(ij=0; ij < DIjAb.params->rowtot[h]; ij++) {
      i = DIjAb.params->roworb[h][ij][0];
      j = DIjAb.params->roworb[h][ij][1];
      isym = DIjAb.params->psym[i];
      jsym = DIjAb.params->qsym[j];
      I = i - occ_off[isym];
      J = j - occ_off[jsym];
      for(ab=0; ab < DIjAb.params->coltot[h]; ab++) {
	a = DIjAb.params->colorb[h][ab][0];
	b = DIjAb.params->colorb[h][ab][1];
	asym = DIjAb.params->rsym[a];
	bsym = DIjAb.params->ssym[b];
	A = a - vir_off[asym];
	B = b - vir_off[bsym];
	tval = FAE.matrix[asym][A][A] + Fae.matrix[bsym][B][B]
	  - FMI.matrix[isym][I][I] - Fmi.matrix[jsym][J][J];
	DIjAb.matrix[h][ij][ab] =
	  ((A >= (virtpi[asym] - openpi[asym])) ||
	   (J >= (occpi[jsym] - openpi[jsym])) ?
	   0.0 : tval);
      }
    }
    dpd_buf4_mat_irrep_wrt(&DIjAb, h);
    dpd_buf4_mat_irrep_close(&DIjAb, h);
  }
  dpd_buf4_close(&DIjAb);

  dpd_file2_close(&FMI);
  dpd_file2_close(&Fmi);
  dpd_file2_close(&Fae);
  dpd_file2_close(&FAE);

  return;
}
