#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void denom_rohf(int L_irr, int e_index);
void denom_uhf(int L_irr, int e_index);

void denom(int L_irr, int e_index) {
  if(params.ref == 0 || params.ref == 1) denom_rohf(L_irr, e_index);
  else if(params.ref == 2) denom_uhf(L_irr, e_index);
}

void denom_uhf(int L_irr, int e_index)
{
  int nirreps, h, i, j, a, b, ij, ab, I, J, A, B, isym, jsym, asym, bsym;
  int *aoccpi, *boccpi, *avirtpi, *bvirtpi; 
  int *aocc_off, *bocc_off, *avir_off, *bvir_off;
  dpdfile2 LFMIt, LFmit, LFaet, LFAEt;
  dpdfile2 dIA, dia;
  dpdfile4 dIJAB, dijab, dIjAb;
  double Fii, Fjj, Faa, Fbb;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi; 
  boccpi = moinfo.boccpi; 
  avirtpi = moinfo.avirtpi;
  bvirtpi = moinfo.bvirtpi;
  aocc_off = moinfo.aocc_off;
  bocc_off = moinfo.bocc_off;
  avir_off = moinfo.avir_off;
  bvir_off = moinfo.bvir_off;

  dpd_file2_init(&LFMIt, CC_OEI, 0, 0, 0, "FMI");
  dpd_file2_mat_init(&LFMIt);
  dpd_file2_mat_rd(&LFMIt);

  dpd_file2_init(&LFmit, CC_OEI, 0, 2, 2, "Fmi");
  dpd_file2_mat_init(&LFmit);
  dpd_file2_mat_rd(&LFmit);

  dpd_file2_init(&LFaet, CC_OEI, 0, 3, 3, "Fae");
  dpd_file2_mat_init(&LFaet);
  dpd_file2_mat_rd(&LFaet);

  dpd_file2_init(&LFAEt, CC_OEI, 0, 1, 1, "FAE");
  dpd_file2_mat_init(&LFAEt);
  dpd_file2_mat_rd(&LFAEt);

  dpd_file2_init(&dIA, CC_OEI, L_irr, 0, 1, "dIA");
  dpd_file2_mat_init(&dIA);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < aoccpi[h]; i++) {
      Fii = LFMIt.matrix[h][i][i];
      for(a=0; a < avirtpi[h^L_irr]; a++) {
        Faa = LFAEt.matrix[h^L_irr][a][a];
        dIA.matrix[h][i][a] = 1.0/(Fii - Faa + params.cceom_energy[L_irr][e_index]);
      }
    }
  }
  dpd_file2_mat_wrt(&dIA);
  dpd_file2_mat_close(&dIA);
  dpd_file2_close(&dIA);

  dpd_file2_init(&dia, CC_OEI, L_irr, 2, 3, "dia");
  dpd_file2_mat_init(&dia);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < boccpi[h]; i++) {
      Fii = LFmit.matrix[h][i][i];
      for(a=0; a < bvirtpi[h^L_irr]; a++) {
        Faa = LFaet.matrix[h^L_irr][a][a];
        dia.matrix[h][i][a] = 1.0/(Fii - Faa + params.cceom_energy[L_irr][e_index]);
      }
    }
  }
  dpd_file2_mat_wrt(&dia);
  dpd_file2_mat_close(&dia);
  dpd_file2_close(&dia);

  dpd_file4_init(&dIJAB, CC_DENOM, L_irr, 1, 6, "dIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&dIJAB, h);
    for(ij=0; ij < dIJAB.params->rowtot[h]; ij++) {
      i = dIJAB.params->roworb[h][ij][0];
      j = dIJAB.params->roworb[h][ij][1];
      isym = dIJAB.params->psym[i];
      jsym = dIJAB.params->qsym[j];
      I = i - aocc_off[isym];
      J = j - aocc_off[jsym];
      Fii = LFMIt.matrix[isym][I][I];
      Fjj = LFMIt.matrix[jsym][J][J];

      for(ab=0; ab < dIJAB.params->coltot[h^L_irr]; ab++) {
	a = dIJAB.params->colorb[h^L_irr][ab][0];
	b = dIJAB.params->colorb[h^L_irr][ab][1];
	asym = dIJAB.params->rsym[a];
	bsym = dIJAB.params->ssym[b];
	A = a - avir_off[asym];
	B = b - avir_off[bsym];
	Faa = LFAEt.matrix[asym][A][A];
	Fbb = LFAEt.matrix[bsym][B][B];

	dIJAB.matrix[h][ij][ab] = 1.0/(Fii + Fjj - Faa - Fbb
            + params.cceom_energy[L_irr][e_index]);
      }
    }
    dpd_file4_mat_irrep_wrt(&dIJAB, h);
    dpd_file4_mat_irrep_close(&dIJAB, h);
  }
  dpd_file4_close(&dIJAB);

  dpd_file4_init(&dijab, CC_DENOM, L_irr, 11, 16, "dijab");

  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&dijab, h);
    for(ij=0; ij < dijab.params->rowtot[h]; ij++) {
      i = dijab.params->roworb[h][ij][0];
      j = dijab.params->roworb[h][ij][1];
      isym = dijab.params->psym[i];
      jsym = dijab.params->qsym[j];
      I = i - bocc_off[isym];
      J = j - bocc_off[jsym];
      Fii = LFmit.matrix[isym][I][I];
      Fjj = LFmit.matrix[jsym][J][J];

      for(ab=0; ab < dijab.params->coltot[h^L_irr]; ab++) {
	a = dijab.params->colorb[h^L_irr][ab][0];
	b = dijab.params->colorb[h^L_irr][ab][1];
	asym = dijab.params->rsym[a];
	bsym = dijab.params->ssym[b];
	A = a - bvir_off[asym];
	B = b - bvir_off[bsym];
	Faa = LFaet.matrix[asym][A][A];
	Fbb = LFaet.matrix[bsym][B][B];

	dijab.matrix[h][ij][ab] = 1.0/(Fii + Fjj - Faa - Fbb
            + params.cceom_energy[L_irr][e_index]);
      }
    }
    dpd_file4_mat_irrep_wrt(&dijab, h);
    dpd_file4_mat_irrep_close(&dijab, h);
  }
  dpd_file4_close(&dijab);

  dpd_file4_init(&dIjAb, CC_DENOM, L_irr, 22, 28, "dIjAb");

  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&dIjAb, h);
    for(ij=0; ij < dIjAb.params->rowtot[h]; ij++) {
      i = dIjAb.params->roworb[h][ij][0];
      j = dIjAb.params->roworb[h][ij][1];
      isym = dIjAb.params->psym[i];
      jsym = dIjAb.params->qsym[j];
      I = i - aocc_off[isym];
      J = j - bocc_off[jsym];
      Fii = LFMIt.matrix[isym][I][I];
      Fjj = LFmit.matrix[jsym][J][J];

      for(ab=0; ab < dIjAb.params->coltot[h^L_irr]; ab++) {
	a = dIjAb.params->colorb[h^L_irr][ab][0];
	b = dIjAb.params->colorb[h^L_irr][ab][1];
	asym = dIjAb.params->rsym[a];
	bsym = dIjAb.params->ssym[b];
	A = a - avir_off[asym];
	B = b - bvir_off[bsym];
	Faa = LFAEt.matrix[asym][A][A];
	Fbb = LFaet.matrix[bsym][B][B];

	dIjAb.matrix[h][ij][ab] = 1.0/(Fii + Fjj - Faa - Fbb
            + params.cceom_energy[L_irr][e_index]);
      }
    }
    dpd_file4_mat_irrep_wrt(&dIjAb, h);
    dpd_file4_mat_irrep_close(&dIjAb, h);
  }
  dpd_file4_close(&dIjAb);

  dpd_file2_mat_close(&LFMIt);
  dpd_file2_mat_close(&LFmit);
  dpd_file2_mat_close(&LFAEt);
  dpd_file2_mat_close(&LFaet);
  dpd_file2_close(&LFMIt);
  dpd_file2_close(&LFmit);
  dpd_file2_close(&LFAEt);
  dpd_file2_close(&LFaet);

  return;
}

void denom_rohf(int L_irr, int e_index)
{
  dpdfile2 LFAEt, LFaet, LFMIt, LFmit;
  dpdfile2 dIA, dia;
  dpdfile4 dIJAB, dijab, dIjAb;
  dpdbuf4 d, bdIJAB, bdijab, bdIjAb;
  double tval;

  int nirreps;
  int h, i, j, a, b, ij, ab;
  int I, J, A, B;
  int isym, jsym, asym, bsym;
  int *occpi, *virtpi;
  int *occ_off, *vir_off;
  int *openpi;
  double Fii, Fjj, Faa, Fbb;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  openpi = moinfo.openpi;
  occ_off = moinfo.occ_off; vir_off = moinfo.vir_off;

  dpd_file2_init(&LFMIt, CC_OEI, 0, 0, 0, "FMI");
  dpd_file2_mat_init(&LFMIt);
  dpd_file2_mat_rd(&LFMIt);

  dpd_file2_init(&LFmit, CC_OEI, 0, 0, 0, "Fmi");
  dpd_file2_mat_init(&LFmit);
  dpd_file2_mat_rd(&LFmit);

  dpd_file2_init(&LFaet, CC_OEI, 0, 1, 1, "Fae");
  dpd_file2_mat_init(&LFaet);
  dpd_file2_mat_rd(&LFaet);

  dpd_file2_init(&LFAEt, CC_OEI, 0, 1, 1, "FAE");
  dpd_file2_mat_init(&LFAEt);
  dpd_file2_mat_rd(&LFAEt);

  /* Alpha one-electron denominator */
  dpd_file2_init(&dIA, CC_OEI, L_irr, 0, 1, "dIA");
  dpd_file2_mat_init(&dIA);
  for(h=0; h < nirreps; h++) { /* irreps of dIA and Fii */
    for(i=0; i < occpi[h]; i++) {
      Fii = LFMIt.matrix[h][i][i];
      for(a=0; a < (virtpi[h^L_irr] - openpi[h^L_irr]); a++) {
        Faa = LFAEt.matrix[h^L_irr][a][a];
        dIA.matrix[h][i][a] = 1.0/(Fii - Faa + params.cceom_energy[L_irr][e_index]);
      }
    }
  }
  dpd_file2_mat_wrt(&dIA);
  dpd_file2_mat_close(&dIA);
  dpd_file2_close(&dIA);

  /* Beta one-electron denominator */
  dpd_file2_init(&dia, CC_OEI, L_irr, 0, 1, "dia");
  dpd_file2_mat_init(&dia);
  for(h=0; h < nirreps; h++) {
    for(i=0; i < (occpi[h] - openpi[h]); i++) {
      Fii = LFmit.matrix[h][i][i];
      for(a=0; a < virtpi[h^L_irr]; a++) {
        Faa = LFaet.matrix[h^L_irr][a][a];
        dia.matrix[h][i][a] = 1.0/(Fii - Faa + params.cceom_energy[L_irr][e_index]);
      }
    }
  }
  dpd_file2_mat_wrt(&dia);
  dpd_file2_mat_close(&dia);
  dpd_file2_close(&dia);

  /* Alpha-alpha two-electron denominator */
  dpd_file4_init(&dIJAB, CC_DENOM, L_irr, 1, 6, "dIJAB");

  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&dIJAB, h);
      /* Loop over the rows */
      for(ij=0; ij < dIJAB.params->rowtot[h]; ij++) {
          i = dIJAB.params->roworb[h][ij][0];
          j = dIJAB.params->roworb[h][ij][1];
          isym = dIJAB.params->psym[i];
          jsym = dIJAB.params->qsym[j];

          /* Convert to relative orbital index */
          I = i - occ_off[isym];
          J = j - occ_off[jsym];

          Fii = LFMIt.matrix[isym][I][I];
          Fjj = LFMIt.matrix[jsym][J][J];

          /* Loop over the columns */
          for(ab=0; ab < dIJAB.params->coltot[h^L_irr]; ab++) {
              a = dIJAB.params->colorb[h^L_irr][ab][0];
              b = dIJAB.params->colorb[h^L_irr][ab][1];
              asym = dIJAB.params->rsym[a];
              bsym = dIJAB.params->ssym[b];

              /* Convert to relative orbital index */
              A = a - vir_off[asym];
              B = b - vir_off[bsym];

              Faa = LFAEt.matrix[asym][A][A];
              Fbb = LFAEt.matrix[bsym][B][B];

              dIJAB.matrix[h][ij][ab] =
                ((A >= (virtpi[asym] - openpi[asym])) ||
                 (B >= (virtpi[bsym] - openpi[bsym])) ?
                 0.0 : 1.0/(Fii + Fjj - Faa - Fbb
                   + params.cceom_energy[L_irr][e_index]));
          }
      }
    dpd_file4_mat_irrep_wrt(&dIJAB, h);
    dpd_file4_mat_irrep_close(&dIJAB, h);
  }
  dpd_file4_close(&dIJAB);

  /* Beta-beta two-electron denominator */
  dpd_file4_init(&dijab, CC_DENOM, L_irr, 1, 6, "dijab");

  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&dijab, h);
    /* Loop over the rows */
    for(ij=0; ij < dijab.params->rowtot[h]; ij++) {
          i = dijab.params->roworb[h][ij][0];
          j = dijab.params->roworb[h][ij][1];
          isym = dijab.params->psym[i];
          jsym = dijab.params->qsym[j];

          /* Convert to relative orbital index */
          I = i - occ_off[isym];
          J = j - occ_off[jsym];

          Fii = LFmit.matrix[isym][I][I];
          Fjj = LFmit.matrix[jsym][J][J];

          /* Loop over the columns */
          for(ab=0; ab < dijab.params->coltot[h^L_irr]; ab++) {
              a = dijab.params->colorb[h^L_irr][ab][0];
              b = dijab.params->colorb[h^L_irr][ab][1];
              asym = dijab.params->rsym[a];
              bsym = dijab.params->ssym[b];

              /* Convert to relative orbital index */
              A = a - vir_off[asym];
              B = b - vir_off[bsym];

              Faa = LFaet.matrix[asym][A][A];
              Fbb = LFaet.matrix[bsym][B][B];

              dijab.matrix[h][ij][ab] =
                ((I >= (occpi[isym] - openpi[isym])) ||
                 (J >= (occpi[jsym] - openpi[jsym])) ?
                 0.0 : 1.0/(Fii + Fjj - Faa - Fbb
                   + params.cceom_energy[L_irr][e_index]));
      }
    }
    dpd_file4_mat_irrep_wrt(&dijab, h);
    dpd_file4_mat_irrep_close(&dijab, h);
  }
  dpd_file4_close(&dijab);


  /* Alpha-beta two-electron denominator */
  dpd_file4_init(&dIjAb, CC_DENOM, L_irr, 0, 5, "dIjAb");

  for(h=0; h < nirreps; h++) {
    dpd_file4_mat_irrep_init(&dIjAb, h);
    /* Loop over the rows */
    for(ij=0; ij < dIjAb.params->rowtot[h]; ij++) {
          i = dIjAb.params->roworb[h][ij][0];
          j = dIjAb.params->roworb[h][ij][1];
          isym = dIjAb.params->psym[i];
          jsym = dIjAb.params->qsym[j];

          /* Convert to relative orbital index */
          I = i - occ_off[isym];
          J = j - occ_off[jsym];
          Fii = LFMIt.matrix[isym][I][I];
          Fjj = LFmit.matrix[jsym][J][J];

          /* Loop over the columns */
          for(ab=0; ab < dIjAb.params->coltot[h^L_irr]; ab++) {
              a = dIjAb.params->colorb[h^L_irr][ab][0];
              b = dIjAb.params->colorb[h^L_irr][ab][1];
              asym = dIjAb.params->rsym[a];
              bsym = dIjAb.params->ssym[b];

              /* Convert to relative orbital index */
              A = a - vir_off[asym];
              B = b - vir_off[bsym];

              Faa = LFAEt.matrix[asym][A][A];
              Fbb = LFaet.matrix[bsym][B][B];

              dIjAb.matrix[h][ij][ab] =
                ((A >= (virtpi[asym] - openpi[asym])) ||
                 (J >= (occpi[jsym] - openpi[jsym])) ?
                 0.0 : 1.0/(Fii + Fjj - Faa - Fbb
                   + params.cceom_energy[L_irr][e_index]));
            }
        }
    dpd_file4_mat_irrep_wrt(&dIjAb, h);
    dpd_file4_mat_irrep_close(&dIjAb, h);
  }
  dpd_file4_close(&dIjAb);

  dpd_file2_mat_close(&LFMIt);
  dpd_file2_mat_close(&LFmit);
  dpd_file2_mat_close(&LFAEt);
  dpd_file2_mat_close(&LFaet);
  dpd_file2_close(&LFMIt);
  dpd_file2_close(&LFmit);
  dpd_file2_close(&LFAEt);
  dpd_file2_close(&LFaet);

  return;
}



