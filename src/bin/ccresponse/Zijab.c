#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void build_ZIjAb(char *cart_x, int irrep_x, double omega_x, char *cart_y, int irrep_y, double omega_y)
{
  int h, ij, ef, i, j, e, f, I, J, E, F;
  int Isym, Jsym, Esym, Fsym;
  int nirreps;
  dpdbuf4 Z1;
  dpdfile2 X1, Y1;
  char lbl[32];

  nirreps = moinfo.nirreps;

  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_y, omega_y);
  dpd_file2_init(&Y1, CC_OEI, irrep_y, 0, 1, lbl);
  dpd_file2_mat_init(&Y1);
  dpd_file2_mat_rd(&Y1);

  sprintf(lbl, "X_%1s_IA (%5.3f)", cart_x, omega_x);
  dpd_file2_init(&X1, CC_OEI, irrep_x, 0, 1, lbl);
  dpd_file2_mat_init(&X1);
  dpd_file2_mat_rd(&X1);

  dpd_buf4_init(&Z1, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab) anti");

  for(h=0; h< nirreps; h++) {

    dpd_buf4_mat_irrep_init(&Z1, h);

    for(ij=0; ij<Z1.params->rowtot[h]; ij++) {
      i = Z1.params->roworb[h][ij][0];
      j = Z1.params->roworb[h][ij][1];
      I = X1.params->rowidx[i];
      J = Y1.params->rowidx[j];
      Isym = X1.params->psym[i];
      Jsym = Y1.params->psym[j];
      for(ef=0; ef < Z1.params->coltot[h]; ef++) {
	e = Z1.params->colorb[h][ef][0];
	f = Z1.params->colorb[h][ef][1];
	E = X1.params->colidx[e];
	F = Y1.params->colidx[f];
	Esym = X1.params->qsym[e];
	Fsym = Y1.params->qsym[f];

	if(((Isym^Esym)==irrep_x) && ((Jsym^Fsym)==irrep_y))
	  Z1.matrix[h][ij][ef] +=
	    (X1.matrix[Isym][I][E] * Y1.matrix[Jsym][J][F]) + 
	    (Y1.matrix[Isym][I][E] * X1.matrix[Jsym][J][F]);
      }
    }
    dpd_buf4_mat_irrep_wrt(&Z1, h);
    dpd_buf4_mat_irrep_close(&Z1, h);
  }
  dpd_buf4_close(&Z1);

  dpd_file2_mat_close(&X1);
  dpd_file2_close(&X1);
  dpd_file2_mat_close(&Y1);
  dpd_file2_close(&Y1);

}
