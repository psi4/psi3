#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Fmi_build(void)
{
  int h,m,i;
  struct oe_dpdfile FMI, Fmi, FMIt, Fmit, fIJ, fij, fIA, fia;
  struct oe_dpdfile tIA, tia, FME, Fme;
  struct dpdbuf E_anti, E, D_anti, D;
  struct dpdbuf tautIJAB, tautijab, tautIjAb;

  dpd_oe_file_init(&fIJ, CC_OEI, 0, 0, "fIJ", 0, outfile);
  dpd_oe_copy(&fIJ, CC_OEI, "FMI", 0, outfile);
  dpd_oe_file_close(&fIJ);
  
  dpd_oe_file_init(&fij, CC_OEI, 0, 0, "fij", 0, outfile);
  dpd_oe_copy(&fij, CC_OEI, "Fmi", 0, outfile);
  dpd_oe_file_close(&fij);

  dpd_oe_file_init(&FMI, CC_OEI, 0, 0, "FMI", 0, outfile);
  dpd_oe_file_init(&Fmi, CC_OEI, 0, 0, "Fmi", 0, outfile);

  dpd_oe_file_mat_init(&FMI);
  dpd_oe_file_mat_rd(&FMI, 0, outfile);
  dpd_oe_file_mat_init(&Fmi);
  dpd_oe_file_mat_rd(&Fmi, 0, outfile);

  for(h=0; h < moinfo.nirreps; h++) {
      for(m=0; m < FMI.params->rowtot[h]; m++) 
	  for(i=0; i < FMI.params->coltot[h]; i++)
	      FMI.matrix[h][m][i] *= (1 - (m==i));

      for(m=0; m < Fmi.params->rowtot[h]; m++) 
	  for(i=0; i < Fmi.params->coltot[h]; i++) 
	      Fmi.matrix[h][m][i] *= (1 - (m==i));
    }

  dpd_oe_file_mat_wrt(&FMI, 0, outfile);
  dpd_oe_file_mat_close(&FMI);
  dpd_oe_file_mat_wrt(&Fmi, 0, outfile);
  dpd_oe_file_mat_close(&Fmi);

  dpd_oe_file_init(&fIA, CC_OEI, 0, 1, "fIA", 0, outfile);
  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract111(&fIA, &tIA, &FMI, 0, 0, 0.5, 1, 0, outfile);
  dpd_oe_file_close(&tIA);
  dpd_oe_file_close(&fIA);
  
  dpd_oe_file_init(&fia, CC_OEI, 0, 1, "fia", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract111(&fia, &tia, &Fmi, 0, 0, 0.5, 1, 0, outfile);
  dpd_oe_file_close(&tia);
  dpd_oe_file_close(&fia);
  
  dpd_buf_init(&E_anti, CC_EINTS, 11, 0, 11, 0, 1, "E <ai|jk>", 0, outfile);
  dpd_buf_init(&E, CC_EINTS, 11, 0, 11, 0, 0, "E <ai|jk>", 0, outfile);
  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);

  dpd_dot13(&tIA, &E_anti, &FMI, 1, 1, 1.0, 1.0, 0, outfile);
  dpd_dot13(&tia, &E, &FMI, 1, 1, 1.0, 1.0, 0, outfile);

  dpd_dot13(&tia, &E_anti, &Fmi, 1, 1, 1.0, 1.0, 0, outfile);
  dpd_dot13(&tIA, &E, &Fmi, 1, 1, 1.0, 1.0, 0, outfile);

  dpd_oe_file_close(&tIA);
  dpd_oe_file_close(&tia);
  dpd_buf_close(&E_anti);
  dpd_buf_close(&E);

  dpd_buf_init(&D_anti, CC_DINTS, 0, 7, 0, 7, 0,
	       "D <ij||ab> (ij,a>b)", 0, outfile);
  dpd_buf_init(&tautIJAB, CC_TAMPS, 0, 7, 2, 7, 0, "tautIJAB", 0, outfile);
  dpd_buf_init(&tautijab, CC_TAMPS, 0, 7, 2, 7, 0, "tautijab", 0, outfile);

  dpd_contract122(&D_anti, &tautIJAB, &FMI, 0, 0, 1, 1, 0, outfile);
  dpd_contract122(&D_anti, &tautijab, &Fmi, 0, 0, 1, 1, 0, outfile);

  dpd_buf_close(&tautIJAB);
  dpd_buf_close(&tautijab);
  dpd_buf_close(&D_anti);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&tautIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tautIjAb", 0, outfile);

  dpd_contract122(&D, &tautIjAb, &FMI, 0, 0, 1, 1, 0, outfile);
  dpd_contract122(&D, &tautIjAb, &Fmi, 1, 1, 1, 1, 0, outfile);

  dpd_buf_close(&tautIjAb);
  dpd_buf_close(&D);

  /* Build the tilde intermediate */
  dpd_oe_copy(&FMI, CC_OEI, "FMIt", 0, outfile);
  dpd_oe_copy(&Fmi, CC_OEI, "Fmit", 0, outfile);

  dpd_oe_file_close(&FMI);
  dpd_oe_file_close(&Fmi);

  dpd_oe_file_init(&FMIt, CC_OEI, 0, 0, "FMIt", 0, outfile);
  dpd_oe_file_init(&Fmit, CC_OEI, 0, 0, "Fmit", 0, outfile);

  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&FME, CC_OEI, 0, 1, "FME", 0, outfile);
  dpd_contract111(&FME, &tIA, &FMIt, 0, 0, 0.5, 1, 0, outfile);
  dpd_oe_file_close(&FME);
  dpd_oe_file_close(&tIA);

  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "Fme", 0, outfile);
  dpd_contract111(&Fme, &tia, &Fmit, 0, 0, 0.5, 1, 0, outfile);
  dpd_oe_file_close(&Fme);
  dpd_oe_file_close(&tia);

  dpd_oe_file_close(&FMIt);
  dpd_oe_file_close(&Fmit);
}
