#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Fae_build(void)
{
  int h,a,e;
  struct oe_dpdfile tIA, tia;
  struct oe_dpdfile FME, Fme;
  struct oe_dpdfile fAB, fab, fIA, fia;
  struct oe_dpdfile FAE, Fae;
  struct oe_dpdfile FAEt, Faet;
  struct dpdbuf F_anti, F, D_anti, D;
  struct dpdbuf tautIJAB, tautijab, tautIjAb;

  dpd_oe_file_init(&fAB, CC_OEI, 1, 1, "fAB", 0, outfile);
  dpd_oe_copy(&fAB, CC_OEI, "FAE", 0, outfile);
  dpd_oe_file_close(&fAB);

  dpd_oe_file_init(&fab, CC_OEI, 1, 1, "fab", 0, outfile);
  dpd_oe_copy(&fab, CC_OEI, "Fae", 0, outfile);
  dpd_oe_file_close(&fab);

  dpd_oe_file_init(&FAE, CC_OEI, 1, 1, "FAE", 0, outfile);
  dpd_oe_file_init(&Fae, CC_OEI, 1, 1, "Fae", 0, outfile);
  
  dpd_oe_file_mat_init(&FAE);
  dpd_oe_file_mat_rd(&FAE, 0, outfile);
  dpd_oe_file_mat_init(&Fae);
  dpd_oe_file_mat_rd(&Fae, 0, outfile);
  
  for(h=0; h < moinfo.nirreps; h++) {
      
      for(a=0; a < FAE.params->rowtot[h]; a++) 
	  for(e=0; e < FAE.params->coltot[h]; e++) 
	      FAE.matrix[h][a][e] *= (1 - (a==e));

      for(a=0; a < Fae.params->rowtot[h]; a++) 
	  for(e=0; e < Fae.params->coltot[h]; e++)
	      Fae.matrix[h][a][e] *= (1 - (a==e));
      
    }

  dpd_oe_file_mat_wrt(&FAE, 0, outfile);
  dpd_oe_file_mat_close(&FAE);
  dpd_oe_file_mat_wrt(&Fae, 0, outfile);
  dpd_oe_file_mat_close(&Fae);

  dpd_oe_file_init(&fIA, CC_OEI, 0, 1, "fIA", 0, outfile);
  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract111(&tIA, &fIA, &FAE, 1, 1, -0.5, 1, 0, outfile);
  dpd_oe_file_close(&tIA);
  dpd_oe_file_close(&fIA);

  dpd_oe_file_init(&fia, CC_OEI, 0, 1, "fia", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract111(&tia, &fia, &Fae, 1, 1, -0.5, 1, 0, outfile);
  dpd_oe_file_close(&tia);
  dpd_oe_file_close(&fia);

  dpd_buf_init(&F_anti, CC_FINTS, 10, 5, 10, 5, 1, "F <ia|bc>", 0, outfile);
  dpd_buf_init(&F, CC_FINTS, 10, 5, 10, 5, 0,"F <ia|bc>", 0, outfile);
  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);

  dpd_dot13(&tIA, &F_anti, &FAE, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_dot13(&tia, &F, &FAE, 0, 0, 1.0, 1.0, 0, outfile);

  dpd_dot13(&tia, &F_anti, &Fae, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_dot13(&tIA, &F, &Fae, 0, 0, 1.0, 1.0, 0, outfile);

  dpd_oe_file_close(&tIA);
  dpd_oe_file_close(&tia);
  dpd_buf_close(&F_anti);
  dpd_buf_close(&F);

  dpd_buf_init(&D_anti, CC_DINTS, 2, 5, 2, 5, 0,
	       "D <ij||ab> (i>j,ab)", 0, outfile);

  dpd_buf_init(&tautIJAB, CC_TAMPS, 2, 5, 2, 7, 0, "tautIJAB", 0, outfile);
  dpd_buf_init(&tautijab, CC_TAMPS, 2, 5, 2, 7, 0, "tautijab", 0, outfile);

  dpd_contract122(&tautIJAB, &D_anti, &FAE, 2, 2, -1, 1, 0, outfile);
  dpd_contract122(&tautijab, &D_anti, &Fae, 2, 2, -1, 1, 0, outfile);

  dpd_buf_close(&D_anti);
  dpd_buf_close(&tautIJAB);
  dpd_buf_close(&tautijab);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_buf_init(&tautIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tautIjAb", 0, outfile);

  dpd_contract122(&tautIjAb, &D, &Fae, 3, 3, -1, 1, 0, outfile);
  dpd_contract122(&tautIjAb, &D, &FAE, 2, 2, -1, 1, 0, outfile);

  dpd_buf_close(&D);
  dpd_buf_close(&tautIjAb);


  /* Build the tilde intermediates */
  dpd_oe_copy(&FAE, CC_OEI, "FAEt", 0, outfile);
  dpd_oe_copy(&Fae, CC_OEI, "Faet", 0, outfile);

  dpd_oe_file_close(&FAE);
  dpd_oe_file_close(&Fae);

  dpd_oe_file_init(&FAEt, CC_OEI, 1, 1, "FAEt", 0, outfile);
  dpd_oe_file_init(&Faet, CC_OEI, 1, 1, "Faet", 0, outfile);

  dpd_oe_file_init(&tIA, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&FME, CC_OEI, 0, 1, "FME", 0, outfile);
  dpd_contract111(&tIA, &FME, &FAEt, 1, 1, -0.5, 1, 0, outfile);
  dpd_oe_file_close(&tIA);
  dpd_oe_file_close(&FME);
  
  dpd_oe_file_init(&tia, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_init(&Fme, CC_OEI, 0, 1, "Fme", 0, outfile);
  dpd_contract111(&tia, &Fme, &Faet, 1, 1, -0.5, 1, 0, outfile);
  dpd_oe_file_close(&tia);
  dpd_oe_file_close(&Fme);

  dpd_oe_file_close(&FAEt);
  dpd_oe_file_close(&Faet); 
}
