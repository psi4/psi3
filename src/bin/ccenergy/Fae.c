#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Fae_build(void)
{
  int h,a,e;
  dpdfile2 tIA, tia;
  dpdfile2 FME, Fme;
  dpdfile2 fAB, fab, fIA, fia;
  dpdfile2 FAE, Fae;
  dpdfile2 FAEt, Faet;
  dpdbuf4 F_anti, F, D_anti, D;
  dpdbuf4 tautIJAB, tautijab, tautIjAb, taut;

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_copy(&fAB, CC_OEI, "FAE");
    dpd_file2_close(&fAB);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_copy(&fAB, CC_OEI, "FAE");
    dpd_file2_close(&fAB);

    dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");
    dpd_file2_copy(&fab, CC_OEI, "Fae");
    dpd_file2_close(&fab);
  }

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
  
    dpd_file2_mat_init(&FAE);
    dpd_file2_mat_rd(&FAE);
  
    for(h=0; h < moinfo.nirreps; h++) {
      for(a=0; a < FAE.params->rowtot[h]; a++) 
	FAE.matrix[h][a][a] = 0;
    }

    dpd_file2_mat_wrt(&FAE);
    dpd_file2_mat_close(&FAE);
    dpd_file2_close(&FAE);
  }

  else if(params.ref == 1) { /** ROHF **/
    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_init(&Fae, CC_OEI, 0, 1, 1, "Fae");
  
    dpd_file2_mat_init(&FAE);
    dpd_file2_mat_rd(&FAE);
    dpd_file2_mat_init(&Fae);
    dpd_file2_mat_rd(&Fae);
  
    for(h=0; h < moinfo.nirreps; h++) {
      
      for(a=0; a < FAE.params->rowtot[h]; a++) 
	for(e=0; e < FAE.params->coltot[h]; e++) 
	  FAE.matrix[h][a][e] *= (1 - (a==e));

      for(a=0; a < Fae.params->rowtot[h]; a++) 
	for(e=0; e < Fae.params->coltot[h]; e++)
	  Fae.matrix[h][a][e] *= (1 - (a==e));
      
    }

    dpd_file2_mat_wrt(&FAE);
    dpd_file2_mat_close(&FAE);
    dpd_file2_mat_wrt(&Fae);
    dpd_file2_mat_close(&Fae);

    dpd_file2_close(&FAE);
    dpd_file2_close(&Fae);
  }

  if(params.ref == 0) { /** RHF **/
    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&fIA);

    dpd_buf4_init(&F_anti, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_dot13(&tIA, &F_anti, &FAE, 0, 0, 1.0, 1.0);
    dpd_dot13(&tIA, &F, &FAE, 0, 0, 1.0, 1.0);
    dpd_file2_close(&tIA);
    dpd_buf4_close(&F_anti);
    dpd_buf4_close(&F);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");
    dpd_contract442(&tautIjAb, &D, &FAE, 2, 2, -1, 1);
    dpd_buf4_close(&D);
    dpd_buf4_close(&tautIjAb);

    /* Build the tilde intermediates */
    dpd_file2_copy(&FAE, CC_OEI, "FAEt");
    dpd_file2_close(&FAE);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&FME);
  
    dpd_file2_close(&FAEt);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_init(&Fae, CC_OEI, 0, 1, 1, "Fae");

    dpd_file2_init(&fIA, CC_OEI, 0, 0, 1, "fIA");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&fIA);

    dpd_file2_init(&fia, CC_OEI, 0, 0, 1, "fia");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&fia);

    dpd_buf4_init(&F_anti, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");

    dpd_dot13(&tIA, &F_anti, &FAE, 0, 0, 1.0, 1.0);
    dpd_dot13(&tia, &F, &FAE, 0, 0, 1.0, 1.0);

    dpd_dot13(&tia, &F_anti, &Fae, 0, 0, 1.0, 1.0);
    dpd_dot13(&tIA, &F, &Fae, 0, 0, 1.0, 1.0);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);
    dpd_buf4_close(&F_anti);
    dpd_buf4_close(&F);

    dpd_buf4_init(&D_anti, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");

    dpd_buf4_init(&tautIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
    dpd_buf4_init(&tautijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautijab");

    dpd_contract442(&tautIJAB, &D_anti, &FAE, 2, 2, -1, 1);
    dpd_contract442(&tautijab, &D_anti, &Fae, 2, 2, -1, 1);

    dpd_buf4_close(&D_anti);
    dpd_buf4_close(&tautIJAB);
    dpd_buf4_close(&tautijab);

    dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_init(&tautIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

    dpd_contract442(&tautIjAb, &D, &Fae, 3, 3, -1, 1);
    dpd_contract442(&tautIjAb, &D, &FAE, 2, 2, -1, 1);

    dpd_buf4_close(&D);
    dpd_buf4_close(&tautIjAb);


    /* Build the tilde intermediates */
    dpd_file2_copy(&FAE, CC_OEI, "FAEt");
    dpd_file2_copy(&Fae, CC_OEI, "Faet");

    dpd_file2_close(&FAE);
    dpd_file2_close(&Fae);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_init(&Faet, CC_OEI, 0, 1, 1, "Faet");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&FME);
  
    dpd_file2_init(&tia, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "Fme");
    dpd_contract222(&tia, &Fme, &Faet, 1, 1, -0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&Fme);

    dpd_file2_close(&FAEt);
    dpd_file2_close(&Faet); 
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_init(&Fae, CC_OEI, 0, 3, 3, "Fae");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");

    dpd_buf4_init(&F, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
    dpd_dot13(&tIA, &F, &FAE, 0, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
    dpd_dot13(&tia, &F, &FAE, 0, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
    dpd_dot13(&tia, &F, &Fae, 0, 0, 1, 0);
    dpd_buf4_close(&F);

    dpd_buf4_init(&F, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
    dpd_dot13(&tIA, &F, &Fae, 0, 0, 1, 1);
    dpd_buf4_close(&F);

    dpd_file2_close(&tIA);
    dpd_file2_close(&tia);

    dpd_buf4_init(&D, CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    dpd_buf4_init(&taut, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
    dpd_contract442(&taut, &D, &FAE, 2, 2, -1, 1);
    dpd_buf4_close(&taut);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_init(&taut, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    dpd_contract442(&taut, &D, &FAE, 2, 2, -1, 1);
    dpd_buf4_close(&taut);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    dpd_buf4_init(&taut, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tautijab");
    dpd_contract442(&taut, &D, &Fae, 2, 2, -1, 1);
    dpd_buf4_close(&taut);
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_init(&taut, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
    dpd_contract442(&taut, &D, &Fae, 3, 3, -1, 1);
    dpd_buf4_close(&taut);
    dpd_buf4_close(&D);

    /* Build the tilde intermediates */
    dpd_file2_copy(&FAE, CC_OEI, "FAEt");
    dpd_file2_copy(&Fae, CC_OEI, "Faet");

    dpd_file2_close(&FAE);
    dpd_file2_close(&Fae);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_init(&Faet, CC_OEI, 0, 3, 3, "Faet");

    dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
    dpd_contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
    dpd_file2_close(&tIA);
    dpd_file2_close(&FME);
  
    dpd_file2_init(&tia, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
    dpd_contract222(&tia, &Fme, &Faet, 1, 1, -0.5, 1);
    dpd_file2_close(&tia);
    dpd_file2_close(&Fme);


    dpd_file2_close(&FAEt);
    dpd_file2_close(&Faet); 
  }
}
