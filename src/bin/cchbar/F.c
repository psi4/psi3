#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

/* FME, Fme, FMI, Fmi, FAE, Fae become real matrix elements */
/* FMEt,Fmet,FMIt,Fmit,FAEt,Faet have zeros on diagonal */

void F_build(void) {
  int h,i,e,a;
  struct oe_dpdfile Faet, FAEt, Fmit, FMIt;
  struct oe_dpdfile Fae, FAE, FMI, Fmi;
  struct oe_dpdfile fab, fAB, fmi, fMI;
  struct oe_dpdfile FME, Fme;
/*
  struct oe_dpdfile LFAEt, LFaet, LFMIt, LFmit;
  struct oe_dpdfile LFAEt2, LFaet2, LFMIt2, LFmit2;
*/

 /* copy all Ft's to F's */
  dpd_oe_file_init(&Faet, CC_OEI, 1, 1, "Faet", 0, outfile);
  dpd_oe_copy(&Faet, CC_OEI, "Fae",0,outfile);
  dpd_oe_file_close(&Faet);

  dpd_oe_file_init(&FAEt, CC_OEI, 1, 1, "FAEt", 0, outfile);
  dpd_oe_copy(&FAEt, CC_OEI, "FAE",0,outfile);
  dpd_oe_file_close(&FAEt);

  dpd_oe_file_init(&Fmit, CC_OEI, 0, 0, "Fmit", 0, outfile);
  dpd_oe_copy(&Fmit, CC_OEI, "Fmi",0,outfile);
  dpd_oe_file_close(&Fmit);
 
  dpd_oe_file_init(&FMIt, CC_OEI, 0, 0, "FMIt", 0, outfile);
  dpd_oe_copy(&FMIt, CC_OEI, "FMI",0,outfile);
  dpd_oe_file_close(&FMIt);


 /* Add fock matrix elements to F's to form real matrix elements */
  dpd_oe_file_init(&fab, CC_OEI, 1, 1, "fab", 0, outfile);
  dpd_oe_file_mat_init(&fab);
  dpd_oe_file_mat_rd(&fab, 0, outfile);

  dpd_oe_file_init(&Fae, CC_OEI, 1, 1, "Fae", 0, outfile);
  dpd_oe_file_mat_init(&Fae);
  dpd_oe_file_mat_rd(&Fae, 0, outfile);

  for(h=0; h < moinfo.nirreps; h++) {
    for(a=0; a < Fae.params->rowtot[h]; a++)
      Fae.matrix[h][a][a] += fab.matrix[h][a][a];
  }

  dpd_oe_file_mat_wrt(&Fae,0,outfile);
  dpd_oe_file_mat_close(&fab);
  dpd_oe_file_mat_close(&Fae);
  dpd_oe_file_close(&fab);
  dpd_oe_file_close(&Fae);


  dpd_oe_file_init(&fAB, CC_OEI, 1, 1, "fAB", 0, outfile);
  dpd_oe_file_mat_init(&fAB);
  dpd_oe_file_mat_rd(&fAB, 0, outfile);

  dpd_oe_file_init(&FAE, CC_OEI, 1, 1, "FAE", 0, outfile);
  dpd_oe_file_mat_init(&FAE);
  dpd_oe_file_mat_rd(&FAE, 0, outfile);

  for(h=0; h < moinfo.nirreps; h++) {
    for(a=0; a < FAE.params->rowtot[h]; a++)
      FAE.matrix[h][a][a] += fAB.matrix[h][a][a];
  }

  dpd_oe_file_mat_wrt(&FAE,0,outfile);
  dpd_oe_file_mat_close(&fAB);
  dpd_oe_file_mat_close(&FAE);
  dpd_oe_file_close(&fAB);
  dpd_oe_file_close(&FAE);


  dpd_oe_file_init(&fmi, CC_OEI, 0, 0, "fij", 0, outfile);
  dpd_oe_file_mat_init(&fmi);
  dpd_oe_file_mat_rd(&fmi, 0, outfile);

  dpd_oe_file_init(&Fmi, CC_OEI, 0, 0, "Fmi", 0, outfile);
  dpd_oe_file_mat_init(&Fmi);
  dpd_oe_file_mat_rd(&Fmi, 0, outfile);

  for(h=0; h < moinfo.nirreps; h++) {
    for(i=0; i < Fmi.params->rowtot[h]; i++)
      Fmi.matrix[h][i][i] += fmi.matrix[h][i][i];
  }

  dpd_oe_file_mat_wrt(&Fmi,0,outfile);
  dpd_oe_file_mat_close(&fmi);
  dpd_oe_file_mat_close(&Fmi);
  dpd_oe_file_close(&fmi);
  dpd_oe_file_close(&Fmi);


  dpd_oe_file_init(&fMI, CC_OEI, 0, 0, "fIJ", 0, outfile);
  dpd_oe_file_mat_init(&fMI);
  dpd_oe_file_mat_rd(&fMI, 0, outfile);

  dpd_oe_file_init(&FMI, CC_OEI, 0, 0, "FMI", 0, outfile);
  dpd_oe_file_mat_init(&FMI);
  dpd_oe_file_mat_rd(&FMI, 0, outfile);

  for(h=0; h < moinfo.nirreps; h++) {
    for(i=0; i < FMI.params->rowtot[h]; i++)
      FMI.matrix[h][i][i] += fMI.matrix[h][i][i];
  }

  dpd_oe_file_mat_wrt(&FMI,0,outfile);
  dpd_oe_file_mat_close(&fMI);
  dpd_oe_file_mat_close(&FMI);
  dpd_oe_file_close(&fMI);
  dpd_oe_file_close(&FMI);



  /* remove diagonal elements from Ft's */
  dpd_oe_file_init(&Faet, CC_OEI, 1, 1, "Faet", 0, outfile);
  dpd_oe_file_mat_init(&Faet);
  dpd_oe_file_mat_rd(&Faet,0,outfile);

  for(h=0; h < moinfo.nirreps; h++)
    for(a=0; a < Faet.params->rowtot[h]; a++)
      Faet.matrix[h][a][a] = 0.0;

  dpd_oe_file_mat_wrt(&Faet,0,outfile);
  dpd_oe_file_mat_close(&Faet);
  dpd_oe_file_close(&Faet);

  dpd_oe_file_init(&FAEt, CC_OEI, 1, 1, "FAEt", 0, outfile);
  dpd_oe_file_mat_init(&FAEt);
  dpd_oe_file_mat_rd(&FAEt,0,outfile);

  for(h=0; h < moinfo.nirreps; h++)
    for(a=0; a < FAEt.params->rowtot[h]; a++)
      FAEt.matrix[h][a][a] = 0.0;

  dpd_oe_file_mat_wrt(&FAEt,0,outfile);
  dpd_oe_file_mat_close(&FAEt);
  dpd_oe_file_close(&FAEt);

  dpd_oe_file_init(&Fmit, CC_OEI, 0, 0, "Fmit", 0, outfile);
  dpd_oe_file_mat_init(&Fmit);
  dpd_oe_file_mat_rd(&Fmit,0,outfile);

  for(h=0; h < moinfo.nirreps; h++)
    for(a=0; a < Fmit.params->rowtot[h]; a++)
      Fmit.matrix[h][a][a] = 0.0;

  dpd_oe_file_mat_wrt(&Fmit,0,outfile);
  dpd_oe_file_mat_close(&Fmit);
  dpd_oe_file_close(&Fmit);

  dpd_oe_file_init(&FMIt, CC_OEI, 0, 0, "FMIt", 0, outfile);
  dpd_oe_file_mat_init(&FMIt);
  dpd_oe_file_mat_rd(&FMIt,0,outfile);

  for(h=0; h < moinfo.nirreps; h++)
    for(a=0; a < FMIt.params->rowtot[h]; a++)
      FMIt.matrix[h][a][a] = 0.0;

  dpd_oe_file_mat_wrt(&FMIt,0,outfile);
  dpd_oe_file_mat_close(&FMIt);
  dpd_oe_file_close(&FMIt);

  return;

}
