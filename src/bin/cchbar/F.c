#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

/* FME, Fme, FMI, Fmi, FAE, Fae become real matrix elements */
/* FMEt,Fmet,FMIt,Fmit,FAEt,Faet have zeros on diagonal */

void F_build(void) {
  int h,i,e,a;
  dpdfile2 Faet, FAEt, Fmit, FMIt;
  dpdfile2 Fae, FAE, FMI, Fmi;
  dpdfile2 fab, fAB, fmi, fMI;
  dpdfile2 FME, Fme;

  if(params.ref == 0 ||  params.ref == 1) { /** RHF or ROHF **/

    /* copy all Ft's to F's */
    dpd_file2_init(&Faet, CC_OEI, 0, 1, 1, "Faet");
    dpd_file2_copy(&Faet, CC_OEI, "Fae");
    dpd_file2_close(&Faet);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_copy(&FAEt, CC_OEI, "FAE");
    dpd_file2_close(&FAEt);

    dpd_file2_init(&Fmit, CC_OEI, 0, 0, 0, "Fmit");
    dpd_file2_copy(&Fmit, CC_OEI, "Fmi");
    dpd_file2_close(&Fmit);
 
    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_copy(&FMIt, CC_OEI, "FMI");
    dpd_file2_close(&FMIt);


    /* Add fock matrix elements to F's to form real matrix elements */
    dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");
    dpd_file2_mat_init(&fab);
    dpd_file2_mat_rd(&fab);

    dpd_file2_init(&Fae, CC_OEI, 0, 1, 1, "Fae");
    dpd_file2_mat_init(&Fae);
    dpd_file2_mat_rd(&Fae);

    for(h=0; h < moinfo.nirreps; h++) {
      for(a=0; a < Fae.params->rowtot[h]; a++)
	Fae.matrix[h][a][a] += fab.matrix[h][a][a];
    }

    dpd_file2_mat_wrt(&Fae);
    dpd_file2_mat_close(&fab);
    dpd_file2_mat_close(&Fae);
    dpd_file2_close(&fab);
    dpd_file2_close(&Fae);


    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_mat_init(&fAB);
    dpd_file2_mat_rd(&fAB);

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_mat_init(&FAE);
    dpd_file2_mat_rd(&FAE);

    for(h=0; h < moinfo.nirreps; h++) {
      for(a=0; a < FAE.params->rowtot[h]; a++)
	FAE.matrix[h][a][a] += fAB.matrix[h][a][a];
    }

    dpd_file2_mat_wrt(&FAE);
    dpd_file2_mat_close(&fAB);
    dpd_file2_mat_close(&FAE);
    dpd_file2_close(&fAB);
    dpd_file2_close(&FAE);


    dpd_file2_init(&fmi, CC_OEI, 0, 0, 0, "fij");
    dpd_file2_mat_init(&fmi);
    dpd_file2_mat_rd(&fmi);

    dpd_file2_init(&Fmi, CC_OEI, 0, 0, 0, "Fmi");
    dpd_file2_mat_init(&Fmi);
    dpd_file2_mat_rd(&Fmi);

    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < Fmi.params->rowtot[h]; i++)
	Fmi.matrix[h][i][i] += fmi.matrix[h][i][i];
    }

    dpd_file2_mat_wrt(&Fmi);
    dpd_file2_mat_close(&fmi);
    dpd_file2_mat_close(&Fmi);
    dpd_file2_close(&fmi);
    dpd_file2_close(&Fmi);


    dpd_file2_init(&fMI, CC_OEI, 0, 0, 0, "fIJ");
    dpd_file2_mat_init(&fMI);
    dpd_file2_mat_rd(&fMI);

    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_mat_init(&FMI);
    dpd_file2_mat_rd(&FMI);

    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < FMI.params->rowtot[h]; i++)
	FMI.matrix[h][i][i] += fMI.matrix[h][i][i];
    }

    dpd_file2_mat_wrt(&FMI);
    dpd_file2_mat_close(&fMI);
    dpd_file2_mat_close(&FMI);
    dpd_file2_close(&fMI);
    dpd_file2_close(&FMI);



    /* remove diagonal elements from Ft's */
    dpd_file2_init(&Faet, CC_OEI, 0, 1, 1, "Faet");
    dpd_file2_mat_init(&Faet);
    dpd_file2_mat_rd(&Faet);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < Faet.params->rowtot[h]; a++)
	Faet.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&Faet);
    dpd_file2_mat_close(&Faet);
    dpd_file2_close(&Faet);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_mat_init(&FAEt);
    dpd_file2_mat_rd(&FAEt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FAEt.params->rowtot[h]; a++)
	FAEt.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&FAEt);
    dpd_file2_mat_close(&FAEt);
    dpd_file2_close(&FAEt);

    dpd_file2_init(&Fmit, CC_OEI, 0, 0, 0, "Fmit");
    dpd_file2_mat_init(&Fmit);
    dpd_file2_mat_rd(&Fmit);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < Fmit.params->rowtot[h]; a++)
	Fmit.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&Fmit);
    dpd_file2_mat_close(&Fmit);
    dpd_file2_close(&Fmit);

    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_mat_init(&FMIt);
    dpd_file2_mat_rd(&FMIt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FMIt.params->rowtot[h]; a++)
	FMIt.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&FMIt);
    dpd_file2_mat_close(&FMIt);
    dpd_file2_close(&FMIt);

  } /** RHF or ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* copy all Ft's to F's */
    dpd_file2_init(&Faet, CC_OEI, 0, 3, 3, "Faet");
    dpd_file2_copy(&Faet, CC_OEI, "Fae");
    dpd_file2_close(&Faet);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_copy(&FAEt, CC_OEI, "FAE");
    dpd_file2_close(&FAEt);

    dpd_file2_init(&Fmit, CC_OEI, 0, 2, 2, "Fmit");
    dpd_file2_copy(&Fmit, CC_OEI, "Fmi");
    dpd_file2_close(&Fmit);
 
    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_copy(&FMIt, CC_OEI, "FMI");
    dpd_file2_close(&FMIt);


    /* Add fock matrix elements to F's to form real matrix elements */
    dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");
    dpd_file2_mat_init(&fab);
    dpd_file2_mat_rd(&fab);

    dpd_file2_init(&Fae, CC_OEI, 0, 3, 3, "Fae");
    dpd_file2_mat_init(&Fae);
    dpd_file2_mat_rd(&Fae);

    for(h=0; h < moinfo.nirreps; h++) {
      for(a=0; a < Fae.params->rowtot[h]; a++)
	Fae.matrix[h][a][a] += fab.matrix[h][a][a];
    }

    dpd_file2_mat_wrt(&Fae);
    dpd_file2_mat_close(&fab);
    dpd_file2_mat_close(&Fae);
    dpd_file2_close(&fab);
    dpd_file2_close(&Fae);


    dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
    dpd_file2_mat_init(&fAB);
    dpd_file2_mat_rd(&fAB);

    dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
    dpd_file2_mat_init(&FAE);
    dpd_file2_mat_rd(&FAE);

    for(h=0; h < moinfo.nirreps; h++) {
      for(a=0; a < FAE.params->rowtot[h]; a++)
	FAE.matrix[h][a][a] += fAB.matrix[h][a][a];
    }

    dpd_file2_mat_wrt(&FAE);
    dpd_file2_mat_close(&fAB);
    dpd_file2_mat_close(&FAE);
    dpd_file2_close(&fAB);
    dpd_file2_close(&FAE);


    dpd_file2_init(&fmi, CC_OEI, 0, 2, 2, "fij");
    dpd_file2_mat_init(&fmi);
    dpd_file2_mat_rd(&fmi);

    dpd_file2_init(&Fmi, CC_OEI, 0, 2, 2, "Fmi");
    dpd_file2_mat_init(&Fmi);
    dpd_file2_mat_rd(&Fmi);

    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < Fmi.params->rowtot[h]; i++)
	Fmi.matrix[h][i][i] += fmi.matrix[h][i][i];
    }

    dpd_file2_mat_wrt(&Fmi);
    dpd_file2_mat_close(&fmi);
    dpd_file2_mat_close(&Fmi);
    dpd_file2_close(&fmi);
    dpd_file2_close(&Fmi);


    dpd_file2_init(&fMI, CC_OEI, 0, 0, 0, "fIJ");
    dpd_file2_mat_init(&fMI);
    dpd_file2_mat_rd(&fMI);

    dpd_file2_init(&FMI, CC_OEI, 0, 0, 0, "FMI");
    dpd_file2_mat_init(&FMI);
    dpd_file2_mat_rd(&FMI);

    for(h=0; h < moinfo.nirreps; h++) {
      for(i=0; i < FMI.params->rowtot[h]; i++)
	FMI.matrix[h][i][i] += fMI.matrix[h][i][i];
    }

    dpd_file2_mat_wrt(&FMI);
    dpd_file2_mat_close(&fMI);
    dpd_file2_mat_close(&FMI);
    dpd_file2_close(&fMI);
    dpd_file2_close(&FMI);


    /* remove diagonal elements from Ft's */
    dpd_file2_init(&Faet, CC_OEI, 0, 3, 3, "Faet");
    dpd_file2_mat_init(&Faet);
    dpd_file2_mat_rd(&Faet);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < Faet.params->rowtot[h]; a++)
	Faet.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&Faet);
    dpd_file2_mat_close(&Faet);
    dpd_file2_close(&Faet);

    dpd_file2_init(&FAEt, CC_OEI, 0, 1, 1, "FAEt");
    dpd_file2_mat_init(&FAEt);
    dpd_file2_mat_rd(&FAEt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FAEt.params->rowtot[h]; a++)
	FAEt.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&FAEt);
    dpd_file2_mat_close(&FAEt);
    dpd_file2_close(&FAEt);

    dpd_file2_init(&Fmit, CC_OEI, 0, 2, 2, "Fmit");
    dpd_file2_mat_init(&Fmit);
    dpd_file2_mat_rd(&Fmit);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < Fmit.params->rowtot[h]; a++)
	Fmit.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&Fmit);
    dpd_file2_mat_close(&Fmit);
    dpd_file2_close(&Fmit);

    dpd_file2_init(&FMIt, CC_OEI, 0, 0, 0, "FMIt");
    dpd_file2_mat_init(&FMIt);
    dpd_file2_mat_rd(&FMIt);

    for(h=0; h < moinfo.nirreps; h++)
      for(a=0; a < FMIt.params->rowtot[h]; a++)
	FMIt.matrix[h][a][a] = 0.0;

    dpd_file2_mat_wrt(&FMIt);
    dpd_file2_mat_close(&FMIt);
    dpd_file2_close(&FMIt);
  } /** UHF **/

  return;

}
