/* cc3_t3z(): Builds the following intermediates, which are required
** for the T3-->L1 contributions in the CC3 model:
**
** Zifln = +1/2 t3(lmn,def) * <im||de>
**
** Zdfan = -1/2 t3(lmn,def) * <lm||ae>
**
** Translating indices in each expression to use t3(ijk,abc):
**
** Zmcik = +1/2 t3(ijk,abc) * <mj||ab>
**
** Zacek = -1/2 t3(ijk,abc) * <ij||eb>
**
** The current version is coded ROHF-like, but works only for RHF.
**
** -TDC, 7/04
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void cc3_t3z_RHF_AAA(void);
void cc3_t3z_RHF_AAB(void);
void cc3_t3z_RHF_BBA(void);

void T3_RHF_AAA(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk, 
		dpdbuf4 *T2, dpdbuf4 *F, dpdbuf4 *E, dpdfile2 *fIJ, dpdfile2 *fAB,
		int *occpi, int *occ_off, int *virtpi, int *vir_off);

void T3_RHF_AAB(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk, 
		dpdbuf4 *T2AA, dpdbuf4 *T2AB, dpdbuf4 *T2BA, dpdbuf4 *FAA, dpdbuf4 *FAB, dpdbuf4 *FBA,
		dpdbuf4 *EAA, dpdbuf4 *EAB, dpdbuf4 *EBA, dpdfile2 *fIJ, dpdfile2 *fij, 
		dpdfile2 *fAB, dpdfile2 *fab, int *aoccpi, int *aocc_off, int *boccpi, int *bocc_off,
		int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off);

void cc3_t3z(void)
{
  if(params.ref == 0) { /** RHF **/
    /*    cc3_t3z_RHF_AAA(); */
    /*    cc3_t3z_RHF_AAB(); */
    cc3_t3z_RHF_BBA();
  }
  else if(params.ref == 2) { /** UHF **/
    /* TBD */
  }
}

void cc3_t3z_RHF_AAA(void)
{
  int nirreps;
  int *occ_off, *occpi;
  int *vir_off, *virtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc, Gab;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int ab;
  double ***W1;
  dpdbuf4 T2, E, F;
  dpdfile2 fIJ, fAB;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  occ_off = moinfo.occ_off;
  virtpi = moinfo.virtpi;
  vir_off = moinfo.vir_off;

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&F, CC3_HET1, 0, 10, 5, 10, 7, 0, "CC3 WABEI (IE,B>A)");
  dpd_buf4_init(&E, CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 WMBIJ (I>J,MB)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  W1[Gab] = dpd_block_matrix(F.params->coltot[Gab], virtpi[Gc]);
	}

	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      T3_RHF_AAA(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2, &F, &E, &fIJ, &fAB, 
			 occpi, occ_off, virtpi, vir_off);

	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk; /* assumes totally symmetric! */

		for(ab=0; ab < F.params->coltot[Gab]; ab++) {
		  A = F.params->colorb[Gab][ab][0];
		  B = F.params->colorb[Gab][ab][1];
		  Ga = F.params->rsym[A];
		  Gb = F.params->ssym[B];
		  a = A - vir_off[Ga];
		  b = B - vir_off[Gb];

		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;

		    if(fabs(W1[Gab][ab][c]) > 1e-6) {
		      fprintf(outfile, "%d %d %d %d %d %d %20.12f\n", I, J, K, A, B, C, W1[Gab][ab][c]);
		    }

		  } /* c */
		} /* ab */
	      } /* Gab */

	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  dpd_free_block(W1[Gab], F.params->coltot[Gab], virtpi[Gc]);
	}

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(W1);

  dpd_buf4_close(&E);
  dpd_buf4_close(&F);
  dpd_buf4_close(&T2);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);
}

void cc3_t3z_RHF_AAB(void)
{
  int nirreps;
  int *occ_off, *occpi;
  int *vir_off, *virtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc, Gab;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int ab;
  double ***W1;
  dpdbuf4 T2AA, T2AB, T2BA, EAA, EAB, EBA, FAA, FAB, FBA;
  dpdfile2 fIJ, fAB, fij, fab;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  occ_off = moinfo.occ_off;
  virtpi = moinfo.virtpi;
  vir_off = moinfo.vir_off;

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij");
  dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");

  dpd_buf4_init(&T2AA, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&T2BA, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_buf4_init(&FAA, CC3_HET1, 0, 10, 5, 10, 7, 0, "CC3 WABEI (IE,B>A)");
  dpd_buf4_init(&FAB, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WaBeI (Ie,Ba)");
  dpd_buf4_init(&FBA, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
  dpd_buf4_init(&EAA, CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 WMBIJ (I>J,MB)");
  dpd_buf4_init(&EAB, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");
  dpd_buf4_init(&EBA, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WmBiJ (iJ,mB)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  W1[Gab] = dpd_block_matrix(FAA.params->coltot[Gab], virtpi[Gc]);
	}

	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      T3_RHF_AAB(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2AA, &T2AB, &T2BA, 
			 &FAA, &FAB, &FBA, &EAA, &EAB, &EBA, &fIJ, &fij, &fAB, &fab,
			 occpi, occ_off, occpi, occ_off, virtpi, vir_off, virtpi, vir_off);

	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk; /* assumes totally symmetric! */

		for(ab=0; ab < FAA.params->coltot[Gab]; ab++) {
		  A = FAA.params->colorb[Gab][ab][0];
		  B = FAA.params->colorb[Gab][ab][1];
		  Ga = FAA.params->rsym[A];
		  Gb = FAA.params->ssym[B];
		  a = A - vir_off[Ga];
		  b = B - vir_off[Gb];

		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;

		    if(fabs(W1[Gab][ab][c]) > 1e-6) {
		      fprintf(outfile, "%d %d %d %d %d %d %20.12f\n", I, J, K, A, B, C, W1[Gab][ab][c]);
		    }

		  } /* c */
		} /* ab */
	      } /* Gab */

	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  dpd_free_block(W1[Gab], FAA.params->coltot[Gab], virtpi[Gc]);
	}

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(W1);

  dpd_buf4_close(&EAA);
  dpd_buf4_close(&EAB);
  dpd_buf4_close(&EBA);
  dpd_buf4_close(&FAA);
  dpd_buf4_close(&FAB);
  dpd_buf4_close(&FBA);
  dpd_buf4_close(&T2AA);
  dpd_buf4_close(&T2AB);
  dpd_buf4_close(&T2BA);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fij);
  dpd_file2_close(&fab);
}

void cc3_t3z_RHF_BBA(void)
{
  int nirreps;
  int *occ_off, *occpi;
  int *vir_off, *virtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc, Gab;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int ab;
  double ***W1;
  dpdbuf4 T2BB, T2AB, T2BA, EBB, EAB, EBA, FBB, FAB, FBA;
  dpdfile2 fIJ, fAB, fij, fab;

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi;
  occ_off = moinfo.occ_off;
  virtpi = moinfo.virtpi;
  vir_off = moinfo.vir_off;

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fij, CC_OEI, 0, 0, 0, "fij");
  dpd_file2_init(&fab, CC_OEI, 0, 1, 1, "fab");

  dpd_buf4_init(&T2BB, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
  dpd_buf4_init(&T2BA, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&FBB, CC3_HET1, 0, 10, 5, 10, 7, 0, "CC3 Wabei (ie,b>a)");
  dpd_buf4_init(&FAB, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WAbEi (iE,bA)");
  dpd_buf4_init(&FBA, CC3_HET1, 0, 10, 5, 10, 5, 0, "CC3 WaBeI (Ie,Ba)");
  dpd_buf4_init(&EBB, CC3_HET1, 0, 0, 10, 2, 10, 0, "CC3 Wmbij (i>j,mb)");
  dpd_buf4_init(&EAB, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WmBiJ (iJ,mB)");
  dpd_buf4_init(&EBA, CC3_HET1, 0, 0, 10, 0, 10, 0, "CC3 WMbIj (Ij,Mb)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  W1[Gab] = dpd_block_matrix(FBB.params->coltot[Gab], virtpi[Gc]);
	}

	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      T3_RHF_AAB(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2BB, &T2AB, &T2BA, 
			 &FBB, &FAB, &FBA, &EBB, &EAB, &EBA, &fIJ, &fij, &fAB, &fab,
			 occpi, occ_off, occpi, occ_off, virtpi, vir_off, virtpi, vir_off);

	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk; /* assumes totally symmetric! */

		for(ab=0; ab < FBB.params->coltot[Gab]; ab++) {
		  A = FBB.params->colorb[Gab][ab][0];
		  B = FBB.params->colorb[Gab][ab][1];
		  Ga = FBB.params->rsym[A];
		  Gb = FBB.params->ssym[B];
		  a = A - vir_off[Ga];
		  b = B - vir_off[Gb];

		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;

		    if(fabs(W1[Gab][ab][c]) > 1e-6) {
		      fprintf(outfile, "%d %d %d %d %d %d %20.12f\n", I, J, K, A, B, C, W1[Gab][ab][c]);
		    }

		  } /* c */
		} /* ab */
	      } /* Gab */

	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  dpd_free_block(W1[Gab], FBB.params->coltot[Gab], virtpi[Gc]);
	}

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(W1);

  dpd_buf4_close(&EBB);
  dpd_buf4_close(&EAB);
  dpd_buf4_close(&EBA);
  dpd_buf4_close(&FBB);
  dpd_buf4_close(&FAB);
  dpd_buf4_close(&FBA);
  dpd_buf4_close(&T2BB);
  dpd_buf4_close(&T2AB);
  dpd_buf4_close(&T2BA);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fij);
  dpd_file2_close(&fab);
}
