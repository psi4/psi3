#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void T3_RHF_AAA(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk, 
		dpdbuf4 *T2, dpdbuf4 *F, dpdbuf4 *E, dpdfile2 *fIJ, dpdfile2 *fAB,
		int *occpi, int *occ_off, int *virtpi, int *vir_off);

void T3_RHF_AAB(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk, 
		dpdbuf4 *T2AA, dpdbuf4 *T2AB, dpdbuf4 *T2BA, dpdbuf4 *FAA, dpdbuf4 *FAB, dpdbuf4 *FBA,
		dpdbuf4 *EAA, dpdbuf4 *EAB, dpdbuf4 *EBA, dpdfile2 *fIJ, dpdfile2 *fij, 
		dpdfile2 *fAB, dpdfile2 *fab, int *aoccpi, int *aocc_off, int *boccpi, int *bocc_off,
		int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off);

void cc3_UHF_AAA(void)
{
  int h, nirreps;
  int *occ_off, *occpi;
  int *vir_off, *virtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int Gij, ij, Gab, ab;
  double ***W1;
  dpdbuf4 T2, E, F;
  dpdfile2 fIJ, fAB;
  dpdfile2 t1new, t1, d1;
  dpdbuf4 D;
  int nrows, ncols;
  dpdbuf4 T2new, D2;
  dpdfile2 Fme;

  nirreps = moinfo.nirreps;
  occpi = moinfo.aoccpi;
  occ_off = moinfo.aocc_off;
  virtpi = moinfo.avirtpi;
  vir_off = moinfo.avir_off;

  dpd_file2_init(&t1new, CC_OEI, 0, 0, 1, "CC3 tIA");
  dpd_file2_mat_init(&t1new);

  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&D, h);
    dpd_buf4_mat_irrep_rd(&D, h);
  }

  dpd_file2_init(&Fme, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);

  dpd_buf4_init(&T2new, CC_MISC, 0, 0, 5, 0, 5, 0, "CC3 tIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2new, h);
  }

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&F, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
  dpd_buf4_init(&E, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gi ^ Gj;
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

	      /*
	      for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;
		for(ab=0; ab < F.params->coltot[Gab]; ab++) {
		  A = F.params->colorb[Gab][ab][0];
		  B = F.params->colorb[Gab][ab][1];
		  Ga = F.params->rsym[A];
		  Gb = F.params->rsym[B];
		  a = vir_off[Ga] - A;
		  b = vir_off[Gb] - B;

		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;
		    if(fabs(W1[Gab][ab][c]) > 1e-10) 
		      fprintf(outfile, "%d %d %d %d %d %d %20.12f\n", I, J, K, A, B, C, W1[Gab][ab][c]);
		  }
		}
	      }
	      */

	      /* t_KC <-- 1/4 t_IJKABC <IJ||AB> */

	      Gc = Gk;    /* assumes T1 is totally symmetric */
	      Gab = Gij;  /* assumes <ij||ab> is totally symmetric */

	      ij = D.params->rowidx[I][J];

	      nrows = D.params->coltot[Gij];
	      ncols = virtpi[Gc];

	      if(nrows && ncols)
		C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, D.matrix[Gij][ij], 1,
			1.0, t1new.matrix[Gk][k], 1);

	      /* t_IJAB <-- t_IJKABC F_KC */
	      Gc = Gk;    /* assumes Fme is totally symmetric */
	      Gab = Gij;  /* Assumes tIJAB is totally symmetric */

	      nrows = T2new.params->coltot[Gij];
	      ncols = virtpi[Gc];
	      ij = T2new.params->rowidx[I][J];

	      if(nrows && ncols)
		C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][0], ncols, Fme.matrix[Gk][k], 1,
		       1.0, T2new.matrix[Gij][ij], 1);

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

  dpd_file2_mat_wrt(&t1new);
  dpd_file2_mat_close(&t1new);
  /* divide T1 contributions by denominators and add to T1 */
  dpd_file2_init(&d1, CC_OEI, 0, 0, 1, "dIA");
  dpd_file2_dirprd(&d1, &t1new);
  dpd_file2_close(&d1);
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "New tIA");
  dpd_file2_axpy(&t1new, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1new);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&D, h);
  }
  dpd_buf4_close(&D);

  dpd_file2_mat_close(&Fme);
  dpd_file2_close(&Fme);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&T2new, h);
    dpd_buf4_mat_irrep_close(&T2new, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 0, 5, 1, 6, 0, "dIJAB");
  dpd_buf4_dirprd(&D2, &T2new);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");
  dpd_buf4_axpy(&T2new, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2new);
}

void cc3_UHF_BBB(void)
{
  int h, nirreps;
  int *occ_off, *occpi;
  int *vir_off, *virtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int Gij, ij, Gab, ab;
  double ***W1;
  dpdbuf4 T2, E, F;
  dpdfile2 fIJ, fAB;
  dpdfile2 t1new, t1, d1;
  dpdbuf4 D;
  int nrows, ncols;
  dpdbuf4 T2new, D2;
  dpdfile2 Fme;

  nirreps = moinfo.nirreps;
  occpi = moinfo.boccpi;
  occ_off = moinfo.bocc_off;
  virtpi = moinfo.bvirtpi;
  vir_off = moinfo.bvir_off;

  dpd_file2_init(&t1new, CC_OEI, 0, 2, 3, "CC3 tia");
  dpd_file2_mat_init(&t1new);

  dpd_buf4_init(&D, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&D, h);
    dpd_buf4_mat_irrep_rd(&D, h);
  }

  dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);

  dpd_buf4_init(&T2new, CC_MISC, 0, 10, 15, 10, 15, 0, "CC3 tijab");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2new, h);
  }

  dpd_file2_init(&fIJ, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_init(&fAB, CC_OEI, 0, 3, 3, "fab");

  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
  dpd_buf4_init(&F, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
  dpd_buf4_init(&E, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gi ^ Gj;
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

	      /* t_kc <-- 1/4 t_ijkabc <ij||ab> */

	      Gc = Gk;    /* assumes T1 is totally symmetric */
	      Gab = Gij;  /* assumes <ij||ab> is totally symmetric */

	      ij = D.params->rowidx[I][J];

	      nrows = D.params->coltot[Gij];
	      ncols = virtpi[Gc];

	      if(nrows && ncols)
		C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, D.matrix[Gij][ij], 1,
			1.0, t1new.matrix[Gk][k], 1);

	      /* t_ijab <-- t_ijkabc F_kc */
	      Gc = Gk;    /* assumes Fme is totally symmetric */
	      Gab = Gij;  /* Assumes tijab is totally symmetric */

	      nrows = T2new.params->coltot[Gij];
	      ncols = virtpi[Gc];
	      ij = T2new.params->rowidx[I][J];

	      if(nrows && ncols)
		C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][0], ncols, Fme.matrix[Gk][k], 1,
		       1.0, T2new.matrix[Gij][ij], 1);

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

  dpd_file2_mat_wrt(&t1new);
  dpd_file2_mat_close(&t1new);
  /* divide T1 contributions by denominators and add to T1 */
  dpd_file2_init(&d1, CC_OEI, 0, 2, 3, "dia");
  dpd_file2_dirprd(&d1, &t1new);
  dpd_file2_close(&d1);
  dpd_file2_init(&t1, CC_OEI, 0, 2, 3, "New tia");
  dpd_file2_axpy(&t1new, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1new);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&D, h);
  }
  dpd_buf4_close(&D);

  dpd_file2_mat_close(&Fme);
  dpd_file2_close(&Fme);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&T2new, h);
    dpd_buf4_mat_irrep_close(&T2new, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 10, 15, 11, 16, 0, "dijab");
  dpd_buf4_dirprd(&D2, &T2new);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "New tijab");
  dpd_buf4_axpy(&T2new, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2new);
}

void cc3_UHF_AAB(void)
{
  int h, nirreps;
  int *aocc_off, *aoccpi;
  int *avir_off, *avirtpi;
  int *bocc_off, *boccpi;
  int *bvir_off, *bvirtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc, Gab;
  int Gij, ij, Gjk, jk, Gbc, bc;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int ab;
  double ***W1;
  dpdbuf4 T2AA, T2AB, T2BA, EAA, EAB, EBA, FAA, FAB, FBA;
  dpdfile2 fIJ, fAB, fij, fab;
  dpdfile2 t1newA, t1newB;
  dpdfile2 t1, d1;
  dpdbuf4 DAA, DAB;
  int nrows, ncols;
  int **W_offset, offset;
  dpdfile2 FME, Fme;
  dpdbuf4 T2AAnew, D2, T2;
  dpdbuf4 T2ABnew;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi;
  aocc_off = moinfo.aocc_off;
  boccpi = moinfo.boccpi;
  bocc_off = moinfo.bocc_off;
  avirtpi = moinfo.avirtpi;
  avir_off = moinfo.avir_off;
  bvirtpi = moinfo.bvirtpi;
  bvir_off = moinfo.bvir_off;

  W_offset = init_int_matrix(nirreps, nirreps);
  for(Gab=0; Gab < nirreps; Gab++) {
    for(Ga=0,offset=0; Ga < nirreps; Ga++) {
      Gb = Ga ^ Gab;
      W_offset[Gab][Ga] = offset;
      offset += avirtpi[Ga] * avirtpi[Gb];
    }
  }

  dpd_file2_init(&t1newA, CC_OEI, 0, 0, 1, "CC3 tIA");
  dpd_file2_mat_init(&t1newA);
  dpd_file2_init(&t1newB, CC_OEI, 0, 2, 3, "CC3 tia");
  dpd_file2_mat_init(&t1newB);

  dpd_buf4_init(&DAA, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
  dpd_buf4_init(&DAB, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&DAA, h);
    dpd_buf4_mat_irrep_rd(&DAA, h);
    dpd_buf4_mat_irrep_init(&DAB, h);
    dpd_buf4_mat_irrep_rd(&DAB, h);
  }

  dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_mat_init(&FME);
  dpd_file2_mat_rd(&FME);

  dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);

  dpd_buf4_init(&T2AAnew, CC_MISC, 0, 0, 5, 0, 5, 0, "CC3 tIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2AAnew, h);
  }

  dpd_buf4_init(&T2ABnew, CC_MISC, 0, 22, 28, 22, 28, 0, "CC3 tIjAb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2ABnew, h);
  }

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");

  dpd_buf4_init(&T2AA, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_init(&T2BA, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
  dpd_buf4_init(&FAA, CC3_HET1, 0, 20, 5, 20, 7, 0, "CC3 WABEI (IE,B>A)");
  dpd_buf4_init(&FAB, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
  dpd_buf4_init(&FBA, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
  dpd_buf4_init(&EAA, CC3_HET1, 0, 0, 20, 2, 20, 0, "CC3 WMBIJ (I>J,MB)");
  dpd_buf4_init(&EAB, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
  dpd_buf4_init(&EBA, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gi ^ Gj;
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;
	Gjk = Gj ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  W1[Gab] = dpd_block_matrix(FAA.params->coltot[Gab], bvirtpi[Gc]);
	}

	for(i=0; i < aoccpi[Gi]; i++) {
	  I = aocc_off[Gi] + i;
	  for(j=0; j < aoccpi[Gj]; j++) {
	    J = aocc_off[Gj] + j;
	    for(k=0; k < boccpi[Gk]; k++) {
	      K = bocc_off[Gk] + k;

	      T3_RHF_AAB(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2AA, &T2AB, &T2BA, 
			 &FAA, &FAB, &FBA, &EAA, &EAB, &EBA, &fIJ, &fij, &fAB, &fab,
			 aoccpi, aocc_off, boccpi, bocc_off, avirtpi, avir_off, bvirtpi, bvir_off);

	      /* t_kc <-- 1/4 t_IJkABc <IJ||AB> */

	      Gc = Gk;   /* assumes T1 is totally symmetric */
	      Gab = Gij; /* assumes <ij||ab> is totally symmetric */

	      ij = DAA.params->rowidx[I][J];

	      nrows = DAA.params->coltot[Gij];
	      ncols = bvirtpi[Gc];

	      if(nrows && ncols)
		C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, DAA.matrix[Gij][ij], 1,
			1.0, t1newB.matrix[Gk][k], 1);

	      /* t_IA <-- t_IJkABc <Jk|Bc> */

	      Ga = Gi;   /* assumes T1 is totally symmetric */
	      Gbc = Gjk; /* assumes <jk|bc> is totally symmetric */

	      jk = DAB.params->rowidx[J][K];

	      for(Gab=0; Gab < nirreps; Gab++) {
		Gb = Ga ^ Gab;
		Gc = Gb ^ Gbc;

		ab = W_offset[Gab][Ga];
		bc = DAB.col_offset[Gjk][Gb];

		nrows = avirtpi[Ga];
		ncols = avirtpi[Gb] * bvirtpi[Gc];

		if(nrows && ncols)
		  C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][ab], ncols, &(DAB.matrix[Gjk][jk][bc]), 1,
			  1.0, t1newA.matrix[Gi][i], 1);
	      }

	      /* t_IJAB <-- t_IJkABc F_kc */
	      Gc = Gk;    /* assumes Fme is totally symmetric */
	      Gab = Gij;  /* Assumes tIJAB is totally symmetric */

	      nrows = T2AAnew.params->coltot[Gij];
	      ncols = bvirtpi[Gc];
	      ij = T2AAnew.params->rowidx[I][J];
 
	      if(nrows && ncols)
		C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][0], ncols, Fme.matrix[Gk][k], 1,
		       1.0, T2AAnew.matrix[Gij][ij], 1);

	      /* t_JkBc <-- t_IJkABc F_IA */
	      Ga = Gi;   /* assumes Fia is totally symmetric */
	      Gbc = Gjk; /* assumes t_jKbC is totally symmetric */

	      jk = T2AB.params->rowidx[J][K];

	      for(Gab=0; Gab < nirreps; Gab++) {
		Gb = Ga ^ Gab;
		Gc = Gb ^ Gbc;

		ab = W_offset[Gab][Ga];
		bc = T2ABnew.col_offset[Gjk][Gb];

		nrows = avirtpi[Ga];
		ncols = avirtpi[Gb] * bvirtpi[Gc];

		if(nrows && ncols)
		  C_DGEMV('t', nrows, ncols, 1.0, W1[Gab][ab], ncols, FME.matrix[Gi][i], 1,
			  1.0, &(T2ABnew.matrix[Gjk][jk][bc]), 1);
	      }

	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  dpd_free_block(W1[Gab], FAA.params->coltot[Gab], bvirtpi[Gc]);
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

  dpd_file2_mat_wrt(&t1newA);
  dpd_file2_init(&d1, CC_OEI, 0, 0, 1, "dIA");
  dpd_file2_dirprd(&d1, &t1newA);
  dpd_file2_close(&d1);
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "New tIA");
  dpd_file2_axpy(&t1newA, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1newA);

  dpd_file2_mat_wrt(&t1newB);
  dpd_file2_init(&d1, CC_OEI, 0, 2, 3, "dia");
  dpd_file2_dirprd(&d1, &t1newB);
  dpd_file2_close(&d1);
  dpd_file2_init(&t1, CC_OEI, 0, 2, 3, "New tia");
  dpd_file2_axpy(&t1newB, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1newB);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&DAA, h);
    dpd_buf4_mat_irrep_close(&DAB, h);
  }
  dpd_buf4_close(&DAA);
  dpd_buf4_close(&DAB);

  free_int_matrix(W_offset, nirreps);

  dpd_file2_mat_close(&FME);
  dpd_file2_close(&FME);

  dpd_file2_mat_close(&Fme);
  dpd_file2_close(&Fme);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&T2AAnew, h);
    dpd_buf4_mat_irrep_close(&T2AAnew, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 0, 5, 1, 6, 0, "dIJAB");
  dpd_buf4_dirprd(&D2, &T2AAnew);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");
  dpd_buf4_axpy(&T2AAnew, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2AAnew);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&T2ABnew, h);
    dpd_buf4_mat_irrep_close(&T2ABnew, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
  dpd_buf4_dirprd(&D2, &T2ABnew);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
  dpd_buf4_axpy(&T2ABnew, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2ABnew);
}

void cc3_UHF_BBA(void)
{
  int h, nirreps;
  int *aocc_off, *aoccpi;
  int *avir_off, *avirtpi;
  int *bocc_off, *boccpi;
  int *bvir_off, *bvirtpi;
  int Gi, Gj, Gk, Gijk;
  int Ga, Gb, Gc, Gab;
  int Gij, ij, Gjk, jk, Gbc, bc;
  int i, j, k, I, J, K;
  int a, b, c, A, B, C;
  int ab;
  double ***W1;
  dpdbuf4 T2BB, T2AB, T2BA, EBB, EAB, EBA, FBB, FAB, FBA;
  dpdfile2 fIJ, fAB, fij, fab;
  dpdfile2 t1newA, t1newB;
  dpdfile2 t1, d1;
  dpdbuf4 DBB, DBA;
  int nrows, ncols;
  int **W_offset, offset;
  dpdfile2 FME, Fme;
  dpdbuf4 T2BBnew, D2, T2;
  dpdbuf4 T2BAnew;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi;
  aocc_off = moinfo.aocc_off;
  boccpi = moinfo.boccpi;
  bocc_off = moinfo.bocc_off;
  avirtpi = moinfo.avirtpi;
  avir_off = moinfo.avir_off;
  bvirtpi = moinfo.bvirtpi;
  bvir_off = moinfo.bvir_off;

  W_offset = init_int_matrix(nirreps, nirreps);
  for(Gab=0; Gab < nirreps; Gab++) {
    for(Ga=0,offset=0; Ga < nirreps; Ga++) {
      Gb = Ga ^ Gab;
      W_offset[Gab][Ga] = offset;
      offset += bvirtpi[Ga] * bvirtpi[Gb];
    }
  }

  dpd_file2_init(&t1newA, CC_OEI, 0, 0, 1, "CC3 tIA");
  dpd_file2_mat_init(&t1newA);

  dpd_file2_init(&t1newB, CC_OEI, 0, 2, 3, "CC3 tia");
  dpd_file2_mat_init(&t1newB);

  dpd_buf4_init(&DBB, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
  dpd_buf4_init(&DBA, CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&DBB, h);
    dpd_buf4_mat_irrep_rd(&DBB, h);
    dpd_buf4_mat_irrep_init(&DBA, h);
    dpd_buf4_mat_irrep_rd(&DBA, h);
  }

  dpd_file2_init(&FME, CC_OEI, 0, 0, 1, "FME");
  dpd_file2_mat_init(&FME);
  dpd_file2_mat_rd(&FME);

  dpd_file2_init(&Fme, CC_OEI, 0, 2, 3, "Fme");
  dpd_file2_mat_init(&Fme);
  dpd_file2_mat_rd(&Fme);

  dpd_buf4_init(&T2BBnew, CC_MISC, 0, 10, 15, 10, 15, 0, "CC3 tijab");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2BBnew, h);
  }

  dpd_buf4_init(&T2BAnew, CC_MISC, 0, 23, 29, 23, 29, 0, "CC3 tiJaB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2BAnew, h);
  }

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");

  dpd_buf4_init(&T2BB, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_init(&T2BA, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
  dpd_buf4_init(&FBB, CC3_HET1, 0, 30, 15, 30, 17, 0, "CC3 Wabei (ie,b>a)");
  dpd_buf4_init(&FAB, CC3_HET1, 0, 24, 28, 24, 28, 0, "CC3 WaBeI (Ie,Ba)");
  dpd_buf4_init(&FBA, CC3_HET1, 0, 27, 29, 27, 29, 0, "CC3 WAbEi (iE,bA)");
  dpd_buf4_init(&EBB, CC3_HET1, 0, 10, 30, 12, 30, 0, "CC3 Wmbij (i>j,mb)");
  dpd_buf4_init(&EAB, CC3_HET1, 0, 22, 24, 22, 24, 0, "CC3 WMbIj (Ij,Mb)");
  dpd_buf4_init(&EBA, CC3_HET1, 0, 23, 27, 23, 27, 0, "CC3 WmBiJ (iJ,mB)");

  /* target T3 amplitudes go in here */
  W1 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      Gij = Gi ^ Gj;
      for(Gk=0; Gk < nirreps; Gk++) {
	Gijk = Gi ^ Gj ^ Gk;
	Gjk = Gj ^ Gk;

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  W1[Gab] = dpd_block_matrix(FBB.params->coltot[Gab], avirtpi[Gc]);
	}

	for(i=0; i < boccpi[Gi]; i++) {
	  I = bocc_off[Gi] + i;
	  for(j=0; j < boccpi[Gj]; j++) {
	    J = bocc_off[Gj] + j;
	    for(k=0; k < aoccpi[Gk]; k++) {
	      K = aocc_off[Gk] + k;

	      T3_RHF_AAB(W1, nirreps, I, Gi, J, Gj, K, Gk, &T2BB, &T2BA, &T2AB, 
			 &FBB, &FBA, &FAB, &EBB, &EBA, &EAB, &fij, &fIJ, &fab, &fAB,
			 boccpi, bocc_off, aoccpi, aocc_off, bvirtpi, bvir_off, avirtpi, avir_off);

	      /* t_KC <-- 1/4 t_ijKabC <ij||ab> */

	      Gc = Gk;  /* assumes T1 is totally symmetric */
	      Gab = Gij; /* assumes <ij||ab> is totally symmetric */

	      ij = DBB.params->rowidx[I][J];

	      nrows = DBB.params->coltot[Gij];
	      ncols = avirtpi[Gc];

	      if(nrows && ncols)
		C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, DBB.matrix[Gij][ij], 1,
			1.0, t1newA.matrix[Gk][k], 1);

	      /* t_ia <-- t_ijKabC <jK|bC> */

	      Ga = Gi;   /* assumes T1 is totally symmetric */
	      Gbc = Gjk; /* assumes <jk|bc> is totally symmetric */

	      jk = DBA.params->rowidx[J][K];

	      for(Gab=0; Gab < nirreps; Gab++) {
		Gb = Ga ^ Gab;
		Gc = Gb ^ Gbc;

		ab = W_offset[Gab][Ga];
		bc = DBA.col_offset[Gjk][Gb];

		nrows = bvirtpi[Ga];
		ncols = bvirtpi[Gb] * avirtpi[Gc];

		if(nrows && ncols)
		  C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][ab], ncols, &(DBA.matrix[Gjk][jk][bc]), 1,
			  1.0, t1newB.matrix[Gi][i], 1);
	      }

	      /* t_ijab <-- t_ijKabC F_KC */
	      Gc = Gk;    /* assumes Fme is totally symmetric */
	      Gab = Gij;  /* Assumes tIJAB is totally symmetric */

	      nrows = T2BBnew.params->coltot[Gij];
	      ncols = avirtpi[Gc];
	      ij = T2BBnew.params->rowidx[I][J];
 
	      if(nrows && ncols)
		C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][0], ncols, FME.matrix[Gk][k], 1,
		       1.0, T2BBnew.matrix[Gij][ij], 1);

	      /* t_jKbC <-- t_ijKabC F_ia */
	      Ga = Gi;   /* assumes Fia is totally symmetric */
	      Gbc = Gjk; /* assumes t_jKbC is totally symmetric */

	      jk = T2BA.params->rowidx[J][K];

	      for(Gab=0; Gab < nirreps; Gab++) {
		Gb = Ga ^ Gab;
		Gc = Gb ^ Gbc;

		ab = W_offset[Gab][Ga];
		bc = T2BAnew.col_offset[Gjk][Gb];

		nrows = bvirtpi[Ga];
		ncols = bvirtpi[Gb] * avirtpi[Gc];

		if(nrows && ncols)
		  C_DGEMV('t', nrows, ncols, 1.0, W1[Gab][ab], ncols, Fme.matrix[Gi][i], 1,
			  1.0, &(T2BAnew.matrix[Gjk][jk][bc]), 1);
	      }

	    } /* k */
	  } /* j */
	} /* i */

	for(Gab=0; Gab < nirreps; Gab++) {
	  /* This will need to change for non-totally-symmetric cases */
	  Gc = Gab ^ Gijk;
	  dpd_free_block(W1[Gab], FBB.params->coltot[Gab], avirtpi[Gc]);
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

  dpd_file2_mat_wrt(&t1newA);
  dpd_file2_init(&d1, CC_OEI, 0, 0, 1, "dIA");
  dpd_file2_dirprd(&d1, &t1newA);
  dpd_file2_close(&d1);
  dpd_file2_init(&t1, CC_OEI, 0, 0, 1, "New tIA");
  dpd_file2_axpy(&t1newA, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1newA);

  dpd_file2_mat_wrt(&t1newB);
  dpd_file2_init(&d1, CC_OEI, 0, 2, 3, "dia");
  dpd_file2_dirprd(&d1, &t1newB);
  dpd_file2_close(&d1);
  dpd_file2_init(&t1, CC_OEI, 0, 2, 3, "New tia");
  dpd_file2_axpy(&t1newB, &t1, 1, 0);
  dpd_file2_close(&t1);
  dpd_file2_close(&t1newB);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&DBB, h);
    dpd_buf4_mat_irrep_close(&DBA, h);
  }
  dpd_buf4_close(&DBB);
  dpd_buf4_close(&DBA);

  free_int_matrix(W_offset, nirreps);

  dpd_file2_mat_close(&FME);
  dpd_file2_close(&FME);

  dpd_file2_mat_close(&Fme);
  dpd_file2_close(&Fme);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&T2BBnew, h);
    dpd_buf4_mat_irrep_close(&T2BBnew, h);
  }
  dpd_buf4_init(&D2, CC_DENOM, 0, 10, 15, 11, 16, 0, "dijab");
  dpd_buf4_dirprd(&D2, &T2BBnew);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "New tijab");
  dpd_buf4_axpy(&T2BBnew, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2BBnew);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_wrt(&T2BAnew, h);
    dpd_buf4_mat_irrep_close(&T2BAnew, h);
  }
  dpd_buf4_sort(&T2BAnew, CC_MISC, qpsr, 22, 28, "CC3 tIjAb");
  dpd_buf4_close(&T2BAnew);
  dpd_buf4_init(&T2BAnew, CC_MISC, 0, 22, 28, 22, 28, 0, "CC3 tIjAb");
  dpd_buf4_init(&D2, CC_DENOM, 0, 22, 28, 22, 28, 0, "dIjAb");
  dpd_buf4_dirprd(&D2, &T2BAnew);
  dpd_buf4_close(&D2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
  dpd_buf4_axpy(&T2BAnew, &T2, 1);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&T2BAnew);
}
