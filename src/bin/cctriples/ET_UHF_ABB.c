#include <stdio.h>
#include <math.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

double ET_UHF_ABB(void)
{
  int cnt;
  int h, nirreps;
  int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
  int Gji, Gij, Gjk, Gkj, Gik, Gki, Gijk;
  int I, J, K, A, B, C;
  int i, j, k, a, b, c;
  int ij, ji, ik, ki, jk, kj;
  int ab, ba, ac, ca, bc, cb;
  int *aoccpi, *avirtpi, *aocc_off, *avir_off;
  int *boccpi, *bvirtpi, *bocc_off, *bvir_off;
  double value_c, value_d, dijk, denom, ET_ABB;
  int numpA, numpB, p, Gp, offset;
  int nrows, ncols, nlinks;
  double t_ia, t_jb, t_jc, t_kb, t_kc;
  double D_jkbc, D_ikac, D_ikab, D_ijac, D_ijab;
  dpdbuf4 T2AB, T2BB, T2BA;
  dpdbuf4 FBBints, FABints, FBAints;
  dpdbuf4 EBBints, EABints, EBAints;
  dpdbuf4 DBBints, DABints;
  int **T2AB_row_start, **T2_BA_row_start, **T2_BB_row_start;
  int **T2_AB_col_start, **T2_BB_col_start;
  int **FAB_row_start, **FBA_row_start, **FBB_row_start;
  int **EAB_col_start, **EBA_col_start, **EBB_col_start;
  dpdfile2 T1A, T1B, fIJ, fij, fAB, fab;
  double ***WAbc, ***VAbc;

  nirreps = moinfo.nirreps;
  aoccpi = moinfo.aoccpi; 
  avirtpi = moinfo.avirtpi;
  aocc_off = moinfo.aocc_off;
  avir_off = moinfo.avir_off;
  boccpi = moinfo.boccpi; 
  bvirtpi = moinfo.bvirtpi;
  bocc_off = moinfo.bocc_off;
  bvir_off = moinfo.bvir_off;

  for(h=0,numpA=0; h < nirreps; h++) numpA += aoccpi[h];
  for(h=0,numpB=0; h < nirreps; h++) numpB += boccpi[h];

  /* Build F integral row offsets */
  FBB_row_start = init_int_matrix(nirreps, numpB);
  FAB_row_start = init_int_matrix(nirreps, numpA);
  FBA_row_start = init_int_matrix(nirreps, numpB);

  for(h=0; h < nirreps; h++) {

    for(p=0; p < numpA; p++) FAB_row_start[h][p] = -1;
    for(p=0; p < numpB; p++) FBB_row_start[h][p] = -1;
    for(p=0; p < numpB; p++) FBA_row_start[h][p] = -1;

    nrows = 0;
    for(Gp=0; Gp < nirreps; Gp++) {
      for(p=0; p < aoccpi[Gp]; p++) {

	if(bvirtpi[Gp^h]) 
	  FAB_row_start[h][aocc_off[Gp] + p] = nrows;

	nrows += bvirtpi[Gp^h];
      }
    }

    nrows = 0;
    for(Gp=0; Gp < nirreps; Gp++) {
      for(p=0; p < boccpi[Gp]; p++) {

	if(bvirtpi[Gp^h]) 
	  FBB_row_start[h][bocc_off[Gp] + p] = nrows;

	nrows += bvirtpi[Gp^h];
      }
    }

    nrows = 0;
    for(Gp=0; Gp < nirreps; Gp++) {
      for(p=0; p < boccpi[Gp]; p++) {

	if(avirtpi[Gp^h]) 
	  FBA_row_start[h][bocc_off[Gp] + p] = nrows;

	nrows += avirtpi[Gp^h];
      }
    }
  }

  /* Build T2 amplitude row offsets */
  T2BB_row_start = init_int_matrix(nirreps, numpB);
  T2AB_row_start = init_int_matrix(nirreps, numpA);
  T2BA_row_start = init_int_matrix(nirreps, numpB);

  for(h=0; h < nirreps; h++) {

    for(p=0; p < numpA; p++) {
      T2AB_row_start[h][p] = -1;
    }

    for(p=0; p < numpB; p++) {
      T2BB_row_start[h][p] = -1;
      T2BA_row_start[h][p] = -1;
    }

    nrows = 0;
    for(Gp=0; Gp < nirreps; Gp++) {
      for(p=0; p < aoccpi[Gp]; p++) {

	if(boccpi[Gp^h]) 
	  T2AB_row_start[h][aocc_off[Gp] + p] = nrows;

	nrows += boccpi[Gp^h];
      }
    }

    nrows = 0;
    for(Gp=0; Gp < nirreps; Gp++) {
      for(p=0; p < boccpi[Gp]; p++) {

	if(aoccpi[Gp^h]) 
	  T2BA_row_start[h][bocc_off[Gp] + p] = nrows;

	nrows += aoccpi[Gp^h];
      }
    }

    nrows = 0;
    for(Gp=0; Gp < nirreps; Gp++) {
      for(p=0; p < boccpi[Gp]; p++) {

	if(boccpi[Gp^h]) 
	  T2BB_row_start[h][bocc_off[Gp] + p] = nrows;

	nrows += boccpi[Gp^h];
      }
    }
  }

  /* Build T2 amplitude column offsets */
  T2BB_col_start = init_int_matrix(nirreps, nirreps);
  T2AB_col_start = init_int_matrix(nirreps, nirreps);

  for(h=0; h < nirreps; h++) {
    for(Gd = 0,offset=0; Gd < nirreps; Gd++) {
      Gc = Gd ^ h;
      T2BB_col_start[h][Gd] = offset;
      offset += bvirtpi[Gd] * bvirtpi[Gc];
    }
  }

  for(h=0; h < nirreps; h++) {
    for(Gd = 0,offset=0; Gd < nirreps; Gd++) {
      Gc = Gd ^ h;
      T2AB_col_start[h][Gd] = offset;
      offset += avirtpi[Gd] * bvirtpi[Gc];
    }
  }

  /* Build E integral column offsets */
  EBB_col_start = init_int_matrix(nirreps, nirreps);
  EAB_col_start = init_int_matrix(nirreps, nirreps);
  EBA_col_start = init_int_matrix(nirreps, nirreps);

  for(h=0; h < nirreps; h++) {
    for(Gl = 0,offset=0; Gl < nirreps; Gl++) {
      Gc = Gl ^ h;
      EBB_col_start[h][Gl] = offset;
      offset += boccpi[Gl] * bvirtpi[Gc];
    }
  }

  for(h=0; h < nirreps; h++) {
    for(Gl = 0,offset=0; Gl < nirreps; Gl++) {
      Gc = Gl ^ h;
      EAB_col_start[h][Gl] = offset;
      offset += aoccpi[Gl] * bvirtpi[Gc];
    }
  }

  for(h=0; h < nirreps; h++) {
    for(Gl = 0,offset=0; Gl < nirreps; Gl++) {
      Gc = Gl ^ h;
      EBA_col_start[h][Gl] = offset;
      offset += boccpi[Gl] * avirtpi[Gc];
    }
  }

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fij, CC_OEI, 0, 2, 2, "fij");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_init(&fab, CC_OEI, 0, 3, 3, "fab");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fij);
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_init(&fab);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_mat_rd(&fij);
  dpd_file2_mat_rd(&fAB);
  dpd_file2_mat_rd(&fab);

  dpd_file2_init(&T1A, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1A);
  dpd_file2_mat_rd(&T1A);
  dpd_file2_init(&T1B, CC_OEI, 0, 2, 3, "tia");
  dpd_file2_mat_init(&T1B);
  dpd_file2_mat_rd(&T1B);

  dpd_buf4_init(&T2BB, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");

  dpd_buf4_init(&FBBints, CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
  dpd_buf4_init(&FABints, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
  dpd_buf4_init(&FBAints, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");

  dpd_buf4_init(&EBBints, CC_EINTS, 0, 10, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&EABints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
  dpd_buf4_init(&EBAints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");

  dpd_buf4_init(&DBBints, CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
  dpd_buf4_init(&DABints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2BB, h);
    dpd_buf4_mat_irrep_rd(&T2BB, h);

    dpd_buf4_mat_irrep_init(&T2AB, h);
    dpd_buf4_mat_irrep_rd(&T2AB, h);

    dpd_buf4_mat_irrep_init(&EBBints, h);
    dpd_buf4_mat_irrep_rd(&EBBints, h);

    dpd_buf4_mat_irrep_init(&EABints, h);
    dpd_buf4_mat_irrep_rd(&EABints, h);

    dpd_buf4_mat_irrep_init(&EBAints, h);
    dpd_buf4_mat_irrep_rd(&EBAints, h);

    dpd_buf4_mat_irrep_init(&DBBints, h);
    dpd_buf4_mat_irrep_rd(&DBBints, h);

    dpd_buf4_mat_irrep_init(&DABints, h);
    dpd_buf4_mat_irrep_rd(&DABints, h);
  }

  cnt = 0;
  ET_ABB = 0.0;

  WAbc = (double ***) malloc(nirreps * sizeof(double **));
  VAbc = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {

	Gij = Gji = Gi ^ Gj;
	Gjk = Gkj = Gj ^ Gk;
	Gik = Gki = Gi ^ Gk;

	Gijk = Gi ^ Gj ^ Gk;

	for(i=0; i < aoccpi[Gi]; i++) {
	  I = aocc_off[Gi] + i;
	  for(j=0; j < boccpi[Gj]; j++) {
	    J = bocc_off[Gj] + j;
	    for(k=0; k < boccpi[Gk]; k++) {
	      K = bocc_off[Gk] + k;

	      if(J >= K) {

		ij = EABints.params->rowidx[I][J];
		ji = EBAints.params->rowidx[J][I];
		jk = EBBints.params->rowidx[J][K];
		kj = EBBints.params->rowidx[K][J];
		ik = EABints.params->rowidx[I][K];
		ki = EBAints.params->rowidx[K][I];

		dijk = 0.0;
		if(fIJ.params->rowtot[Gi])
		  dijk += fIJ.matrix[Gi][i][i];
		if(fij.params->rowtot[Gj])
		  dijk += fij.matrix[Gj][j][j];
		if(fij.params->rowtot[Gk])
		  dijk += fij.matrix[Gk][k][k];

		/* Begin connected triples */

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  WAbc[Gab] = dpd_block_matrix(FABints.params->coltot[Gab], bvirtpi[Gc]);
		}

		/* Add disconnected triples and finish W and V */
		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  VAbc[Gab] = dpd_block_matrix(FABints.params->coltot[Gab], bvirtpi[Gc]);

		  for(ab=0; ab < FABints.params->coltot[Gab]; ab++) {
		    A = FABints.params->colorb[Gab][ab][0];
		    Ga = FABints.params->rsym[A];
		    a = A - avir_off[Ga];
		    B = FABints.params->colorb[Gab][ab][1];
		    Gb = FABints.params->ssym[B];
		    b = B - bvir_off[Gb];

		    Gbc = Gb ^ Gc;
		    Gac = Ga ^ Gc;

		    for(c=0; c < bvirtpi[Gc]; c++) {
		      C = bvir_off[Gc] + c;

		      bc = DBBints.params->colidx[B][C];
		      ac = DABints.params->colidx[A][C];

		      /* +t_IA * D_jkbc */
		      if(Gi == Ga && Gjk == Gbc) {
			t_ia = D_jkbc = 0.0;

			if(T1A.params->rowtot[Gi] && T1A.params->coltot[Gi])
			  t_ia = T1A.matrix[Gi][i][a];

			if(DBBints.params->rowtot[Gjk] && DBBints.params->coltot[Gjk])
			  D_jkbc = DBBints.matrix[Gjk][jk][bc];

			VAbc[Gab][ab][c] += t_ia * D_jkbc;
		      }

		      /* +t_jb * D_IkAc */
		      if(Gj == Gb && Gik == Gac) {
			t_jb = D_ikac = 0.0;

			if(T1B.params->rowtot[Gj] && T1B.params->coltot[Gj])
			  t_jb = T1B.matrix[Gj][j][b];

			if(DABints.params->rowtot[Gik] && DABints.params->coltot[Gik])
			  D_ikac = DABints.matrix[Gik][ik][ac];

			VAbc[Gab][ab][c] += t_jb * D_ikac;
		      }

		      /* -t_jc * D_IkAb */
		      if(Gj == Gc && Gik == Gba) {
			t_jc = D_ikab = 0.0;

			if(T1B.params->rowtot[Gj] && T1B.params->coltot[Gj])
			  t_jc = T1B.matrix[Gj][j][c];

			if(DABints.params->rowtot[Gik] && DABints.params->coltot[Gik])
			  D_ikab = DABints.matrix[Gik][ik][ab];

			GAbc[Gab][ab][c] -= t_jc * D_ikab;
		      }

		      /* -t_kb * D_IjAc */
		      if(Gk == Gb && Gji == Gac) {
			t_kb = D_ijac = 0.0;

			if(T1B.params->rowtot[Gk] && T1B.params->coltot[Gk])
			  t_kb = T1B.matrix[Gk][k][b];

			if(DABints.params->rowtot[Gji] && DABints.params->coltot[Gji])
			  D_ijac = DABints.matrix[Gji][ij][ac];

			VAbc[Gab][ab][c] -= t_kb * D_ijac;
		      }

		      /* +t_kc * D_IjAb */
		      if(Gk == Gc && Gji == Gba) {
			t_kc = D_ijab = 0.0;

			if(T1B.params->rowtot[Gk] && T1B.params->coltot[Gk])
			  t_kc = T1B.matrix[Gk][k][c];

			if(DABints.params->rowtot[Gji] && DABints.params->coltot[Gji])
			  D_ijab = DABints.matrix[Gji][ij][ab];

			VAbc[Gab][ab][c] += t_kc * D_ijab;
		      }

		      /* Sum V and W into V */
		      VAbc[Gab][ab][c] += WAbc[Gab][ab][c];

		      /* Build the rest of the denominator and divide it into W */
		      denom = dijk;
		      if(fAB.params->rowtot[Ga])
			denom -= fAB.matrix[Ga][a][a];
		      if(fab.params->rowtot[Gb])
			denom -= fab.matrix[Gb][b][b];
		      if(fab.params->rowtot[Gc])
			denom -= fab.matrix[Gc][c][c];

		      WAbc[Gab][ab][c] /= denom;

		    } /* c */
		  } /* ab */
		} /* Gab */

		/* 1/2 Dot product of final V and W is the energy for this ijk triple */
		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;
		  ET_ABB += dot_block(WAbc[Gab], VAbc[Gab], FABints.params->coltot[Gab], bvirtpi[Gc], 0.5);
		  dpd_free_block(WAbc[Gab], FABints.params->coltot[Gab], bvirtpi[Gc]);
		  dpd_free_block(VAbc[Gab], FABints.params->coltot[Gab], bvirtpi[Gc]);
		}

	      } /* J >= K */

	    } /* k */
	  } /* j */
	} /* i */

      } /* Gk */
    } /* Gj */
  } /* Gi */

  free(WAbc);
  free(VAbc);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&T2BB, h);
    dpd_buf4_mat_irrep_close(&T2AB, h);
    dpd_buf4_mat_irrep_close(&EBBints, h);
    dpd_buf4_mat_irrep_close(&EABints, h);
    dpd_buf4_mat_irrep_close(&EBAints, h);
    dpd_buf4_mat_irrep_close(&DBBints, h);
    dpd_buf4_mat_irrep_close(&DABints, h);
  }

  dpd_buf4_close(&T2BB);
  dpd_buf4_close(&T2AB);
  dpd_buf4_close(&FBBints);
  dpd_buf4_close(&FABints);
  dpd_buf4_close(&FBAints);
  dpd_buf4_close(&EBBints);
  dpd_buf4_close(&EABints);
  dpd_buf4_close(&EBAints);
  dpd_buf4_close(&DBBints);
  dpd_buf4_close(&DABints);

  dpd_file2_mat_close(&T1A);
  dpd_file2_close(&T1A);
  dpd_file2_mat_close(&T1B);
  dpd_file2_close(&T1B);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fij);
  dpd_file2_mat_close(&fAB);
  dpd_file2_mat_close(&fab);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fij);
  dpd_file2_close(&fAB);
  dpd_file2_close(&fab);

  free_int_matrix(FBB_row_start, nirreps);
  free_int_matrix(FAB_row_start, nirreps);
  free_int_matrix(FBA_row_start, nirreps);
  free_int_matrix(T2BB_row_start, nirreps);
  free_int_matrix(T2AB_row_start, nirreps);
  free_int_matrix(T2BA_row_start, nirreps);

  free_int_matrix(T2BB_col_start, nirreps);
  free_int_matrix(T2AB_col_start, nirreps);
  free_int_matrix(EAA_col_start, nirreps);
  free_int_matrix(EAB_col_start, nirreps);
  free_int_matrix(EBA_col_start, nirreps);

  return ET_ABB;
}
