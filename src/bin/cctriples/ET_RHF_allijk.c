#include <stdio.h>
#include <math.h>
#include <libciomr.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

enum pattern {abc, acb, cab, cba, bca, bac};

void W_sort(double ***Win, double ***Wout, int nirreps, int h, int *coltot, int **colidx, 
	    int ***colorb, int *asym, int *bsym, int *aoff, int *boff,
	    int *cpi, int *coff, enum pattern index);

double ET_RHF_allijk(void)
{
  int h, nirreps;
  int Gp, p, nump;
  int nrows, ncols, nlinks;
  int Gijk, Gid, Gkd, Gjd, Gil, Gkl, Gjl; 
  int Gab, Gba, Gbc, Gcb, Gac, Gca;
  int ab, ba, bc, cb, ac, ca;
  int cd, bd, ad, lc, lb, la;
  int il, jl, kl;
  int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
  int Gij, Gji, Gjk, Gkj, Gik, Gki;
  int I, J, K, A, B, C, D, L;
  int i, j, k, a, b, c, d, l;
  int ij, ji, ik, ki, jk, kj;
  int *occpi, *virtpi, *occ_off, *vir_off;
  int **F_row_start, **T2_col_start, **E_col_start, **T2_row_start, offset;
  double t_ia, t_jb, t_kc, D_jkbc, D_ikac, D_ijab;
  double value, dijk, value1, value2, denom, ET;
  double ***W0, ***W1, ***V;
  dpdbuf4 T2, Fints, Eints, Dints;
  dpdfile2 fIJ, fAB, T1;
  int nijk, mijk;
  FILE *ijkfile;

  timer_on("ET_RHF_all");

  nirreps = moinfo.nirreps;
  occpi = moinfo.occpi; virtpi = moinfo.virtpi;
  occ_off = moinfo.occ_off;
  vir_off = moinfo.vir_off;

  /* Compute starting row for index p in submatrix h */
  for(h=0,nump=0; h < nirreps; h++) nump += occpi[h];
  F_row_start = init_int_matrix(nirreps, nump);
  for(h=0; h < nirreps; h++) {

    /* Initialize the array for error checking */
    for(p=0; p < nump; p++) F_row_start[h][p] = -1;

    nrows = 0;
    for(Gp=0; Gp < nirreps; Gp++) {
      for(p=0; p < occpi[Gp]; p++) {

	if(virtpi[Gp^h]) 
	  F_row_start[h][occ_off[Gp] + p] = nrows;

	nrows += virtpi[Gp^h];
      }
    }
  }

  /* Compute starting row for index p in submatrix h */
  for(h=0,nump=0; h < nirreps; h++) nump += occpi[h];
  T2_row_start = init_int_matrix(nirreps, nump);
  for(h=0; h < nirreps; h++) {

    /* Initialize the array for error checking */
    for(p=0; p < nump; p++) T2_row_start[h][p] = -1;

    nrows = 0;
    for(Gp=0; Gp < nirreps; Gp++) {
      for(p=0; p < occpi[Gp]; p++) {

	if(occpi[Gp^h]) 
	  T2_row_start[h][occ_off[Gp] + p] = nrows;

	nrows += occpi[Gp^h];
      }
    }
  }

  T2_col_start = init_int_matrix(nirreps, nirreps);
  for(h=0; h < nirreps; h++) {
    for(Gd = 0,offset=0; Gd < nirreps; Gd++) {
      Gc = Gd ^ h;
      T2_col_start[h][Gd] = offset;
      offset += virtpi[Gd] * virtpi[Gc];
    }
  }

  E_col_start = init_int_matrix(nirreps, nirreps);
  for(h=0; h < nirreps; h++) {
    for(Gl = 0,offset=0; Gl < nirreps; Gl++) {
      Gc = Gl ^ h;
      E_col_start[h][Gl] = offset;
      offset += occpi[Gl] * virtpi[Gc];
    }
  }

  dpd_file2_init(&fIJ, CC_OEI, 0, 0, 0, "fIJ");
  dpd_file2_init(&fAB, CC_OEI, 0, 1, 1, "fAB");
  dpd_file2_mat_init(&fIJ);
  dpd_file2_mat_init(&fAB);
  dpd_file2_mat_rd(&fIJ);
  dpd_file2_mat_rd(&fAB);

  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1);
  dpd_file2_mat_rd(&T1);

  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&Fints, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_buf4_init(&Eints, CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
  dpd_buf4_init(&Dints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2, h);
    dpd_buf4_mat_irrep_rd(&T2, h);

    dpd_buf4_mat_irrep_init(&Eints, h);
    dpd_buf4_mat_irrep_rd(&Eints, h);

    dpd_buf4_mat_irrep_init(&Dints, h);
    dpd_buf4_mat_irrep_rd(&Dints, h);
  }

  /* Compute the number of IJK combinations */
  nijk = 0;
  for(Gi=0; Gi < nirreps; Gi++)
    for(Gj=0; Gj < nirreps; Gj++)
      for(Gk=0; Gk < nirreps; Gk++)
	for(i=0; i < occpi[Gi]; i++)
	  for(j=0; j < occpi[Gj]; j++)
	    for(k=0; k < occpi[Gk]; k++) nijk++;

  ffile(&ijkfile, "ijk.dat", 0);
  fprintf(ijkfile, "Number of IJK combintions: %d\n", nijk);
  fprintf(ijkfile, "\nCurrent IJK Combination\n");

  ET = 0.0;

  W0 = (double ***) malloc(nirreps * sizeof(double **));
  W1 = (double ***) malloc(nirreps * sizeof(double **));
  V = (double ***) malloc(nirreps * sizeof(double **));

  mijk = 0;
  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {

	Gkj = Gjk = Gk ^ Gj;
	Gji = Gij = Gi ^ Gj;
	Gik = Gki = Gi ^ Gk;

	Gijk = Gi ^ Gj ^ Gk;

	for(i=0; i < occpi[Gi]; i++) {
	  I = occ_off[Gi] + i;
	  for(j=0; j < occpi[Gj]; j++) {
	    J = occ_off[Gj] + j;
	    for(k=0; k < occpi[Gk]; k++) {
	      K = occ_off[Gk] + k;

	      mijk++;
	      fprintf(ijkfile, "%d of %d\n", mijk, nijk);
	      fflush(ijkfile);

	      ij = T2.params->rowidx[I][J];
	      ji = T2.params->rowidx[J][I];
	      ik = T2.params->rowidx[I][K];
	      ki = T2.params->rowidx[K][I];
	      jk = T2.params->rowidx[J][K];
	      kj = T2.params->rowidx[K][J];

	      dijk = 0.0;
	      if(fIJ.params->rowtot[Gi])
		dijk += fIJ.matrix[Gi][i][i];
	      if(fIJ.params->rowtot[Gj])
		dijk += fIJ.matrix[Gj][j][j];
	      if(fIJ.params->rowtot[Gk])
		dijk += fIJ.matrix[Gk][k][k];

              /* Malloc space for the W intermediate */

	      timer_on("malloc");
              for(Gab=0; Gab < nirreps; Gab++) {
                Gc = Gab ^ Gijk;

                W0[Gab] = dpd_block_matrix(Fints.params->coltot[Gab],virtpi[Gc]);
                W1[Gab] = dpd_block_matrix(Fints.params->coltot[Gab],virtpi[Gc]);
                V[Gab] = dpd_block_matrix(Fints.params->coltot[Gab],virtpi[Gc]);
              }
	      timer_off("malloc");

	      timer_on("N7 Terms");

	      /* +F_idab * t_kjcd */
	      for(Gd=0; Gd < nirreps; Gd++) {

		Gab = Gid = Gi ^ Gd;
		Gc = Gkj ^ Gd;

		/* Set up F integrals */
		Fints.matrix[Gid] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gid]);
		dpd_buf4_mat_irrep_rd_block(&Fints, Gid, F_row_start[Gid][I], virtpi[Gd]);

		/* Set up T2 amplitudes */
		cd = T2_col_start[Gkj][Gc];

		/* Set up multiplication parameters */
		nrows = Fints.params->coltot[Gid];
		ncols = virtpi[Gc];
		nlinks = virtpi[Gd];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			  &(Fints.matrix[Gid][0][0]), nrows, 
			  &(T2.matrix[Gkj][kj][cd]), nlinks, 0.0,
			  &(W0[Gab][0][0]), ncols);

		dpd_free_block(Fints.matrix[Gid], virtpi[Gd], Fints.params->coltot[Gid]);
	      }

	      /* -E_jklc * t_ilab */
	      for(Gl=0; Gl < nirreps; Gl++) {

		Gab = Gil = Gi ^ Gl;
		Gc = Gjk ^ Gl;

		/* Set up E integrals */
		lc = E_col_start[Gjk][Gl];

		/* Set up T2 amplitudes */
		il = T2_row_start[Gil][I];

		/* Set up multiplication parameters */
		nrows = T2.params->coltot[Gil];
		ncols = virtpi[Gc];
		nlinks = occpi[Gl];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			  &(T2.matrix[Gil][il][0]), nrows,
			  &(Eints.matrix[Gjk][jk][lc]), ncols, 1.0,
			  &(W0[Gab][0][0]), ncols);
	      }

	      /* Sort W[ab][c] --> W[ac][b] */
	      W_sort(W0, W1, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx, 
		     Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
		     vir_off, vir_off, virtpi, vir_off, acb);

	      /* +F_idac * t_jkbd */
	      for(Gd=0; Gd < nirreps; Gd++) {

		Gac = Gid = Gi ^ Gd;
		Gb = Gjk ^ Gd;

		Fints.matrix[Gid] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gid]);
		dpd_buf4_mat_irrep_rd_block(&Fints, Gid, F_row_start[Gid][I], virtpi[Gd]);

		bd = T2_col_start[Gjk][Gb];

		nrows = Fints.params->coltot[Gid];
		ncols = virtpi[Gb];
		nlinks = virtpi[Gd];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			  &(Fints.matrix[Gid][0][0]), nrows, 
			  &(T2.matrix[Gjk][jk][bd]), nlinks, 1.0,
			  &(W1[Gac][0][0]), ncols);

		dpd_free_block(Fints.matrix[Gid], virtpi[Gd], Fints.params->coltot[Gid]);
	      }

	      /* -E_kjlb * t_ilac */
	      for(Gl=0; Gl < nirreps; Gl++) {

		Gac = Gil = Gi ^ Gl;
		Gb = Gkj ^ Gl;

		lb = E_col_start[Gkj][Gl];

		il = T2_row_start[Gil][I];

		nrows = T2.params->coltot[Gil];
		ncols = virtpi[Gb];
		nlinks = occpi[Gl];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			  &(T2.matrix[Gil][il][0]), nrows,
			  &(Eints.matrix[Gkj][kj][lb]), ncols, 1.0,
			  &(W1[Gac][0][0]), ncols);
	      }

	      /* Sort W[ac][b] --> W[ca][b] */
	      W_sort(W1, W0, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx, 
		     Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
		     vir_off, vir_off, virtpi, vir_off, bac);

	      /* +F_kdca * t_jibd */
	      for(Gd=0; Gd < nirreps; Gd++) {

		Gca = Gkd = Gk ^ Gd;
		Gb = Gji ^ Gd;

		Fints.matrix[Gkd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gkd]);
		dpd_buf4_mat_irrep_rd_block(&Fints, Gkd, F_row_start[Gkd][K], virtpi[Gd]);

		bd = T2_col_start[Gji][Gb];

		nrows = Fints.params->coltot[Gkd];
		ncols = virtpi[Gb];
		nlinks = virtpi[Gd];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			  &(Fints.matrix[Gkd][0][0]), nrows, 
			  &(T2.matrix[Gji][ji][bd]), nlinks, 1.0,
			  &(W0[Gca][0][0]), ncols);

		dpd_free_block(Fints.matrix[Gkd], virtpi[Gd], Fints.params->coltot[Gkd]);
	      }

	      /* -E_ijlb * t_klca */
	      for(Gl=0; Gl < nirreps; Gl++) {

		Gca = Gkl = Gk ^ Gl;
		Gb = Gij ^ Gl;

		lb = E_col_start[Gij][Gl];

		kl = T2_row_start[Gkl][K];

		nrows = T2.params->coltot[Gkl];
		ncols = virtpi[Gb];
		nlinks = occpi[Gl];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			  &(T2.matrix[Gkl][kl][0]), nrows,
			  &(Eints.matrix[Gij][ij][lb]), ncols, 1.0,
			  &(W0[Gca][0][0]), ncols);
	      }

	      /* Sort W[ca][b] --> W[cb][a] */
	      W_sort(W0, W1, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx, 
		     Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
		     vir_off, vir_off, virtpi, vir_off, acb);

	      /* +F_kdcb * t_ijad */
	      for(Gd=0; Gd < nirreps; Gd++) {

		Gcb = Gkd = Gk ^ Gd;
		Ga = Gij ^ Gd;

		Fints.matrix[Gkd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gkd]);
		dpd_buf4_mat_irrep_rd_block(&Fints, Gkd, F_row_start[Gkd][K], virtpi[Gd]);

		ad = T2_col_start[Gij][Ga];

		nrows = Fints.params->coltot[Gkd];
		ncols = virtpi[Ga];
		nlinks = virtpi[Gd];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			  &(Fints.matrix[Gkd][0][0]), nrows, 
			  &(T2.matrix[Gij][ij][ad]), nlinks, 1.0,
			  &(W1[Gcb][0][0]), ncols);

		dpd_free_block(Fints.matrix[Gkd], virtpi[Gd], Fints.params->coltot[Gkd]);
	      }

	      /* -E_jila * t_klcb */
	      for(Gl=0; Gl < nirreps; Gl++) {

		Gcb = Gkl = Gk ^ Gl;
		Ga = Gji ^ Gl;

		la = E_col_start[Gji][Gl];

		kl = T2_row_start[Gkl][K];

		nrows = T2.params->coltot[Gkl];
		ncols = virtpi[Ga];
		nlinks = occpi[Gl];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			  &(T2.matrix[Gkl][kl][0]), nrows,
			  &(Eints.matrix[Gji][ji][la]), ncols, 1.0,
			  &(W1[Gcb][0][0]), ncols);
	      }

	      /* Sort W[cb][a] --> W[bc][a] */
	      W_sort(W1, W0, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx, 
		     Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
		     vir_off, vir_off, virtpi, vir_off, bac);

	      /* +F_jdbc * t_ikad */
	      for(Gd=0; Gd < nirreps; Gd++) {

		Gbc = Gjd = Gj ^ Gd;
		Ga = Gik ^ Gd;

		Fints.matrix[Gjd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gjd]);
		dpd_buf4_mat_irrep_rd_block(&Fints, Gjd, F_row_start[Gjd][J], virtpi[Gd]);

		ad = T2_col_start[Gik][Ga];

		nrows = Fints.params->coltot[Gjd];
		ncols = virtpi[Ga];
		nlinks = virtpi[Gd];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			  &(Fints.matrix[Gjd][0][0]), nrows, 
			  &(T2.matrix[Gik][ik][ad]), nlinks, 1.0,
			  &(W0[Gbc][0][0]), ncols);

		dpd_free_block(Fints.matrix[Gjd], virtpi[Gd], Fints.params->coltot[Gjd]);
	      }

	      /* -E_kila * t_jlbc */
	      for(Gl=0; Gl < nirreps; Gl++) {

		Gbc = Gjl = Gj ^ Gl;
		Ga = Gki ^ Gl;

		la = E_col_start[Gki][Gl];

		jl = T2_row_start[Gjl][J];

		nrows = T2.params->coltot[Gjl];
		ncols = virtpi[Ga];
		nlinks = occpi[Gl];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			  &(T2.matrix[Gjl][jl][0]), nrows,
			  &(Eints.matrix[Gki][ki][la]), ncols, 1.0,
			  &(W0[Gbc][0][0]), ncols);
	      }

	      /* Sort W[bc][a] --> W[ba][c] */
	      W_sort(W0, W1, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx, 
		     Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
		     vir_off, vir_off, virtpi, vir_off, acb);

	      /* +F_jdba * t_kicd */
	      for(Gd=0; Gd < nirreps; Gd++) {

		Gba = Gjd = Gj ^ Gd;
		Gc = Gki ^ Gd;

		Fints.matrix[Gjd] = dpd_block_matrix(virtpi[Gd], Fints.params->coltot[Gjd]);
		dpd_buf4_mat_irrep_rd_block(&Fints, Gjd, F_row_start[Gjd][J], virtpi[Gd]);

		cd = T2_col_start[Gki][Gc];

		nrows = Fints.params->coltot[Gjd];
		ncols = virtpi[Gc];
		nlinks = virtpi[Gd];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			  &(Fints.matrix[Gjd][0][0]), nrows, 
			  &(T2.matrix[Gki][ki][cd]), nlinks, 1.0,
			  &(W1[Gba][0][0]), ncols);

		dpd_free_block(Fints.matrix[Gjd], virtpi[Gd], Fints.params->coltot[Gjd]);
	      }

	      /* -E_iklc * t_jlba */
	      for(Gl=0; Gl < nirreps; Gl++) {

		Gba = Gjl = Gj ^ Gl;
		Gc = Gik ^ Gl;

		lc = E_col_start[Gik][Gl];

		jl = T2_row_start[Gjl][J];

		nrows = T2.params->coltot[Gjl];
		ncols = virtpi[Gc];
		nlinks = occpi[Gl];

		if(nrows && ncols && nlinks)
		  C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			  &(T2.matrix[Gjl][jl][0]), nrows,
			  &(Eints.matrix[Gik][ik][lc]), ncols, 1.0,
			  &(W1[Gba][0][0]), ncols);
	      }

	      /* Sort W[ba][c] --> W[ab][c] */
	      W_sort(W1, W0, nirreps, Gijk, Fints.params->coltot, Fints.params->colidx, 
		     Fints.params->colorb, Fints.params->rsym, Fints.params->ssym, 
		     vir_off, vir_off, virtpi, vir_off, bac);

	      timer_off("N7 Terms");

	      /* Copy W intermediate into V */
              for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;

		for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {
		  for(c=0; c < virtpi[Gc]; c++) {

		    V[Gab][ab][c] = W0[Gab][ab][c];
		  }
		}
	      }

	      timer_on("EST Terms");

	      /* Add EST terms to V */

	      for(Gab=0; Gab < nirreps; Gab++) {

		Gc = Gab ^ Gijk;

		for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {

		  A = Fints.params->colorb[Gab][ab][0];
		  Ga = Fints.params->rsym[A];
		  a = A - vir_off[Ga];
		  B = Fints.params->colorb[Gab][ab][1];
		  Gb = Fints.params->ssym[B];
		  b = B - vir_off[Gb];

		  Gbc = Gb ^ Gc;
		  Gac = Ga ^ Gc;

		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;

		    bc = Dints.params->colidx[B][C];
		    ac = Dints.params->colidx[A][C];

		    /* +t_ia * D_jkbc */
		    if(Gi == Ga && Gjk == Gbc) {
		      t_ia = D_jkbc = 0.0;

		      if(T1.params->rowtot[Gi] && T1.params->coltot[Gi])
			t_ia = T1.matrix[Gi][i][a];

		      if(Dints.params->rowtot[Gjk] && Dints.params->coltot[Gjk])
			D_jkbc = Dints.matrix[Gjk][jk][bc];

		      V[Gab][ab][c] += t_ia * D_jkbc;

		    }

		    /* +t_jb * D_ikac */
		    if(Gj == Gb && Gik == Gac) {
		      t_jb = D_ikac = 0.0;

		      if(T1.params->rowtot[Gj] && T1.params->coltot[Gj])
			t_jb = T1.matrix[Gj][j][b];

		      if(Dints.params->rowtot[Gik] && Dints.params->coltot[Gik])
			D_ikac = Dints.matrix[Gik][ik][ac];

		      V[Gab][ab][c] += t_jb * D_ikac;
		    }

		    /* +t_kc * D_ijab */
		    if(Gk == Gc && Gij == Gab) {
		      t_kc = D_ijab = 0.0;

		      if(T1.params->rowtot[Gk] && T1.params->coltot[Gk])
			t_kc = T1.matrix[Gk][k][c];

		      if(Dints.params->rowtot[Gij] && Dints.params->coltot[Gij])
			D_ijab = Dints.matrix[Gij][ij][ab];

		      V[Gab][ab][c] += t_kc * D_ijab;
		    }
		  }
		}
	      }

	      timer_off("EST Terms");

	      timer_on("Energy");

	      for(Gab=0; Gab < nirreps; Gab++) {

		Gc = Gab ^ Gijk;

		for(ab=0; ab < Fints.params->coltot[Gab]; ab++) {

		  A = Fints.params->colorb[Gab][ab][0];
		  Ga = Fints.params->rsym[A];
		  a = A - vir_off[Ga];
		  B = Fints.params->colorb[Gab][ab][1];
		  Gb = Fints.params->ssym[B];
		  b = B - vir_off[Gb];

		  Gbc = Gcb = Gb ^ Gc;
		  Gca = Gc ^ Ga;

		  for(c=0; c < virtpi[Gc]; c++) {
		    C = vir_off[Gc] + c;

		    bc = Fints.params->colidx[B][C];
		    ca = Fints.params->colidx[C][A];
		    cb = Fints.params->colidx[C][B];

		    value1 = 4.0 * W0[Gab][ab][c];
		    value1 += W0[Gbc][bc][a];
		    value1 += W0[Gca][ca][b];

		    value2 = V[Gab][ab][c] - V[Gcb][cb][a];

		    denom = dijk;
		    if(fAB.params->rowtot[Ga])
		      denom -= fAB.matrix[Ga][a][a];
		    if(fAB.params->rowtot[Gb])
		      denom -= fAB.matrix[Gb][b][b];
		    if(fAB.params->rowtot[Gc])
		      denom -= fAB.matrix[Gc][c][c];

		    ET += value1 * value2 / denom;

		  }
		}
	      }

	      timer_off("Energy");

              /* Free the W and V intermediates */
	      timer_on("malloc");
              for(Gab=0; Gab < nirreps; Gab++) {
		Gc = Gab ^ Gijk;

                dpd_free_block(W0[Gab],Fints.params->coltot[Gab],virtpi[Gc]);
                dpd_free_block(W1[Gab],Fints.params->coltot[Gab],virtpi[Gc]);
                dpd_free_block(V[Gab],Fints.params->coltot[Gab],virtpi[Gc]);
              }
	      timer_off("malloc");

	    } /* k */
	  } /* j */
	} /* i */

      } /* Gk */
    } /* Gj */
  } /* Gi */

  fclose(ijkfile);

  ET /= 3.0;

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&T2, h);
    dpd_buf4_mat_irrep_close(&Eints, h);
    dpd_buf4_mat_irrep_close(&Dints, h);
  }

  dpd_buf4_close(&T2);
  dpd_buf4_close(&Fints);
  dpd_buf4_close(&Eints);
  dpd_buf4_close(&Dints);

  dpd_file2_mat_close(&T1);
  dpd_file2_close(&T1);

  dpd_file2_mat_close(&fIJ);
  dpd_file2_mat_close(&fAB);
  dpd_file2_close(&fIJ);
  dpd_file2_close(&fAB);

  timer_off("ET_RHF_all");

  return ET;
}
