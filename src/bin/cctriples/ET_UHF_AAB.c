#include <stdio.h>
#include <math.h>
#include <libciomr.h>
#include <qt.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

enum pattern {abc, acb, cab, cba, bca, bac};

void W_sort(double ***Win, double ***Wout, int nirreps, int h, int *coltot, int **colidx, 
	    int ***colorb, int *asym, int *bsym, int *aoff, int *boff,
	    int *cpi, int *coff, int **colidx_out, enum pattern index, int sum);

double ET_UHF_AAB(void)
{
  int cnt;
  int h, nirreps;
  int Gi, Gj, Gk, Ga, Gb, Gc, Gd, Gl;
  int Gji, Gij, Gjk, Gkj, Gik, Gki, Gijk;
  int Gab, Gbc, Gac, Gcb, Gca;
  int Gid, Gjd, Gkd;
  int Gil, Gjl, Gkl;
  int I, J, K, A, B, C;
  int i, j, k, a, b, c;
  int ij, ji, ik, ki, jk, kj;
  int ab, ba, ac, ca, bc, cb;
  int dc, ad, bd;
  int lc, la, lb;
  int id, jd, kd;
  int il, jl, kl;
  int *aoccpi, *avirtpi, *aocc_off, *avir_off;
  int *boccpi, *bvirtpi, *bocc_off, *bvir_off;
  double value_c, value_d, dijk, denom, ET_AAB;
  double t_ia, t_ib, t_ja, t_jb, t_kc;
  double D_jkbc, D_jkac, D_ikbc, D_ikac, D_jiab;
  int numpA, numpB, p, Gp,offset;
  int nrows, ncols, nlinks;
  dpdbuf4 T2AB, T2AA, T2BA;
  dpdbuf4 FAAints, FABints, FBAints;
  dpdbuf4 EAAints, EABints, EBAints;
  dpdbuf4 DAAints, DABints;
  int **T2AB_row_start, **T2AA_row_start, **T2BA_row_start;
  int **T2AB_col_start, **T2AA_col_start;
  int **FAB_row_start, **FBA_row_start, **FAA_row_start;
  int **EAB_col_start, **EBA_col_start, **EAA_col_start;
  dpdfile2 T1A, T1B, fIJ, fij, fAB, fab;
  double ***WABc, ***WBcA, ***WAcB, ***WcAB, ***WcBA, ***VABc;
  int nijk, mijk;
  FILE *ijkfile;

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
  FAA_row_start = init_int_matrix(nirreps, numpA);
  FAB_row_start = init_int_matrix(nirreps, numpA);
  FBA_row_start = init_int_matrix(nirreps, numpB);

  for(h=0; h < nirreps; h++) {

    for(p=0; p < numpA; p++) FAB_row_start[h][p] = -1;
    for(p=0; p < numpA; p++) FAA_row_start[h][p] = -1;
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
      for(p=0; p < aoccpi[Gp]; p++) {

	if(avirtpi[Gp^h]) 
	  FAA_row_start[h][aocc_off[Gp] + p] = nrows;

	nrows += avirtpi[Gp^h];
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
  T2AA_row_start = init_int_matrix(nirreps, numpA);
  T2AB_row_start = init_int_matrix(nirreps, numpA);
  T2BA_row_start = init_int_matrix(nirreps, numpB);

  for(h=0; h < nirreps; h++) {

    for(p=0; p < numpA; p++) {
      T2AA_row_start[h][p] = -1;
      T2AB_row_start[h][p] = -1;
    }

    for(p=0; p < numpB; p++) {
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
      for(p=0; p < aoccpi[Gp]; p++) {

	if(aoccpi[Gp^h]) 
	  T2AA_row_start[h][aocc_off[Gp] + p] = nrows;

	nrows += aoccpi[Gp^h];
      }
    }
  }

  /* Build T2 amplitude column offsets */
  T2AA_col_start = init_int_matrix(nirreps, nirreps);
  T2AB_col_start = init_int_matrix(nirreps, nirreps);

  for(h=0; h < nirreps; h++) {
    for(Gd = 0,offset=0; Gd < nirreps; Gd++) {
      Gc = Gd ^ h;
      T2AA_col_start[h][Gd] = offset;
      offset += avirtpi[Gd] * avirtpi[Gc];
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
  EAA_col_start = init_int_matrix(nirreps, nirreps);
  EAB_col_start = init_int_matrix(nirreps, nirreps);
  EBA_col_start = init_int_matrix(nirreps, nirreps);

  for(h=0; h < nirreps; h++) {
    for(Gl = 0,offset=0; Gl < nirreps; Gl++) {
      Gc = Gl ^ h;
      EAA_col_start[h][Gl] = offset;
      offset += aoccpi[Gl] * avirtpi[Gc];
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

  dpd_buf4_init(&T2AA, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&T2AB, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_init(&T2BA, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");

  dpd_buf4_init(&FAAints, CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
  dpd_buf4_init(&FABints, CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
  dpd_buf4_init(&FBAints, CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");

  dpd_buf4_init(&EAAints, CC_EINTS, 0, 0, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
  dpd_buf4_init(&EABints, CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
  dpd_buf4_init(&EBAints, CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");

  dpd_buf4_init(&DAAints, CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
  dpd_buf4_init(&DABints, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2AA, h);
    dpd_buf4_mat_irrep_rd(&T2AA, h);

    dpd_buf4_mat_irrep_init(&T2AB, h);
    dpd_buf4_mat_irrep_rd(&T2AB, h);

    dpd_buf4_mat_irrep_init(&T2BA, h);
    dpd_buf4_mat_irrep_rd(&T2BA, h);

    dpd_buf4_mat_irrep_init(&EAAints, h);
    dpd_buf4_mat_irrep_rd(&EAAints, h);

    dpd_buf4_mat_irrep_init(&EABints, h);
    dpd_buf4_mat_irrep_rd(&EABints, h);

    dpd_buf4_mat_irrep_init(&EBAints, h);
    dpd_buf4_mat_irrep_rd(&EBAints, h);

    dpd_buf4_mat_irrep_init(&DAAints, h);
    dpd_buf4_mat_irrep_rd(&DAAints, h);

    dpd_buf4_mat_irrep_init(&DABints, h);
    dpd_buf4_mat_irrep_rd(&DABints, h);
  }

  /* Compute the number of IJK combinations in this spin case */
  nijk = 0;
  for(Gi=0; Gi < nirreps; Gi++)
    for(Gj=0; Gj < nirreps; Gj++)
      for(Gk=0; Gk < nirreps; Gk++)
	for(i=0; i < aoccpi[Gi]; i++) {
	  I = aocc_off[Gi] + i;
	  for(j=0; j < aoccpi[Gj]; j++) {
	    J = aocc_off[Gj] + j;
	    for(k=0; k < boccpi[Gk]; k++) {
	      K = bocc_off[Gk] + k;

	      if(I > J) nijk++;
	    }
	  }
	}

  ffile(&ijkfile, "ijk.dat", 0);
  fprintf(ijkfile, "Spin Case: AAB\n");
  fprintf(ijkfile, "Number of IJK combintions: %d\n", nijk);
  fprintf(ijkfile, "\nCurrent IJK Combination:\n");

  mijk = 0;
  ET_AAB = 0.0;

  WABc = (double ***) malloc(nirreps * sizeof(double **));
  WBcA = (double ***) malloc(nirreps * sizeof(double **));
  WAcB = (double ***) malloc(nirreps * sizeof(double **));
  WcAB = (double ***) malloc(nirreps * sizeof(double **));
  WcBA = (double ***) malloc(nirreps * sizeof(double **));
  VABc = (double ***) malloc(nirreps * sizeof(double **));

  for(Gi=0; Gi < nirreps; Gi++) {
    for(Gj=0; Gj < nirreps; Gj++) {
      for(Gk=0; Gk < nirreps; Gk++) {

	Gij = Gji = Gi ^ Gj;
	Gjk = Gkj = Gj ^ Gk;
	Gik = Gki = Gi ^ Gk;

	Gijk = Gi ^ Gj ^ Gk;

	for(i=0; i < aoccpi[Gi]; i++) {
	  I = aocc_off[Gi] + i;
	  for(j=0; j < aoccpi[Gj]; j++) {
	    J = aocc_off[Gj] + j;
	    for(k=0; k < boccpi[Gk]; k++) {
	      K = bocc_off[Gk] + k;

	      if(I > J) {

		mijk++;
		fprintf(ijkfile, "%d\n", mijk);
		fflush(ijkfile);

		ij = EAAints.params->rowidx[I][J];
		ji = EAAints.params->rowidx[J][I];
		jk = EABints.params->rowidx[J][K];
		kj = EBAints.params->rowidx[K][J];
		ik = EABints.params->rowidx[I][K];
		ki = EBAints.params->rowidx[K][I];

		dijk = 0.0;
		if(fIJ.params->rowtot[Gi])
		  dijk += fIJ.matrix[Gi][i][i];
		if(fIJ.params->rowtot[Gj])
		  dijk += fIJ.matrix[Gj][j][j];
		if(fij.params->rowtot[Gk])
		  dijk += fij.matrix[Gk][k][k];

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  WABc[Gab] = dpd_block_matrix(FAAints.params->coltot[Gab], bvirtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {
		  /* +t_JkDc * F_IDAB */
		  Gab = Gid = Gi ^ Gd;
		  Gc = Gjk ^ Gd;

		  dc = T2AB_col_start[Gjk][Gd];
		  id = FAA_row_start[Gid][I];

		  FAAints.matrix[Gid] = dpd_block_matrix(avirtpi[Gd], FAAints.params->coltot[Gid]);
		  dpd_buf4_mat_irrep_rd_block(&FAAints, Gid, id, avirtpi[Gd]);

 		  nrows = FAAints.params->coltot[Gid];
		  ncols = bvirtpi[Gc];
		  nlinks = avirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(FAAints.matrix[Gid][0][0]), nrows,
			    &(T2AB.matrix[Gjk][jk][dc]), ncols, 1.0,
			    &(WABc[Gab][0][0]), ncols);

		  dpd_free_block(FAAints.matrix[Gid], avirtpi[Gd], FAAints.params->coltot[Gid]);

		  /* -t_IkDc * F_JDAB */
		  Gab = Gjd = Gj ^ Gd;
		  Gc = Gik ^ Gd;

		  dc = T2AB_col_start[Gik][Gd];
		  jd = FAA_row_start[Gjd][J];

		  FAAints.matrix[Gjd] = dpd_block_matrix(avirtpi[Gd], FAAints.params->coltot[Gjd]);
		  dpd_buf4_mat_irrep_rd_block(&FAAints, Gjd, jd, avirtpi[Gd]);

 		  nrows = FAAints.params->coltot[Gjd];
		  ncols = bvirtpi[Gc];
		  nlinks = avirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(FAAints.matrix[Gjd][0][0]), nrows,
			    &(T2AB.matrix[Gik][ik][dc]), ncols, 1.0,
			    &(WABc[Gab][0][0]), ncols);

		  dpd_free_block(FAAints.matrix[Gjd], avirtpi[Gd], FAAints.params->coltot[Gjd]);

		}

		for(Gl=0; Gl < nirreps; Gl++) {
		  /* -t_ILAB * E_JkLc */
		  Gab = Gil = Gi ^ Gl;
		  Gc = Gjk ^ Gl;

		  lc = EAB_col_start[Gjk][Gl];
		  il = T2AA_row_start[Gil][I];

		  nrows = T2AA.params->coltot[Gil];
		  ncols = bvirtpi[Gc];
		  nlinks = aoccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(T2AA.matrix[Gil][il][0]), nrows,
			    &(EABints.matrix[Gjk][jk][lc]), ncols, 1.0,
			    &(WABc[Gab][0][0]), ncols);


		  /* +t_JLAB * E_IkLc */
		  Gab = Gjl = Gj ^ Gl;
		  Gc = Gik ^ Gl;

		  lc = EAB_col_start[Gik][Gl];
		  jl = T2AA_row_start[Gjl][J];

		  nrows = T2AA.params->coltot[Gjl];
		  ncols = bvirtpi[Gc];
		  nlinks = aoccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2AA.matrix[Gjl][jl][0]), nrows,
			    &(EABints.matrix[Gik][ik][lc]), ncols, 1.0,
			    &(WABc[Gab][0][0]), ncols);
		}

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  WBcA[Gab] = dpd_block_matrix(FABints.params->coltot[Gab], avirtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {

		  /* -t_JkAd * F_IdBc */
		  Gbc = Gid = Gi ^ Gd;
		  Ga = Gjk ^ Gd;

		  ad = T2AB_col_start[Gjk][Ga];
		  id = FAB_row_start[Gid][I];

		  FABints.matrix[Gid] = dpd_block_matrix(bvirtpi[Gd], FABints.params->coltot[Gid]);
		  dpd_buf4_mat_irrep_rd_block(&FABints, Gid, id, bvirtpi[Gd]);

 		  nrows = FABints.params->coltot[Gid];
		  ncols = avirtpi[Ga];
		  nlinks = bvirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
			    &(FABints.matrix[Gid][0][0]), nrows,
			    &(T2AB.matrix[Gjk][jk][ad]), nlinks, 1.0,
			    &(WBcA[Gbc][0][0]), ncols);

		  dpd_free_block(FABints.matrix[Gid], bvirtpi[Gd], FABints.params->coltot[Gid]);


		  /* +t_IkAd * F_JdBc */
		  Gbc = Gjd = Gj ^ Gd;
		  Ga = Gik ^ Gd;

		  ad = T2AB_col_start[Gik][Ga];
		  jd = FAB_row_start[Gjd][J];

		  FABints.matrix[Gjd] = dpd_block_matrix(bvirtpi[Gd], FABints.params->coltot[Gjd]);
		  dpd_buf4_mat_irrep_rd_block(&FABints, Gjd, jd, bvirtpi[Gd]);

 		  nrows = FABints.params->coltot[Gjd];
		  ncols = avirtpi[Ga];
		  nlinks = bvirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			    &(FABints.matrix[Gjd][0][0]), nrows,
			    &(T2AB.matrix[Gik][ik][ad]), nlinks, 1.0,
			    &(WBcA[Gbc][0][0]), ncols);

		  dpd_free_block(FABints.matrix[Gjd], bvirtpi[Gd], FABints.params->coltot[Gjd]);
		}

		for(Gl=0; Gl < nirreps; Gl++) {

		  /* +t_IlBc * E_kJlA */
		  Gbc = Gil = Gi ^ Gl;
		  Ga = Gkj ^ Gl;

		  la = EBA_col_start[Gkj][Gl];
		  il = T2AB_row_start[Gil][I];

		  nrows = T2AB.params->coltot[Gil];
		  ncols = avirtpi[Ga];
		  nlinks = boccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2AB.matrix[Gil][il][0]), nrows,
			    &(EBAints.matrix[Gkj][kj][la]), ncols, 1.0,
			    &(WBcA[Gbc][0][0]), ncols);


		  /* -t_JlBc * E_kIlA */
		  Gbc = Gjl = Gj ^ Gl;
		  Ga = Gki ^ Gl;

		  la = EBA_col_start[Gki][Gl];
		  jl = T2AB_row_start[Gjl][J];

		  nrows = T2AB.params->coltot[Gjl];
		  ncols = avirtpi[Ga];
		  nlinks = boccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(T2AB.matrix[Gjl][jl][0]), nrows,
			    &(EBAints.matrix[Gki][ki][la]), ncols, 1.0,
			    &(WBcA[Gbc][0][0]), ncols);

		}

		W_sort(WBcA, WABc, nirreps, Gijk, FABints.params->coltot, FABints.params->colidx,
		       FABints.params->colorb, FABints.params->rsym, FABints.params->ssym,
		       avir_off, bvir_off, avirtpi, avir_off, FAAints.params->colidx, cab, 1);

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  dpd_free_block(WBcA[Gab], FABints.params->coltot[Gab], avirtpi[Gc]);

		  WAcB[Gab] = dpd_block_matrix(FABints.params->coltot[Gab], avirtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {

		  /* +t_JkBd * F_IdAc */
		  Gac = Gid = Gi ^ Gd;
		  Gb = Gjk ^ Gd;

		  bd = T2AB_col_start[Gjk][Gb];
		  id = FAB_row_start[Gid][I];

		  FABints.matrix[Gid] = dpd_block_matrix(bvirtpi[Gd], FABints.params->coltot[Gid]);
		  dpd_buf4_mat_irrep_rd_block(&FABints, Gid, id, bvirtpi[Gd]);

 		  nrows = FABints.params->coltot[Gid];
		  ncols = avirtpi[Gb];
		  nlinks = bvirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			    &(FABints.matrix[Gid][0][0]), nrows,
			    &(T2AB.matrix[Gjk][jk][bd]), nlinks, 1.0,
			    &(WAcB[Gac][0][0]), ncols);

		  dpd_free_block(FABints.matrix[Gid], bvirtpi[Gd], FABints.params->coltot[Gid]);


		  /* -t_IkBd * F_JdAc */
		  Gac = Gjd = Gj ^ Gd;
		  Gb = Gik ^ Gd;

		  bd = T2AB_col_start[Gik][Gb];
		  jd = FAB_row_start[Gjd][J];

		  FABints.matrix[Gjd] = dpd_block_matrix(bvirtpi[Gd], FABints.params->coltot[Gjd]);
		  dpd_buf4_mat_irrep_rd_block(&FABints, Gjd, jd, bvirtpi[Gd]);

 		  nrows = FABints.params->coltot[Gjd];
		  ncols = avirtpi[Gb];
		  nlinks = bvirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
			    &(FABints.matrix[Gjd][0][0]), nrows,
			    &(T2AB.matrix[Gik][ik][bd]), nlinks, 1.0,
			    &(WAcB[Gac][0][0]), ncols);

		  dpd_free_block(FABints.matrix[Gjd], bvirtpi[Gd], FABints.params->coltot[Gjd]);

		}

		for(Gl=0; Gl < nirreps; Gl++) {

		  /* -t_IlAc * E_kJlB */
		  Gac = Gil = Gi ^ Gl;
		  Gb = Gkj ^ Gl;

		  lb = EBA_col_start[Gkj][Gl];
		  il = T2AB_row_start[Gil][I];

		  nrows = T2AB.params->coltot[Gil];
		  ncols = avirtpi[Gb];
		  nlinks = boccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(T2AB.matrix[Gil][il][0]), nrows,
			    &(EBAints.matrix[Gkj][kj][lb]), ncols, 1.0,
			    &(WAcB[Gac][0][0]), ncols);


		  /* +t_JlAc * E_kIlB */
		  Gac = Gjl = Gj ^ Gl;
		  Gb = Gki ^ Gl;

		  lb = EBA_col_start[Gki][Gl];
		  jl = T2AB_row_start[Gjl][J];

		  nrows = T2AB.params->coltot[Gjl];
		  ncols = avirtpi[Gb];
		  nlinks = boccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2AB.matrix[Gjl][jl][0]), nrows,
			    &(EBAints.matrix[Gki][ki][lb]), ncols, 1.0,
			    &(WAcB[Gac][0][0]), ncols);
		}

		W_sort(WAcB, WABc, nirreps, Gijk, FABints.params->coltot, FABints.params->colidx,
		       FABints.params->colorb, FABints.params->rsym, FABints.params->ssym,
		       avir_off, bvir_off, avirtpi, avir_off, FAAints.params->colidx, acb, 1);

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  dpd_free_block(WAcB[Gab], FABints.params->coltot[Gab], avirtpi[Gc]);

		  WcBA[Gab] = dpd_block_matrix(FBAints.params->coltot[Gab], avirtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {

		  /* -t_JIAD * F_kDcB */
		  Gcb = Gkd = Gk ^ Gd;
		  Ga = Gji ^ Gd;

		  ad = T2AA_col_start[Gji][Ga];
		  kd = FBA_row_start[Gkd][K];

		  FBAints.matrix[Gkd] = dpd_block_matrix(avirtpi[Gd], FBAints.params->coltot[Gkd]);
		  dpd_buf4_mat_irrep_rd_block(&FBAints, Gkd, kd, avirtpi[Gd]);

 		  nrows = FBAints.params->coltot[Gkd];
		  ncols = avirtpi[Ga];
		  nlinks = avirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0,
			    &(FBAints.matrix[Gkd][0][0]), nrows,
			    &(T2AA.matrix[Gji][ji][ad]), nlinks, 1.0,
			    &(WcBA[Gcb][0][0]), ncols);

		  dpd_free_block(FBAints.matrix[Gkd], avirtpi[Gd], FBAints.params->coltot[Gkd]);
		}

		for(Gl=0; Gl < nirreps; Gl++) {

		  /* -t_kLcB * E_JILA */
		  Gcb = Gkl = Gk ^ Gl;
		  Ga = Gji ^ Gl;

		  la = EAA_col_start[Gji][Gl];
		  kl = T2BA_row_start[Gkl][K];

		  nrows = T2BA.params->coltot[Gkl];
		  ncols = avirtpi[Ga];
		  nlinks = aoccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0,
			    &(T2BA.matrix[Gkl][kl][0]), nrows,
			    &(EAAints.matrix[Gji][ji][la]), ncols, 1.0,
			    &(WcBA[Gcb][0][0]), ncols);
		}

		W_sort(WcBA, WABc, nirreps, Gijk, FBAints.params->coltot, FBAints.params->colidx,
		       FBAints.params->colorb, FBAints.params->rsym, FBAints.params->ssym,
		       bvir_off, avir_off, avirtpi, avir_off, FAAints.params->colidx, cba, 1);

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  dpd_free_block(WcBA[Gab], FBAints.params->coltot[Gab], avirtpi[Gc]);

		  WcAB[Gab] = dpd_block_matrix(FBAints.params->coltot[Gab], avirtpi[Gc]);
		}

		for(Gd=0; Gd < nirreps; Gd++) {

		  /* +t_JIBD * F_kDcA */
		  Gca = Gkd = Gk ^ Gd;
		  Gb = Gji ^ Gd;

		  bd = T2AA_col_start[Gji][Gb];
		  kd = FBA_row_start[Gkd][K];

		  FBAints.matrix[Gkd] = dpd_block_matrix(avirtpi[Gd], FBAints.params->coltot[Gkd]);
		  dpd_buf4_mat_irrep_rd_block(&FBAints, Gkd, kd, avirtpi[Gd]);

 		  nrows = FBAints.params->coltot[Gkd];
		  ncols = avirtpi[Gb];
		  nlinks = avirtpi[Gd];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0,
			    &(FBAints.matrix[Gkd][0][0]), nrows,
			    &(T2AA.matrix[Gji][ji][bd]), nlinks, 1.0,
			    &(WcAB[Gca][0][0]), ncols);

		  dpd_free_block(FBAints.matrix[Gkd], avirtpi[Gd], FBAints.params->coltot[Gkd]);
		}

		for(Gl=0; Gl < nirreps; Gl++) {

		  /* t_kLcA * E_JILB */
		  Gca = Gkl = Gk ^ Gl;
		  Gb = Gji ^ Gl;

		  lb = EAA_col_start[Gji][Gl];
		  kl = T2BA_row_start[Gkl][K];

		  nrows = T2BA.params->coltot[Gkl];
		  ncols = avirtpi[Gb];
		  nlinks = aoccpi[Gl];

		  if(nrows && ncols && nlinks)
		    C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0,
			    &(T2BA.matrix[Gkl][kl][0]), nrows,
			    &(EAAints.matrix[Gji][ji][lb]), ncols, 1.0,
			    &(WcAB[Gca][0][0]), ncols);
		}

		W_sort(WcAB, WABc, nirreps, Gijk, FBAints.params->coltot, FBAints.params->colidx,
		       FBAints.params->colorb, FBAints.params->rsym, FBAints.params->ssym,
		       bvir_off, avir_off, avirtpi, avir_off, FAAints.params->colidx, bca, 1);


		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  dpd_free_block(WcAB[Gab], FBAints.params->coltot[Gab], avirtpi[Gc]);

		  VABc[Gab] = dpd_block_matrix(FAAints.params->coltot[Gab], bvirtpi[Gc]);
		}

		/* Add disconnected triples and finish W and V arrays */
		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;

		  for(ab=0; ab < FAAints.params->coltot[Gab]; ab++) {
		    A = FAAints.params->colorb[Gab][ab][0];
		    Ga = FAAints.params->rsym[A];
		    a = A - avir_off[Ga];
		    B = FAAints.params->colorb[Gab][ab][1];
		    Gb = FAAints.params->ssym[B];
		    b = B - avir_off[Gb];

		    Gbc = Gb ^ Gc;
		    Gac = Ga ^ Gc;

		    for(c=0; c < bvirtpi[Gc]; c++) {
		      C = bvir_off[Gc] + c;

		      bc = DABints.params->colidx[B][C];
		      ac = DABints.params->colidx[A][C];

		      /* +t_IA * D_JkBc */
		      if(Gi == Ga && Gjk == Gbc) {
			t_ia = D_jkbc = 0.0;

			if(T1A.params->rowtot[Gi] && T1A.params->coltot[Gi])
			  t_ia = T1A.matrix[Gi][i][a];

			if(DABints.params->rowtot[Gjk] && DABints.params->coltot[Gjk])
			  D_jkbc = DABints.matrix[Gjk][jk][bc];

			VABc[Gab][ab][c] += t_ia * D_jkbc;
		      }

		      /* -t_IB * D_JkAc */
		      if(Gi == Gb && Gjk == Gac) {
			t_ib = D_jkac = 0.0;

			if(T1A.params->rowtot[Gi] && T1A.params->coltot[Gi])
			  t_ib = T1A.matrix[Gi][i][b];

			if(DABints.params->rowtot[Gjk] && DABints.params->coltot[Gjk])
			  D_jkac = DABints.matrix[Gjk][jk][ac];

			VABc[Gab][ab][c] -= t_ib * D_jkac;
		      }

		      /* -t_JA * D_IkBc */
		      if(Gj == Ga && Gik == Gbc) {
			t_ja = D_ikbc = 0.0;

			if(T1A.params->rowtot[Gj] && T1A.params->coltot[Gj])
			  t_ja = T1A.matrix[Gj][j][a];

			if(DABints.params->rowtot[Gik] && DABints.params->coltot[Gik])
			  D_ikbc = DABints.matrix[Gik][ik][bc];

			VABc[Gab][ab][c] -= t_ja * D_ikbc;
		      }

		      /* +t_JB * D_IkAc */
		      if(Gj == Gb && Gik == Gac) {
			t_jb = D_ikac = 0.0;

			if(T1A.params->rowtot[Gj] && T1A.params->coltot[Gj])
			  t_jb = T1A.matrix[Gj][j][b];

			if(DABints.params->rowtot[Gik] && DABints.params->coltot[Gik])
			  D_ikac = DABints.matrix[Gik][ik][ac];

			VABc[Gab][ab][c] += t_jb * D_ikac;
		      }

		      /* -t_kc * D_JIAB */
		      if(Gk == Gc && Gji == Gab) {
			t_kc = D_jiab = 0.0;

			if(T1B.params->rowtot[Gk] && T1B.params->coltot[Gk])
			  t_kc = T1B.matrix[Gk][k][c];

			if(DAAints.params->rowtot[Gji] && DAAints.params->coltot[Gji])
			  D_jiab = DAAints.matrix[Gji][ji][ab];

			VABc[Gab][ab][c] -= t_kc * D_jiab;
		      }

		      /* Sum V and W into V */
		      VABc[Gab][ab][c] += WABc[Gab][ab][c];

		      /* Build the rest of the denominator and divide it into W */
		      denom = dijk;
		      if(fAB.params->rowtot[Ga])
			denom -= fAB.matrix[Ga][a][a];
		      if(fAB.params->rowtot[Gb])
			denom -= fAB.matrix[Gb][b][b];
		      if(fab.params->rowtot[Gc])
			denom -= fab.matrix[Gc][c][c];

		      WABc[Gab][ab][c] /= denom;

		    } /* c */
		  } /* ab */
		} /* Gab */

		/* 1/2 Dot product of final V and W is the energy for this ijk triple */
		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;
		  ET_AAB += dot_block(WABc[Gab], VABc[Gab], FAAints.params->coltot[Gab], bvirtpi[Gc], 0.5);
		}

		for(Gab=0; Gab < nirreps; Gab++) {
		  Gc = Gab ^ Gijk;
		  dpd_free_block(WABc[Gab], FAAints.params->coltot[Gab], bvirtpi[Gc]);
		  dpd_free_block(VABc[Gab], FAAints.params->coltot[Gab], bvirtpi[Gc]);
		}

	      } /* I >= J */

	    } /* k */
	  } /* j */
	} /* i */

      } /* Gk */
    } /* Gj */
  } /* Gi */

  /*  fprintf(outfile, "cnt = %d\n", cnt); */
  /*  fprintf(outfile, "ET_AAB = %20.14f\n", ET_AAB); */

  free(WABc);
  free(WBcA);
  free(WAcB);
  free(WcAB);
  free(WcBA);
  free(VABc);

  fclose(ijkfile);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_close(&T2AA, h);
    dpd_buf4_mat_irrep_close(&T2AB, h);
    dpd_buf4_mat_irrep_close(&T2BA, h);
    dpd_buf4_mat_irrep_close(&EAAints, h);
    dpd_buf4_mat_irrep_close(&EABints, h);
    dpd_buf4_mat_irrep_close(&EBAints, h);
    dpd_buf4_mat_irrep_close(&DAAints, h);
    dpd_buf4_mat_irrep_close(&DABints, h);
  }

  dpd_buf4_close(&T2AA);
  dpd_buf4_close(&T2AB);
  dpd_buf4_close(&T2BA);
  dpd_buf4_close(&FAAints);
  dpd_buf4_close(&FABints);
  dpd_buf4_close(&FBAints);
  dpd_buf4_close(&EAAints);
  dpd_buf4_close(&EABints);
  dpd_buf4_close(&EBAints);
  dpd_buf4_close(&DAAints);
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

  free_int_matrix(FAA_row_start, nirreps);
  free_int_matrix(FAB_row_start, nirreps);
  free_int_matrix(FBA_row_start, nirreps);
  free_int_matrix(T2AA_row_start, nirreps);
  free_int_matrix(T2AB_row_start, nirreps);
  free_int_matrix(T2BA_row_start, nirreps);

  free_int_matrix(T2AA_col_start, nirreps);
  free_int_matrix(T2AB_col_start, nirreps);
  free_int_matrix(EAA_col_start, nirreps);
  free_int_matrix(EAB_col_start, nirreps);
  free_int_matrix(EBA_col_start, nirreps);

  return ET_AAB;
}
