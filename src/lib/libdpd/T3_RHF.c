
/* T3_RHF(): Computes all T3(IJK,ABC) amplitudes for a given I, J,
** K combination for input T2, F, and E intermediates.  This function
** is specifically for AAB spin cases with RHF orbitals.
**
** Arguments: 
**
**   double ***W1: The target triples amplitudes in an nirreps x AB x
**   C array.  The memory for this must be allocated externally.
**
**   int nirreps: Number of irreps.
**
**   int I: Absolute index of orbital I.
**
**   int Gi: Irrep of I.
**
**   int J: Absolute index of orbital J.
**
**   int Gj: Irrep of J.
**
**   int K: Absolute index of orbital K.
**
**   int Gk: Irrep of K.
**
**   dpdbuf4 *T2: Pointer to dpd buffer for double excitation amps,
**   ordered (IJ,AB).
**
**   dpdbuf4 *F: Pointer to dpd buffer for three-virtual-index
**   intermediate, ordered (IA,BC).
**
**   dpdbuf4 *E: Pointer to dpd buffer for three-occupied-index
**   intermediate, ordered (IJ,KA).
**
**   dpdfile2 *fIJ: Pointer to the dpd file2 for the occ-occ block of
**   the Fock matrix (or other appropriate one-electron operator).
**   
**   dpdfile2 *fAB: Pointer to the dpd file2 for the vir-vir block of
**   the Fock matrix (or other appropriate one-electron operator).
**   
**   int *occpi: Number of occupied orbitals per irrep lookup array.
**
**   int *occ_off: Offset lookup for translating between absolute and
**   relative orbital indices for occupied space.
**
**   int *virtpi: Number of virtual orbitals per irrep lookup array.
**
**   int *vir_off: Offset lookup for translating between absolute and
**   relative orbital indices for virtual space.
**
** TDC, July 2004
** modified for just RHF RAK
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libdpd/dpd.h>
#include <ccfiles.h>

void T3_RHF(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk, 
		dpdbuf4 *T2, dpdbuf4 *F, dpdbuf4 *E, dpdfile2 *fIJ, dpdfile2 *fAB, 
		int *occpi, int *occ_off, int *virtpi, int *vir_off)
{
  int h;
  int i, j, k;
  int ij, ji, ik, ki, jk, kj;
  int Gij, Gji, Gik, Gki, Gjk, Gkj, Gijk;
  int Ga, Gb, Gc;
  int Gd, Gl;
  int Gid, Gjd, Gkd;
  int Gab, Gcb, Gca;
  int Gla, Glb, Glc;
  int Gil, Gjl, Gkl;
  int a, b, c, A, B, C;
  int ab;
  int cd, bd, ad,dc;
  int id, jd, kd;
  int la, lb, lc;
  int il, jl, kl;
  int nrows, ncols, nlinks;
  double dijk, denom;
  double ***W2;

  dpd_file2_mat_init(fIJ);
  dpd_file2_mat_init(fAB);
  dpd_file2_mat_rd(fIJ);
  dpd_file2_mat_rd(fAB);

  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(T2, h);
    dpd_buf4_mat_irrep_rd(T2, h);

    dpd_buf4_mat_irrep_init(E, h);
    dpd_buf4_mat_irrep_rd(E, h);
  }

  i = I - occ_off[Gi];
  j = J - occ_off[Gj];
  k = K - occ_off[Gk];

  Gij = Gji = Gi ^ Gj;
  Gik = Gki = Gi ^ Gk;
  Gjk = Gkj = Gj ^ Gk;
  Gijk = Gi ^ Gj ^ Gk;

  ij = T2->params->rowidx[I][J];
  ji = T2->params->rowidx[J][I];
  jk = T2->params->rowidx[J][K];
  kj = T2->params->rowidx[K][J];
  ik = T2->params->rowidx[I][K];
  ki = T2->params->rowidx[K][I];

  dijk = 0.0;
  if(fIJ->params->rowtot[Gi]) dijk += fIJ->matrix[Gi][i][i];
  if(fIJ->params->rowtot[Gj]) dijk += fIJ->matrix[Gj][j][j];
  if(fIJ->params->rowtot[Gk]) dijk += fIJ->matrix[Gk][k][k];

  W2 = (double ***) malloc(nirreps * sizeof(double **));

  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk; /* assumes totally symmetric! */

    W2[Gab] = dpd_block_matrix(F->params->coltot[Gab], virtpi[Gc]);

    if(F->params->coltot[Gab] && virtpi[Gc]) {
      memset(W1[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
    }
  }

  /***** do terms that are (ab,c) and don't need sorted *****/
  for(Gd=0; Gd < nirreps; Gd++) {

    /* term 1F */
    /* +t_JkDc * F_IDAB */
    Gab = Gid = Gi ^ Gd; /* assumes totally symmetric! */
    Gc = Gjk ^ Gd;       /* assumes totally symmetric! */

    cd = T2->col_offset[Gjk][Gc];
    id = F->row_offset[Gid][I];

    F->matrix[Gid] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gid]);
    dpd_buf4_mat_irrep_rd_block(F, Gid, id, virtpi[Gd]);

    nrows = F->params->coltot[Gid];
    ncols = virtpi[Gc];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gid][0], nrows,
	      &(T2->matrix[Gjk][kj][cd]), nlinks, 1.0, W1[Gab][0], ncols);

    dpd_free_block(F->matrix[Gid], virtpi[Gd], F->params->coltot[Gid]);
  }

  for(Gl=0; Gl < nirreps; Gl++) {

    /* term 1E */
    /* -t_ILAB * E_JkLc */
    Gab = Gil = Gi ^ Gl; /* assumes totally symmetric! */
    Gc = Gjk ^ Gl;       /* assumes totally symmetric! */

    lc = E->col_offset[Gjk][Gl];
    il = T2->row_offset[Gil][I];
 
    nrows = T2->params->coltot[Gil];
    ncols = virtpi[Gc];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2->matrix[Gil][il], nrows, 
	      &(E->matrix[Gjk][jk][lc]), ncols, 1.0, W1[Gab][0], ncols);

  }

  /***** Do terms (ba,c) that need bac sorts *****/

  /* zero out W2 memory */
  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk; /* assumes totally symmetric! */
    if(F->params->coltot[Gab] && virtpi[Gc]) {
      memset(W2[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
    }
  }

  for(Gd=0; Gd < nirreps; Gd++) {
    /* term 6F */
    /* + F_JdBa * t_IkDc */
    Gab = Gjd = Gj ^ Gd; /* assumes totally symmetric! */
    Gc = Gik ^ Gd;       /* assumes totally symmetric! */

    cd = T2->col_offset[Gik][Gc];
    jd = F->row_offset[Gjd][J];

    F->matrix[Gjd] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gjd]);
    dpd_buf4_mat_irrep_rd_block(F, Gjd, jd, virtpi[Gd]);

    nrows = F->params->coltot[Gjd];
    ncols = virtpi[Gc];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gjd][0], nrows,
	      &(T2->matrix[Gik][ki][cd]), nlinks, 1.0, W2[Gab][0], ncols);

    dpd_free_block(F->matrix[Gjd], virtpi[Gd], F->params->coltot[Gjd]);
  }

  for(Gl=0; Gl < nirreps; Gl++) {
    /* term 6E */
    /* +t_JlAb * E_IkLc */
    Gab = Gjl = Gj ^ Gl; /* assumes totally symmetric! */
    Gc = Gik ^ Gl;       /* assumes totally symmetric! */

    lc = E->col_offset[Gik][Gl];
    jl = T2->row_offset[Gjl][J];

    nrows = T2->params->coltot[Gjl];
    ncols = virtpi[Gc];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2->matrix[Gjl][jl], nrows,
              &(E->matrix[Gik][ik][lc]), ncols, 1.0, W2[Gab][0], ncols);
  }

  dpd_3d_sort(W2, W1, nirreps, Gijk, F->params->coltot, F->params->colidx,
              F->params->colorb, F->params->rsym, F->params->ssym, vir_off,
              vir_off, virtpi, vir_off, F->params->colidx, bac, 1);


  /***** do (ac,b) terms that require acb sort *****/

  /* zero out W2 memory */
  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk; /* assumes totally symmetric! */
    if(F->params->coltot[Gab] && virtpi[Gc]) {
      memset(W2[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
    }
  }

  for(Gd=0; Gd < nirreps; Gd++) {
    /* term 2F */
    /* + F_IdAc * t_JkBd */
    Gca = Gid = Gi ^ Gd; /* assumes totally symmetric! */
    Gb = Gjk ^ Gd;       /* assumes totally symmetric! */

    bd = T2->col_offset[Gjk][Gb];
    id = F->row_offset[Gid][I];

    F->matrix[Gid] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gid]);
    dpd_buf4_mat_irrep_rd_block(F, Gid, id, virtpi[Gd]);

    nrows = F->params->coltot[Gid];
    ncols = virtpi[Gb];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gid][0], nrows,
	      &(T2->matrix[Gjk][jk][bd]), nlinks, 1.0, W2[Gca][0], ncols);

    dpd_free_block(F->matrix[Gid], virtpi[Gd], F->params->coltot[Gid]);
  }

  for(Gl=0; Gl < nirreps; Gl++) {
    /* term 2E */
    /* -t_IlAc * E_KjLb */
    Gca = Gil = Gi ^ Gl; /* assumes totally symmetric! */
    Gb = Gjk ^ Gl;       /* assumes totally symmetric! */

    lb = E->col_offset[Gjk][Gl];
    il = T2->row_offset[Gil][I];

    nrows = T2->params->coltot[Gil];
    ncols = virtpi[Gb];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2->matrix[Gil][il], nrows, 
	      &(E->matrix[Gjk][kj][lb]), ncols, 1.0, W2[Gca][0], ncols);
  }


  dpd_3d_sort(W2, W1, nirreps, Gijk, F->params->coltot, F->params->colidx,
              F->params->colorb, F->params->rsym, F->params->ssym, vir_off,
              vir_off, virtpi, vir_off, F->params->colidx, acb, 1);


  /***** do (ca,b) terms that need a bca sort *****/

  /* zero out W2 memory */
  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk; /* assumes totally symmetric! */
    if(F->params->coltot[Gab] && virtpi[Gc]) {
      memset(W2[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
    }
  }

  for(Gd=0; Gd < nirreps; Gd++) {
    /* term 3F */
    /* + F_KdCa * t_JiBd */
    Gca = Gkd = Gk ^ Gd; /* assumes totally symmetric! */
    Gb = Gij ^ Gd;       /* assumes totally symmetric! */

    bd = T2->col_offset[Gij][Gb];
    kd = F->row_offset[Gkd][K];

    F->matrix[Gkd] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gkd]);
    dpd_buf4_mat_irrep_rd_block(F, Gkd, kd, virtpi[Gd]);

    nrows = F->params->coltot[Gkd];
    ncols = virtpi[Gb];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gkd][0], nrows,
	      &(T2->matrix[Gij][ji][bd]), nlinks, 1.0, W2[Gca][0], ncols);

    dpd_free_block(F->matrix[Gkd], virtpi[Gd], F->params->coltot[Gkd]);
  }

  for(Gl=0; Gl < nirreps; Gl++) {
    /* term 3E */
    /* - t_KlCa * E_IjLb */
    Gca = Gkl = Gk ^ Gl; /* assumes totally symmetric! */
    Gb = Gji ^ Gl;       /* assumes totally symmetric! */

    lb = E->col_offset[Gji][Gl];
    kl = T2->row_offset[Gkl][K];

    nrows = T2->params->coltot[Gkl];
    ncols = virtpi[Gb];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2->matrix[Gkl][kl], nrows, 
	      &(E->matrix[Gji][ij][lb]), ncols, 1.0, W2[Gca][0], ncols);
  }

  dpd_3d_sort(W2, W1, nirreps, Gijk, F->params->coltot, F->params->colidx,
              F->params->colorb, F->params->rsym, F->params->ssym, vir_off,
              vir_off, virtpi, vir_off, F->params->colidx, bca, 1);


  /***** do (cb,a) terms that require a cba sort *****/
  /* zero out W2 memory */
  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk; /* assumes totally symmetric! */
    if(F->params->coltot[Gab] && virtpi[Gc]) {
      memset(W2[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
    }
  }

  for(Gd=0; Gd < nirreps; Gd++) {
    /* term 4F */
    /* + F_KdCb * t_IjAd */
    Gcb = Gkd = Gk ^ Gd; /* assumes totally symmetric! */
    Ga = Gij ^ Gd;       /* assumes totally symmetric! */

    ad = T2->col_offset[Gij][Ga];
    kd = F->row_offset[Gkd][K];

    F->matrix[Gkd] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gkd]);
    dpd_buf4_mat_irrep_rd_block(F, Gkd, kd, virtpi[Gd]);

    nrows = F->params->coltot[Gkd];
    ncols = virtpi[Ga];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gkd][0], nrows,
	      &(T2->matrix[Gij][ij][ad]), nlinks, 1.0, W2[Gcb][0], ncols);

    dpd_free_block(F->matrix[Gkd], virtpi[Gd], F->params->coltot[Gkd]);
  }

  for(Gl=0; Gl < nirreps; Gl++) {
    /* term 4E */
    /*  -t_KlCb * E_JiLa */
    Gcb  = Gkl = Gk ^ Gl; /* assumes totally symmetric! */
    Ga = Gji ^ Gl;       /* assumes totally symmetric! */

    la = E->col_offset[Gji][Gl];
    kl = T2->row_offset[Gkl][K];

    nrows = T2->params->coltot[Gkl];
    ncols = virtpi[Ga];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2->matrix[Gkl][kl], nrows, 
	      &(E->matrix[Gji][ji][la]), ncols, 1.0, W2[Gcb][0], ncols);
  }

  dpd_3d_sort(W2, W1, nirreps, Gijk, F->params->coltot, F->params->colidx,
	      F->params->colorb, F->params->rsym, F->params->ssym, vir_off, 
	      vir_off, virtpi, vir_off, F->params->colidx, cba, 1);

  /***** do (bc,a) terms that require a cab sort *****/
  /* zero out W2 memory */
  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk; /* assumes totally symmetric! */
    if(F->params->coltot[Gab] && virtpi[Gc]) {
      memset(W2[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
    }
  }

  for(Gd=0; Gd < nirreps; Gd++) {
    /* term 5F */
    /* + t_IkAd * F_JdBc */
    Gcb = Gjd = Gj ^ Gd; /* assumes totally symmetric! */
    Ga = Gik ^ Gd;       /* assumes totally symmetric! */

    ad = T2->col_offset[Gik][Ga];
    jd = F->row_offset[Gjd][J];

    F->matrix[Gjd] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gjd]);
    dpd_buf4_mat_irrep_rd_block(F, Gjd, jd, virtpi[Gd]);

    nrows = F->params->coltot[Gjd];
    ncols = virtpi[Ga];
    nlinks = virtpi[Gd];

    if(nrows && ncols && nlinks)
      C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gjd][0], nrows,
	      &(T2->matrix[Gik][ik][ad]), nlinks, 1.0, W2[Gcb][0], ncols);

    dpd_free_block(F->matrix[Gjd], virtpi[Gd], F->params->coltot[Gjd]);
  }

  for(Gl=0; Gl < nirreps; Gl++) {
    /* term 5E */
    /* - t_JlBc * E_Kila */
    Gcb  = Gjl = Gj ^ Gl; /* assumes totally symmetric! */
    Ga = Gik ^ Gl;       /* assumes totally symmetric! */

    la = E->col_offset[Gik][Gl];
    jl = T2->row_offset[Gjl][J];

    nrows = T2->params->coltot[Gjl];
    ncols = virtpi[Ga];
    nlinks = occpi[Gl];

    if(nrows && ncols && nlinks)
      C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, T2->matrix[Gjl][jl], nrows, 
	      &(E->matrix[Gik][ki][la]), ncols, 1.0, W2[Gcb][0], ncols);
  }

  dpd_3d_sort(W2, W1, nirreps, Gijk, F->params->coltot, F->params->colidx,
	      F->params->colorb, F->params->rsym, F->params->ssym, vir_off, 
	      vir_off, virtpi, vir_off, F->params->colidx, cab, 1);

  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk; /* assumes totally symmetric! */

    for(ab=0; ab < F->params->coltot[Gab]; ab++) {
      A = F->params->colorb[Gab][ab][0];
      B = F->params->colorb[Gab][ab][1];
      Ga = F->params->rsym[A];
      Gb = F->params->ssym[B];
      a = A - vir_off[Ga];
      b = B - vir_off[Gb];

      for(c=0; c < virtpi[Gc]; c++) {
	C = vir_off[Gc] + c;

	denom = dijk;
	if(fAB->params->rowtot[Ga]) denom -= fAB->matrix[Ga][a][a];
	if(fAB->params->rowtot[Gb]) denom -= fAB->matrix[Gb][b][b];
	if(fAB->params->rowtot[Gc]) denom -= fAB->matrix[Gc][c][c];

	W1[Gab][ab][c] /= denom;

      } /* c */
    } /* ab */
  } /* Gab */


  for(Gab=0; Gab < nirreps; Gab++) {
    Gc = Gab ^ Gijk; /* assumes totally symmetric! */
    dpd_free_block(W2[Gab], F->params->coltot[Gab], virtpi[Gc]);
  }
  free(W2);

  dpd_file2_mat_close(fIJ);
  dpd_file2_mat_close(fAB);
}
