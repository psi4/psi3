#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Gijab(void)
{
  int h, nirreps, i, a, m, e, I, A, M, E, Isym, Asym, Msym, Esym, row, col;
  int j, b, J, B, Jsym, Bsym;
  double value;
  struct oe_dpdfile T1, L1, g, ZZ, ZZ2, T1A, T1B;
  struct dpdbuf G, L, T, V, Z, Z1, Z2;

  nirreps = moinfo.nirreps;

  /* ( g(I,M) + L(M,E) T(I,E) ) --> Z(I,M)(TMP0)  */
  dpd_oe_file_init(&g, CC_OEI, 0, 0, "GMI", 0, outfile);
  dpd_oe_copy(&g, CC_TMP0, "Z(I,M)", 0, outfile);
  dpd_oe_file_close(&g);
  dpd_oe_file_init(&ZZ, CC_TMP0, 0, 0, "Z(I,M)", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract111(&T1, &L1, &ZZ, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&ZZ);

  /* ( g(i,m) + L(m,e) T(i,e) ) --> Z(i,m)(TMP1)  */
  dpd_oe_file_init(&g, CC_OEI, 0, 0, "Gmi", 0, outfile);
  dpd_oe_copy(&g, CC_TMP1, "Z(i,m)", 0, outfile);
  dpd_oe_file_close(&g);
  dpd_oe_file_init(&ZZ, CC_TMP1, 0, 0, "Z(i,m)", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract111(&T1, &L1, &ZZ, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&ZZ);

  /* ( g(E,A) - L(M,E) T(M,A) ) --> Z(E,A)(TMP2) */
  dpd_oe_file_init(&g, CC_OEI, 1, 1, "GAE", 0, outfile);
  dpd_oe_copy(&g, CC_TMP2, "Z(E,A)", 0, outfile);
  dpd_oe_file_close(&g);
  dpd_oe_file_init(&ZZ, CC_TMP2, 1, 1, "Z(E,A)", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract111(&L1, &T1, &ZZ, 1, 1, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&ZZ);

  /* ( g(e,a) - L(m,e) T(m,a) ) --> Z(e,a)(TMP3) */
  dpd_oe_file_init(&g, CC_OEI, 1, 1, "Gae", 0, outfile);
  dpd_oe_copy(&g, CC_TMP3, "Z(e,a)", 0, outfile);
  dpd_oe_file_close(&g);
  dpd_oe_file_init(&ZZ, CC_TMP3, 1, 1, "Z(e,a)", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract111(&L1, &T1, &ZZ, 1, 1, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&ZZ);

  dpd_oe_file_init(&T1A, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1A);
  dpd_oe_file_mat_rd(&T1A, 0, outfile);
  dpd_oe_file_init(&T1B, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&T1B);
  dpd_oe_file_mat_rd(&T1B, 0, outfile);

  /* ( - T(IA,ME) + 2 * T(I,E) T(M,A) ) --> Z(IA,ME) */
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_copy(&T, CC_TMP4, "Z(IA,ME)", 0, outfile);
  dpd_buf_close(&T);
  dpd_buf_init(&Z, CC_TMP4, 10, 10, 10, 10, 0, "Z(IA,ME)", 0, outfile);
  dpd_scm(&Z, -1.0, 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&Z, h);
      dpd_buf_mat_irrep_rd(&Z, h, 0, outfile);
      for(row=0; row < Z.params->rowtot[h]; row++) {
	  i = Z.params->roworb[h][row][0];
	  a = Z.params->roworb[h][row][1];
	  I = T1A.params->rowidx[i];  Isym = T1A.params->psym[i];
	  A = T1A.params->colidx[a];  Asym = T1A.params->qsym[a];
	  
	  for(col=0; col < Z.params->coltot[h]; col++) {
	      m = Z.params->colorb[h][col][0];
	      e = Z.params->colorb[h][col][1];
	      M = T1A.params->rowidx[m];  Msym = T1A.params->psym[m];
	      E = T1A.params->colidx[e];  Esym = T1A.params->qsym[e];

	      if((Isym==Esym) && (Msym==Asym))
		  Z.matrix[h][row][col] += (2* T1A.matrix[Isym][I][E] *
					    T1A.matrix[Msym][M][A]);
	    }
	}
      dpd_buf_mat_irrep_wrt(&Z, h, 0, outfile);
      dpd_buf_mat_irrep_close(&Z, h);
    }
  dpd_buf_close(&Z);

  /* ( - T(ia,me) + 2 * T(i,e) T(m,a) ) --> Z(ia,me) */
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_copy(&T, CC_TMP5, "Z(ia,me)", 0, outfile);
  dpd_buf_close(&T);
  dpd_buf_init(&Z, CC_TMP5, 10, 10, 10, 10, 0, "Z(ia,me)", 0, outfile);
  dpd_scm(&Z, -1.0, 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&Z, h);
      dpd_buf_mat_irrep_rd(&Z, h, 0, outfile);
      for(row=0; row < Z.params->rowtot[h]; row++) {
	  i = Z.params->roworb[h][row][0];
	  a = Z.params->roworb[h][row][1];
	  I = T1B.params->rowidx[i];  Isym = T1B.params->psym[i];
	  A = T1B.params->colidx[a];  Asym = T1B.params->qsym[a];
	  
	  for(col=0; col < Z.params->coltot[h]; col++) {
	      m = Z.params->colorb[h][col][0];
	      e = Z.params->colorb[h][col][1];
	      M = T1B.params->rowidx[m];  Msym = T1B.params->psym[m];
	      E = T1B.params->colidx[e];  Esym = T1B.params->qsym[e];

	      if((Isym==Esym) && (Msym==Asym))
		  Z.matrix[h][row][col] += (2* T1B.matrix[Isym][I][E] *
					    T1B.matrix[Msym][M][A]);
	    }
	}
      dpd_buf_mat_irrep_wrt(&Z, h, 0, outfile);
      dpd_buf_mat_irrep_close(&Z, h);
    }
  dpd_buf_close(&Z);

  /* ( - T(iA,Me) + 2 * T(i,e) T(M,A) ) --> Z(iA,Me) */
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tjAIb", 0, outfile);
  dpd_copy(&T, CC_TMP6, "Z(iA,Me)", 0, outfile);
  dpd_buf_close(&T);
  dpd_buf_init(&Z, CC_TMP6, 10, 10, 10, 10, 0, "Z(iA,Me)", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&Z, h);
      dpd_buf_mat_irrep_rd(&Z, h, 0, outfile);
      for(row=0; row < Z.params->rowtot[h]; row++) {
	  i = Z.params->roworb[h][row][0];
	  a = Z.params->roworb[h][row][1];
	  I = T1B.params->rowidx[i];  Isym = T1B.params->psym[i];
	  A = T1A.params->colidx[a];  Asym = T1A.params->qsym[a];
	  
	  for(col=0; col < Z.params->coltot[h]; col++) {
	      m = Z.params->colorb[h][col][0];
	      e = Z.params->colorb[h][col][1];
	      M = T1A.params->rowidx[m];  Msym = T1A.params->psym[m];
	      E = T1B.params->colidx[e];  Esym = T1B.params->qsym[e];

	      if((Isym==Esym) && (Msym==Asym))
		  Z.matrix[h][row][col] += (2* T1B.matrix[Isym][I][E] *
					    T1A.matrix[Msym][M][A]);
	    }
	}
      dpd_buf_mat_irrep_wrt(&Z, h, 0, outfile);
      dpd_buf_mat_irrep_close(&Z, h);
    }
  dpd_buf_close(&Z);

  /* ( - T(Ia,mE) + 2 * T(I,E) T(m,a) ) --> Z(Ia,mE) */
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tIbjA", 0, outfile);
  dpd_copy(&T, CC_TMP7, "Z(Ia,mE)", 0, outfile);
  dpd_buf_close(&T);
  dpd_buf_init(&Z, CC_TMP7, 10, 10, 10, 10, 0, "Z(Ia,mE)", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&Z, h);
      dpd_buf_mat_irrep_rd(&Z, h, 0, outfile);
      for(row=0; row < Z.params->rowtot[h]; row++) {
	  i = Z.params->roworb[h][row][0];
	  a = Z.params->roworb[h][row][1];
	  I = T1A.params->rowidx[i];  Isym = T1A.params->psym[i];
	  A = T1B.params->colidx[a];  Asym = T1B.params->qsym[a];
	  
	  for(col=0; col < Z.params->coltot[h]; col++) {
	      m = Z.params->colorb[h][col][0];
	      e = Z.params->colorb[h][col][1];
	      M = T1B.params->rowidx[m];  Msym = T1B.params->psym[m];
	      E = T1A.params->colidx[e];  Esym = T1A.params->qsym[e];

	      if((Isym==Esym) && (Msym==Asym))
		  Z.matrix[h][row][col] += (2* T1A.matrix[Isym][I][E] *
					    T1B.matrix[Msym][M][A]);
	    }
	}
      dpd_buf_mat_irrep_wrt(&Z, h, 0, outfile);
      dpd_buf_mat_irrep_close(&Z, h);
    }
  dpd_buf_close(&Z);

  dpd_oe_file_mat_close(&T1A);
  dpd_oe_file_close(&T1A);
  dpd_oe_file_mat_close(&T1B);
  dpd_oe_file_close(&T1B);

  /* L(IJ,AB) */
  dpd_buf_init(&L, CC_LAMPS, 2, 7, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_copy(&L, CC_GAMMA, "GIJAB", 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_init(&G, CC_GAMMA, 2, 7, 2, 7, 0, "GIJAB", 0, outfile);
  /* Tau(IJ,AB) */
  dpd_buf_init(&T, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_axpy(&T, &G, 1.0, 0, outfile);
  dpd_buf_close(&T);
  /* V(IJ,MN) Tau(MN,AB) */
  dpd_buf_init(&T, CC_TAMPS, 2, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 2, 2, 2, 2, 0, "VMNIJ", 0, outfile);
  dpd_contract222(&V, &T, &G, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  dpd_buf_close(&G);
  /* - ( Z(I,M) Tau(MJ,AB) - Z(J,M) Tau(MI,AB) ) */
  dpd_buf_init(&Z1, CC_TMP8, 0, 7, 0, 7, 0, "Z1(IJ,AB)", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 0, 7, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_oe_file_init(&ZZ, CC_TMP0, 0, 0, "Z(I,M)", 0, outfile);
  dpd_contract212(&ZZ, &T, &Z1, 1, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&ZZ);
  dpd_buf_close(&T);
  dpd_swap12(&Z1, CC_TMP9, 0, 7, "Z2(JI,AB)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP9, 0, 7, 0, 7, 0, "Z2(JI,AB)", 0, outfile);
  dpd_axpy(&Z2, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&G, CC_GAMMA, 0, 7, 2, 7, 0, "GIJAB", 0, outfile);
  dpd_axpy(&Z1, &G, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&G);
  /* - ( Z(E,A) Tau(IJ,BE) - Z(E,B) Tau(IJ,AE) ) */
  dpd_buf_init(&Z1, CC_TMP8, 2, 5, 2, 5, 0, "Z1(IJ,AB)", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 2, 5, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_oe_file_init(&ZZ, CC_TMP2, 1, 1, "Z(E,A)", 0, outfile);
  dpd_contract221(&T, &ZZ, &Z1, 3, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&ZZ);
  dpd_buf_close(&T);
  dpd_swap34(&Z1, CC_TMP9, 2, 5, "Z2(IJ,BA)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP9, 2, 5, 2, 5, 0, "Z2(IJ,BA)", 0, outfile);
  dpd_axpy(&Z2, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&G, CC_GAMMA, 2, 5, 2, 7, 0, "GIJAB", 0, outfile);
  dpd_axpy(&Z1, &G, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&G);
  /* - P(IJ) P(AB) ( T'(IA,ME) (TMP4) V(JB,ME) + T(IA,me) (T2) V(JB,me) ) */
  dpd_buf_init(&Z, CC_TMP8, 10, 10, 10, 10, 0, "Z(IA,JB)", 0, outfile);
  dpd_buf_init(&T, CC_TMP4, 10, 10, 10, 10, 0, "Z(IA,ME)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIAJB", 0, outfile);
  dpd_contract222(&T, &V, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIAjb", 0, outfile);
  dpd_contract222(&T, &V, &Z, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  dpd_swap13(&Z, CC_TMP9, 10, 10, "Z(JA,IB)", 0, outfile);
  dpd_swap24(&Z, CC_TMP10, 10, 10, "Z(IB,JA)", 0, outfile);
  dpd_swapbk(&Z, CC_TMP11, 10, 10, "Z(JB,IA)", 0, outfile);
  dpd_buf_init(&Z1, CC_TMP9, 10, 10, 10, 10, 0, "Z(JA,IB)", 0, outfile);
  dpd_axpy(&Z1, &Z, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP10, 10, 10, 10, 10, 0, "Z(IB,JA)", 0, outfile);
  dpd_axpy(&Z1, &Z, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP11, 10, 10, 10, 10, 0, "Z(JB,IA)", 0, outfile);
  dpd_axpy(&Z1, &Z, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_swap23(&Z, CC_TMP9, 0, 5, "Z(IJ,AB)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 2, 7, 0, "GIJAB", 0, outfile);
  dpd_buf_init(&Z, CC_TMP9, 0, 5, 0, 5, 0, "Z(IJ,AB)", 0, outfile);
  /* I don't understand this factor of 1/2 that shows up here */
  dpd_axpy(&Z, &G, -0.5, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&G);
  /* T'(IA,ME) (TMP4) L(M,E) + T'(IA,me) (T2) L(m,e) --> ZZ(I,A) */
  dpd_oe_file_init(&ZZ, CC_TMP8, 0, 1, "ZZ(I,A)", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_buf_init(&T, CC_TMP4, 10, 10, 10, 10, 0, "Z(IA,ME)", 0, outfile);
  dpd_contract121(&T, &L1, &ZZ, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_oe_file_close(&L1);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_contract121(&T, &L1, &ZZ, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&ZZ);
  /* - P(IJ) P(AB) ZZ(I,A) T(J,B) --> G(IJ,AB) */
  dpd_oe_file_init(&ZZ, CC_TMP8, 0, 1, "ZZ(I,A)", 0, outfile);
  dpd_oe_file_mat_init(&ZZ);
  dpd_oe_file_mat_rd(&ZZ, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 2, 7, 2, 7, 0, "GIJAB", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);
      for(row=0; row < G.params->rowtot[h]; row++) {
	  i = G.params->roworb[h][row][0];
	  j = G.params->roworb[h][row][1];
	  
	  for(col=0; col < G.params->coltot[h]; col++) {
	      a = G.params->colorb[h][col][0];
	      b = G.params->colorb[h][col][1];

	      value = 0.0;

	      I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	      J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];
	      A = ZZ.params->colidx[a]; Asym = ZZ.params->qsym[a];
	      B = T1.params->colidx[b]; Bsym = T1.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  value += ZZ.matrix[Isym][I][A] * T1.matrix[Jsym][J][B];

	      I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	      J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];

	      if((Jsym==Asym) && (Isym==Bsym))
		  value -= ZZ.matrix[Jsym][J][A] * T1.matrix[Isym][I][B];

	      I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	      J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];
	      A = T1.params->colidx[a]; Asym = T1.params->qsym[a];
	      B = ZZ.params->colidx[b]; Bsym = ZZ.params->qsym[b];

	      if((Isym==Bsym) && (Jsym==Asym))
		  value -= ZZ.matrix[Isym][I][B] * T1.matrix[Jsym][J][A];

	      I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	      J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];

	      if((Isym==Asym) && (Jsym==Bsym))
		  value += T1.matrix[Isym][I][A] * ZZ.matrix[Jsym][J][B];

	      G.matrix[h][row][col] -= value;
	      
	    }
	}
      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  dpd_buf_close(&G);

  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);
  dpd_oe_file_mat_close(&ZZ);
  dpd_oe_file_close(&ZZ);

  /* T(J,E) L(M,E) --> ZZ(J,M) */
  dpd_oe_file_init(&ZZ, CC_TMP8, 0, 0, "Z(J,M)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_contract111(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&T1);
  /* ZZ(J,M) T(M,B) --> ZZ2(J,B) */
  dpd_oe_file_init(&ZZ2, CC_TMP9, 0, 1, "ZZ2(J,B)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract111(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_close(&ZZ);
  dpd_oe_file_close(&ZZ2);

  /* 3 P(IJ) P(AB) T(I,A) ZZ(J,B) --> G(IJ,AB) */
  dpd_oe_file_init(&ZZ, CC_TMP9, 0, 1, "ZZ2(J,B)", 0, outfile);
  dpd_oe_file_mat_init(&ZZ);
  dpd_oe_file_mat_rd(&ZZ, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 2, 7, 2, 7, 0, "GIJAB", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);
      for(row=0; row < G.params->rowtot[h]; row++) {
	  i = G.params->roworb[h][row][0];
	  j = G.params->roworb[h][row][1];
	  
	  for(col=0; col < G.params->coltot[h]; col++) {
	      a = G.params->colorb[h][col][0];
	      b = G.params->colorb[h][col][1];

	      value = 0.0;

	      I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	      J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];
	      A = T1.params->colidx[a]; Asym = T1.params->qsym[a];
	      B = ZZ.params->colidx[b]; Bsym = ZZ.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  value += T1.matrix[Isym][I][A] * ZZ.matrix[Jsym][J][B];

	      I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	      J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];

	      if((Jsym==Asym) && (Isym==Bsym))
		  value -= T1.matrix[Jsym][J][A] * ZZ.matrix[Isym][I][B];

	      I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	      J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];
	      A = ZZ.params->colidx[a]; Asym = ZZ.params->qsym[a];
	      B = T1.params->colidx[b]; Bsym = T1.params->qsym[b];

	      if((Isym==Bsym) && (Jsym==Asym))
		  value -= T1.matrix[Isym][I][B] * ZZ.matrix[Jsym][J][A];

	      I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	      J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];

	      if((Isym==Asym) && (Jsym==Bsym))
		  value += ZZ.matrix[Isym][I][A] * T1.matrix[Jsym][J][B];

	      G.matrix[h][row][col] += 3.0 * value;
	      
	    }
	}
      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  dpd_scm(&G, 0.5, 0, outfile);
  dpd_buf_close(&G);

  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);
  dpd_oe_file_mat_close(&ZZ);
  dpd_oe_file_close(&ZZ);



  /* L(ij,ab) */
  dpd_buf_init(&L, CC_LAMPS, 2, 7, 2, 7, 0, "Lijab", 0, outfile);
  dpd_copy(&L, CC_GAMMA, "Gijab", 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_init(&G, CC_GAMMA, 2, 7, 2, 7, 0, "Gijab", 0, outfile);
  /* Tau(ij,ab) */
  dpd_buf_init(&T, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_axpy(&T, &G, 1.0, 0, outfile);
  dpd_buf_close(&T);
  /* V(ij,mn) Tau(mn,ab) */
  dpd_buf_init(&T, CC_TAMPS, 2, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 2, 2, 2, 2, 0, "Vmnij", 0, outfile);
  dpd_contract222(&V, &T, &G, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  dpd_buf_close(&G);
  /* - ( Z(i,m) Tau(mj,ab) - Z(j,m) Tau(mi,ab) ) */
  dpd_buf_init(&Z1, CC_TMP8, 0, 7, 0, 7, 0, "Z1(ij,ab)", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 0, 7, 2, 7, 0, "tauijab", 0, outfile);
  dpd_oe_file_init(&ZZ, CC_TMP1, 0, 0, "Z(i,m)", 0, outfile);
  dpd_contract212(&ZZ, &T, &Z1, 1, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&ZZ);
  dpd_buf_close(&T);
  dpd_swap12(&Z1, CC_TMP9, 0, 7, "Z2(ji,ab)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP9, 0, 7, 0, 7, 0, "Z2(ji,ab)", 0, outfile);
  dpd_axpy(&Z2, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&G, CC_GAMMA, 0, 7, 2, 7, 0, "Gijab", 0, outfile);
  dpd_axpy(&Z1, &G, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&G);
  /* - ( Z(e,a) Tau(ij,be) - Z(e,b) Tau(ij,ae) ) */
  dpd_buf_init(&Z1, CC_TMP8, 2, 5, 2, 5, 0, "Z1(ij,ab)", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 2, 5, 2, 7, 0, "tauijab", 0, outfile);
  dpd_oe_file_init(&ZZ, CC_TMP3, 1, 1, "Z(e,a)", 0, outfile);
  dpd_contract221(&T, &ZZ, &Z1, 3, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&ZZ);
  dpd_buf_close(&T);
  dpd_swap34(&Z1, CC_TMP9, 2, 5, "Z2(ij,ba)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP9, 2, 5, 2, 5, 0, "Z2(ij,ba)", 0, outfile);
  dpd_axpy(&Z2, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&G, CC_GAMMA, 2, 5, 2, 7, 0, "Gijab", 0, outfile);
  dpd_axpy(&Z1, &G, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&G);
  /* - P(ij) P(ab) ( T'(ia,me) (TMP5) V(jb,me) + T(ia,ME) (T2) V(jb,ME) ) */
  dpd_buf_init(&Z, CC_TMP8, 10, 10, 10, 10, 0, "Z(ia,jb)", 0, outfile);
  dpd_buf_init(&T, CC_TMP5, 10, 10, 10, 10, 0, "Z(ia,me)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "Viajb", 0, outfile);
  dpd_contract222(&T, &V, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "ViaJB", 0, outfile);
  dpd_contract222(&T, &V, &Z, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  dpd_swap13(&Z, CC_TMP9, 10, 10, "Z(ja,ib)", 0, outfile);
  dpd_swap24(&Z, CC_TMP10, 10, 10, "Z(ib,ja)", 0, outfile);
  dpd_swapbk(&Z, CC_TMP11, 10, 10, "Z(jb,ia)", 0, outfile);
  dpd_buf_init(&Z1, CC_TMP9, 10, 10, 10, 10, 0, "Z(ja,ib)", 0, outfile);
  dpd_axpy(&Z1, &Z, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP10, 10, 10, 10, 10, 0, "Z(ib,ja)", 0, outfile);
  dpd_axpy(&Z1, &Z, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP11, 10, 10, 10, 10, 0, "Z(jb,ia)", 0, outfile);
  dpd_axpy(&Z1, &Z, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_swap23(&Z, CC_TMP9, 0, 5, "Z(ij,ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 2, 7, 0, "Gijab", 0, outfile);
  dpd_buf_init(&Z, CC_TMP9, 0, 5, 0, 5, 0, "Z(ij,ab)", 0, outfile);
  /* I don't understand this factor of 1/2 that shows up here */
  dpd_axpy(&Z, &G, -0.5, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&G);
  /* T'(ia,me) (TMP5) L(m,e) + T'(ia,ME) (T2) L(M,E) --> ZZ(i,a) */
  dpd_oe_file_init(&ZZ, CC_TMP8, 0, 1, "ZZ(i,a)", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_buf_init(&T, CC_TMP5, 10, 10, 10, 10, 0, "Z(ia,me)", 0, outfile);
  dpd_contract121(&T, &L1, &ZZ, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_oe_file_close(&L1);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_contract121(&T, &L1, &ZZ, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&ZZ);
  /* - P(ij) P(ab) ZZ(i,a) T(j,b) --> G(ij,ab) */
  dpd_oe_file_init(&ZZ, CC_TMP8, 0, 1, "ZZ(i,a)", 0, outfile);
  dpd_oe_file_mat_init(&ZZ);
  dpd_oe_file_mat_rd(&ZZ, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 2, 7, 2, 7, 0, "Gijab", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);
      for(row=0; row < G.params->rowtot[h]; row++) {
	  i = G.params->roworb[h][row][0];
	  j = G.params->roworb[h][row][1];
	  
	  for(col=0; col < G.params->coltot[h]; col++) {
	      a = G.params->colorb[h][col][0];
	      b = G.params->colorb[h][col][1];

	      value = 0.0;

	      I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	      J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];
	      A = ZZ.params->colidx[a]; Asym = ZZ.params->qsym[a];
	      B = T1.params->colidx[b]; Bsym = T1.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  value += ZZ.matrix[Isym][I][A] * T1.matrix[Jsym][J][B];

	      I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	      J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];

	      if((Jsym==Asym) && (Isym==Bsym))
		  value -= ZZ.matrix[Jsym][J][A] * T1.matrix[Isym][I][B];

	      I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	      J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];
	      A = T1.params->colidx[a]; Asym = T1.params->qsym[a];
	      B = ZZ.params->colidx[b]; Bsym = ZZ.params->qsym[b];

	      if((Isym==Bsym) && (Jsym==Asym))
		  value -= ZZ.matrix[Isym][I][B] * T1.matrix[Jsym][J][A];

	      I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	      J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];

	      if((Isym==Asym) && (Jsym==Bsym))
		  value += T1.matrix[Isym][I][A] * ZZ.matrix[Jsym][J][B];

	      G.matrix[h][row][col] -= value;
	      
	    }
	}
      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  dpd_buf_close(&G);

  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);
  dpd_oe_file_mat_close(&ZZ);
  dpd_oe_file_close(&ZZ);

  /* T(j,e) L(m,e) --> ZZ(j,m) */
  dpd_oe_file_init(&ZZ, CC_TMP8, 0, 0, "Z(j,m)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_contract111(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&T1);
  /* ZZ(j,m) T(m,b) --> ZZ2(j,b) */
  dpd_oe_file_init(&ZZ2, CC_TMP9, 0, 1, "ZZ2(j,b)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract111(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_close(&ZZ);
  dpd_oe_file_close(&ZZ2);

  /* 3 P(ij) P(ab) T(i,a) ZZ(j,b) --> G(ij,ab) */
  dpd_oe_file_init(&ZZ, CC_TMP9, 0, 1, "ZZ2(j,b)", 0, outfile);
  dpd_oe_file_mat_init(&ZZ);
  dpd_oe_file_mat_rd(&ZZ, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 2, 7, 2, 7, 0, "Gijab", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);
      for(row=0; row < G.params->rowtot[h]; row++) {
	  i = G.params->roworb[h][row][0];
	  j = G.params->roworb[h][row][1];
	  
	  for(col=0; col < G.params->coltot[h]; col++) {
	      a = G.params->colorb[h][col][0];
	      b = G.params->colorb[h][col][1];

	      value = 0.0;

	      I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	      J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];
	      A = T1.params->colidx[a]; Asym = T1.params->qsym[a];
	      B = ZZ.params->colidx[b]; Bsym = ZZ.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  value += T1.matrix[Isym][I][A] * ZZ.matrix[Jsym][J][B];

	      I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	      J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];

	      if((Jsym==Asym) && (Isym==Bsym))
		  value -= T1.matrix[Jsym][J][A] * ZZ.matrix[Isym][I][B];

	      I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	      J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];
	      A = ZZ.params->colidx[a]; Asym = ZZ.params->qsym[a];
	      B = T1.params->colidx[b]; Bsym = T1.params->qsym[b];

	      if((Isym==Bsym) && (Jsym==Asym))
		  value -= T1.matrix[Isym][I][B] * ZZ.matrix[Jsym][J][A];

	      I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	      J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];

	      if((Isym==Asym) && (Jsym==Bsym))
		  value += ZZ.matrix[Isym][I][A] * T1.matrix[Jsym][J][B];

	      G.matrix[h][row][col] += 3.0 * value;
	      
	    }
	}
      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  dpd_scm(&G, 0.5, 0, outfile);
  dpd_buf_close(&G);

  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);
  dpd_oe_file_mat_close(&ZZ);
  dpd_oe_file_close(&ZZ);



  /* L(Ij,Ab) */
  dpd_buf_init(&L, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_copy(&L, CC_GAMMA, "GIjAb", 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  /* Tau(Ij,Ab) */
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  dpd_axpy(&T, &G, 1.0, 0, outfile);
  dpd_buf_close(&T);
  /* V(Ij,Mn) Tau(Mn,Ab) */
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 0, 0, 0, 0, 0, "VMnIj", 0, outfile);
  dpd_contract222(&V, &T, &G, 0, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  dpd_buf_close(&G);
  /* - ( Z(I,M) Tau(Mj,Ab) - Z(j,m) Tau(mI,bA) ) */
  dpd_buf_init(&Z1, CC_TMP8, 0, 5, 0, 5, 0, "Z1(Ij,Ab)", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  dpd_oe_file_init(&ZZ, CC_TMP0, 0, 0, "Z(I,M)", 0, outfile);
  dpd_contract212(&ZZ, &T, &Z1, 1, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&ZZ);
  dpd_buf_close(&T);
  dpd_buf_init(&Z2, CC_TMP9, 0, 5, 0, 5, 0, "Z2(jI,bA)", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauiJaB", 0, outfile);
  dpd_oe_file_init(&ZZ, CC_TMP1, 0, 0, "Z(i,m)", 0, outfile);
  dpd_contract212(&ZZ, &T, &Z2, 1, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&ZZ);
  dpd_buf_close(&T);
  dpd_swap12(&Z2, CC_TMP10, 0, 5, "Z2(Ij,bA)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP10, 0, 5, 0, 5, 0, "Z2(Ij,bA)", 0, outfile);
  dpd_swap34(&Z2, CC_TMP9, 0, 5, "Z2(Ij,Ab)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP9, 0, 5, 0, 5, 0, "Z2(Ij,Ab)", 0, outfile);
  dpd_axpy(&Z2, &Z1, 1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  dpd_axpy(&Z1, &G, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&G);
  /* - ( Z(E,A) Tau(Ij,bE) - Z(e,b) Tau(Ij,Ae) ) */
  dpd_buf_init(&Z1, CC_TMP8, 0, 5, 0, 5, 0, "Z1(Ij,Ab)", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  dpd_oe_file_init(&ZZ, CC_TMP3, 1, 1, "Z(e,a)", 0, outfile);
  dpd_contract221(&T, &ZZ, &Z1, 3, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&ZZ);
  dpd_buf_close(&T);
  dpd_buf_init(&Z2, CC_TMP9, 0, 5, 0, 5, 0, "Z2(jI,bA)", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauiJaB", 0, outfile);
  dpd_oe_file_init(&ZZ, CC_TMP2, 1, 1, "Z(E,A)", 0, outfile);
  dpd_contract221(&T, &ZZ, &Z2, 3, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&ZZ);
  dpd_buf_close(&T);
  dpd_swap12(&Z2, CC_TMP10, 0, 5, "Z2(Ij,bA)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP10, 0, 5, 0, 5, 0, "Z2(Ij,bA)", 0, outfile);
  dpd_swap34(&Z2, CC_TMP9, 0, 5, "Z2(Ij,Ab)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z2, CC_TMP9, 0, 5, 0, 5, 0, "Z2(Ij,Ab)", 0, outfile);
  dpd_axpy(&Z2, &Z1, 1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  dpd_axpy(&Z1, &G, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&G);
  /* - P(Ij) P(Ab) ( T'(IA,me) (T2) V(jb,me) + T'(IA,ME) (TMP4) V(jb,ME) ) */
  dpd_buf_init(&Z, CC_TMP8, 10, 10, 10, 10, 0, "Z(IA,jb)", 0, outfile);
  dpd_buf_init(&T, CC_TMP4, 10, 10, 10, 10, 0, "Z(IA,ME)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "ViaJB", 0, outfile);
  dpd_contract222(&T, &V, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "Viajb", 0, outfile);
  dpd_contract222(&T, &V, &Z, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  /* T'(jA,Me) V(Ib,Me) */
  dpd_buf_init(&Z1, CC_TMP9, 10, 10, 10, 10, 0, "Z(jA,Ib)", 0, outfile);
  dpd_buf_init(&T, CC_TMP6, 10, 10, 10, 10, 0, "Z(iA,Me)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIaJb", 0, outfile);
  dpd_contract222(&T, &V, &Z1, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  dpd_swap13(&Z1, CC_TMP10, 10, 10, "Z(IA,jb)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP10, 10, 10, 10, 10, 0, "Z(IA,jb)", 0, outfile);
  dpd_axpy(&Z1, &Z, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  /* T'(Ib,mE) V(jA,mE) */
  dpd_buf_init(&Z1, CC_TMP9, 10, 10, 10, 10, 0, "Z(Ib,jA)", 0, outfile);
  dpd_buf_init(&T, CC_TMP7, 10, 10, 10, 10, 0, "Z(Ia,mE)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "ViAjB", 0, outfile);
  dpd_contract222(&T, &V, &Z1, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  dpd_swap24(&Z1, CC_TMP10, 10, 10, "Z(IA,jb)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP10, 10, 10, 10, 10, 0, "Z(IA,jb)", 0, outfile);
  dpd_axpy(&Z1, &Z, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  /* T'(jb,ME) (T2) V(IA,ME) + T'(jb,me) (TMP5) V(IA,me) */
  dpd_buf_init(&Z1, CC_TMP9, 10, 10, 10, 10, 0, "Z(IA,jb)", 0, outfile);
  dpd_buf_init(&T, CC_TMP5, 10, 10, 10, 10, 0, "Z(ia,me)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIAjb", 0, outfile);
  dpd_contract222(&V, &T, &Z1, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIAJB", 0, outfile);
  dpd_contract222(&V, &T, &Z1, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_close(&T);
  dpd_axpy(&Z1, &Z, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  /* - Z(IA,jb) --> G(Ij,Ab) */
  dpd_swap23(&Z, CC_TMP9, 0, 5, "Z(Ij,Ab)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  dpd_buf_init(&Z, CC_TMP9, 0, 5, 0, 5, 0, "Z(Ij,Ab)", 0, outfile);
  dpd_axpy(&Z, &G, -0.5, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&G);
  /* T'(IA,me) (T2) L(m,e) + T'(IA,ME) (TMP4) L(M,E) --> ZZ(I,A) */
  dpd_oe_file_init(&ZZ, CC_TMP8, 0, 1, "ZZ(I,A)", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_buf_init(&T, CC_TMP4, 10, 10, 10, 10, 0, "Z(IA,ME)", 0, outfile);
  dpd_contract121(&T, &L1, &ZZ, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_oe_file_close(&L1);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_contract121(&T, &L1, &ZZ, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&ZZ);
  /* - ZZ(I,A) T(j,b) --> G(Ij,Ab) */
  dpd_oe_file_init(&ZZ, CC_TMP8, 0, 1, "ZZ(I,A)", 0, outfile);
  dpd_oe_file_mat_init(&ZZ);
  dpd_oe_file_mat_rd(&ZZ, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);
      for(row=0; row < G.params->rowtot[h]; row++) {
	  i = G.params->roworb[h][row][0];
	  j = G.params->roworb[h][row][1];
	  
	  for(col=0; col < G.params->coltot[h]; col++) {
	      a = G.params->colorb[h][col][0];
	      b = G.params->colorb[h][col][1];

	      value = 0.0;

	      I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	      J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];
	      A = ZZ.params->colidx[a]; Asym = ZZ.params->qsym[a];
	      B = T1.params->colidx[b]; Bsym = T1.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  value += ZZ.matrix[Isym][I][A] * T1.matrix[Jsym][J][B];

	      G.matrix[h][row][col] -= value;
	      
	    }
	}
      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  dpd_buf_close(&G);

  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);
  dpd_oe_file_mat_close(&ZZ);
  dpd_oe_file_close(&ZZ);

  /* T'(jb,ME) (T2) L(M,E) + T'(jb,me) (TMP5) L(m,e) --> ZZ(j,b) */
  dpd_oe_file_init(&ZZ, CC_TMP8, 0, 1, "ZZ(j,b)", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_buf_init(&T, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_contract121(&T, &L1, &ZZ, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_oe_file_close(&L1);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_buf_init(&T, CC_TMP5, 10, 10, 10, 10, 0, "Z(ia,me)", 0, outfile);
  dpd_contract121(&T, &L1, &ZZ, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&T);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&ZZ);
  /* - ZZ(j,b) T(I,A) --> G(Ij,Ab) */
  dpd_oe_file_init(&ZZ, CC_TMP8, 0, 1, "ZZ(j,b)", 0, outfile);
  dpd_oe_file_mat_init(&ZZ);
  dpd_oe_file_mat_rd(&ZZ, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);
      for(row=0; row < G.params->rowtot[h]; row++) {
	  i = G.params->roworb[h][row][0];
	  j = G.params->roworb[h][row][1];
	  
	  for(col=0; col < G.params->coltot[h]; col++) {
	      a = G.params->colorb[h][col][0];
	      b = G.params->colorb[h][col][1];

	      value = 0.0;

	      I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	      J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];
	      A = T1.params->colidx[a]; Asym = T1.params->qsym[a];
	      B = ZZ.params->colidx[b]; Bsym = ZZ.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  value += T1.matrix[Isym][I][A] * ZZ.matrix[Jsym][J][B];

	      G.matrix[h][row][col] -= value;
	      
	    }
	}
      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  dpd_buf_close(&G);

  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);
  dpd_oe_file_mat_close(&ZZ);
  dpd_oe_file_close(&ZZ);

  /* T(j,e) L(m,e) --> ZZ(j,m) */
  dpd_oe_file_init(&ZZ, CC_TMP8, 0, 0, "Z(j,m)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_contract111(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&T1);
  /* ZZ(j,m) T(m,b) --> ZZ2(j,b) */
  dpd_oe_file_init(&ZZ2, CC_TMP9, 0, 1, "ZZ2(j,b)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract111(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_close(&ZZ);
  dpd_oe_file_close(&ZZ2);

  /* 3 T(I,A) ZZ(j,b) --> G(Ij,Ab) */
  dpd_oe_file_init(&ZZ, CC_TMP9, 0, 1, "ZZ2(j,b)", 0, outfile);
  dpd_oe_file_mat_init(&ZZ);
  dpd_oe_file_mat_rd(&ZZ, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);
      for(row=0; row < G.params->rowtot[h]; row++) {
	  i = G.params->roworb[h][row][0];
	  j = G.params->roworb[h][row][1];
	  
	  for(col=0; col < G.params->coltot[h]; col++) {
	      a = G.params->colorb[h][col][0];
	      b = G.params->colorb[h][col][1];

	      value = 0.0;

	      I = T1.params->rowidx[i]; Isym = T1.params->psym[i];
	      J = ZZ.params->rowidx[j]; Jsym = ZZ.params->psym[j];
	      A = T1.params->colidx[a]; Asym = T1.params->qsym[a];
	      B = ZZ.params->colidx[b]; Bsym = ZZ.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  value += T1.matrix[Isym][I][A] * ZZ.matrix[Jsym][J][B];

	      G.matrix[h][row][col] += 3.0 * value;
	      
	    }
	}
      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  dpd_buf_close(&G);

  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);
  dpd_oe_file_mat_close(&ZZ);
  dpd_oe_file_close(&ZZ);

  /* T(I,E) L(M,E) --> ZZ(I,M) */
  dpd_oe_file_init(&ZZ, CC_TMP8, 0, 0, "Z(I,M)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_contract111(&T1, &L1, &ZZ, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&T1);

  /* ZZ(I,M) T(M,A) --> ZZ2(I,A) */
  dpd_oe_file_init(&ZZ2, CC_TMP9, 0, 1, "ZZ2(I,A)", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract111(&ZZ, &T1, &ZZ2, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_close(&ZZ);
  dpd_oe_file_close(&ZZ2);

  /* 3 T(j,b) ZZ(I,A) --> G(Ij,Ab) */
  dpd_oe_file_init(&ZZ, CC_TMP9, 0, 1, "ZZ2(I,A)", 0, outfile);
  dpd_oe_file_mat_init(&ZZ);
  dpd_oe_file_mat_rd(&ZZ, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 0, 5, 0, 5, 0, "GIjAb", 0, outfile);
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);
      for(row=0; row < G.params->rowtot[h]; row++) {
	  i = G.params->roworb[h][row][0];
	  j = G.params->roworb[h][row][1];
	  
	  for(col=0; col < G.params->coltot[h]; col++) {
	      a = G.params->colorb[h][col][0];
	      b = G.params->colorb[h][col][1];

	      value = 0.0;

	      I = ZZ.params->rowidx[i]; Isym = ZZ.params->psym[i];
	      J = T1.params->rowidx[j]; Jsym = T1.params->psym[j];
	      A = ZZ.params->colidx[a]; Asym = ZZ.params->qsym[a];
	      B = T1.params->colidx[b]; Bsym = T1.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
		  value += ZZ.matrix[Isym][I][A] * T1.matrix[Jsym][J][B];

	      G.matrix[h][row][col] += 3.0 * value;
	      
	    }
	}
      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  dpd_scm(&G, 0.5, 0, outfile);
  dpd_buf_close(&G);

  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);
  dpd_oe_file_mat_close(&ZZ);
  dpd_oe_file_close(&ZZ);
  
  
}
