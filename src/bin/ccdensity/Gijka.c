#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void Gijka(void)
{
  int h, nirreps, i, j, k, a, I, J, K, A, Isym, Jsym, Ksym, Asym, row, col;
  double value;
  struct oe_dpdfile L1, T1, g;
  struct dpdbuf G, V, T, L, Z, Z1, Z2;

  nirreps = moinfo.nirreps;

  dpd_buf_init(&G, CC_GAMMA, 2, 10, 2, 10, 0, "GIJKA", 0, outfile);
  /* - tau(IJ,EA) l(K,E) */
  dpd_buf_init(&T, CC_TAMPS, 2, 5, 2, 7, 0, "tauIJAB", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_contract212(&L1, &T, &G, 1, 2, 1, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_buf_close(&T);
  /* - L(IJ,EA) t(K,E) */
  dpd_buf_init(&L, CC_LAMPS, 2, 5, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract212(&T1, &L, &G, 1, 2, 1, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&L);
  /* V(IJ,KM) t(M,A) */
  dpd_buf_init(&V, CC_MISC, 2, 0, 2, 2, 0, "VMNIJ", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&V, &T1, &G, 3, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&V);
  dpd_buf_close(&G);
  /* V(IA,KF) T(J,F) --> Z(KA,IJ) */
  dpd_buf_init(&Z, CC_TMP0, 10, 0, 10, 0, 0, "Z(IA,KJ)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIAJB", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&V, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&V);
  dpd_swap13(&Z, CC_TMP1, 10, 0, "Z(KA,IJ)", 0, outfile);
  dpd_buf_close(&Z);
  /* Z(KA,IJ) - Z(KA,JI) --> G(IJ,KA) */
  dpd_buf_init(&Z1, CC_TMP1, 10, 0, 10, 0, 0, "Z(KA,IJ)", 0, outfile);
  dpd_swap34(&Z1, CC_TMP0, 10, 0, "Z(KA,JI)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 10, 0, 10, 0, 0, "Z(KA,JI)", 0, outfile);
  dpd_axpy(&Z2, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_swapbk(&Z1, CC_TMP0, 0, 10, "Z(IJ,KA)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 2, 10, 0, "GIJKA", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, "Z(IJ,KA)", 0, outfile);
  dpd_axpy(&Z, &G, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&G);

  /* - ( g(I,K) T(J,A) - g(J,K) T(I,A) ) --> G(IJ,KA) */
  dpd_oe_file_init(&g, CC_OEI, 0, 0, "GMI", 0, outfile);
  dpd_oe_file_mat_init(&g);
  dpd_oe_file_mat_rd(&g, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 2, 10, 2, 10, 0, "GIJKA", 0, outfile);
  
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);

      for(row=0; row < G.params->rowtot[h]; row++) {
	  i = G.params->roworb[h][row][0];
	  j = G.params->roworb[h][row][1];
	  for(col=0; col < G.params->coltot[h]; col++) {
	      k = G.params->colorb[h][col][0];
	      a = G.params->colorb[h][col][1];

	      value = 0.0;

	      I = g.params->rowidx[i];  J = T1.params->rowidx[j];
	      Isym = g.params->psym[i]; Jsym = T1.params->psym[j];
	      K = g.params->colidx[k];  A = T1.params->colidx[a];
	      Ksym = g.params->qsym[k];  Asym = T1.params->qsym[a];
	      
	      if((Isym==Ksym) && (Jsym==Asym))
		  value += g.matrix[Isym][I][K] * T1.matrix[Jsym][J][A];

	      J = g.params->rowidx[j];  I = T1.params->rowidx[i];
	      Jsym = g.params->psym[j]; Isym = T1.params->psym[i];
	      
	      if((Jsym==Ksym) && (Isym==Asym))
		  value -= g.matrix[Jsym][J][K] * T1.matrix[Isym][I][A];

	      G.matrix[h][row][col] -= value;
	    }
	}

      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }

  dpd_scm(&G, 0.5, 0, outfile);
  dpd_buf_close(&G);
  
  dpd_oe_file_mat_close(&g);
  dpd_oe_file_close(&g);
  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);


  dpd_buf_init(&G, CC_GAMMA, 2, 10, 2, 10, 0, "Gijka", 0, outfile);
  /* - tau(ij,ea) l(k,e) */
  dpd_buf_init(&T, CC_TAMPS, 2, 5, 2, 7, 0, "tauijab", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_contract212(&L1, &T, &G, 1, 2, 1, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_buf_close(&T);
  /* -L(ij,ea) t(k,e) */
  dpd_buf_init(&L, CC_LAMPS, 2, 5, 2, 7, 0, "Lijab", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract212(&T1, &L, &G, 1, 2, 1, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&L);
  /* V(ij,km) t(m,a) */
  dpd_buf_init(&V, CC_MISC, 2, 0, 2, 2, 0, "Vmnij", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&V, &T1, &G, 3, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&V);
  dpd_buf_close(&G);
  /* V(ia,kf) T(j,f) --> Z(ka,ij) */
  dpd_buf_init(&Z, CC_TMP0, 10, 0, 10, 0, 0, "Z(ia,kj)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "Viajb", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&V, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&V);
  dpd_swap13(&Z, CC_TMP1, 10, 0, "Z(ka,ij)", 0, outfile);
  dpd_buf_close(&Z);
  /* Z(ka,ij) - Z(ka,ji) --> G(ij,ka) */
  dpd_buf_init(&Z1, CC_TMP1, 10, 0, 10, 0, 0, "Z(ka,ij)", 0, outfile);
  dpd_swap34(&Z1, CC_TMP0, 10, 0, "Z(ka,ji)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 10, 0, 10, 0, 0, "Z(ka,ji)", 0, outfile);
  dpd_axpy(&Z2, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_swapbk(&Z1, CC_TMP0, 0, 10, "Z(ij,ka)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 2, 10, 0, "Gijka", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, "Z(ij,ka)", 0, outfile);
  dpd_axpy(&Z, &G, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&G);

  /* - ( g(i,k) T(j,a) - g(j,k) T(i,a) ) --> G(ij,ka) */
  dpd_oe_file_init(&g, CC_OEI, 0, 0, "Gmi", 0, outfile);
  dpd_oe_file_mat_init(&g);
  dpd_oe_file_mat_rd(&g, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 2, 10, 2, 10, 0, "Gijka", 0, outfile);
  
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);

      for(row=0; row < G.params->rowtot[h]; row++) {
	  i = G.params->roworb[h][row][0];
	  j = G.params->roworb[h][row][1];
	  for(col=0; col < G.params->coltot[h]; col++) {
	      k = G.params->colorb[h][col][0];
	      a = G.params->colorb[h][col][1];

	      value = 0.0;

	      I = g.params->rowidx[i];  J = T1.params->rowidx[j];
	      Isym = g.params->psym[i]; Jsym = T1.params->psym[j];
	      K = g.params->colidx[k];  A = T1.params->colidx[a];
	      Ksym = g.params->qsym[k];  Asym = T1.params->qsym[a];
	      
	      if((Isym==Ksym) && (Jsym==Asym))
		  value += g.matrix[Isym][I][K] * T1.matrix[Jsym][J][A];

	      J = g.params->rowidx[j];  I = T1.params->rowidx[i];
	      Jsym = g.params->psym[j]; Isym = T1.params->psym[i];
	      
	      if((Jsym==Ksym) && (Isym==Asym))
		  value -= g.matrix[Jsym][J][K] * T1.matrix[Isym][I][A];

	      G.matrix[h][row][col] -= value;
	    }
	}

      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }

  dpd_scm(&G, 0.5, 0, outfile);
  dpd_buf_close(&G);
  
  dpd_oe_file_mat_close(&g);
  dpd_oe_file_close(&g);
  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);


  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GIjKa", 0, outfile);
  /* - tau(Ij,Ea) l(K,E) */
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauIjAb", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_contract212(&L1, &T, &G, 1, 2, 1, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_buf_close(&T);
  /* -L(Ij,Ea) t(K,E) */
  dpd_buf_init(&L, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract212(&T1, &L, &G, 1, 2, 1, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&L);
  /* V(Ij,Km) t(m,a) */
  dpd_buf_init(&V, CC_MISC, 0, 0, 0, 0, 0, "VMnIj", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&V, &T1, &G, 3, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&V);
  dpd_buf_close(&G);
  /* V(Ia,Kf) T(j,f) --> Z(Ka,Ij) */
  dpd_buf_init(&Z, CC_TMP0, 10, 0, 10, 0, 0, "Z(Ia,Kj)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIaJb", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&V, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&V);
  dpd_swap13(&Z, CC_TMP1, 10, 0, "Z(Ka,Ij)", 0, outfile);
  dpd_buf_close(&Z);
  /* V(ja,KF) T(I,F) --> Z(Ka,jI) */
  dpd_buf_init(&Z, CC_TMP0, 10, 0, 10, 0, 0, "Z(ja,KI)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "ViaJB", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&V, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&V);
  dpd_swap13(&Z, CC_TMP2, 10, 0, "Z(Ka,jI)", 0, outfile);
  dpd_buf_close(&Z);
  /* Z(Ka,Ij) - Z(Ka,jI) --> G(Ij,Ka) */
  dpd_buf_init(&Z2, CC_TMP2, 10, 0, 10, 0, 0, "Z(Ka,jI)", 0, outfile);
  dpd_swap34(&Z2, CC_TMP0, 10, 0, "Z(Ka,Ij)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z1, CC_TMP1, 10, 0, 10, 0, 0, "Z(Ka,Ij)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 10, 0, 10, 0, 0, "Z(Ka,Ij)", 0, outfile);
  dpd_axpy(&Z2, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_swapbk(&Z1, CC_TMP0, 0, 10, "Z(Ij,Ka)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GIjKa", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, "Z(Ij,Ka)", 0, outfile);
  dpd_axpy(&Z, &G, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&G);

  /* - g(I,K) T(j,a) --> G(Ij,Ka) */
  dpd_oe_file_init(&g, CC_OEI, 0, 0, "GMI", 0, outfile);
  dpd_oe_file_mat_init(&g);
  dpd_oe_file_mat_rd(&g, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GIjKa", 0, outfile);
  
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);

      for(row=0; row < G.params->rowtot[h]; row++) {
	  i = G.params->roworb[h][row][0];
	  j = G.params->roworb[h][row][1];
	  for(col=0; col < G.params->coltot[h]; col++) {
	      k = G.params->colorb[h][col][0];
	      a = G.params->colorb[h][col][1];

	      value = 0.0;

	      I = g.params->rowidx[i];  J = T1.params->rowidx[j];
	      Isym = g.params->psym[i]; Jsym = T1.params->psym[j];
	      K = g.params->colidx[k];  A = T1.params->colidx[a];
	      Ksym = g.params->qsym[k];  Asym = T1.params->qsym[a];
	      
	      if((Isym==Ksym) && (Jsym==Asym))
		  value += g.matrix[Isym][I][K] * T1.matrix[Jsym][J][A];

	      G.matrix[h][row][col] -= value;
	    }
	}

      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }

  dpd_scm(&G, 0.5, 0, outfile);
  dpd_buf_close(&G);
  
  dpd_oe_file_mat_close(&g);
  dpd_oe_file_close(&g);
  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);


  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GiJkA", 0, outfile);
  /* - tau(iJ,eA) l(k,e) */
  dpd_buf_init(&T, CC_TAMPS, 0, 5, 0, 5, 0, "tauiJaB", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_contract212(&L1, &T, &G, 1, 2, 1, -1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_buf_close(&T);
  /* -L(iJ,eA) t(k,e) */
  dpd_buf_init(&L, CC_LAMPS, 0, 5, 0, 5, 0, "LiJaB", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract212(&T1, &L, &G, 1, 2, 1, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&L);
  dpd_buf_close(&G);
  /* V(iJ,kM) t(M,A) */
  dpd_buf_init(&Z, CC_TMP0, 0, 11, 0, 11, 0, "Z(Ji,Ak)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 0, 0, 0, 0, 0, "VMnIj", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract212(&T1, &V, &Z, 0, 2, 1, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&V);
  dpd_swap12(&Z, CC_TMP1, 0, 11, "Z(iJ,Ak)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&Z, CC_TMP1, 0, 11, 0, 11, 0, "Z(iJ,Ak)", 0, outfile);
  dpd_swap34(&Z, CC_TMP0, 0, 10, "Z(iJ,kA)", 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GiJkA", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, "Z(iJ,kA)", 0, outfile);
  dpd_axpy(&Z, &G, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&G);
  /* V(iA,kF) T(J,F) --> Z(kA,iJ) */
  dpd_buf_init(&Z, CC_TMP0, 10, 0, 10, 0, 0, "Z(iA,kJ)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "ViAjB", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&V, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&V);
  dpd_swap13(&Z, CC_TMP1, 10, 0, "Z(kA,iJ)", 0, outfile);
  dpd_buf_close(&Z);
  /* V(iA,kf) T(i,f) --> Z(kA,Ji) */
  dpd_buf_init(&Z, CC_TMP0, 10, 0, 10, 0, 0, "Z(JA,ki)", 0, outfile);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIAjb", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&V, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_buf_close(&V);
  dpd_swap13(&Z, CC_TMP2, 10, 0, "Z(kA,Ji)", 0, outfile);
  dpd_buf_close(&Z);
  /* Z(kA,iJ) - Z(kA,Ji) --> G(iJ,kA) */
  dpd_buf_init(&Z2, CC_TMP2, 10, 0, 10, 0, 0, "Z(kA,Ji)", 0, outfile);
  dpd_swap34(&Z2, CC_TMP0, 10, 0, "Z(kA,iJ)", 0, outfile);
  dpd_buf_close(&Z2);
  dpd_buf_init(&Z1, CC_TMP1, 10, 0, 10, 0, 0, "Z(kA,iJ)", 0, outfile);
  dpd_buf_init(&Z2, CC_TMP0, 10, 0, 10, 0, 0, "Z(kA,iJ)", 0, outfile);
  dpd_axpy(&Z2, &Z1, -1.0, 0, outfile);
  dpd_buf_close(&Z2);
  dpd_swapbk(&Z1, CC_TMP0, 0, 10, "Z(iJ,kA)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GiJkA", 0, outfile);
  dpd_buf_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, "Z(iJ,kA)", 0, outfile);
  dpd_axpy(&Z, &G, 1.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_buf_close(&G);

  /* - g(i,k) T(J,A) --> G(iJ,kA) */
  dpd_oe_file_init(&g, CC_OEI, 0, 0, "Gmi", 0, outfile);
  dpd_oe_file_mat_init(&g);
  dpd_oe_file_mat_rd(&g, 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1);
  dpd_oe_file_mat_rd(&T1, 0, outfile);

  dpd_buf_init(&G, CC_GAMMA, 0, 10, 0, 10, 0, "GiJkA", 0, outfile);
  
  for(h=0; h < nirreps; h++) {
      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);

      for(row=0; row < G.params->rowtot[h]; row++) {
	  i = G.params->roworb[h][row][0];
	  j = G.params->roworb[h][row][1];
	  for(col=0; col < G.params->coltot[h]; col++) {
	      k = G.params->colorb[h][col][0];
	      a = G.params->colorb[h][col][1];

	      value = 0.0;

	      I = g.params->rowidx[i];  J = T1.params->rowidx[j];
	      Isym = g.params->psym[i]; Jsym = T1.params->psym[j];
	      K = g.params->colidx[k];  A = T1.params->colidx[a];
	      Ksym = g.params->qsym[k];  Asym = T1.params->qsym[a];
	      
	      if((Isym==Ksym) && (Jsym==Asym))
		  value += g.matrix[Isym][I][K] * T1.matrix[Jsym][J][A];

	      G.matrix[h][row][col] -= value;
	    }
	}

      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }

  dpd_scm(&G, 0.5, 0, outfile);
  dpd_buf_close(&G);
  
  dpd_oe_file_mat_close(&g);
  dpd_oe_file_close(&g);
  dpd_oe_file_mat_close(&T1);
  dpd_oe_file_close(&T1);
  
}

