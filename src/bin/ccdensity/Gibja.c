#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

/*
** GIBJA(): Constructs the ibja block of the two-electron density
** matrix, defined based on the energy contribution:
**
**  E(TWO) <-- sum_ibja G(ib,ja) <ja||ib>
**
** or by spin-case:
**
** E(AA) <-- sum_IBJA G(IB,JA) <JA||IB>
** E(BB) <-- sum_ibja G(ib,ja) <ja||ib>
** E(AB) <-- sum_IbJa ( G(Ib,Ja) <Ja|Ib> + G(iB,jA) <jA|iB> -
**                      G(Ib,jA) <jA|bI> - G(iB,Ja) <Ja|Bi> )
**
** The (spin-orbital) equation for G(ib,ja) is:
**
** G(ib,ja) = -L(i,a) T(j,b) - L(im,ae) [ T(jm,be) - T(j,e) T(m,b) ]
**
** I actually build the negative of all the terms above for each spin
** case and then multiply each of these by -1.  Hence the dpd_scm()
** calls you see in the code.
**
** Each Gibja spin-case is built here with the storage G(ia,jb) in
** CC_MISC, but finally resorted to a "proper" G(ib,ja) ordering at
** the end of this routine and placed in CC_GAMMA.  In addition, all
** blocks are "bra-ket symmetrized" (as required by the
** backtransformation) after all contributions have been accounted
** for.  */

void Gibja(void)
{
  int h, nirreps, row, col;
  int i, j, a, b, I, J, A, B, Isym, Jsym, Asym, Bsym;
  struct oe_dpdfile T1, L1, T1A, T1B, L1A, L1B;
  struct dpdbuf G, L, T, Z, Z1, V, G1, G2;

  nirreps = moinfo.nirreps;

  /* G(ia,jb) <-- L(im,ae) T(jm,be) */
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIAJB", 0, outfile);
  dpd_swapbk(&V, CC_MISC, 10, 10, "GIAJB", 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "Viajb", 0, outfile);
  dpd_swapbk(&V, CC_MISC, 10, 10, "Giajb", 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "ViaJB", 0, outfile);
  dpd_swapbk(&V, CC_MISC, 10, 10, "GIAjb", 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIAjb", 0, outfile);
  dpd_swapbk(&V, CC_MISC, 10, 10, "GiaJB", 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "VIaJb", 0, outfile);
  dpd_swapbk(&V, CC_MISC, 10, 10, "GIaJb", 0, outfile);
  dpd_buf_close(&V);
  dpd_buf_init(&V, CC_MISC, 10, 10, 10, 10, 0, "ViAjB", 0, outfile);
  dpd_swapbk(&V, CC_MISC, 10, 10, "GiAjB", 0, outfile);
  dpd_buf_close(&V);

  /* G(IA,JB) <-- - L(IM,AE) T(J,E) T(M,B) */
  dpd_buf_init(&Z, CC_TMP0, 0, 11, 0, 11, 0, "Z(IM,AJ)", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 0, 5, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_init(&Z1, CC_TMP1, 10, 11, 10, 11, 0, "Z(IB,AJ)", 0, outfile);
  dpd_contract221(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_oe_file_close(&T1);
  dpd_swap23(&Z1, CC_TMP0, 10, 11, "Z(IA,BJ)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP0, 10, 11, 10, 11, 0, "Z(IA,BJ)", 0, outfile);
  dpd_swap34(&Z1, CC_TMP1, 10, 10, "Z(IA,JB)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP1, 10, 10, 10, 10, 0, "Z(IA,JB)", 0, outfile);
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GIAJB", 0, outfile);
  dpd_axpy(&Z1, &G, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&G);

  /* G(ia,jb) <-- - L(im,ae) T(j,e) T(m,b) */
  dpd_buf_init(&Z, CC_TMP0, 0, 11, 0, 11, 0, "Z(im,aj)", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 0, 5, 2, 7, 0, "Lijab", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_init(&Z1, CC_TMP1, 10, 11, 10, 11, 0, "Z(ib,aj)", 0, outfile);
  dpd_contract221(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_oe_file_close(&T1);
  dpd_swap23(&Z1, CC_TMP0, 10, 11, "Z(ia,bj)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP0, 10, 11, 10, 11, 0, "Z(ia,bj)", 0, outfile);
  dpd_swap34(&Z1, CC_TMP1, 10, 10, "Z(ia,jb)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP1, 10, 10, 10, 10, 0, "Z(ia,jb)", 0, outfile);
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "Giajb", 0, outfile);
  dpd_axpy(&Z1, &G, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&G);

  /* G(IA,jb) <-- - L(Im,Ae) T(j,e) T(m,b) */
  dpd_buf_init(&Z, CC_TMP0, 0, 11, 0, 11, 0, "Z(Im,Aj)", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract221(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_init(&Z1, CC_TMP1, 10, 11, 10, 11, 0, "Z(Ib,Aj)", 0, outfile);
  dpd_contract221(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_oe_file_close(&T1);
  dpd_swap23(&Z1, CC_TMP0, 10, 11, "Z(IA,bj)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP0, 10, 11, 10, 11, 0, "Z(IA,bj)", 0, outfile);
  dpd_swap34(&Z1, CC_TMP1, 10, 10, "Z(IA,jb)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP1, 10, 10, 10, 10, 0, "Z(IA,jb)", 0, outfile);
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GIAjb", 0, outfile);
  dpd_axpy(&Z1, &G, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&G);

  /* G(ia,JB) <-- - L(iM,aE) T(J,E) T(M,B) */
  dpd_buf_init(&Z, CC_TMP0, 0, 11, 0, 11, 0, "Z(iM,aJ)", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 0, 5, 0, 5, 0, "LiJaB", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract221(&L, &T1, &Z, 3, 1, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_buf_init(&Z1, CC_TMP1, 10, 11, 10, 11, 0, "Z(iB,aJ)", 0, outfile);
  dpd_contract221(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_oe_file_close(&T1);
  dpd_swap23(&Z1, CC_TMP0, 10, 11, "Z(ia,BJ)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP0, 10, 11, 10, 11, 0, "Z(ia,BJ)", 0, outfile);
  dpd_swap34(&Z1, CC_TMP1, 10, 10, "Z(ia,JB)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP1, 10, 10, 10, 10, 0, "Z(ia,JB)", 0, outfile);
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GiaJB", 0, outfile);
  dpd_axpy(&Z1, &G, -1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_close(&G);

  /* G(Ia,Jb) <-- - L(Im,Ea) T(J,E) T(m,b) */
  dpd_buf_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, "Z(Im,Ja)", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract212(&T1, &L, &Z, 1, 2, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_oe_file_close(&T1);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_buf_init(&Z1, CC_TMP1, 10, 10, 10, 10, 0, "Z(Ib,Ja)", 0, outfile);
  dpd_contract221(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_oe_file_close(&T1);
  dpd_swap24(&Z1, CC_TMP0, 10, 10, "Z(Ia,Jb)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP0, 10, 10, 10, 10, 0, "Z(Ia,Jb)", 0, outfile);
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GIaJb", 0, outfile);
  dpd_axpy(&Z1, &G, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_scm(&G, -1.0, 0, outfile);
  dpd_buf_close(&G);

  /* G(iA,jB) <-- - L(iM,eA) T(j,e) T(M,B) */
  dpd_buf_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, "Z(iM,jA)", 0, outfile);
  dpd_buf_init(&L, CC_LAMPS, 0, 5, 0, 5, 0, "LiJaB", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract212(&T1, &L, &Z, 1, 2, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L);
  dpd_oe_file_close(&T1);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_buf_init(&Z1, CC_TMP1, 10, 10, 10, 10, 0, "Z(iB,jA)", 0, outfile);
  dpd_contract221(&Z, &T1, &Z1, 1, 0, 1, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&Z);
  dpd_oe_file_close(&T1);
  dpd_swap24(&Z1, CC_TMP0, 10, 10, "Z(iA,jB)", 0, outfile);
  dpd_buf_close(&Z1);
  dpd_buf_init(&Z1, CC_TMP0, 10, 10, 10, 10, 0, "Z(iA,jB)", 0, outfile);
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GiAjB", 0, outfile);
  dpd_axpy(&Z1, &G, 1.0, 0, outfile);
  dpd_buf_close(&Z1);
  dpd_scm(&G, -1.0, 0, outfile);
  dpd_buf_close(&G);

  dpd_oe_file_init(&T1A, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_mat_init(&T1A);
  dpd_oe_file_mat_rd(&T1A, 0, outfile);
  dpd_oe_file_init(&T1B, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_mat_init(&T1B);
  dpd_oe_file_mat_rd(&T1B, 0, outfile);
  dpd_oe_file_init(&L1A, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_oe_file_mat_init(&L1A);
  dpd_oe_file_mat_rd(&L1A, 0, outfile);
  dpd_oe_file_init(&L1B, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_oe_file_mat_init(&L1B);
  dpd_oe_file_mat_rd(&L1B, 0, outfile);

  /* G(IA,JB) <-- L(I,A) T(J,B) */
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GIAJB", 0, outfile);
  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);

      for(row=0; row < G.params->rowtot[h]; row++) {
          i = G.params->roworb[h][row][0];
          I = L1A.params->rowidx[i]; Isym = L1A.params->psym[i];
          a = G.params->roworb[h][row][1];
          A = L1A.params->colidx[a]; Asym = L1A.params->qsym[a];

          for(col=0; col < G.params->coltot[h]; col++) {
              j = G.params->colorb[h][col][0];
              J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
              b = G.params->colorb[h][col][1];
              B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
              G.matrix[h][row][col] += L1A.matrix[Isym][I][A] * 
                                       T1A.matrix[Jsym][J][B];

            }
        }

      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  dpd_scm(&G, -1.0, 0, outfile);
  dpd_buf_close(&G);


  /* G(ia,jb) <-- L(i,a) T(j,b) */
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "Giajb", 0, outfile);
  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);

      for(row=0; row < G.params->rowtot[h]; row++) {
          i = G.params->roworb[h][row][0];
          I = L1B.params->rowidx[i]; Isym = L1B.params->psym[i];
          a = G.params->roworb[h][row][1];
          A = L1B.params->colidx[a]; Asym = L1B.params->qsym[a];

          for(col=0; col < G.params->coltot[h]; col++) {
              j = G.params->colorb[h][col][0];
              J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
              b = G.params->colorb[h][col][1];
              B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
              G.matrix[h][row][col] += L1B.matrix[Isym][I][A] * 
                                       T1B.matrix[Jsym][J][B];
            }
        }

      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  dpd_scm(&G, -1.0, 0, outfile);
  dpd_buf_close(&G);

  /* G(IA,jb) <-- L(I,A) T(j,b) */
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GIAjb", 0, outfile);
  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);

      for(row=0; row < G.params->rowtot[h]; row++) {
          i = G.params->roworb[h][row][0];
          I = L1A.params->rowidx[i]; Isym = L1A.params->psym[i];
          a = G.params->roworb[h][row][1];
          A = L1A.params->colidx[a]; Asym = L1A.params->qsym[a];

          for(col=0; col < G.params->coltot[h]; col++) {
              j = G.params->colorb[h][col][0];
              J = T1B.params->rowidx[j]; Jsym = T1B.params->psym[j];
              b = G.params->colorb[h][col][1];
              B = T1B.params->colidx[b]; Bsym = T1B.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
              G.matrix[h][row][col] += L1A.matrix[Isym][I][A] * 
                                       T1B.matrix[Jsym][J][B];
            }
        }

      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  dpd_scm(&G, -1.0, 0, outfile);
  dpd_buf_close(&G);

  /* G(ia,JB) <-- L(i,a) T(J,B) */
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GiaJB", 0, outfile);
  for(h=0; h < nirreps; h++) {

      dpd_buf_mat_irrep_init(&G, h);
      dpd_buf_mat_irrep_rd(&G, h, 0, outfile);

      for(row=0; row < G.params->rowtot[h]; row++) {
          i = G.params->roworb[h][row][0];
          I = L1B.params->rowidx[i]; Isym = L1B.params->psym[i];
          a = G.params->roworb[h][row][1];
          A = L1B.params->colidx[a]; Asym = L1B.params->qsym[a];

          for(col=0; col < G.params->coltot[h]; col++) {
              j = G.params->colorb[h][col][0];
              J = T1A.params->rowidx[j]; Jsym = T1A.params->psym[j];
              b = G.params->colorb[h][col][1];
              B = T1A.params->colidx[b]; Bsym = T1A.params->qsym[b];

	      if((Isym==Asym) && (Jsym==Bsym))
              G.matrix[h][row][col] += L1B.matrix[Isym][I][A] * 
                                       T1A.matrix[Jsym][J][B];
            }
        }

      dpd_buf_mat_irrep_wrt(&G, h, 0, outfile);
      dpd_buf_mat_irrep_close(&G, h);
    }
  dpd_scm(&G, -1.0, 0, outfile);
  dpd_buf_close(&G);

  dpd_oe_file_mat_close(&L1A);
  dpd_oe_file_close(&L1A);
  dpd_oe_file_mat_close(&L1B);
  dpd_oe_file_close(&L1B);
  dpd_oe_file_mat_close(&T1A);
  dpd_oe_file_close(&T1A);
  dpd_oe_file_mat_close(&T1B);
  dpd_oe_file_close(&T1B);

  /* Sort all spin cases to correct ordering */
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GIAJB", 0, outfile);
  dpd_buf_sort(&G, CC_GAMMA, psrq, 10, 10, "GIBJA", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "Giajb", 0, outfile);
  dpd_buf_sort(&G, CC_GAMMA, psrq, 10, 10, "Gibja", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GIAjb", 0, outfile);
  dpd_buf_sort(&G, CC_GAMMA, psrq, 10, 10, "GIbjA", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GiaJB", 0, outfile);
  dpd_buf_sort(&G, CC_GAMMA, psrq, 10, 10, "GiBJa", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GIaJb", 0, outfile);
  dpd_buf_sort(&G, CC_GAMMA, psrq, 10, 10, "GIbJa", 0, outfile);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_MISC, 10, 10, 10, 10, 0, "GiAjB", 0, outfile);
  dpd_buf_sort(&G, CC_GAMMA, psrq, 10, 10, "GiBjA", 0, outfile);
  dpd_buf_close(&G);

  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIBJA", 0, outfile);
  dpd_buf_symm(&G);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "Gibja", 0, outfile);
  dpd_buf_symm(&G);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GIbJa", 0, outfile);
  dpd_buf_symm(&G);
  dpd_buf_close(&G);
  dpd_buf_init(&G, CC_GAMMA, 10, 10, 10, 10, 0, "GiBjA", 0, outfile);
  dpd_buf_symm(&G);
  dpd_buf_close(&G);
  dpd_buf_init(&G1, CC_GAMMA, 10, 10, 10, 10, 0, "GIbjA", 0, outfile);
  dpd_buf_init(&G2, CC_GAMMA, 10, 10, 10, 10, 0, "GiBJa", 0, outfile);
  dpd_buf_symm2(&G1, &G2);
  dpd_buf_close(&G2);
  dpd_buf_sort(&G1, CC_GAMMA, rspq, 10, 10, "GiBJa", 0, outfile);
  dpd_buf_close(&G1);
}
