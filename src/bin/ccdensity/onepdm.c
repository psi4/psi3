#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

/* onepdm(): Computes the one-particle density matrix for CC wave functions.
** 
** The spin-orbital expressions for the onepdm components are:
**
** D_ij = -1/2 t_im^ef L^jm_ef - t_i^e L^j_e
**
** D_ab = 1/2 L^mn_ae t_mn^be + L^m_a t_m^b
**
** D_ia = t_i^a + (t_im^ae - t_i^e t_m^a) L^m_e 
**        - 1/2 L^mn_ef (t_in^ef t_m^a + t_i^e t_mn^af)
** 
** D_ai = L^i_a
**
** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).]
**
** TDC, July 2002
*/

void onepdm(void)
{
  dpdfile2 D, T1, L1, Z;
  dpdbuf4 T2, L2;
  double trace=0.0;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/

    dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D); 
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D); 
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D); 
    dpd_file2_close(&D);

    fprintf(outfile, "\n\tTrace of onepdm = %20.15f\n", trace);

    /* this term does not include l, so we must multiply by R0 explicitly */
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_copy(&T1, CC_OEI, "DIA");
    dpd_file2_close(&T1);

    dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
    if (!params.ground) {
      //dpd_file2_scm(&D, params.R0);
      dpd_file2_scm(&D, 0.0);
    }

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */
    dpd_file2_init(&Z, CC_TMP0, 0, 1, 1, "Z(A,E)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_copy(&T1, CC_OEI, "Dia");
    dpd_file2_close(&T1);

    dpd_file2_init(&D, CC_OEI, 0, 0, 1, "Dia");
    if (!params.ground) {
      //dpd_file2_scm(&D, params.R0);
      dpd_file2_scm(&D, 0.0);
    }

    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(i,m)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(i,m)");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    dpd_file2_init(&Z, CC_TMP0, 0, 1, 1, "Z(a,e)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "Lijab");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 5, 0, 5, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    /* Note that these blocks are still stored occ/vir */
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_OEI, "DAI");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "Lia");
    dpd_file2_copy(&L1, CC_OEI, "Dai");
    dpd_file2_close(&L1);

  }
  else if(params.ref == 2) { /** UHF **/

    dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2); 
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D);
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 2, 2, "Dij");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_buf4_init(&L2, CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D); 
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D); 
    dpd_file2_close(&D);

    dpd_file2_init(&D, CC_OEI, 0, 3, 3, "Dab");
    dpd_buf4_init(&L2, CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_file2_close(&T1);
    trace += dpd_file2_trace(&D); 
    dpd_file2_close(&D);

    fprintf(outfile, "\n\tTrace of onepdm = %20.15f\n", trace);

    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_copy(&T1, CC_OEI, "DIA");
    dpd_file2_close(&T1);
    dpd_file2_init(&D, CC_OEI, 0, 0, 1, "DIA");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    /*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */
    dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
    dpd_buf4_init(&L2, CC_GLG, 0, 0, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    /* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */
    dpd_file2_init(&Z, CC_TMP0, 0, 1, 1, "Z(A,E)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_init(&L2, CC_GLG, 0, 2, 5, 2, 7, 0, "LIJAB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_init(&L2, CC_GLG, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_copy(&T1, CC_OEI, "Dia");
    dpd_file2_close(&T1);
    dpd_file2_init(&D, CC_OEI, 0, 2, 3, "Dia");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
    dpd_file2_close(&L1);
    dpd_buf4_close(&T2);
    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_file2_init(&Z, CC_TMP0, 0, 2, 2, "Z(i,m)");
    dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
    dpd_file2_close(&L1);
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_init(&Z, CC_TMP0, 0, 2, 2, "Z(i,m)");
    dpd_buf4_init(&L2, CC_GLG, 0, 10, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 1.0);
    dpd_buf4_close(&T2);
    dpd_buf4_close(&L2);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
    dpd_file2_close(&Z);
    dpd_file2_close(&T1);
    dpd_file2_init(&Z, CC_TMP0, 0, 3, 3, "Z(a,e)");
    dpd_buf4_init(&T2, CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    dpd_buf4_init(&L2, CC_GLG, 0, 12, 15, 12, 17, 0, "Lijab");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_buf4_init(&T2, CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    dpd_buf4_init(&L2, CC_GLG, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
    dpd_buf4_close(&L2);
    dpd_buf4_close(&T2);
    dpd_file2_init(&T1, CC_OEI, 0, 2, 3, "tia");
    dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
    dpd_file2_close(&T1);
    dpd_file2_close(&Z);
    dpd_file2_close(&D);

    /* Note that these blocks are still stored occ/vir */
    dpd_file2_init(&L1, CC_GLG, 0, 0, 1, "LIA");
    dpd_file2_copy(&L1, CC_OEI, "DAI");
    dpd_file2_close(&L1);

    dpd_file2_init(&L1, CC_GLG, 0, 2, 3, "Lia");
    dpd_file2_copy(&L1, CC_OEI, "Dai");
    dpd_file2_close(&L1);

  }
}
