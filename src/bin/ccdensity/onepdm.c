#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void onepdm(void)
{
  dpdfile2 D, T1, L1, Z;
  dpdbuf4 T2, L2;
  double trace=0.0;

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 7, 2, 7, 0, "LIJAB");
  dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&T2); 
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&T2); 
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "LIA");
  dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
  dpd_file2_close(&L1);
  dpd_file2_close(&T1);
  trace += dpd_file2_trace(&D);
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "Dij");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 7, 2, 7, 0, "Lijab");
  dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&T2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LiJaB");
  dpd_contract442(&T2, &L2, &D, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&T2);
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "Lia");
  dpd_contract222(&T1, &L1, &D, 0, 0, -1.0, 1.0);
  dpd_file2_close(&L1);
  dpd_file2_close(&T1);
  trace += dpd_file2_trace(&D); 
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 5, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
  dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LiJaB");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&T2);
  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "LIA");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_contract222(&L1, &T1, &D, 1, 1, 1.0, 1.0);
  dpd_file2_close(&L1);
  dpd_file2_close(&T1);
  trace += dpd_file2_trace(&D); 
  dpd_file2_close(&D);

  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "Dab");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 5, 2, 7, 0, "Lijab");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
  dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract442(&L2, &T2, &D, 3, 3, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&T2);
  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "Lia");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
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
  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "LIA");
  dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
  dpd_file2_close(&L1);
  dpd_buf4_close(&T2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "Lia");
  dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
  dpd_file2_close(&L1);
  dpd_buf4_close(&T2);
  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "LIA");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
  dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
  dpd_file2_close(&L1);
  dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_file2_close(&Z);
/*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */
  dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(I,M)");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 7, 2, 7, 0, "LIJAB");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
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
  dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 5, 2, 7, 0, "LIJAB");
  dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&T2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
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
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "Lia");
  dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
  dpd_file2_close(&L1);
  dpd_buf4_close(&T2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "LIA");
  dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0);
  dpd_file2_close(&L1);
  dpd_buf4_close(&T2);
  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "Lia");
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(i,m)");
  dpd_contract222(&T1, &L1, &Z, 0, 0, 1.0, 0.0);
  dpd_file2_close(&L1);
  dpd_contract222(&Z, &T1, &D, 0, 1, -1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_file2_close(&Z);
  dpd_file2_init(&Z, CC_TMP0, 0, 0, 0, "Z(i,m)");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 7, 2, 7, 0, "Lijab");
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
  dpd_contract442(&T2, &L2, &Z, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&T2);
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LiJaB");
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
  dpd_buf4_init(&L2, CC_LAMPS, 0, 2, 5, 2, 7, 0, "Lijab");
  dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 0.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&T2);
  dpd_buf4_init(&T2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LiJaB");
  dpd_contract442(&T2, &L2, &Z, 2, 2, 1.0, 1.0);
  dpd_buf4_close(&L2);
  dpd_buf4_close(&T2);
  dpd_file2_init(&T1, CC_OEI, 0, 0, 1, "tia");
  dpd_contract222(&T1, &Z, &D, 0, 0, -1.0, 1.0);
  dpd_file2_close(&T1);
  dpd_file2_close(&Z);
  dpd_file2_close(&D);


  /* Note that these blocks are still stored occ/vir */
  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "LIA");
  dpd_file2_copy(&L1, CC_OEI, "DAI");
  dpd_file2_close(&L1);

  dpd_file2_init(&L1, CC_OEI, 0, 0, 1, "Lia");
  dpd_file2_copy(&L1, CC_OEI, "Dai");
  dpd_file2_close(&L1);
}
