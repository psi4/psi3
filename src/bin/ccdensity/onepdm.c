#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void onepdm(void)
{
  struct oe_dpdfile D, T1, L1, Z;
  struct dpdbuf T2, L2;
  double trace=0.0;

  dpd_oe_file_init(&D, CC_OEI, 0, 0, "DIJ", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 7, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 7, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_contract122(&T2, &L2, &D, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&T2); 
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_contract122(&T2, &L2, &D, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&T2); 
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_contract111(&T1, &L1, &D, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&T1);
  trace += dpd_oe_trace(&D, 0, outfile);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 0, 0, "Dij", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 7, 2, 7, 0, "tijab", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 7, 2, 7, 0, "Lijab", 0, outfile);
  dpd_contract122(&T2, &L2, &D, 0, 0, -1.0, 0.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LiJaB", 0, outfile);
  dpd_contract122(&T2, &L2, &D, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&T2);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_contract111(&T1, &L1, &D, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&T1);
  trace += dpd_oe_trace(&D, 0, outfile); 
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 1, 1, "DAB", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 2, 5, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 2, 5, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_contract122(&L2, &T2, &D, 3, 3, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&L2);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LiJaB", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_contract122(&L2, &T2, &D, 3, 3, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&T2);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract111(&L1, &T1, &D, 1, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&T1);
  trace += dpd_oe_trace(&D, 0, outfile); 
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&D, CC_OEI, 1, 1, "Dab", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 2, 5, 2, 7, 0, "Lijab", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 2, 5, 2, 7, 0, "tijab", 0, outfile);
  dpd_contract122(&L2, &T2, &D, 3, 3, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&L2);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_contract122(&L2, &T2, &D, 3, 3, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&T2);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract111(&L1, &T1, &D, 1, 1, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_oe_file_close(&T1);
  trace += dpd_oe_trace(&D, 0, outfile); 
  dpd_oe_file_close(&D);

  fprintf(outfile, "\n\tTrace of onepdm = %20.15f\n", trace);

  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_copy(&T1, CC_OEI, "DIA", 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "DIA", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_buf_close(&T2);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_oe_file_init(&Z, CC_TMP0, 0, 0, "Z(I,M)", 0, outfile);
  dpd_contract111(&T1, &L1, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_contract111(&Z, &T1, &D, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_close(&Z);
/*  D(I,A) << L2(MN,EF) T2(IN,EF) T(M,A) + L2(Mn,Ef) T2(In,Ef) T(M,A) */
  dpd_oe_file_init(&Z, CC_TMP0, 0, 0, "Z(I,M)", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 7, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 7, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_contract122(&T2, &L2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&L2);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_contract122(&T2, &L2, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&L2);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract111(&Z, &T1, &D, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&Z);
  dpd_oe_file_close(&T1);
/* T2(MN,AF) L2(MN,EF) T(I,E) + T2(Mn,Af) L2(Mn,Ef) T(I,E) */
  dpd_oe_file_init(&Z, CC_TMP0, 1, 1, "Z(A,E)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 2, 5, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 2, 5, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_contract122(&T2, &L2, &Z, 2, 2, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_contract122(&T2, &L2, &Z, 2, 2, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&T2);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tIA", 0, outfile);
  dpd_contract111(&T1, &Z, &D, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_close(&Z);
  dpd_oe_file_close(&D);

  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_copy(&T1, CC_OEI, "Dia", 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_init(&D, CC_OEI, 0, 1, "Dia", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 2, 7, 0, "tijab", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_dot24(&L1, &T2, &D, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_buf_close(&T2);
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_oe_file_init(&Z, CC_TMP0, 0, 0, "Z(i,m)", 0, outfile);
  dpd_contract111(&T1, &L1, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_oe_file_close(&L1);
  dpd_contract111(&Z, &T1, &D, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_close(&Z);
  dpd_oe_file_init(&Z, CC_TMP0, 0, 0, "Z(i,m)", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 7, 2, 7, 0, "Lijab", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 7, 2, 7, 0, "tijab", 0, outfile);
  dpd_contract122(&T2, &L2, &Z, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&L2);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LiJaB", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_contract122(&T2, &L2, &Z, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&T2);
  dpd_buf_close(&L2);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract111(&Z, &T1, &D, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&Z);
  dpd_oe_file_close(&T1);
  dpd_oe_file_init(&Z, CC_TMP0, 1, 1, "Z(a,e)", 0, outfile);
  dpd_buf_init(&T2, CC_TAMPS, 2, 5, 2, 7, 0, "tijab", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 2, 5, 2, 7, 0, "Lijab", 0, outfile);
  dpd_contract122(&T2, &L2, &Z, 2, 2, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&T2);
  dpd_buf_init(&T2, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LiJaB", 0, outfile);
  dpd_contract122(&T2, &L2, &Z, 2, 2, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_close(&T2);
  dpd_oe_file_init(&T1, CC_OEI, 0, 1, "tia", 0, outfile);
  dpd_contract111(&T1, &Z, &D, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_oe_file_close(&T1);
  dpd_oe_file_close(&Z);
  dpd_oe_file_close(&D);


  /* Note that these blocks are still stored occ/vir */
  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "LIA", 0, outfile);
  dpd_oe_copy(&L1, CC_OEI, "DAI", 0, outfile);
  dpd_oe_file_close(&L1);

  dpd_oe_file_init(&L1, CC_OEI, 0, 1, "Lia", 0, outfile);
  dpd_oe_copy(&L1, CC_OEI, "Dai", 0, outfile);
  dpd_oe_file_close(&L1);
}
