#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void G_build(void) {
  struct dpdbuf LIJAB, Lijab, LiJaB, LIjAb, LijAB, LIJab;
  struct dpdbuf tIJAB, tijab, tiJaB, tIjAb, tijAB, tIJab;
  struct oe_dpdfile GAE, Gae, GMI, Gmi;

  dpd_oe_file_init(&GMI, CC_OEI, 0, 0, "GMI", 0, outfile);
  dpd_oe_file_init(&Gmi, CC_OEI, 0, 0, "Gmi", 0, outfile);

  /* T2(MJ,AB) * L2(IJ,AB) --> G(M,I) */
  dpd_buf_init(&tIJAB, CC_TAMPS, 0, 7, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_buf_init(&LIJAB, CC_LAMPS, 0, 7, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_contract122(&tIJAB, &LIJAB, &GMI, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&tIJAB);
  dpd_buf_close(&LIJAB);  

  /* T2(Mj,Ab) * L2(Ij,Ab) --> G(m,i) */
  dpd_buf_init(&tIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_buf_init(&LIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_contract122(&tIjAb, &LIjAb, &GMI, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&tIjAb);
  dpd_buf_close(&LIjAb);

  /* T2(mj,ab) * L2(ij,ab) --> G(m,i) */
  dpd_buf_init(&tijab, CC_TAMPS, 0, 7, 2, 7, 0, "tijab", 0, outfile);
  dpd_buf_init(&Lijab, CC_LAMPS, 0, 7, 2, 7, 0, "Lijab", 0, outfile);
  dpd_contract122(&tijab, &Lijab, &Gmi, 0, 0, 1.0, 0.0, 0, outfile);
  dpd_buf_close(&tijab);
  dpd_buf_close(&Lijab); 

  /* T2(mJ,aB) * L2(iJ,aB) --> G(m,i) */
  dpd_buf_init(&tiJaB, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_buf_init(&LiJaB, CC_LAMPS, 0, 5, 0, 5, 0, "LiJaB", 0, outfile);
  dpd_contract122(&tiJaB, &LiJaB, &Gmi, 0, 0, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&tiJaB);
  dpd_buf_close(&LiJaB);

  dpd_oe_file_close(&Gmi);
  dpd_oe_file_close(&GMI);


  
  dpd_oe_file_init(&GAE, CC_OEI, 1, 1, "GAE", 0, outfile);
  dpd_oe_file_init(&Gae, CC_OEI, 1, 1, "Gae", 0, outfile);

  /* T2(IJ,AB) * L2(IJ,EB) --> G(A,E) */
  dpd_buf_init(&tIJAB, CC_TAMPS, 2, 5, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_buf_init(&LIJAB, CC_LAMPS, 2, 5, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_contract122(&LIJAB, &tIJAB, &GAE, 2, 2, -1.0, 0.0, 0, outfile);
  dpd_buf_close(&tIJAB);
  dpd_buf_close(&LIJAB);

  /* T2(Ij,Ab) * L2(Ij,Eb) --> G(A,E) */
  dpd_buf_init(&tIjAb, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_buf_init(&LIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_contract122(&LIjAb, &tIjAb, &GAE, 2, 2, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&tIjAb);
  dpd_buf_close(&LIjAb);

  /* T2(ij,ab) * L2(ij,eb) --> G(a,e) */
  dpd_buf_init(&tijab, CC_TAMPS, 2, 5, 2, 7, 0, "tijab", 0, outfile);
  dpd_buf_init(&Lijab, CC_LAMPS, 2, 5, 2, 7, 0, "Lijab", 0, outfile);
  dpd_contract122(&Lijab, &tijab, &Gae, 2, 2, -1.0, 0.0, 0, outfile);
  dpd_buf_close(&tijab);
  dpd_buf_close(&Lijab);

  /* T2(iJ,aB) * L2(iJ,eB) --> G(a,e) */
  dpd_buf_init(&tiJaB, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_buf_init(&LiJaB, CC_LAMPS, 0, 5, 0, 5, 0, "LiJaB", 0, outfile);
  dpd_contract122(&LiJaB, &tiJaB, &Gae, 2, 2, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&tiJaB);
  dpd_buf_close(&LiJaB);

  dpd_oe_file_close(&GAE);
  dpd_oe_file_close(&Gae);

  return;
}

