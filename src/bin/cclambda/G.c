#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void G_build(void) {
  dpdbuf4 LIJAB, Lijab, LiJaB, LIjAb, LijAB, LIJab;
  dpdbuf4 tIJAB, tijab, tiJaB, tIjAb, tijAB, tIJab;
  dpdfile2 GAE, Gae, GMI, Gmi;

  dpd_file2_init(&GMI, CC_OEI, 0, 0, 0, "GMI");
  dpd_file2_init(&Gmi, CC_OEI, 0, 0, 0, "Gmi");

  /* T2(MJ,AB) * L2(IJ,AB) --> G(M,I) */
  dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 0, 7, 2, 7, 0, "LIJAB");
  dpd_contract442(&tIJAB, &LIJAB, &GMI, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&tIJAB);
  dpd_buf4_close(&LIJAB);  

  /* T2(Mj,Ab) * L2(Ij,Ab) --> G(m,i) */
  dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&LIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract442(&tIjAb, &LIjAb, &GMI, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&tIjAb);
  dpd_buf4_close(&LIjAb);

  /* T2(mj,ab) * L2(ij,ab) --> G(m,i) */
  dpd_buf4_init(&tijab, CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
  dpd_buf4_init(&Lijab, CC_LAMPS, 0, 0, 7, 2, 7, 0, "Lijab");
  dpd_contract442(&tijab, &Lijab, &Gmi, 0, 0, 1.0, 0.0);
  dpd_buf4_close(&tijab);
  dpd_buf4_close(&Lijab); 

  /* T2(mJ,aB) * L2(iJ,aB) --> G(m,i) */
  dpd_buf4_init(&tiJaB, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_buf4_init(&LiJaB, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LiJaB");
  dpd_contract442(&tiJaB, &LiJaB, &Gmi, 0, 0, 1.0, 1.0);
  dpd_buf4_close(&tiJaB);
  dpd_buf4_close(&LiJaB);

  dpd_file2_close(&Gmi);
  dpd_file2_close(&GMI);


  
  dpd_file2_init(&GAE, CC_OEI, 0, 1, 1, "GAE");
  dpd_file2_init(&Gae, CC_OEI, 0, 1, 1, "Gae");

  /* T2(IJ,AB) * L2(IJ,EB) --> G(A,E) */
  dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&LIJAB, CC_LAMPS, 0, 2, 5, 2, 7, 0, "LIJAB");
  dpd_contract442(&LIJAB, &tIJAB, &GAE, 2, 2, -1.0, 0.0);
  dpd_buf4_close(&tIJAB);
  dpd_buf4_close(&LIJAB);

  /* T2(Ij,Ab) * L2(Ij,Eb) --> G(A,E) */
  dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&LIjAb, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_contract442(&LIjAb, &tIjAb, &GAE, 2, 2, -1.0, 1.0);
  dpd_buf4_close(&tIjAb);
  dpd_buf4_close(&LIjAb);

  /* T2(ij,ab) * L2(ij,eb) --> G(a,e) */
  dpd_buf4_init(&tijab, CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
  dpd_buf4_init(&Lijab, CC_LAMPS, 0, 2, 5, 2, 7, 0, "Lijab");
  dpd_contract442(&Lijab, &tijab, &Gae, 2, 2, -1.0, 0.0);
  dpd_buf4_close(&tijab);
  dpd_buf4_close(&Lijab);

  /* T2(iJ,aB) * L2(iJ,eB) --> G(a,e) */
  dpd_buf4_init(&tiJaB, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
  dpd_buf4_init(&LiJaB, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LiJaB");
  dpd_contract442(&LiJaB, &tiJaB, &Gae, 2, 2, -1.0, 1.0);
  dpd_buf4_close(&tiJaB);
  dpd_buf4_close(&LiJaB);

  dpd_file2_close(&GAE);
  dpd_file2_close(&Gae);

  return;
}

