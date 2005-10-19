#include <math.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void rhf_lag(void);
void uhf_lag(void);

void lag(void)
{
  if(params.ref == 0) rhf_lag();
  else if(params.ref == 2) uhf_lag();
}

void rhf_lag(void)
{
  dpdfile2 D; 
  dpdfile2 L; 
  dpdbuf4 I;  
  dpdbuf4 T;  

  dpd_file2_init(&L, CC_OEI, 0, 1, 0, "LAI");
  dpd_buf4_init(&I, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
  dpd_file2_init(&D, CC_OEI, 0, 0, 0, "DIJ");
  dpd_dot24(&D, &I, &L, 0, 0, 2.0, 0.0); 
  dpd_dot23(&D, &I, &L, 0, 0, -1.0, 1.0);
  dpd_buf4_close(&I);
  dpd_file2_close(&D);
  dpd_buf4_init(&I, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
  dpd_file2_init(&D, CC_OEI, 0, 1, 1, "DAB");
  dpd_dot24(&D, &I, &L, 0, 1, 2.0, 1.0);
  dpd_dot23(&D, &I, &L, 0, 1, -1.0, 1.0);
  dpd_buf4_close(&I);
  dpd_file2_close(&D);
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
  dpd_buf4_init(&I, CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
  dpd_contract442(&T, &I, &L, 2, 0, -2.0, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_init(&I, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
  dpd_contract442(&I, &T, &L, 0, 0, 2.0, 1.0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_file2_close(&L);
}

void uhf_lag(void)
{

}
