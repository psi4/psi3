#include <math.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void rhf_lag(void);
void uhf_lag(void);

void lag(void)
{
  if(params.ref == 0) return(rhf_lag());
  else if(params.ref == 2) return(uhf_lag());
}

void rhf_lag(void)
{
  dpdfile2 D; /* One-particle density matrix */
  dpdfile2 L; /* Orbital Lagrangian */
  dpdbuf4 I;  /* Two-electron integral */
  dpdbuf4 T;  /* T2 amplitude */

  dpd_file2_init(&L, CC_OEI, 0, 1, 0, "LAI");
  dpd_buf4_init(&I, CC_EINTS, 0, 0, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
  dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_contract442(&T, &I, &L, 2, 2, -2, 0);
  dpd_buf4_close(&I);
  dpd_buf4_close(&T);
  dpd_file2_close(&L);
}

void uhf_lag(void)
{
  dpdfile2 D; /* One-particle density matrix */
  dpdbuf4 I;  /* Two-electron integral */
  dpdbuf4 T;  /* T2 amplitude */
  
}
