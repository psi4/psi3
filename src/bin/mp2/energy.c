#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

double rhf_energy(void);
double rohf_energy(void);
double uhf_energy(void);

double energy(void)
{
  if(params.ref == 0) return(rhf_energy());
  else if(params.ref == 1) return(rohf_energy());
  else if(params.ref == 2) return(uhf_energy());
}

double rhf_energy(void) 
{
  double E = 0;
  dpdbuf4 tIjAb;
  dpdbuf4 D;
  
  dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  dpd_buf4_init(&D, CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  E = dpd_buf4_dot(&D, &tIjAb);
  dpd_buf4_close(&tIjAb);
  dpd_buf4_close(&D);

  return(E);
}

double rohf_energy(void)
{
}

double uhf_energy(void)
{
  double E2AA = 0;
  double E2BB = 0;
  double E2AB = 0;
  dpdbuf4 tIJAB;
  dpdbuf4 tijab;
  dpdbuf4 tIjAb;
  dpdbuf4 D;
  
  dpd_buf4_init(&tIJAB, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&D, CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
  E2AA = dpd_buf4_dot(&D, &tIJAB);
  dpd_buf4_close(&D);
  dpd_buf4_close(&tIJAB);

  dpd_buf4_init(&tijab, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
  dpd_buf4_init(&D, CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
  E2BB = dpd_buf4_dot(&D, &tijab);
  dpd_buf4_close(&D);
  dpd_buf4_close(&tijab);
  
  dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  dpd_buf4_init(&D, CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
  E2AB = dpd_buf4_dot(&D, &tIjAb);
  dpd_buf4_close(&D);
  dpd_buf4_close(&tIjAb);

  return(E2AA+E2BB+E2AB);
}

