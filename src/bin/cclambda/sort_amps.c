#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void sort_amps(void)
{
  dpdbuf4 L2;

  /* Build L2iJaB list */
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_buf4_sort(&L2, CC_TMP0, qprs, 0, 5, "LjIAb");
  dpd_buf4_close(&L2);
  dpd_buf4_init(&L2, CC_TMP0, 0, 0, 5, 0, 5, 0, "LjIAb");
  dpd_buf4_sort(&L2, CC_LAMPS, pqsr,0, 5, "LiJaB");
  dpd_buf4_close(&L2);

  /* Build L2IAJB List */
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 2, 7, 0, "LIJAB");
  dpd_buf4_sort(&L2, CC_LAMPS, prqs, 10, 10, "LIAJB");
  dpd_buf4_close(&L2);

  /* Build L2iajb List */
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 2, 7, 0, "Lijab");
  dpd_buf4_sort(&L2, CC_LAMPS, prqs, 10, 10, "Liajb");
  dpd_buf4_close(&L2);

  /* Build L2IAjb List */
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
  dpd_buf4_sort(&L2, CC_LAMPS, prqs, 10, 10, "LIAjb");
  dpd_buf4_close(&L2);

  /* Build L2iaJB List */
  dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LiJaB");
  dpd_buf4_sort(&L2, CC_LAMPS, prqs, 10, 10, "LiaJB");
  dpd_buf4_close(&L2);

  /* Build L2IbjA and L2 jAIb Lists */
  dpd_buf4_init(&L2, CC_LAMPS, 0, 10, 10, 10, 10, 0, "LIAjb");
  dpd_buf4_sort(&L2, CC_LAMPS, psrq, 10, 10, "LIbjA");
  dpd_buf4_sort(&L2, CC_LAMPS, rqps, 10, 10, "LjAIb");
  dpd_buf4_close(&L2);

}

