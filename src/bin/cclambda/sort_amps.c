#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void sort_amps(void)
{
  dpdbuf4 L2;

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    /* Build L2iJaB list */
    dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_LAMPS, qpsr, 0, 5, "LiJaB");
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
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&L2, CC_LAMPS, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_LAMPS, qpsr, 23, 29, "LiJaB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 0, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_sort(&L2, CC_LAMPS, prqs, 20, 20, "LIAJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 10, 15, 12, 17, 0, "Lijab");
    dpd_buf4_sort(&L2, CC_LAMPS, prqs, 30, 30, "Liajb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_LAMPS, prqs, 20, 30, "LIAjb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_sort(&L2, CC_LAMPS, prqs, 30, 20, "LiaJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, 0, 20, 30, 20, 30, 0, "LIAjb");
    dpd_buf4_sort(&L2, CC_LAMPS, psrq, 24, 27, "LIbjA");
    dpd_buf4_sort(&L2, CC_LAMPS, rqps, 27, 24, "LjAIb");
    dpd_buf4_close(&L2);
  }

}

