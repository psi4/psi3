/*! \file sort_Ls.c
    \ingroup (CCDENSITY)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void sort_Ls(void)
{
  dpdbuf4 L2;
  int L_irr = 0; /* only ground state for now at least */

  if(params.ref == 0 || params.ref == 1) { /** RHF/ROHF **/
    /* Build L2iJaB list */
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_LAMPS, qpsr, 0, 5, "LiJaB");
    dpd_buf4_close(&L2);

    /* Build L2IAJB List */
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_sort(&L2, CC_LAMPS, prqs, 10, 10, "LIAJB");
    dpd_buf4_close(&L2);

    /* Build L2iajb List */
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 2, 7, 0, "Lijab");
    dpd_buf4_sort(&L2, CC_LAMPS, prqs, 10, 10, "Liajb");
    dpd_buf4_close(&L2);

    /* Build L2IAjb List */
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_LAMPS, prqs, 10, 10, "LIAjb");
    dpd_buf4_close(&L2);

    /* Build L2iaJB List */
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    dpd_buf4_sort(&L2, CC_LAMPS, prqs, 10, 10, "LiaJB");
    dpd_buf4_close(&L2);

    /* Build L2IbjA and L2 jAIb Lists */
    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 10, 10, 10, 10, 0, "LIAjb");
    dpd_buf4_sort(&L2, CC_LAMPS, psrq, 10, 10, "LIbjA");
    dpd_buf4_sort(&L2, CC_LAMPS, rqps, 10, 10, "LjAIb");
    dpd_buf4_close(&L2);
  }
  else if(params.ref == 2) { /** UHF **/

    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_LAMPS, qpsr, 23, 29, "LiJaB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 0, 5, 2, 7, 0, "LIJAB");
    dpd_buf4_sort(&L2, CC_LAMPS, prqs, 20, 20, "LIAJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 10, 15, 12, 17, 0, "Lijab");
    dpd_buf4_sort(&L2, CC_LAMPS, prqs, 30, 30, "Liajb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    dpd_buf4_sort(&L2, CC_LAMPS, prqs, 20, 30, "LIAjb");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    dpd_buf4_sort(&L2, CC_LAMPS, prqs, 30, 20, "LiaJB");
    dpd_buf4_close(&L2);

    dpd_buf4_init(&L2, CC_LAMPS, L_irr, 20, 30, 20, 30, 0, "LIAjb");
    dpd_buf4_sort(&L2, CC_LAMPS, psrq, 24, 27, "LIbjA");
    dpd_buf4_sort(&L2, CC_LAMPS, rqps, 27, 24, "LjAIb");
    dpd_buf4_close(&L2);
  }

}

