#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void sort_amps(void)
{
  struct dpdbuf L2;

  /* Build L2iJaB list */
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_swap12(&L2, CC_TMP0, 0, 5, "LjIAb", 0, outfile);
  dpd_buf_close(&L2);
  dpd_buf_init(&L2, CC_TMP0, 0, 5, 0, 5, 0, "LjIAb", 0, outfile);
  dpd_swap34(&L2, CC_LAMPS, 0, 5, "LiJaB", 0, outfile);
  dpd_buf_close(&L2);

  /* Build L2IAJB List */
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 2, 7, 0, "LIJAB", 0, outfile);
  dpd_swap23(&L2, CC_LAMPS, 10, 10, "LIAJB", 0, outfile);
  dpd_buf_close(&L2);

  /* Build L2iajb List */
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 2, 7, 0, "Lijab", 0, outfile);
  dpd_swap23(&L2, CC_LAMPS, 10, 10, "Liajb", 0, outfile);
  dpd_buf_close(&L2);

  /* Build L2IAjb List */
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LIjAb", 0, outfile);
  dpd_swap23(&L2, CC_LAMPS, 10, 10, "LIAjb", 0, outfile);
  dpd_buf_close(&L2);

  /* Build L2iaJB List */
  dpd_buf_init(&L2, CC_LAMPS, 0, 5, 0, 5, 0, "LiJaB", 0, outfile);
  dpd_swap23(&L2, CC_LAMPS, 10, 10, "LiaJB", 0, outfile);
  dpd_buf_close(&L2);

  /* Build L2IbjA and L2 jAIb Lists */
  dpd_buf_init(&L2, CC_LAMPS, 10, 10, 10, 10, 0, "LIAjb", 0, outfile);
  dpd_swap24(&L2, CC_LAMPS, 10, 10, "LIbjA", 0, outfile);
  dpd_swap13(&L2, CC_LAMPS, 10, 10, "LjAIb", 0, outfile);
  dpd_buf_close(&L2);

}

