#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void sort_amps(void)
{
  struct dpdbuf t2;

/*  timer_on("sort amps", outfile); */

  /* Build t2iJaB list */
  dpd_buf_init(&t2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_swap12(&t2, CC_TMP0, 0, 5, "tjIAb", 0, outfile);
  dpd_buf_close(&t2);
  dpd_buf_init(&t2, CC_TMP0, 0, 5, 0, 5, 0, "tjIAb", 0, outfile);
  dpd_swap34(&t2, CC_TAMPS, 0, 5, "tiJaB", 0, outfile);
  dpd_buf_close(&t2);

  /* Build t2IAJB List */
  dpd_buf_init(&t2, CC_TAMPS, 0, 5, 2, 7, 0, "tIJAB", 0, outfile);
  dpd_swap23(&t2, CC_TAMPS, 10, 10, "tIAJB", 0, outfile);
  dpd_buf_close(&t2);

  /* Build t2iajb List */
  dpd_buf_init(&t2, CC_TAMPS, 0, 5, 2, 7, 0, "tijab", 0, outfile);
  dpd_swap23(&t2, CC_TAMPS, 10, 10, "tiajb", 0, outfile);
  dpd_buf_close(&t2);

  /* Build t2IAjb List */
  dpd_buf_init(&t2, CC_TAMPS, 0, 5, 0, 5, 0, "tIjAb", 0, outfile);
  dpd_swap23(&t2, CC_TAMPS, 10, 10, "tIAjb", 0, outfile);
  dpd_buf_close(&t2);

  /* Build t2iaJB List */
  dpd_buf_init(&t2, CC_TAMPS, 0, 5, 0, 5, 0, "tiJaB", 0, outfile);
  dpd_swap23(&t2, CC_TAMPS, 10, 10, "tiaJB", 0, outfile);
  dpd_buf_close(&t2);

  /* Build t2IbjA and t2 jAIb Lists */
  dpd_buf_init(&t2, CC_TAMPS, 10, 10, 10, 10, 0, "tIAjb", 0, outfile);
  dpd_swap24(&t2, CC_TAMPS, 10, 10, "tIbjA", 0, outfile);
  dpd_swap13(&t2, CC_TAMPS, 10, 10, "tjAIb", 0, outfile);
  dpd_buf_close(&t2);

/*  timer_off("sort amps", outfile); */
}
