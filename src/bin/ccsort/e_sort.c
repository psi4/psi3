#include <stdio.h>
#include <dpd.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

void e_sort(void)
{
  struct dpdbuf E;

  dpd_buf_init(&E, CC_EINTS, 11, 0, 11, 0, 0, "E <ai|jk>", 0, outfile);
  dpd_swap13(&E, CC_TMP0, 0, 11, "E <ai|jk> (ji,ak)", 0, outfile);
  dpd_buf_close(&E);
  dpd_buf_init(&E, CC_TMP0, 0, 11, 0, 11, 0, "E <ai|jk> (ji,ak)", 0, outfile);
  dpd_swap24(&E, CC_TMP1, 0, 11, "E <ai|jk> (jk,ai)", 0, outfile);
  dpd_buf_close(&E);
  dpd_buf_init(&E, CC_TMP1, 0, 11, 0, 11, 0, "E <ai|jk> (jk,ai)", 0, outfile);
  dpd_swap12(&E, CC_TMP0, 0, 11, "E <ai|jk> (kj,ai)", 0, outfile);
  dpd_buf_close(&E);
  dpd_buf_init(&E, CC_TMP0, 0, 11, 0, 11, 0, "E <ai|jk> (kj,ai)", 0, outfile);
  dpd_swap34(&E, CC_EINTS, 0, 10, "E <ij|ka>", 0, outfile);
  dpd_buf_close(&E);

  dpd_buf_init(&E, CC_EINTS, 11, 0, 11, 0, 1, "E <ai|jk>", 0, outfile);
  dpd_swap13(&E, CC_TMP0, 0, 11, "E <ai||jk> (ji,ak)", 0, outfile);
  dpd_buf_close(&E);
  dpd_buf_init(&E, CC_TMP0, 0, 11, 0, 11, 0, "E <ai||jk> (ji,ak)", 0, outfile);
  dpd_swap24(&E, CC_TMP1, 0, 11, "E <ai||jk> (jk,ai)", 0, outfile);
  dpd_buf_close(&E);
  dpd_buf_init(&E, CC_TMP1, 0, 11, 0, 11, 0, "E <ai||jk> (jk,ai)", 0, outfile);
  dpd_swap12(&E, CC_TMP0, 0, 11, "E <ai||jk> (kj,ai)", 0, outfile);
  dpd_buf_close(&E);
  dpd_buf_init(&E, CC_TMP0, 0, 11, 0, 11, 0, "E <ai||jk> (kj,ai)", 0, outfile);
  dpd_swap34(&E, CC_EINTS, 2, 10, "E <ij||ka> (i>j,ka)", 0, outfile);
  dpd_buf_close(&E);
}
