#include <stdio.h>
#include <dpd.h>
#include <qt.h>
#define EXTERN
#include "globals.h"

void X_build(void)
{
  struct dpdbuf tIAJB, tiajb, X, Y, tIbjA, tiaJB;

/*  timer_on("X_build", outfile); */

  dpd_buf_init(&tIAJB, CC_TAMPS, 10, 10, 10, 10, 0, "tIAJB", 0, outfile);
  dpd_copy(&tIAJB, CC_MISC, "XIAJB", 0, outfile);
  dpd_buf_close(&tIAJB);
  dpd_buf_init(&X, CC_MISC, 10, 10, 10, 10, 0, "XIAJB", 0, outfile);
  dpd_scm(&X, 0.5, 0, outfile);
  dpd_buf_init(&Y, CC_MISC, 10, 10, 10, 10, 0, "YIAJB", 0, outfile);
  dpd_axpy(&Y, &X, -1, 0, outfile);
  dpd_buf_close(&Y);  dpd_buf_close(&X);

  dpd_buf_init(&tiajb, CC_TAMPS, 10, 10, 10, 10, 0, "tiajb", 0, outfile);
  dpd_copy(&tiajb, CC_MISC, "Xiajb", 0, outfile);
  dpd_buf_close(&tiajb);
  dpd_buf_init(&X, CC_MISC, 10, 10, 10, 10, 0, "Xiajb", 0, outfile);
  dpd_scm(&X, 0.5, 0, outfile);
  dpd_buf_init(&Y, CC_MISC, 10, 10, 10, 10, 0, "Yiajb", 0, outfile);
  dpd_axpy(&Y, &X, -1, 0, outfile);
  dpd_buf_close(&Y);  dpd_buf_close(&X);

  dpd_buf_init(&tIbjA, CC_TAMPS, 10, 10, 10, 10, 0, "tIbjA", 0, outfile);
  dpd_copy(&tIbjA, CC_MISC, "XIajB", 0, outfile);
  dpd_buf_init(&X, CC_MISC, 10, 10, 10, 10, 0, "XIajB", 0, outfile);
  dpd_scm(&X, 0.5, 0, outfile);
  dpd_buf_init(&Y, CC_MISC, 10, 10, 10, 10, 0, "YIajB", 0, outfile);
  dpd_axpy(&Y, &X, 1, 0, outfile);
  dpd_buf_close(&Y);  dpd_buf_close(&X);

  dpd_buf_init(&tiaJB, CC_TAMPS, 10, 10, 10, 10, 0, "tiaJB", 0, outfile);
  dpd_swap24(&tiaJB, CC_MISC, 10, 10, "XiAJb", 0, outfile);
  dpd_buf_close(&tiaJB);
  dpd_buf_init(&X, CC_MISC, 10, 10, 10, 10, 0, "XiAJb", 0, outfile);
  dpd_scm(&X, 0.5, 0, outfile);
  dpd_buf_init(&Y, CC_MISC, 10, 10, 10, 10, 0, "YiAJb", 0, outfile);
  dpd_axpy(&Y, &X, 1, 0, outfile);
  dpd_buf_close(&Y); dpd_buf_close(&X);

/*  timer_off("X_build", outfile); */
}
