#include <dpd.h>
#define EXTERN
#include "globals.h"

/* BUILD_X(): Construct the orbital rotation gradient, XAI, for CC
** gradient calculations:
**
**  Xai = I'ia - I'ai
** */

void build_X(void)
{
  struct oe_dpdfile X, I, X2;

  dpd_oe_file_init(&I, CC_OEI, 1, 0, "I'AI", 0, outfile);
  dpd_oe_copy(&I, CC_OEI, "XAI", 0, outfile);
  dpd_oe_file_close(&I);

  dpd_oe_file_init(&X, CC_OEI, 1, 0, "XAI", 0, outfile);
  dpd_oe_scm(&X, -1.0, 0, outfile);
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'IA", 0, outfile);
  dpd_oe_axpy(&I, &X, 1.0, 1, 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_close(&X);

  dpd_oe_file_init(&I, CC_OEI, 1, 0, "I'ai", 0, outfile);
  dpd_oe_copy(&I, CC_OEI, "Xai", 0, outfile);
  dpd_oe_file_close(&I);

  dpd_oe_file_init(&X, CC_OEI, 1, 0, "Xai", 0, outfile);
  dpd_oe_scm(&X, -1.0, 0, outfile);
  dpd_oe_file_init(&I, CC_OEI, 0, 1, "I'ia", 0, outfile);
  dpd_oe_axpy(&I, &X, 1.0, 1, 0, outfile);
  dpd_oe_file_close(&I);
  dpd_oe_file_close(&X);

  /* Build spatial orbital version of X for Zvector:
     X(A,I) = 1/2[X(A,I)+X(a,i)]  */
  dpd_oe_file_init(&X, CC_OEI, 1, 0, "XAI", 0, outfile);
  dpd_oe_copy(&X, CC_MISC, "X(A,I)", 0, outfile);
  dpd_oe_file_close(&X);
  dpd_oe_file_init(&X, CC_MISC, 1, 0, "X(A,I)", 0, outfile);
  dpd_oe_file_init(&X2, CC_OEI, 1, 0, "Xai", 0, outfile);
  dpd_oe_axpy(&X2, &X, 1.0, 0, 0, outfile);
  dpd_oe_file_close(&X2);
  dpd_oe_scm(&X, 0.5, 0, outfile);
  dpd_oe_file_close(&X);
}
