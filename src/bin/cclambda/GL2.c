#include <stdio.h>
#include <dpd.h>
#define EXTERN
#include "globals.h"

void GaeL2(void)
{
  struct dpdbuf newLijab, newLIJAB, newLIjAb;
  struct dpdbuf D;
  struct oe_dpdfile GAE, Gae;

  /* RHS += P(ab)<ij||ae>Gbe */
  dpd_oe_file_init(&GAE, CC_OEI, 1, 1, "GAE", 0, outfile);
  dpd_oe_file_init(&Gae, CC_OEI, 1, 1, "Gae", 0, outfile);

  dpd_buf_init(&newLIJAB, CC_LAMPS, 2, 5, 2, 7, 0, "New LIJAB",
           0, outfile);

  dpd_buf_init(&D, CC_DINTS, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)",
          0, outfile);
  dpd_contract221(&D, &GAE, &newLIJAB, 3, 1, 0, 1.0, 1.0, 0, outfile);
  dpd_contract212(&GAE, &D, &newLIJAB, 1, 2, 1, 1.0, 1.0, 0, outfile);

  dpd_buf_close(&newLIJAB);


  dpd_buf_init(&newLijab, CC_LAMPS, 2, 5, 2, 7, 0, "New Lijab",
           0, outfile);

  dpd_contract221(&D, &Gae, &newLijab, 3, 1, 0, 1.0, 1.0, 0, outfile);
  dpd_contract212(&Gae, &D, &newLijab, 1, 2, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&D);

  dpd_buf_close(&newLijab);


  dpd_buf_init(&newLIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb",
          0, outfile);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_contract221(&D, &Gae, &newLIjAb, 3, 1, 0, 1.0, 1.0, 0, outfile);
  dpd_contract212(&GAE, &D, &newLIjAb, 1, 2, 1, 1.0, 1.0, 0, outfile);
  dpd_buf_close(&D);

  dpd_buf_close(&newLIjAb);

  dpd_oe_file_close(&GAE);
  dpd_oe_file_close(&Gae);

}

void GmiL2(void)
{

  struct dpdbuf newLijab, newLIJAB, newLIjAb;
  struct dpdbuf D;
  struct oe_dpdfile GMI, Gmi;

  /* RHS -= P(ij) * <im||ab> * Gmj */

  dpd_oe_file_init(&GMI, CC_OEI, 0, 0, "GMI", 0, outfile);
  dpd_oe_file_init(&Gmi, CC_OEI, 0, 0, "Gmi", 0, outfile);

  dpd_buf_init(&newLIJAB, CC_LAMPS, 0, 7, 2, 7, 0, "New LIJAB",
           0, outfile);

  dpd_buf_init(&D, CC_DINTS, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)",
          0, outfile);
  dpd_contract221(&D, &GMI, &newLIJAB, 1, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_contract212(&GMI, &D, &newLIJAB, 0, 0, 0, -1.0, 1.0, 0, outfile);

  dpd_buf_close(&newLIJAB);


  dpd_buf_init(&newLijab, CC_LAMPS, 0, 7, 2, 7, 0, "New Lijab",
           0, outfile);

  dpd_contract221(&D, &Gmi, &newLijab, 1, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_contract212(&Gmi, &D, &newLijab, 0, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&D);

  dpd_buf_close(&newLijab);


  dpd_buf_init(&newLIjAb, CC_LAMPS, 0, 5, 0, 5, 0, "New LIjAb",
          0, outfile);

  dpd_buf_init(&D, CC_DINTS, 0, 5, 0, 5, 0, "D <ij|ab>", 0, outfile);
  dpd_contract221(&D, &Gmi, &newLIjAb, 1, 0, 1, -1.0, 1.0, 0, outfile);
  dpd_contract212(&GMI, &D, &newLIjAb, 0, 0, 0, -1.0, 1.0, 0, outfile);
  dpd_buf_close(&D);

  dpd_buf_close(&newLIjAb);

  dpd_oe_file_close(&Gmi);
  dpd_oe_file_close(&GMI);

}
