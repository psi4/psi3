#include <stdio.h>
#include "dpd.h"

int dpd_oe_copy(struct oe_dpdfile *InFile, int outfilenum, 
		char *label, int print_flag, FILE *outfile)
{
  int h;
  double ***matrix;
  struct oe_dpdfile OutFile;

  dpd_oe_file_init(&OutFile, outfilenum, InFile->params->pnum,
		   InFile->params->qnum, label, print_flag, outfile);

  /* Save the OutFile's matrix pointer */
  matrix = OutFile.matrix;

  /* Temporarily assign the OutFile's matrix pointer to InBuf */
  OutFile.matrix = InFile->matrix;

  dpd_oe_file_mat_init(InFile);
  dpd_oe_file_mat_rd(InFile, print_flag, outfile);

  dpd_oe_file_mat_wrt(&OutFile, print_flag, outfile);

  dpd_oe_file_mat_close(InFile);

  /* Clean up the OutFile's matrix pointer */
  OutFile.matrix = matrix;

  dpd_oe_file_close(&OutFile);

  return 0;
}
