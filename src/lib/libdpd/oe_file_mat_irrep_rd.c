#include <stdio.h>
#include <psio.h>
#include "dpd.h"

int dpd_oe_file_mat_irrep_rd(struct oe_dpdfile *File, int irrep, int print_flag,
			     FILE *outfile)
{
  int rowtot, coltot;
  psio_address irrep_ptr, next_address;

  irrep_ptr = File->lfiles[irrep];
  rowtot = File->params->rowtot[irrep];
  coltot = File->params->coltot[irrep];

  if(rowtot && coltot)
     psio_read(File->filenum, File->label, (char *) File->matrix[irrep][0],
	       rowtot*coltot*sizeof(double), irrep_ptr, &next_address);

  if(print_flag) {
      fprintf(outfile, "\n\tOE DPD File: %s\n", File->label);
      fprintf(outfile,   "\tMatrix for Irrep %1d\n", irrep);
      fprintf(outfile,   "\t-----------------------------------------\n");
      dpd_oe_mat_irrep_print(File->matrix[irrep], File->params,
                             irrep, outfile);
    }

  return 0;

}
