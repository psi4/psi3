#include <stdio.h>
#include <libciomr.h>
#include "dpd.h"

int dpd_file_mat_irrep_wrt_block(struct dpdfile *File, int irrep, int start_pq,
				 int num_pq, int print_flag, FILE *outfile)
{
  int rowtot, coltot;
  psio_address irrep_ptr, next_address;

  irrep_ptr = File->lfiles[irrep];
  rowtot = num_pq;
  coltot = File->params->coltot[irrep];

  /* Advance file pointer to current row */
  irrep_ptr = psio_get_address(irrep_ptr, start_pq*coltot*sizeof(double));

  if(rowtot && coltot)
     psio_write(File->filenum, File->label, (char *) File->matrix[irrep][0],
		rowtot*coltot*sizeof(double), irrep_ptr, &next_address);

  /* No printing for this code yet */
  if(0) {
      fprintf(outfile, "\n\tDPD File: %s\n", File->label);
      fprintf(outfile,   "\tMatrix for Irrep %1d\n", irrep);
      fprintf(outfile,   "\t-----------------------------------------\n");
      dpd_mat_irrep_print(File->matrix[irrep], File->params,
                          irrep, outfile);
    }

  return 0;

}
