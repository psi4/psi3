#include <stdio.h>
#include <qt.h>
#include "dpd.h"

/* oe_dirprd(): Computes the direct product between two one-electron dpd
** files.
**
** Arguments:
**   struct oe_dpdfile *FileA: A pointer to one of the one-electron dpd
**                             files.
**   struct oe_dpdfile *FileB: A pointer to the other one-electron dpd
**                             file.
**   int print_flag: A value for the print routines.
**   FILE *outfile: A pointer to the output file stream.
*/

int dpd_oe_dirprd(struct oe_dpdfile *FileA, struct oe_dpdfile *FileB, 
		  int print_flag, FILE *outfile)
{
  int h, nirreps;

  nirreps = FileA->params->nirreps;

  for(h=0; h < nirreps; h++) {

      dpd_oe_file_mat_irrep_init(FileA, h);
      dpd_oe_file_mat_irrep_init(FileB, h);
      dpd_oe_file_mat_irrep_rd(FileA, h, print_flag, outfile);
      dpd_oe_file_mat_irrep_rd(FileB, h, print_flag, outfile);

      dirprd_block(FileA->matrix[h], FileB->matrix[h],
		   FileA->params->rowtot[h], FileA->params->coltot[h]);

      dpd_oe_file_mat_irrep_wrt(FileB, h, print_flag, outfile);
      dpd_oe_file_mat_irrep_close(FileA, h);
      dpd_oe_file_mat_irrep_close(FileB, h);
    }

  return 0;
}
      
