#include <stdio.h>
#include <qt.h>
#include "dpd.h"

/* oe_axpy(): Evaluates the standard operation a * X + Y -> Y for dpd
** one-electron files.
**
** Arguments:
**   struct oe_dpdfile *FileA: A pointer to the leftmost oe_dpdfile.
**   struct oe_dpdfile *FileB: A pointer to the rightmost (and target)
**                             oe_dpdfile.
**   double alpha: The scalar prefactor in the multiplication.
**   int transA: A boolean indicating that we should use the transpose of
**               FileA
**   int print_flag: A booelan for the print routines.
**   FILE *outfile: A pointer to the formatted output stream.
*/

int dpd_oe_axpy(struct oe_dpdfile *FileA, struct oe_dpdfile *FileB, 
		double alpha, int transA, int print_flag, FILE *outfile)
{
  int h, nirreps;
  int row, col;

  nirreps = FileA->params->nirreps;

  for(h=0; h < nirreps; h++) {

      dpd_oe_file_mat_irrep_init(FileA, h);
      dpd_oe_file_mat_irrep_init(FileB, h);
      dpd_oe_file_mat_irrep_rd(FileA, h, print_flag, outfile);
      dpd_oe_file_mat_irrep_rd(FileB, h, print_flag, outfile);

      if(!transA) {

	  for(row=0; row < FileA->params->rowtot[h]; row++)
	      for(col=0; col < FileA->params->coltot[h]; col++)
		  FileB->matrix[h][row][col] += alpha*FileA->matrix[h][row][col]; 

	}
      else {
	  for(row=0; row < FileB->params->rowtot[h]; row++)
	      for(col=0; col < FileB->params->coltot[h]; col++)
		  FileB->matrix[h][row][col] += alpha*FileA->matrix[h][col][row];
	}

      dpd_oe_file_mat_irrep_wrt(FileB, h, print_flag, outfile);
      dpd_oe_file_mat_irrep_close(FileA, h);
      dpd_oe_file_mat_irrep_close(FileB, h);
    }

  return 0;
}
      
