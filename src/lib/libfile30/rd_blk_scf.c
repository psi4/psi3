#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_blk_scf():  Reads in a symmetry block of the SCF eigenvector 
**   (or whatever is stored in its place).
**  
**   arguments: int irrep   designates the particular irrep to which the block
**     belongs.
**
**   STB(10/29/99)-This is now a wrapper function for rd_alpha_blk_scf();
**  
**   returns: double **scf_vector    This rectangular matrix has orbspi[irrep] 
**     rows.
*/


double **file30_rd_blk_scf(int irrep)
{
  return file30_rd_alpha_blk_scf(irrep);
}
