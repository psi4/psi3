#include "file30.h"
#include "file30.gbl"

/* file30_rd_scf(): Shawn T. Brown(10/28/99) Reads in just the alpha eigenvector matrix.  
Modified for UHF but to keep compatibility so that we didn't have to redo all of the code.
Since the alpha matrix is treated as the closed shell matrix when writing the file30, now 
all that is needed is to use this wrapper function.

returns the same matrix described in file30_rd_alpha_scf() */

double **file30_rd_scf(void)
{
    return file30_rd_alpha_scf();
}
