/*!
  \file rd_evals.c
*/

#include <stdio.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_rd_evals():  Reads in the SCF eigenvalues: the SCF orbital energies.
**
**  This is now a wrapper function due to the fact that we have a pointer 
**  structure for the SCF data in file30 and possibility of UHF
**
**  returns: double *evals   an array of _all_ of the SCF eigenvalues,
**      ordered by irrep, and by increasing energy within each irrep.  
**      (i.e. for sto water, the four a1 eigenvalues all come first, and 
**      those four are ordered from lowest energy to highest energy,
**      followed by the single b1 eigenvalue, etc.)
*/

double *file30_rd_evals(void)
{
  return file30_rd_alpha_evals();
}
