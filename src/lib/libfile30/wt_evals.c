/*!
  \file wt_evals.c
*/

#include <stdio.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_wt_evals():  Writes the SCF eigenvalues: the SCF orbital energies.
**
**  This is now a wrapper function for wt_alpha_evals.
**
**  arguments:
** \param double *evals   an array of _all_ of the SCF eigenvalues,
**      ordered by irrep, and by increasing energy within each irrep.  
**      (i.e. for sto water, the four a1 eigenvalues all come first, and 
**      those four are ordered from lowest energy to highest energy,
**      followed by the single b1 eigenvalue, etc.)
**
**  returns nothing
*/

void file30_wt_evals(double *evals)
{
    file30_wt_alpha_evals(evals);
}

