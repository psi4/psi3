/*!
  \file wt_blk_scf.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/*!
** file30_wt_blk_scf():  Writes in the SCF eigenvector (or whatever is to 
**     be stored in its place).
**
**   This is now a wrapper function for wt_alpha_blk_scf.
**
**   arguments: 
**
** \param int irrep   the number of the irrep to which the symmetry block 
**         belongs (this includes irreps with orbspi[irrep] == 0)
**         n.b. this routine assumes that the first irrep will have
**         irrep == 0.
**
** \param double **scf_vector    This should be a single symmetry
**         block of the SCF eigenvector.  Its dimension should be 
**         sopi[irrep]*orbspi[irrep];
**
**   returns: nothing.
*/

void file30_wt_blk_scf(double **scf_vector, int irrep)
{
    file30_wt_alpha_blk_scf(scf_vector,irrep);
} 
  
  

