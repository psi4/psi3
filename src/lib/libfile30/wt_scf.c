#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_wt_scf():  Writes out the SCF eigenvector (or whatever is stored in
**     its place).
**
**  This is now a wrapper function for wt_alpha_scf.
**
**   arguments: double **scf_vector    This rectangular matrix has dimentions nso
**     by nmo (see: rd_nmo()). For STO water, scf_vector 
**     should look something like the following:
**
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         0.0 0.0 0.0 0.0 *** 0.0 0.0
**         0.0 0.0 0.0 0.0 0.0 *** ***
**         0.0 0.0 0.0 0.0 0.0 *** ***
**
**     where the *** represent the non-zero values, and the 0.0 entries
**     represent (double)0.
**
**   returns: nothing.
*/

void file30_wt_scf(double **scf_vector)
{
  file30_wt_alpha_scf(scf_vector);
} 
  
  

