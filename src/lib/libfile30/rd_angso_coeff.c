#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_angso_coeff():     Read in coefficients the array of AO to SO 
** coefficients for every angular momentum present in the basis. 
**
**   takes no arguments.
**
**   returns: double **angso_coeff  an array of mx_angso_coeff double
**
*/


double *file30_rd_angso_coeff(void)
{
  double *angso_coeff;
  int mx_angso_coeff;
  PSI_FPTR junk;
  PSI_FPTR angso_coeff_ptr;

  mx_angso_coeff = file30_rd_mx_angso_coeff();
  angso_coeff_ptr = (PSI_FPTR) (info30_.mpoint[16]-1)*sizeof(int);
  angso_coeff = init_array(mx_angso_coeff);

  wreadw(info30_.filenum, (char *) angso_coeff, 
         (int) mx_angso_coeff*sizeof(double), angso_coeff_ptr, &junk);

  return angso_coeff;
}
