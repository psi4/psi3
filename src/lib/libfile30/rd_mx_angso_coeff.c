#include "file30.h"
#include "file30.gbl"

/*
** file30_rd_mx_angso_coeff():  Reads in the total number of elements in
** the array of AO to SO coefficients for every angular momentum present
** in the basis.
**
**   takes no arguments.
**
**   returns: int mx_angso_coeff
*/


int file30_rd_mx_angso_coeff(void)
{
  return info30_.mconst[34];
}
