#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_angso_labels():     Read in labels (see libfile30.txt) for
** the array of AO to SO coefficients for every angular momentum
** present in the basis. angso_labels is stored in packed form
** (see pack_4int.c)
**
**   takes no arguments.
**
**   returns: int **angso_labels   a matrix of mx_angso_coeff by 4 integers
**
*/


int **file30_rd_angso_labels(void)
{
  int **angso_labels, *tmp;
  int mx_angso_coeff;
  PSI_FPTR junk;
  PSI_FPTR angso_labels_ptr;

  mx_angso_coeff = file30_rd_mx_angso_coeff();
  angso_labels_ptr = (PSI_FPTR) (info30_.mpoint[17]-1)*sizeof(int);
  tmp = init_int_array(mx_angso_coeff);
  angso_labels = init_int_matrix(mx_angso_coeff,4);
  wreadw(info30_.filenum, (char *) tmp, (int) mx_angso_coeff*sizeof(int),
	 angso_labels_ptr, &junk);

  unpack_4int(tmp, angso_labels, mx_angso_coeff, 4);
  free(tmp);

  return angso_labels;
}
