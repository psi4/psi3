/*!
  \file cartrep.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** chkpt_rd_cartrep():  Reads the point group representation in the basis of
**     cartesian unit vectors.
**
**   takes no arguments.
**
**   returns: double **cartrep  a vector of block matrices of doubles. Each 
**     row corresponds to a particular symmetry operation, each column is 
**     a 3x3 block matrix.
**  \ingroup (CHKPT)
*/

double **chkpt_rd_cartrep(void)
{
  int i, nirrep;
  double **cartrep;
  psio_address ptr;
  char *keyword;
  keyword = chkpt_build_keyword("Cart. Repr. Matrices");

  nirrep = chkpt_rd_nirreps();
  ptr = PSIO_ZERO;
  cartrep = block_matrix(nirrep,9);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) cartrep[0], 
      9*nirrep*sizeof(double));

  free(keyword);
  return cartrep;
}


/*!
** chkpt_wt_cartrep():  Writes the point group representation in the basis of
**     cartesian unit vectors.
**
** \param cartrep = a vector of block matrices of doubles. Each row 
**                  corresponds to a particular symmetry operation, each 
**                  column is a 3x3 block matrix.
**
** returns nothing.
** \ingroup (CHKPT)
*/

void chkpt_wt_cartrep(double **cartrep)
{
  int i, nirrep;
  psio_address ptr;
  char *keyword;
  keyword = chkpt_build_keyword("Cart. Repr. Matrices");

  nirrep = chkpt_rd_nirreps();

  ptr = PSIO_ZERO;
  for(i=0; i < nirrep; i++)
    psio_write(PSIF_CHKPT, keyword, (char *) cartrep[i], 
               9*sizeof(double), ptr, &ptr);

  free(keyword);
}
