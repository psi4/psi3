/*!
  \file symoper.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_symoper()
** Reads in the mapping array between "canonical" ordering of symmetry
** operations in the point group and the one defined in symmetry.h
**
**  returns: int *symoper    Array nirrep long
*/


int *chkpt_rd_symoper(void)
{
  int *symoper;
  int nirreps;

  nirreps = chkpt_rd_nirreps();
  symoper = init_int_array(nirreps);

  psio_read_entry(PSIF_CHKPT, "::Cotton -> local map", (char *) symoper, nirreps*sizeof(int));

  return symoper;
}

/*!
** void chkpt_wt_symoper(int *)
** Writes out the mapping array between "canonical" ordering of symmetry
** operations in the point group and the one defined in symmetry.h
**
**  arguments: 
**   \param int *symoper    Array nirrep long
**
**  returns: none
*/


void chkpt_wt_symoper(int *symoper)
{
  int nirreps;

  nirreps = chkpt_rd_nirreps();

  psio_write_entry(PSIF_CHKPT, "::Cotton -> local map", (char *) symoper, nirreps*sizeof(int));
}
