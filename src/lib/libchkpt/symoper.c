/*!
  \file symoper.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_symoper()
** Reads in the mapping array between "canonical" ordering of symmetry
** operations in the point group and the one defined in symmetry.h
**
** returns: symoper = Array nirrep long
**
** \ingroup (CHKPT)
*/

int *chkpt_rd_symoper(void)
{
  int *symoper;
  int nirreps;
  char *keyword;
  keyword = chkpt_build_keyword("Cotton -> local map");

  nirreps = chkpt_rd_nirreps();
  symoper = init_int_array(nirreps);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) symoper, nirreps*sizeof(int));

  free(keyword);
  return symoper;
}


/*!
** void chkpt_wt_symoper(int *)
** Writes out the mapping array between "canonical" ordering of symmetry
** operations in the point group and the one defined in symmetry.h
**
** \param symoper = Array nirrep long
**
** returns: none
**
** \ingroup (CHKPT)
*/

void chkpt_wt_symoper(int *symoper)
{
  int nirreps;
  char *keyword;
  keyword = chkpt_build_keyword("Cotton -> local map");

  nirreps = chkpt_rd_nirreps();

  psio_write_entry(PSIF_CHKPT, keyword, (char *) symoper, nirreps*sizeof(int));

  free(keyword);
}
