/*!
  \file ua2a.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_ua2a()
** Reads in a mapping array from the symmetry-unique atom
** list to the full atom list
**
** returns: ua2a = Read in an array num_unique_atom long
** 
** \ingroup (CHKPT)
*/

int *chkpt_rd_ua2a(void)
{
  int *ua2a;
  int num_unique_atoms;
  char *keyword;
  keyword = chkpt_build_keyword("Unique atom -> full atom map");

  num_unique_atoms = chkpt_rd_num_unique_atom();
  ua2a = init_int_array(num_unique_atoms);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) ua2a, num_unique_atoms*sizeof(int));

  free(keyword);
  return ua2a;
}


/*!
** void chkpt_wt_ua2a(int *)
** Writes out a mapping array from the symmetry-unique atom
** list to the full atom list
**
** \param ua2a = An array num_unique_atom long
**
** returns: none
*/

void chkpt_wt_ua2a(int *ua2a)
{
  int num_unique_atoms;
  char *keyword;
  keyword = chkpt_build_keyword("Unique atom -> full atom map");

  num_unique_atoms = chkpt_rd_num_unique_atom();

  psio_write_entry(PSIF_CHKPT, keyword, (char *) ua2a, num_unique_atoms*sizeof(int));

  free(keyword);
}
