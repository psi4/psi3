/*!
  \file atom_dummy.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_atom_dummy()  
**
** Reads the array of flags which indicate whether the atom in full_geom is dummy
**
** takes no arguments.
**
** returns: atom_dummy = array of integers nallatom long.
** \ingroup (CHKPT)
*/

int* chkpt_rd_atom_dummy(void)
{
  int num_allatoms;
  int *atom_dummy;
  char *keyword;
  keyword = chkpt_build_keyword("Dummy atom flags");

  num_allatoms = chkpt_rd_nallatom();
  atom_dummy = (int *) malloc(sizeof(int)*num_allatoms);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) atom_dummy, 
		  num_allatoms*sizeof(int));

  free(keyword);
  return atom_dummy;
}


/*!
** chkpt_wt_atom_dummy()
**
** Writes the array of flags which indicate whether the atom in full_geom is dummy
**
** \param atom_dummy = array of integers nallatom long.
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_atom_dummy(int* atom_dummy)
{
  int num_allatoms = chkpt_rd_nallatom();
  char *keyword;
  keyword = chkpt_build_keyword("Dummy atom flags");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) atom_dummy, 
		   num_allatoms*sizeof(int));

  free(keyword);
}
