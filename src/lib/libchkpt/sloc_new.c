/*!
  \file sloc_new.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_sloc_new()	
** Read in an array of the numbers of the first basis 
** functions (not AOs as rd_sloc does)  from the shells.
**
** returns: 
**   sloc = Read in an array nshell long of the numbers of 
**          the first basis functions from the shells.
*/

int *chkpt_rd_sloc_new(void)
{
  int *sloc_new;
  int nshell;

  nshell = chkpt_rd_nshell();
  sloc_new = init_int_array(nshell);

  psio_read_entry(PSIF_CHKPT, "::First BF per shell", (char *) sloc_new, 
                  nshell*sizeof(int));

  return sloc_new;
}


/*!
** void chkpt_wt_sloc_new(int *)	
** Writes out an array of the numbers of the first basis 
** functions (not AOs as rd_sloc does)  from the shells.
**
** \param sloc = An array nshell long of the numbers of 
**               the first basis functions from the shells.
**
** returns: none
*/

void chkpt_wt_sloc_new(int *sloc_new)
{
  int nshell;

  nshell = chkpt_rd_nshell();

  psio_write_entry(PSIF_CHKPT, "::First BF per shell", 
		   (char *) sloc_new, nshell*sizeof(int));
}
