/*!
  \file nentry.c
  \ingroup (CHKPT)
*/

#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_nentry()  
**
** If zmatrix, reads total number of entries (including dummy atoms)
**
** takes no arguments.
**
** returns: nentry = total number of entries.
** \ingroup (CHKPT)
*/

int chkpt_rd_nentry(void)
{
  int num_entries;

  psio_read_entry(PSIF_CHKPT, "::Num. entries", (char *) &num_entries, 
		  sizeof(int));
  return num_entries;
}


/*!
** chkpt_wt_nentry()
**
** If zmatrix, writes the total number of entries (including dummy atoms)
**
** \param nentry  = total number of entries.
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_nentry(int num_entries)
{
  psio_write_entry(PSIF_CHKPT, "::Num. entries", (char *) &num_entries, 
		   sizeof(int));
}
