/*!
  \file us2c.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>


/*!
** int *chkpt_rd_us2s()
** Read in a mapping array betwen unique shell and 
** full shell lists
**
** returns: us2s = Read in an array num_unique_shell
** 
** \ingroup (CHKPT)
*/

int *chkpt_rd_us2s(void)
{
  int *us2s;
  int num_unique_shells;
  char *keyword;
  keyword = chkpt_build_keyword("Unique shell -> full shell map");

  num_unique_shells = chkpt_rd_num_unique_shell();
  us2s = init_int_array(num_unique_shells);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) us2s, num_unique_shells*sizeof(int));

  free(keyword);
  return us2s;
}


/*!
** void chkpt_wt_us2s(int *)
** Writes out a mapping array betwen unique shell and 
** full shell lists.
**
** \param us2s = An array num_unique_shell
**
** returns: none
** 
** \ingroup (CHKPT)
*/

void chkpt_wt_us2s(int *us2s)
{
  int num_unique_shells;
  char *keyword;
  keyword = chkpt_build_keyword("Unique shell -> full shell map");

  num_unique_shells = chkpt_rd_num_unique_shell();

  psio_write_entry(PSIF_CHKPT, keyword, (char *) us2s, num_unique_shells*sizeof(int));

  free(keyword);
}
