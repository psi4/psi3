/*! \defgroup CHKPT libchkpt: The Checkpoint Interface library */

/*!
  \file am2canon_shell_order.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** int *chkpt_rd_am2canon_shell_order() 
** Reads in the mapping array from the am-ordered
** to the canonical (in the order of appearance)
** list of shells.
**
** returns: int *am2can_shell_order
** \ingroup(CHKPT)
*/


int *chkpt_rd_am2canon_shell_order(void)
{
  int *am2can_sh_ord, nshell;
  char *keyword;
  keyword = chkpt_build_keyword("AM -> canonical shell map"); 

  nshell = chkpt_rd_nshell();
  am2can_sh_ord = init_int_array(nshell);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) am2can_sh_ord,
      nshell*sizeof(int));

  free(keyword);
  return am2can_sh_ord;
}


/*!
** void chkpt_wt_am2canon_shell_order(int *) 
** Writes out the mapping array from the am-ordered
** to the canonical (in the order of appearance)
** list of shells.
**
**  \param am2can_shell_order = array to store the mapping array
**
** returns: none
** \ingroup(CHKPT)
*/

void chkpt_wt_am2canon_shell_order(int *am2can_sh_ord)
{
  int nshell;
  char *keyword;
  keyword = chkpt_build_keyword("AM -> canonical shell map");

  nshell = chkpt_rd_nshell();

  psio_write_entry(PSIF_CHKPT, keyword, (char *) am2can_sh_ord,
      nshell*sizeof(int));

  free(keyword);
}
