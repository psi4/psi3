/*!
  \file label.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_label():  Reads the main chkpt label.
**
**   takes no arguments.
**
**   returns: pointer to the checkpoint label
** \ingroup (CHKPT)
*/

char *chkpt_rd_label(void)
{
  char *label;
  char *keyword;
  keyword = chkpt_build_keyword("Label");

  label = (char *) malloc(80 * sizeof(char));

  psio_read_entry(PSIF_CHKPT, keyword, (char *) label, 80*sizeof(char));

  free(keyword);
  return label;
}

/*!
** chkpt_wt_label():  Writes the main chkpt label.
**
**  arguments:
**  \param label = The calculation label.
**
**   returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_label(char *label)
{
  char *keyword;
  keyword = chkpt_build_keyword("Label");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) label, 80*sizeof(char));

  free(keyword);
}
