/*!
  \file label.c
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
**   returns: char *
*/

char *chkpt_rd_label(void)
{
  char *label;

  label = (char *) malloc(80 * sizeof(char));
  psio_read_entry(PSIF_CHKPT, "::Label", (char *) label, 80*sizeof(char));

  return label;
}

/*!
** chkpt_wt_label():  Writes the main chkpt label.
**
**  arguments:
**  \param char *label  The calculation label.
**
**   returns: nothing
*/

void chkpt_wt_label(char *label)
{
  psio_write_entry(PSIF_CHKPT, "::Label", (char *) label, 80*sizeof(char));
}
