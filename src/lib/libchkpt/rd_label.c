/*!
  \file rd_label.c
*/

#include <string.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

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
