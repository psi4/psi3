/*!
  \file sym_label.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_sym_label():  Reads in the symmetry label.
**
**   takes no arguments.
**
**   returns: symmetry = symmetry label.
**
** \ingroup (CHKPT)
*/

char *chkpt_rd_sym_label(void)
{
  char *sym_label;
  char *key;

  sym_label = (char *) malloc(4*sizeof(char));

  key = chkpt_build_keyword("Symmetry label");
  psio_read_entry(PSIF_CHKPT, key, (char *) sym_label, 4*sizeof(char));
  sym_label[3] = '\0';
  free(key);

  return sym_label;  
}


/*!
** chkpt_wt_sym_label():  Writes out the symmetry label.
**
** \param symmetry = symmetry label.
**
** returns none
**
** \ingroup (CHKPT)
*/

void chkpt_wt_sym_label(char *sym_label)
{
  char *key;
  key = chkpt_build_keyword("Symmetry label");
  psio_write_entry(PSIF_CHKPT, "::Symmetry label", (char *) sym_label, 
                   4*sizeof(char));
  free(key);
}

