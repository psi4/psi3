/*!
  \file init.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/* first definition of chkpt_prefix */
char chkpt_prefix[CHKPT_PREFIX_LEN];

/*!
**  chkpt_init()  Initializes the checkpoint file for other chkpt_
**    functions to perform their duties.
**
**  arguments: 
**    int status: boolean indicating if the chkpt file should be
**                initialized (PSIO_OPEN_NEW) or the old chkpt 
**                file should be used (PSIO_OPEN_OLD).
**
**  returns: zero.  Perhaps this will change some day.
** \ingroup (CHKPT)
*/

int chkpt_init(int status)
{
  char *prefix;
  psio_tocentry *this_entry;

  psio_open(PSIF_CHKPT, status);

  if(psio_tocscan(PSIF_CHKPT, "Default prefix") != NULL) {
    prefix = chkpt_rd_prefix();
    chkpt_set_prefix(prefix);
    free(prefix);
  }
  else {
    chkpt_set_prefix("");
    chkpt_commit_prefix();  /* we assume that no default prefix existed in PSIF_CHKPT */
  }

  return 0;
}

