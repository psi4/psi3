/*!
  \file close.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
**  chkpt_close()  closes up the checkpoint file.
** 
**  arguments: none, but chkpt_init must already have been called for 
**    this to work.  
**
**  returns: zero.  Perhaps this, too, will change one day.
*/


int chkpt_close(void)
{
  psio_close(PSIF_CHKPT, 1);
  return 0;
}
