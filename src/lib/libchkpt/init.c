/*!
  \file init.c
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
**  chkpt_init()  Initializes the checkpoint file for other chkpt_
**    functions to perform their duties.
**
**  arguments: none, but it requires that the input parser be initialized
**    so that it can open the file.  
**
**  returns: zero.  Perhaps this will change some day.
*/

int chkpt_init(void)
{
  psio_open(PSIF_CHKPT, PSIO_OPEN_OLD);

  return 0;
}








