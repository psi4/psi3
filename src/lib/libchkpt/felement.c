/*!
  \file felement.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** chkpt_rd_felement():  Reads in element labels including dummy atoms
**
**   takes no arguments.
**
**   returns: char **label element label matrix
** \ingroup (CHKPT)
*/

char **chkpt_rd_felement(void)
{
  char **label;
  int nentry, i;
  psio_address ptr;

  nentry = chkpt_rd_nentry();

  label = (char **)malloc(nentry*sizeof(char*));
  for(i=0; i < nentry; i++) 
    label[i] = (char *) malloc(MAX_ELEMNAME*sizeof(char));

  ptr = PSIO_ZERO;
  for(i=0; i < nentry; i++)
    psio_read(PSIF_CHKPT, "::Full atom labels", (char *) label[i], 
              MAX_ELEMNAME*sizeof(char), ptr, &ptr);

  return label;  
}


/*!
** chkpt_wt_felement():  Writes out element labels including dummy atoms
**
** arguments: 
**   \param label = element label matrix.
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_felement(char **label)
{
  int nentry, i;
  psio_address ptr;

  nentry = chkpt_rd_nentry();

  ptr = PSIO_ZERO;
  for(i=0; i < nentry; i++)
    psio_write(PSIF_CHKPT, "::Full atom labels", (char *) label[i], 
               MAX_ELEMNAME*sizeof(char), ptr, &ptr);
}
