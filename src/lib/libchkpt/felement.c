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
  int nallatom, i;
  psio_address ptr;
  char *keyword;
  keyword = chkpt_build_keyword("Full atom labels");

  nallatom = chkpt_rd_nallatom();

  label = (char **)malloc(nallatom*sizeof(char*));
  for(i=0; i < nallatom; i++) 
    label[i] = (char *) malloc(MAX_ELEMNAME*sizeof(char));

  ptr = PSIO_ZERO;
  for(i=0; i < nallatom; i++)
    psio_read(PSIF_CHKPT, keyword, (char *) label[i], 
              MAX_ELEMNAME*sizeof(char), ptr, &ptr);

  free(keyword);
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
  int nallatom, i;
  psio_address ptr;
  char *keyword;
  keyword = chkpt_build_keyword("Full atom labels");

  nallatom = chkpt_rd_nallatom();

  ptr = PSIO_ZERO;
  for(i=0; i < nallatom; i++)
    psio_write(PSIF_CHKPT, keyword, (char *) label[i], 
               MAX_ELEMNAME*sizeof(char), ptr, &ptr);

  free(keyword);
}
