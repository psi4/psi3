/*!
  \file override_occ.c
  \ingroup (CHKPT)
*/

#include <stdlib.h>
#include "chkpt.h"
#include <libpsio/psio.h>
#include <psifiles.h>

/*!
** chkpt_rd_override_occ(): Reads flag which tells cscf to ignore docc/socc
** vectors and use occupations in chkpt file instead
**  
** takes no arguments.
**  
** returns: 1 if chkpt occupations should be forced; 0 otherwise
**        
** \ingroup (CHKPT)
*/

int chkpt_rd_override_occ(void)
{
  int override=0, i, keylen;
  char *keyword;
  /*
  keyword = chkpt_build_keyword("Override Occupations");
  if ( chkpt_exist(keyword) )
    psio_read_entry(PSIF_CHKPT, keyword, (char *) &override, sizeof(int));
  free(keyword);
  return override;
  */

  /* new way - keep override in root area always */
  keylen = 22;
  keyword = (char *) malloc((keylen+1)*sizeof(char));
  sprintf(keyword, "::%s", "Override Occupations");
  keyword[keylen] = '\0';

  if ( chkpt_exist(keyword) )
    psio_read_entry(PSIF_CHKPT, keyword, (char *) &override, sizeof(int));

  free(keyword);
  return override;
}

/*!
** chkpt_wt_override_occ(): Writes flag which tells cscf to ignore docc/socc
** vectors and use occupations in chkpt file instead
**  
** arguments: (int) 1 to set override; 0 otherwise
**  
** returns: none
**        
** \ingroup (CHKPT)
*/

void chkpt_wt_override_occ(int override)
{
  char *keyword;
  int keylen;
  /*
  keyword = chkpt_build_keyword("Override Occupations");
  psio_write_entry(PSIF_CHKPT, keyword, (char *) &override, sizeof(int));
  */

  /* new way - keep override in root area always */
  keylen = 22;
  keyword = (char *) malloc((keylen+1)*sizeof(char));
  sprintf(keyword, "::%s", "Override Occupations");
  keyword[keylen] = '\0';

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &override, sizeof(int));

  free(keyword);
  return;
}
