/*!
  \file natom.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <stdlib.h>
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
**  int chkpt_rd_natom() 
**  Reads in the total number of atoms.
**
**  returns: natom = total number of atoms.
** \ingroup (CHKPT)
*/

int chkpt_rd_natom(void)
{  
  int natom;
  char *keyword;
  keyword = chkpt_build_keyword("Num. atoms");

  psio_read_entry(PSIF_CHKPT, keyword, (char *) &natom, sizeof(int));

  free(keyword);
  return natom;
}


/*!
**  void chkpt_wt_natom(int) 
**  Writes out the total number of atoms.
**
**  \param natom = total number of atoms.
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_natom(int natom)
{  
  char *keyword;
  keyword = chkpt_build_keyword("Num. atoms");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) &natom, sizeof(int));

  free(keyword);
}

