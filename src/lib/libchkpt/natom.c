/*!
  \file natom.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
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

  psio_read_entry(PSIF_CHKPT, "::Num. atoms", (char *) &natom, sizeof(int));
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
  psio_write_entry(PSIF_CHKPT, "::Num. atoms", (char *) &natom, sizeof(int));
}

