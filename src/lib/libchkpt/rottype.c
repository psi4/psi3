/*!
  \file rottype.c
  \ingroup (CHKPT)
*/

#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>

/*!
** int chkpt_rd_rottype()  
** Reads in type of the rigid rotor molecule represents.
**
** returns: rottype = type of rigid rotor. Allowed values are:
**            0 - asymmetric top
**            1 - symmetric top
**            2 - spherical top
**            3 - linear molecule
**            6 - atom
** \ingroup (CHKPT)
*/

int chkpt_rd_rottype(void)
{
  int rottype;

  psio_read_entry(PSIF_CHKPT, "::Rotor type", (char *) &rottype, 
                  sizeof(int) );
  return rottype;
}


/*!
** void chkpt_wt_rottype(int)  
** Reads in type of the rigid rotor molecule represents.
**
** \param rottype = type of rigid rotor. Allowed values are:
**            0 - asymmetric top
**            1 - symmetric top
**            2 - spherical top
**            3 - linear molecule
**            6 - atom
**
** returns: none
** \ingroup (CHKPT)
*/

void chkpt_wt_rottype(int rottype)
{
  psio_write_entry(PSIF_CHKPT, "::Rotor type", (char *) &rottype, sizeof(int));
}
