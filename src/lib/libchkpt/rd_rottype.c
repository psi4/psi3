/*!
  \file rd_rottype.c
*/

#include "chkpt.h"
#include <libpsio/psio.h>

/*!
** int chkpt_rd_rottype()  
** Reads in type of the rigid rotor molecule represents.
**
** returns: int rottype   type of rigid rotor. Allowed values are:
**            0 - asymmetric top
**            1 - symmetric top
**            2 - spherical top
**            3 - linear molecule
**            6 - atom
*/


int chkpt_rd_rottype(void)
{
  int rottype;

  psio_read_entry(PSIF_CHKPT, "::Rottype", (char *) &rottype, 
                  sizeof(int) );
  return rottype;
}
