#include "file30.h"
#include "file30.gbl"

/*
** rd_rottype():  Reads in type of the rigid rotor molecule represents.
**
**   takes no arguments.
**
**   returns: int rottype   type of rigid rotor. Allowed values are:
**            0 - asymmetric top
**            1 - symmetric top
**            2 - spherical top
**            3 - linear molecule
**            6 - atom
*/


int file30_rd_rottype(void)
{
  return info30_.mconst[6];
}
