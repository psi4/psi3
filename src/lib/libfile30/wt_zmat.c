/*!
  \file wt_zmat.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>
 
/*!
** file30_wt_zmat(struct z_entry *z_geom, int num_atoms)
**   : Write the z-matrix geometry to file30
**
**  arguments:
** \param struct z_entry *z_geom -- the z-matrix as an array of z_entry
**                                       structs
** \param int num_atoms          -- number of atoms
**
**  returns: none
*/
 
void file30_wt_zmat(struct z_entry *z_geom, int num_atoms) {

  PSI_FPTR junk;
  PSI_FPTR z_ptr;

  z_ptr = (PSI_FPTR) (info30_.mpoint[46] - 1) * sizeof(int);

  wwritw(info30_.filenum,(char *) z_geom, num_atoms*(sizeof(struct z_entry)), z_ptr, &junk);

  return;
}
