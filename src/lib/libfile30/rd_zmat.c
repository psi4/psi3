#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_zmat():  Reads in the z_matrix.
**
**   takes no arguments.

**   returns: double *z_geom   An array natom long which contains 
**     a z_entry struct for each atom
*/

 struct z_entry *file30_rd_zmat(void)
{
  int natom;
  PSI_FPTR z_geom_ptr, junk;
  struct z_entry *z_geom;

  natom = file30_rd_natom();

  z_geom = (struct z_entry *) malloc(natom*(sizeof(struct z_entry)));

  z_geom_ptr = (PSI_FPTR) (info30_.mpoint[46] - 1)*sizeof(int);

  wreadw(info30_.filenum, (char *) z_geom, (int) sizeof(struct z_entry)*natom,
         z_geom_ptr, &junk);

  return  z_geom;
}
