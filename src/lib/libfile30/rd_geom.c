#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>

/* file30_rd_geom(): Reads in the cartesian geometry from file30
**
**  takes no arguments.
**
**  returns: double **geom   The cartesian geometry is returned as a matrix
**     of doubles.  The row index is the atomic index, and the column is the
**     cartesian direction index (x=0, y=1, z=2).  Therefore, geom[2][0] 
**     would be the x-coordinate of the third atom.
*/


double **file30_rd_geom(void)
{
  double **geom;
  double *temp_geom;
  int natom, i;
  PSI_FPTR junk;
  PSI_FPTR geom_ptr;

  natom = file30_rd_natom();
  geom_ptr = (PSI_FPTR) (info30_.mcalcs[0] + 60 + 20 -1)*sizeof(int);

  temp_geom = init_array(3*natom);
  geom = block_matrix(natom, 3);

  wreadw(info30_.filenum, (char *) temp_geom, (int) 3*natom*sizeof(double),
	 geom_ptr, &junk);

  for(i=0; i < natom; i++) {
      geom[i][0] = temp_geom[3*i];
      geom[i][1] = temp_geom[3*i+1];
      geom[i][2] = temp_geom[3*i+2];
    }

  free(temp_geom);

  return geom;
}
