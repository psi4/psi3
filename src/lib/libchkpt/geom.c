/*!
  \file geom.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>

/* chkpt_rd_geom(): Reads in the cartesian geometry from chkpt
**
**  takes no arguments.
**
**  returns: double **geom   The cartesian geometry is returned as a matrix
**     of doubles.  The row index is the atomic index, and the column is the
**     cartesian direction index (x=0, y=1, z=2).  Therefore, geom[2][0] 
**     would be the x-coordinate of the third atom.
** \ingroup (CHKPT)
*/


double **chkpt_rd_geom(void)
{
  double **geom;
  int natom;

  natom = chkpt_rd_natom();

  geom = block_matrix(natom, 3);

  psio_read_entry(PSIF_CHKPT, "::Geometry", (char *) geom[0], 
                  (int) 3*natom*sizeof(double));

  return geom;
}


/* chkpt_wt_geom(): Writes out the cartesian geometry to chkpt
**
** arguments: 
**  \param geom =  The cartesian geometry is supplied as a matrix
**     of doubles.  The row index is the atomic index, and the column is the
**     cartesian direction index (x=0, y=1, z=2).  Therefore, geom[2][0] 
**     would be the x-coordinate of the third atom.
** \ingroup (CHKPT)
*/

void chkpt_wt_geom(double **geom)
{
  int natom;

  natom = chkpt_rd_natom();

  psio_write_entry(PSIF_CHKPT, "::Geometry", (char *) geom[0], 
                   (int) 3*natom*sizeof(double));
}
