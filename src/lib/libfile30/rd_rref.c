/*!
  \file rd_rref.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr/libciomr.h>

/*!
** file30_rd_rref():  Reads in a 3x3 matrix used to rotate back to the reference frame.
**
**   takes no arguments.
**
**   returns: double *z_geom          A 3x3 matrix describing the rotation back to the reference frame,
** 			              Reference frame is a coordinate system defined by the "raw"
**				      geometry specification (either Z-matrix or geometry array
**				      in input.dat or file30). Can be used to transform quantities
**				      corresponding to different but similar calculations
**				      (gradients at displaced geometries) to a common
**				      frame
*/

double **file30_rd_rref(void)
{
  PSI_FPTR Rref_ptr;
  double **Rref;

  Rref = block_matrix(3,3);
  Rref_ptr = (PSI_FPTR) (info30_.mpoint[49] - 1)*sizeof(int);

  wreadw(info30_.filenum, (char *) Rref[0], sizeof(double)*9, Rref_ptr, &Rref_ptr);

  return Rref;
}
