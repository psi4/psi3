/*!
  \file rref.c
  \ingroup (CHKPT)
*/

#include <stdio.h>
#include <stdlib.h>
#include "chkpt.h"
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/*!
** chkpt_rd_rref()  
** Reads in a 3x3 matrix used to rotate back to the reference frame.
**
**   takes no arguments.
**
**   returns: rref = A 3x3 matrix describing the rotation back to the 
**            reference frame.  The reference frame is a coordinate system 
**            defined by the "raw" geometry specification (either Z-matrix 
**            or geometry array in input.dat or chkpt). Can be used to 
**            transform quantities corresponding to different but similar 
**            calculations (gradients at displaced geometries) to a 
**            common frame.
**
** \ingroup (CHKPT)
*/

double **chkpt_rd_rref(void)
{
  char *key;
  double **Rref;
  char *keyword;
  keyword = chkpt_build_keyword("Transmat to reference frame");

  Rref = block_matrix(3,3);

  psio_read_entry(PSIF_CHKPT, keyword, (char *) Rref[0], sizeof(double)*9);

  free(keyword);
  return Rref;
}


/*!
** chkpt_wt_rref()
** Writes out a 3x3 matrix used to rotate back to the reference frame.
**
** \params rref = A 3x3 matrix describing the rotation back to the reference 
**                frame.  The reference frame is a coordinate system defined 
**                by the "raw" geometry specification (either Z-matrix or 
**                geometry array in input.dat or chkpt). Can be used to 
**                transform quantities corresponding to different but 
**                similar calculations (gradients at displaced geometries) 
**                to a common frame.
**
** returns: none
**
** \ingroup (CHKPT)
*/

void chkpt_wt_rref(double **Rref)
{
  char *keyword;
  keyword = chkpt_build_keyword("Transmat to reference frame");

  psio_write_entry(PSIF_CHKPT, keyword, (char *) Rref[0], sizeof(double)*9);

  free(keyword);
}
