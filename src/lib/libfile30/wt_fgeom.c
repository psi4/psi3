/*!
  \file wt_fgeom.c
*/

#include <stdio.h>
#include <stdlib.h>
#include "file30.h"
#include "file30.gbl"
#include <libciomr.h>
 
/*!
** file30_wt_fgeom(double **fgeom, int num_atoms)
**   : Write the full cartesian geometry including dummy atoms to file30
**
**  arguments:
** \param double **fgeom, the full geometry matrix 
** \param int num_atoms
**
**  returns: none
*/
 
void file30_wt_fgeom(double **fgeom) {

  double *temp_geom;
  int num_entries,i;
  PSI_FPTR junk;
  PSI_FPTR f_ptr;

  num_entries = file30_rd_nentry();
  f_ptr = (PSI_FPTR) (info30_.mpoint[50] - 1) * sizeof(int);

  temp_geom = init_array(3*num_entries);

  for(i=0;i<num_entries;++i) {
	temp_geom[3*i] = fgeom[i][0];
	temp_geom[3*i+1] = fgeom[i][1];
	temp_geom[3*i+2] = fgeom[i][2];
  }

  wwritw(info30_.filenum,(char *) temp_geom, 	
	 (int) num_entries*sizeof(double)*3, f_ptr, &junk);

  free(temp_geom);

  return;
}
