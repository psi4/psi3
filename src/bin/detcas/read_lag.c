/*
** READ_LAG.C
**
** Read the lagrangian
**
** C. David Sherrill
** University of California, Berkeley
**
** April 1998
*/

#include <stdio.h>
#include <iwl.h>
#include <libciomr.h>
#include "globaldefs.h"
#include "globals.h"


void read_lagrangian(void)
{
  int nmo;
  PSI_FPTR lag_fptr=0;

  nmo = CalcInfo.nbfso;
  
  CalcInfo.lag = block_matrix(nmo, nmo);
  rfile(Params.lag_file);
  wreadw(Params.lag_file, (char *) CalcInfo.lag[0], nmo*nmo*sizeof(double),
         lag_fptr, &lag_fptr);
  if (Params.print_lvl > 3) {
    fprintf(outfile, "Lagrangian matrix\n");
    print_mat(CalcInfo.lag, nmo, nmo, outfile);
  }
  rclose(Params.lag_file, Params.lag_erase ? 4 : 3);

} 


