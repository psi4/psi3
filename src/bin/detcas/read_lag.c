/*
** READ_LAG.C
**
** Read the lagrangian
**
** C. David Sherrill
** University of California, Berkeley
**
** April 1998
** Updated to new libpsio libraries 8/03
*/

#include <stdlib.h>
#include <stdio.h>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include "globaldefs.h"
#include "globals.h"


void read_lagrangian(void)
{
  int nmo;
  PSI_FPTR lag_fptr=0;

  nmo = CalcInfo.nbfso;
  
  CalcInfo.lag = block_matrix(nmo, nmo);

  /*
  rfile(Params.lag_file);
  wreadw(Params.lag_file, (char *) CalcInfo.lag[0], nmo*nmo*sizeof(double),
         lag_fptr, &lag_fptr);
  */

  psio_open(Params.lag_file, PSIO_OPEN_OLD);  
  psio_read_entry(Params.lag_file, "MO-basis Lagrangian", 
    (char *) CalcInfo.lag[0], nmo*nmo*sizeof(double));

  if (Params.print_lvl > 3) {
    fprintf(outfile, "Lagrangian matrix\n");
    print_mat(CalcInfo.lag, nmo, nmo, outfile);
  }

  /* rclose(Params.lag_file, Params.lag_erase ? 4 : 3); */
  psio_close(Params.lag_file, Params.lag_erase ? 0 : 1);

} 


