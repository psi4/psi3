#include <stdio.h>
#include <stdlib.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"
#include <physconst.h>

#define _hartree2nm 0.02194746313710

void td_print(void)
{
  int i;

  fprintf(outfile,"\n\t                   Excitation Energy                OS       RS\n");
  fprintf(outfile,"\tState    (eV)    (cm^-1)    (nm)       (au)                 (au)\n");
  for(i=0; i<params.nstates; i++) {
    fprintf(outfile,"\t %d%3s %8.3lf %9.1lf %7.1lf %14.10lf %8.4lf %8.4lf\n",
            td_params[i].root+1,moinfo.labels[td_params[i].irrep],
            td_params[i].cceom_energy*_hartree2ev,
            td_params[i].cceom_energy*_hartree2wavenumbers,
            1/(td_params[i].cceom_energy*_hartree2nm),
            td_params[i].cceom_energy,td_params[i].OS,
            td_params[i].RS);
  }
  fprintf(outfile,"\n");
}
