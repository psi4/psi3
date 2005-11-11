#include <stdio.h>
#include <stdlib.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"
#include <physconst.h>

#define _hartree2nm 0.02194746313710

void td_print(struct RHO_Params *rho_params)
{
  int i;

  fprintf(outfile,"\n\t                   Excitation Energy                OS       RS\n");
  fprintf(outfile,"\tState    (eV)    (cm^-1)    (nm)       (au)                 (au)\n");
  for(i=0; i<params.nstates; i++) {
    fprintf(outfile,"\t %d%3s %8.3lf %9.1lf %7.1lf %14.10lf %8.4lf %8.4lf\n",
            rho_params[i].L_root+1,moinfo.labels[rho_params[i].L_irr],
            rho_params[i].cceom_energy*_hartree2ev,
            rho_params[i].cceom_energy*_hartree2wavenumbers,
            1/(rho_params[i].cceom_energy*_hartree2nm),
            rho_params[i].cceom_energy, params.OS[i],
            params.RS[i]);
  }
  fprintf(outfile,"\n");
}
