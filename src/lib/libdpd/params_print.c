#include <stdio.h>
#include "dpd.h"

int dpd_params_print(struct dpdparams *Params, FILE *outfile)
{
  int i;

  fprintf(outfile, "\n\tDPD Parameters:\n");
  fprintf(outfile,   "\t---------------\n");
  fprintf(outfile,   "\tpqnum = %d   rsnum = %d\n",
	  Params->pqnum, Params->rsnum);
  fprintf(outfile, "\t   Row and column dimensions for DPD Block:\n");
  fprintf(outfile, "\t   ----------------------------------------\n");
  for(i=0; i < Params->nirreps; i++)
      fprintf(outfile,   "\t   Irrep: %1d row = %5d\tcol = %5d\n", i,
	      Params->rowtot[i], Params->coltot[i]);
  fflush(outfile);
  
  return 0;
}
