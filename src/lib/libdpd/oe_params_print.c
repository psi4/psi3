#include <stdio.h>
#include "dpd.h"

int dpd_oe_params_print(struct oe_dpdparams *Params, FILE *outfile)
{
  int i;

  fprintf(outfile,   "\tOE DPD Parameters:\n");
  fprintf(outfile,   "\t------------------\n");
  fprintf(outfile,   "\tpnum = %d   qnum = %d\n",Params->pnum, Params->qnum);
  fprintf(outfile,   "\tIrreps = %1d\n\n", Params->nirreps);
  fprintf(outfile, "\t   Row and column dimensions for DPD Block:\n");
  fprintf(outfile, "\t   ----------------------------------------\n");
  for(i=0; i < Params->nirreps; i++)
      fprintf(outfile,   "\t   Irrep: %1d row = %5d\tcol = %5d\n", i,
	      Params->rowtot[i], Params->coltot[i]);
  fflush(outfile);

  return 0;
}
