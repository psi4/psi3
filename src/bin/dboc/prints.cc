#include <stdio.h>
#include "params.h"

#define PRINT_INTRO 1
#define PRINT_PARAMS 1

extern FILE *outfile;
extern Params_t Params;

void print_intro()
{
  if (Params.print_lvl >= PRINT_INTRO) {
     fprintf(outfile,"                  ----------------------------------------------\n");
     fprintf(outfile,"                    DBOC: diagonal Born-Oppenheimer correction\n");
     fprintf(outfile,"                      evaluation by numerical differentiation\n");
     fprintf(outfile,"                  ----------------------------------------------\n\n");
  }
}

void print_params()
{
  if (Params.print_lvl >= PRINT_PARAMS) {
    fprintf(outfile,"\n  -OPTIONS:\n");
    fprintf(outfile,"    Print level                 = %d\n",Params.print_lvl);
    fprintf(outfile,"    Displacement size           = %lf a.u.\n",Params.delta);
    if (strcmp(Params.wfn,"SCF"))
      fprintf(outfile,"    Wave function               = %s\n",Params.wfn);
    else if (Params.reftype == rhf)
      fprintf(outfile,"    Wave function               = RHF SCF\n");
    else if (Params.reftype == rohf)
      fprintf(outfile,"    Wave function               = ROHF SCF\n");
    else if (Params.reftype == uhf)
      fprintf(outfile,"    Wave function               = UHF SCF\n");
  }
}
  
