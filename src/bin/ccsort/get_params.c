#include <stdio.h>
#include <math.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <psifiles.h>
#include "Params.h"
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

void get_params()
{
  int i, errcod, tol;

  params.print_lvl = 1;
  errcod = ip_data("PRINT_LVL","%d",&(params.print_lvl),0);
  params.TEIFile = PSIF_MO_TEI;
  errcod = ip_data("TEI_FILE","%d",&(params.TEIFile),0);
  params.keep_TEIFile = 1;
  errcod = ip_boolean("KEEP_TEIFILE",&(params.keep_TEIFile),0);
  params.OEIFile = PSIF_MO_OEI;
  errcod = ip_data("OEI_FILE","%d",&(params.OEIFile),0);
  params.keep_OEIFile = 1;
  errcod = ip_boolean("KEEP_OEIFILE",&(params.keep_OEIFile),0);
  params.FZCFile = PSIF_MO_FZC;
  errcod = ip_data("FZC_FILE","%d",&(params.FZCFile),0);

  params.tolerance = 1e-14;
  errcod = ip_data("TOLERANCE", "%d", &(tol),0);
  if(errcod == IPE_OK) params.tolerance = 1.0*pow(10.0,(double) -tol);

  fndcor(&(params.memory), infile, outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);
}

