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
  char *junk;

  errcod = ip_string("REFERENCE", &(junk),0);
  if(!strcmp(junk, "RHF")) params.ref = 0;
  else if(!strcmp(junk, "ROHF")) params.ref = 1;
  else if(!strcmp(junk, "UHF")) params.ref = 2;
  else { 
    printf("Invalid value of input keyword REFERENCE: %s\n", junk);
    exit(2); 
  }
  free(junk);

  errcod = ip_string("DERTYPE", &(junk),0);
  if(!strcmp(junk,"NONE")) params.dertype = 0;
  else if(!strcmp(junk,"FIRST")) params.dertype = 1;
  else {
    printf("Invalid value of input keyword DERTYPE: %s\n", junk);
    exit(2); 
  }
  free(junk);

  params.print_lvl = 1;
  errcod = ip_data("PRINT_LVL","%d",&(params.print_lvl),0);

  params.keep_TEIFile = 1;
  errcod = ip_boolean("KEEP_TEIFILE",&(params.keep_TEIFile),0);

  params.keep_OEIFile = 1;
  errcod = ip_boolean("KEEP_OEIFILE",&(params.keep_OEIFile),0);

  params.tolerance = 1e-14;
  errcod = ip_data("TOLERANCE", "%d", &(tol),0);
  if(errcod == IPE_OK) params.tolerance = 1.0*pow(10.0,(double) -tol);

  fndcor(&(params.memory), infile, outfile);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);
}
