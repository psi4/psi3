#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <ccfiles.h>
#include "moinfo.h"
#include "params.h"
#define EXTERN
#include "globals.h"

void get_params()
{
  int errcod;
  char *cachetype = NULL;
  
  errcod = ip_string("WFN", &(params.wfn), 0);

  errcod = ip_string("REFERENCE", &(params.ref),0);
  if (strcmp(params.ref,"RHF")) {
    fprintf(outfile, "\nIncorrect Reference: RHF only\n");
    abort();
  }
  
  params.print = 0;
  errcod = ip_data("PRINT", "%d", &(params.print),0);
  
  errcod = ip_boolean("OPDM", &(params.opdm),0);
  
  if (params.opdm) {
    params.opdm_write = 1;
  }  
  else {
    params.opdm_write = 0;
  }
  errcod = ip_boolean("OPDM_WRITE", &(params.opdm_write),0);
  
  params.opdm_print = 0;
  errcod = ip_boolean("OPDM_PRINT", &(params.opdm_print),0);

  params.cachelev = 2;
  errcod = ip_data("CACHELEV", "%d", &(params.cachelev),0);
  
  params.cachetype = 1;
  errcod = ip_string("CACHETYPE", &(cachetype),0);
  if (cachetype != NULL && strlen(cachetype)) {
    if (!strcmp(cachetype,"LOW")) 
      params.cachetype = 1;
    else if (!strcmp(cachetype,"LRU")) 
      params.cachetype = 0;
    else {
      fprintf(outfile, "Invalide CACHETYPE = %s\n",cachetype);
      abort();
    }
    free(cachetype);
  }
  
  fndcor(&(params.memory),infile,outfile);
 
  fprintf(outfile, "\n");
  fprintf(outfile, "\tInput parameters:\n");
  fprintf(outfile, "\t-----------------\n");
  fprintf(outfile, "\tWave function \t=\t%s\n", params.wfn);
  fprintf(outfile, "\tReference WFN \t=\t%s\n", params.ref);
  fprintf(outfile, "\tCache Level   \t=\t%d\n", params.cachelev);
  fprintf(outfile, "\tCache Type    \t=\t%s\n", params.cachetype ? "LOW":"LRU");
  fprintf(outfile, "\tMemory (MB)   \t=\t%.1f\n",params.memory/1e6);
  fprintf(outfile, "\tPrint Level   \t=\t%d\n", params.print);
  fprintf(outfile, "\tOPDM          \t=\t%s\n", params.opdm ? "YES":"NO");
  fprintf(outfile, "\tWrite OPDM    \t=\t%s\n", params.opdm_write ? "YES":"NO");
  fprintf(outfile, "\tPrint OPDM    \t=\t%s\n", params.opdm_print ? "YES":"NO");
}
