#define EXTERN
#include "includes.h"
#include "globals.h"
#include "prototypes.h"
#include <libqt/qt.h>
#include <psifiles.h>
#include <ccfiles.h>

void get_opdm_lbl(void) { 
  int i;

  nrho = 1;

  if (cc_wfn(wfn)) {
    psio_open(CC_INFO,1);
    psio_read_entry(CC_INFO, "Num. of CC densities", (char *) &(nrho), sizeof(int));
    fprintf(outfile,"\tNumber of CC densities to analyze is: %d\n", nrho);
    psio_close(CC_INFO,1);
  } 

  if ( !strcmp(ref,"RHF") || !strcmp(ref,"ROHF") ) {
    opdm_lbl = (char **) malloc(sizeof(char *) * nrho);
    for (i=0;i<nrho;++i) {
      opdm_lbl[i] = (char *) malloc(32*sizeof(char));
      if (i==0) sprintf(opdm_lbl[i], "MO-basis OPDM");
      else sprintf(opdm_lbl[i], "MO-basis OPDM %d", i);
    }
  }
  else {
    opdm_a_lbl = (char **) malloc(sizeof(char *) * nrho);
    opdm_b_lbl = (char **) malloc(sizeof(char *) * nrho);
    for (i=0;i<nrho;++i) {
      opdm_a_lbl[i] = (char *) malloc(32*sizeof(char));
      opdm_b_lbl[i] = (char *) malloc(32*sizeof(char));
      if (i==0) {
        sprintf(opdm_a_lbl[i], "MO-basis Alpha OPDM");
        sprintf(opdm_b_lbl[i], "MO-basis Beta OPDM");
      }
      else {
        sprintf(opdm_a_lbl[i], "MO-basis Alpha OPDM %d", i);
        sprintf(opdm_b_lbl[i], "MO-basis Beta OPDM %d", i);
      }
    }
  }
	return;
}
