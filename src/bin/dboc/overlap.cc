#include <stdio.h>
#include <stdlib.h>
#include <string.h>
extern "C" {
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <psifiles.h>
}
#include "params.h"
#include "mo_overlap.h"

extern Params_t Params;
extern FILE *outfile;

extern void done(const char *);
extern double eval_rhf_derwfn_overlap();
extern double eval_rohf_derwfn_overlap();
extern double eval_uhf_derwfn_overlap();
extern double eval_rccsd_derwfn_overlap();
extern double eval_rci_derwfn_overlap();
extern double eval_roci_derwfn_overlap();

double eval_derwfn_overlap()
{
  double S;

  if (!strcmp(Params.wfn,"SCF")) {
    if (Params.reftype == Params_t::rhf) {
      S = eval_rhf_derwfn_overlap();
    }
    else if (Params.reftype == Params_t::rohf) {
      S = eval_rohf_derwfn_overlap();
    }
    else if (Params.reftype == Params_t::uhf) {
      S = eval_uhf_derwfn_overlap();
    }
    else
      done("This HF SCF method is not supported at the moment");
  }
  else if (!strcmp(Params.wfn,"CCSD")) {
    if (Params.reftype == Params_t::rhf) {
      //      S = eval_rccsd_derwfn_overlap();
    }
    else
      done("CCSD method with this reference is not supported at the moment");
  }
  else if (!strcmp(Params.wfn,"DETCI") || !strcmp(Params.wfn,"DETCAS")) {
    if (Params.reftype == Params_t::rhf) {
      S = eval_rci_derwfn_overlap();
    }
    else if (Params.reftype == Params_t::rohf) {
      S = eval_roci_derwfn_overlap();
    }
    else
      done("CI method with this reference is not supported at the moment");
  }

  return S;
}  
