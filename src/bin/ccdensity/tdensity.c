#include <stdio.h>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

void ltdensity_rohf(struct RHO_Params S);
void ltdensity_uhf(struct RHO_Params S);
void ltdensity_intermediates(struct RHO_Params S);
void sort_ltd_rohf(struct RHO_Params S);
void sort_ltd_uhf(struct RHO_Params S);
void rtdensity(struct RHO_Params S);
void sort_rtd_rohf(struct RHO_Params S);
void sort_rtd_uhf(struct RHO_Params S);

void tdensity(struct RHO_Params S) {

  if(params.ref == 0 || params.ref == 1) {
    ltdensity_rohf(S);
    sort_ltd_rohf(S);
    rtdensity(S);
    sort_rtd_rohf(S);
  }
  else if(params.ref == 2) {
    ltdensity_uhf(S);
    sort_ltd_uhf(S);
    rtdensity(S);
    sort_rtd_uhf(S);
  }

  return;
}

