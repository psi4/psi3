/*##################################################################################
#
#  opt_step()
#
#  computes optimization steps
#
#  no parameters or returns
#
##################################################################################*/

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include <physconst.h>
}

#define EXTERN
#include "opt.h"


void opt_step() {

  int i, j;
  double *s;

  s = init_array(num_coords);

  for(i=0;i<num_coords;++i) {
      for(j=0;j<num_coords;++j) {
	  s[i] += -H[i][j] * grad_vec[j];
	}
    }

  for(i=0;i<num_coords;++i) {
      if( (fabs(s[i]) > 0.1) && (s[i] > 0.0) )
	  s[i]=0.1;
      if( (fabs(s[i]) > 0.1) && (s[i] < 0.0) )
	  s[i]=-0.1;
  }

  fprintf(outfile,"\nNew coordinate vector:\n");
  for(i=0;i<num_coords;++i) {
      coord_vec[i] += s[i];
      fprintf(outfile,"%lf\n",coord_vec[i]);
    }

  free(s);
  return;
}
