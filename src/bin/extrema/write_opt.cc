/*#######################################################################

  WRITE_OPT

  purpose: writes hessian, gradients and coordinate values to
     opt.dat so they can be used in next iteration

  no parameters or returns
  ######################################################################*/
  
#include <stdio.h>
#include <stdlib.h>

extern "C" {
#include <libciomr.h>
#include <ip_libv1.h>
}

#define EXTERN
#include "opt.h"

void write_opt() {

  int place, i, r;
  
  FILE *opt_ptr;
  
  ffile(&opt_ptr,"opt.dat",0);
  ip_set_uppercase(1);
  ip_initialize(opt_ptr,outfile);
  ip_cwk_add("OPT_INFO");

  fprintf(opt_ptr,"opt_info: (\n\n");

  /*write iteration number*/
  fprintf(opt_ptr,"  iteration = %d\n\n",iteration);

  /*write coordinate vector*/
  place = 0;
  fprintf(opt_ptr,"  coord = ( ");
  for (i=0;i<num_coords;++i) {
      if( place==8 ) {
	  place = 0;
	  fprintf(opt_ptr,"\n            ");
	}
      fprintf(opt_ptr,"%lf  ",coord_vec[i]);
      ++place;
    }
  fprintf(opt_ptr,")\n\n");

  /*write gradient vector*/
  place = 0;
  fprintf(opt_ptr,"  grad = ( ");
  for (i=0;i<num_coords;++i) {
      if( place==8 ) {
	  place = 0;
	  fprintf(opt_ptr,"\n            ");
	}
      fprintf(opt_ptr,"%lf  ",grad_vec[i]);
      ++place;
    }
  fprintf(opt_ptr,")\n\n");

  /*write Hessian one row at a time*/
  place=0;
  fprintf(opt_ptr,"  hmat = ( ");
  for(r=0;r<num_coords;++r) {
      if( place==0 )
	  fprintf(opt_ptr,"\n         ( ");
      for (i=0;i<num_coords;++i) {
          if( place==8 ) {
	      place = 0;
	      fprintf(opt_ptr,"\n           ");
	    }
          fprintf(opt_ptr,"%lf  ",H[r][i]);
          ++place;
        }
      fprintf(opt_ptr,")\n         ");
      place = 0;
    }
  fprintf(opt_ptr,")\n\n");

  fprintf(opt_ptr,"          )\n");

  fclose(opt_ptr);
  ip_done();
  return;
}




