/*#######################################################################

  READ_OPT

  purpose: reads old hessian, gradients and coordinate values from
     opt.dat so they can be used in update

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

int read_opt() {

  FILE *opt_ptr;

  int i, j, error, num_elems, num_iter;
  
  ffile(&opt_ptr,"opt.dat",0);
  ip_set_uppercase(1);
  ip_initialize(opt_ptr,outfile);
  ip_cwk_add("OPT_INFO");

  if( ip_exist("ITERATION",0) ) {
 
      ip_data("ITERATION","%d",&num_iter,1,0);
      ++num_iter;
      
      /*read old coordinate vector*/
      error = 0;
      error += ( !ip_exist("COORD",0) );
      ip_count("COORD",&num_elems,0);
      error += ( num_elems != 1 );
      ip_count("COORD",&num_elems,1,0);
      error += ( num_elems != num_coords );

      for (i=0;i<num_coords;++i) {
	  error += ip_data("COORD","%lf",&coord_old[i],2,0,i);
	}
      
      if(error != 0)
	  punt("Problem reading old coordinate values from opt.dat");
  
      /*read old gradient vector*/
      error += ( !ip_exist("GRAD",0) );
      ip_count("GRAD",&num_elems,0);
      error += ( num_elems != 1 );
      ip_count("GRAD",&num_elems,1,0);
      error += ( num_elems != num_coords );

      for (i=0;i<num_coords;++i) {
	  error += ip_data("GRAD","%lf",&grad_old[i],2,0,i);
	}

      if(error != 0)
	  punt("Problem reading old gradient from opt.dat");
  
      /*read hmat, the inverse of the hessian*/
      error += (!ip_exist("HMAT",0));
      ip_count("HMAT",&num_elems,0);
      error += ( num_elems != num_coords );
      ip_count("HMAT",&num_elems,1,0);
      error += ( num_elems != num_coords );
      
      for (i=0;i<num_coords;++i) {
	  for (j=0;j<num_coords;++j) {
	      error += ip_data("HMAT","%lf", &H_old[i][j], 2, i, j);
	    }
	}
      
      if(error != 0)
	  punt("Problem reading old hessian from opt.dat");
    }
  else { num_iter = 1; }

  return num_iter;
}




