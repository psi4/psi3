/*###################################################################
#
#  extrema.cc
#
#  main function for extrema
#                                                     J.Kenny 7-22-00
####################################################################*/						      
#include "extrema.h"

void set_up(coord_base* base);

void main() {

  parsing();


/*------------
  CARTESIANS
  ----------*/
  if(coord_type==1) {
      
      carts cart_geom;
      set_up(&cart_geom);

  }
  

/*---------
  ZMATRIX
  -------*/
  else if(coord_type==2) {

      zmat zmat_geom;
      zmat_geom.read_carts();
      zmat_geom.print_carts(1.0);
      zmat_geom.read_file11();
      zmat_geom.print_c_grads();
      zmat_geom.read_opt();
      zmat_geom.optimize_internals((internals*) &zmat_geom);

  }


/*-------------
  DELOCALIZED
  -----------*/
  else if(coord_type==3) {

      // deloc deloc_geom;
      // coord_base* b = &deloc_geom;
      // set_up(b);

      //internals* obj = &deloc_geom;
      //optimize_internals(obj);
  }

  if(converged)
      fprintf(outfile,"\n  Optimization completed\n");
  file30_close();
  ip_done();
  tstop(outfile);
  fclose(infile);
  fclose(outfile);
  if(converged)
      exit(1);
  if(!converged)
      exit(0);
}




void set_up(coord_base* base) {
    (*base).read_carts();
    (*base).print_carts(1.0);
    (*base).read_file11();
    (*base).print_c_grads();
    (*base).read_opt();
    return;
}






















