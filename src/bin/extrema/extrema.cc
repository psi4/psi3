/*###################################################################
#
#  extrema.cc
#
#  main function for extrema
#                                                     J.Kenny 7-22-00
####################################################################*/						      
#include "extrema.h"

void set_up(coord_base* base);
void optimize_internals(internals* iobj);

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
      set_up(&zmat_geom);
      optimize_internals(&zmat_geom);

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




void optimize_internals(internals* iobj) {
    (*iobj).print_internals();
    (*iobj).compute_B();
    (*iobj).print_B();
    (*iobj).grad_trans();
    (*iobj).print_grads();
    switch(iteration) {
    case 1: (*iobj).initial_H(); break;
    default: (*iobj).update_bfgs(); break; }
    (*iobj).opt_step();
    (*iobj).back_transform();
    (*iobj).write_opt();
    (*iobj).write_file30();
    return;
}

















