/*###################################################################
#
#  opt.cc
#
#  main function for extrema
#                                                     J.Kenny 7-22-00
####################################################################*/						      

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern "C" {
#include <libciomr.h>
#include <ip_libv1.h>
#include <physconst.h>
} 

#include "opt.h"

void input();
void read_file11(coord_base* coord);

void main() {

  input(); 

  typedef z_class type;
  if(coord_type!=0) 
      punt("Extrema currently capable of z-matrix optimization only");
  type coord_set;

  /*read previous iteration info from opt.dat*/
  iteration = coord_set.read_opt();
  fprintf(outfile,"\nStarting Optimization Iteration %d\n",iteration); 

  /*get cartesian coordinates and gradients*/
  read_file11(dynamic_cast<coord_base*>(&coord_set));
  coord_set.print_carts();
  coord_set.print_c_grads();

  /*compute B*/
  coord_set.compute_B();
  coord_set.print_B();

  /*transform gradients from cartesian to internal coordinates*/
  coord_set.grad_trans();

  /*form or update H*/
  coord_set.update_H();

  coord_set.opt_step();
  coord_set.back_transform(); 

  coord_set.write_opt();
  coord_set.write_file30();

  printf("\nNormal termination\n");
  exit(0);
}























