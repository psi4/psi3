/*###################################################################
#
#  opt.cc
#
#  main function for extrema
#                                                     J.Kenny 7-22-00
####################################################################*/						      

#include <stdio.h>
#include <stdlib.h>

extern "C" {
#include <libciomr.h>
#include <ip_libv1.h>
}

#include "opt.h"
#include "simple_internal.h"
#include "coord_base.h"
#include "z_class.h"  

void main() {

  int i;

  /*set up i/o stuff*/
  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);
 
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_clear();
  ip_cwk_add(":INPUT");

  /*read in cartesian coordinates and gradients from file11*/
  read_file11();

  /*if using z-matrix construct z_class*/
  if(ip_exist("ZMAT",0)){

    if( num_atoms==2 )
	num_coords = 1;
    else if( num_atoms>2 )
	num_coords = 3*num_atoms - 6;
    else punt("Bad number of atoms or trouble reading from file30"); 

    /*allocate memory*/
    global_allocate();

    fprintf(outfile,"\nH: %d", H);
    fprintf(outfile,"\nH[0]: %d", H[0]);
    fprintf(outfile,"\n&H[0][1]: %d", &H[0][1]);
    
    z_class z_coord(num_coords);
    ip_done();

    int iteration;
    iteration = read_opt();

      if(iteration == 1) {
 
       for(i=0;i<num_coords;++i) {
           if( z_coord.coord_arr[i].get_type() == 0 )
               H[i][i] = 1.00;
           if( z_coord.coord_arr[i].get_type() == 1 || z_coord.coord_arr[i].get_type() == 2 )
               H[i][i] = 0.25;
         }
 
       fprintf(outfile,"\nForming empirical Hessian:\n");
       print_mat(H,num_coords,num_coords,outfile);
         }

      }
  else punt("z-matrix not found in input ... can't do whatever it is you want to do yet");
  
  printf("\nNormal termination\n");
  exit(0);
  
}





