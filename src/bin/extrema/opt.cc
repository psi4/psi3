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
#include <file30.h>
}

#include "opt.h"
#include "simple_internal.h"
#include "coord_base.h"
#include "z_class.h"

void main() {

  int i, j;
  FILE *temp;

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
    else punt("Bad number of atoms or trouble reading from input.dat"); 

    /*allocate memory*/
    global_allocate();

    z_class z_coord(num_coords);
    ip_done();

    iteration = read_opt();

    fprintf(outfile,"\niteration: %d\n",iteration);
    
      if(iteration == 1) {
 
           for(i=0;i<num_coords;++i) {
               if( z_coord.coord_arr[i].get_type() == 0 )
                   H[i][i] = 1.0;
               if( z_coord.coord_arr[i].get_type() == 1 || z_coord.coord_arr[i].get_type() == 2 )
                   H[i][i] = 4.0;
               coord_vec[i] = z_coord.coord_arr[i].get_val();
	     }
 
           fprintf(outfile,"\nForming empirical Hessian:\n");
           print_mat(H,num_coords,num_coords,outfile);
       }

      else {
          for(i=0;i<num_coords;++i) {
             coord_vec[i] = z_coord.coord_arr[i].get_val();
           }
          update_H();
	}

      /*compute optimization step*/
      opt_step();

      fprintf(outfile,"\nnew coordinate values:\n");
      for(i=0;i<num_coords;++i) {
	  fprintf(outfile,"%lf\n",coord_vec[i]);
	}

      /*update z_geom with new coordinates*/
      j=0;
      for(i=1;i<num_atoms;++i) {
	  if(i==1) {
	      z_coord.z_geom[i].bond_val = coord_vec[j];
	      ++j; fprintf(outfile,"\nj: %d",j);
	    }
	  if(i==2) {
	      z_coord.z_geom[i].bond_val = coord_vec[j];
	      ++j; fprintf(outfile,"\nj: %d",j);
	      z_coord.z_geom[i].angle_val = coord_vec[j];
	      ++j; fprintf(outfile,"\nj: %d",j);
	    }
	  else if(i>2) {
	      z_coord.z_geom[i].bond_val = coord_vec[j];
	      ++j;
	      z_coord.z_geom[i].angle_val = coord_vec[j];
	      ++j;
	      z_coord.z_geom[i].tors_val - coord_vec[j];
	      ++j;
	    }
	}

      /*test z_to_cart*/
      double **cart_geom;
      cart_geom = file30_z_to_cart(z_coord.z_geom,num_atoms);
      print_mat(cart_geom,num_atoms,3,outfile);
      fflush(outfile);

      file30_wt_geom(cart_geom);
       
      write_opt();
      
    }
  else punt("z-matrix not found in input ... can't do whatever it is you want to do yet");

  global_free();
  printf("\nNormal termination\n");
  exit(0);
  
}





