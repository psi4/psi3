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
#include <physconst.h>
}

#include "simple_internal.h"
#include "coord_base.h"
#include "z_class.h"
#include "opt.h"

void form_B(struct z_class& z);
void back_transform(struct z_class& z);

void main() {

  int i, j;
  double *cgrad_vec,*cart_vec, coord;
  FILE *temp;
  
  cart_vec = init_array(3*num_atoms);

  /*set up i/o stuff*/
  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);
 
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_clear();
  ip_cwk_add(":INPUT");
  
  /*read in cartesian coordinates and gradients from file11*/
  read_file11();
  for(i=0;i<num_atoms;++i){
      for(j=0;j<3;++j){
	  cart_geom[i][j] = cart_geom[i][j] /  _bohr2angstroms;
      }
  }

  int spot=0;
  for(i=0;i<num_atoms;++i) {
      for(j=0;j<3;++j) {
	  cart_vec[spot]=cart_geom[i][j];
	  ++spot;
      }
  }

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

    fflush(outfile);
    B_mat = (double **) malloc(num_coords*sizeof(double *));
    form_B(z_coord);
    fprintf(outfile,"\nB_mat:\n");
    print_mat(B_mat,num_coords,3*num_atoms,outfile);

    /*transform gradients to internal coordinates*/
    cgrad_vec = init_array(3*num_atoms);
 
    int pos=0;
    for (i=0;i<num_atoms;++i) {
        for(j=0;j<3;++j) {
            cgrad_vec[pos] = cart_grad[i][j];
            ++pos;
          }
      }
 
    for (i=0;i<num_coords;++i) {
        for (j=0;j<(3*num_atoms);++j) {
            grad_vec[i] += B_mat[i][j]*cgrad_vec[j];
          }
      }
 
    fprintf(outfile,"\ngradient vector in cartesian coordinates:\n");
    for (i=0;i<3*num_atoms;++i) {
        fprintf(outfile,"%lf\n",cgrad_vec[i]);
      }
 
    fprintf(outfile,"\ngradient vector in internal coordinates:\n");
    for (i=0;i<num_coords;++i) {
        fprintf(outfile,"%lf\n",grad_vec[i]);
      }

    free(cgrad_vec);

    iteration = read_opt();

    fprintf(outfile,"\niteration: %d\n",iteration);
    old_coord_vec = init_array(num_coords);    

    if(iteration == 1) {
 
        for(i=0;i<num_coords;++i) {
	    if( z_coord.coord_arr[i].get_type() == 0 ){
	        H[i][i] = 1.0;
		coord_vec[i] = old_coord_vec[i] = z_coord.coord_arr[i].get_val() * _bohr2angstroms;
	    }
	    if( z_coord.coord_arr[i].get_type() == 1 || z_coord.coord_arr[i].get_type() == 2 ){
	        H[i][i] = 4.0;
		coord_vec[i] = old_coord_vec[i] = z_coord.coord_arr[i].get_val() * _pi/180.0;
	    }
        }
 
        fprintf(outfile,"\nForming empirical Hessian:\n");
        print_mat(H,num_coords,num_coords,outfile);
    }

    else {
        for(i=0;i<num_coords;++i) {
	    coord_vec[i] = old_coord_vec[i] = z_coord.coord_arr[i].get_val()*_bohr2angstroms;
        }
        update_H();
    }

    write_opt();

    /*compute optimization step*/
    opt_step();

    fprintf(outfile,"\nnew coordinate values:\n");
    fprintf(outfile,"%lf\n",coord_vec[0]);
    fprintf(outfile,"%lf\n",coord_vec[1]);
    fprintf(outfile,"%lf\n",coord_vec[2]);


    /*update z_geom with new coordinates*/
    j=0;
    for(i=1;i<num_atoms;++i) {
        if(i==1) {
	    z_coord.z_geom[i].bond_val = coord_vec[j];
	    ++j;
        }
        if(i==2) {
	    z_coord.z_geom[i].bond_val = coord_vec[j];
	    ++j; 
	    z_coord.z_geom[i].angle_val = coord_vec[j];
	    ++j; 
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
    //double **cart_geom;
    //cart_geom = file30_z_to_cart(z_coord.z_geom,num_atoms);
    //print_mat(cart_geom,num_atoms,3,outfile);

    /*need to use iterative transform for now*/
    back_transform(z_coord);
    fprintf(outfile,"\nnew cartesian coordinates\n");
    print_mat(cart_geom,num_atoms,3,outfile);

    file30_wt_geom(cart_geom);
    file30_wt_zmat(z_coord.z_geom,num_atoms);
      
  }
  else punt("z-matrix not found in input ... can't do whatever it is you want to do yet");

  global_free();
  printf("\nNormal termination\n");
  exit(0);
  
}





