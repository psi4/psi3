/*##################################################################################
#
#  z_class.cc
#
#  definitions for the z-matrix derived class
#
##################################################################################*/

extern "C" {
#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include <file30.h>
#include <ip_libv1.h>
}

#define EXTERN
#include "simple_internal.h"
#include "coord_base.h" 
#include "z_class.h"
#include "opt.h"



/***********************************************************************************
*
*  z_class :: z_class
*
*  purpose: constructor for z-matrix derived class
*
*  parameters: none
*
*  returns: nothing
***********************************************************************************/

z_class :: z_class(int num_coord)
  : coord_base<simple_internal>(num_coord)
   {

  int i, j, pos;            /*counter variables, pos is position of next row of B matrix*/
  double *B_row0, *B_row1,  /*arrays to rows of the B matrix*/ 
         *B_row2,
         *cgrad_vec;

  /*read z_mat from file30*/
  file30_init();
  z_geom = file30_rd_zmat();

  /*print z-mat to output*/
  fprintf(outfile,"\nz-matrix:\n");
  for(i=0;i<num_atoms;++i) {
    if( i==0 ) 
      fprintf(outfile,"1\n");
    else if( i==1) 
      fprintf(outfile,"%i  %i  %lf\n",i+1,z_geom[i].bond_atom,z_geom[i].bond_val);
    else if( i==2 ) 
      fprintf(outfile,"%i  %i  %lf  %i  %lf\n",i+1,z_geom[i].bond_atom,
              z_geom[i].bond_val,z_geom[i].angle_atom,z_geom[i].angle_val); 
    else 
      fprintf(outfile,"%i  %i  %lf  %i  %lf  %i  %lf\n",i+1,z_geom[i].bond_atom,
              z_geom[i].bond_val,z_geom[i].angle_atom,z_geom[i].angle_val,
              z_geom[i].tors_atom,z_geom[i].tors_val);  
    }

  /*print opt flags*/
  fprintf(outfile,"\nopt flags:\n");
  for(i=1;i<num_atoms;++i) {
    if( i==1) 
      fprintf(outfile,"%i\n",z_geom[i].bond_opt);
    else if( i==2 ) 
      fprintf(outfile,"%i  %i\n",z_geom[i].bond_opt,z_geom[i].angle_opt); 
    else 
      fprintf(outfile,"%i  %i  %i\n",z_geom[i].bond_opt,z_geom[i].angle_opt,z_geom[i].tors_opt);
    }

  /*write z_mat to the array of simple_internal objects in coord_base class*/
  for(i=1;i<num_atoms;++i) {
      if( i==1 ) {
	  coord_arr[0].set(0,z_geom[1].bond_val);
	}
      else if( i==2 ) {
	  coord_arr[1].set(0,z_geom[2].bond_val);
	  coord_arr[2].set(1,z_geom[2].angle_val);
	  j=3;
	}
      else if( i>2 ) {
	  coord_arr[j].set(0,z_geom[i].bond_val);
	  coord_arr[j+1].set(1,z_geom[i].angle_val);
	  coord_arr[j+2].set(2,z_geom[i].tors_val);
	  j+=3;
	}
    }

  /*PRINT LEVEL?  check internals array*/
  fprintf(outfile,"\ninternal coordinate info:\n");
  for(i=0;i<num_coord;++i) {
      fprintf(outfile,"coord %i  type: %i  val: %lf\n",i,coord_arr[i].get_type(),coord_arr[i].get_val());
    }

  return;
 }


z_class::~z_class() {
  free(z_geom);
}












