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

#include <stdio.h>
#include <stdlib.h>

extern "C" {
  #include <libciomr.h>
  #include <ip_libv1.h>
  #include <physconst.h>
}

#define EXTERN
#include "opt.h"

z_class :: z_class()
    : coord_base()
   {

  int i, j, pos;        

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
	  simple_arr[0].set_simple(0,z_geom[1].bond_val,2,z_geom[1].bond_atom,-1,-1);
	}
      else if( i==2 ) {
	  simple_arr[1].set_simple(0,z_geom[2].bond_val,3,z_geom[2].bond_atom,-1,-1);
	  simple_arr[2].set_simple(1,z_geom[2].angle_val *_pi/180.0,3,z_geom[2].bond_atom,z_geom[2].angle_atom,-1);
	  j=3;
	}
      else if( i>2 ) {
	  simple_arr[j].set_simple(0,z_geom[i].bond_val,i+1,z_geom[i].bond_atom,-1,-1);
	  simple_arr[j+1].set_simple(1,z_geom[i].angle_val*_pi/180.0,i+1,z_geom[i].bond_atom,z_geom[i].angle_atom,-1);
	  simple_arr[j+2].set_simple(2,z_geom[i].tors_val*_pi/180.0,i+1,z_geom[i].bond_atom,z_geom[i].angle_atom,z_geom[i].tors_atom);
	  j+=3;
	}
    }

  /*PRINT LEVEL?  check internals array*/
  fprintf(outfile,"\ninternal coordinate info:\n");
  for(i=0;i<num_simples;++i) {
      fprintf(outfile,"coord %i  type: %i  val: %lf\n",i,simple_arr[i].get_type(),simple_arr[i].get_val());
    }
  
  for(i=0;i<num_coords;++i) 
      coord_arr[i] = simple_arr[i].get_val();

  return;
}


/***********************************************************************************
*
*  z_class :: form_B
*
*  forms B matrix for simple internal coordinates
***********************************************************************************/
  double *B_row_bond(double* c_arr, int atom1, int atom2 );
  double *B_row_angle(double* c_arr, int atom1, int atom3, int atom2 );
  double *B_row_tors(double* c_arr, int atom1, int atom2, int atom3, int atom4 );

void z_class::compute_B() {  

    int i, j, pos=0;
    double *B_row0, *B_row1, *B_row2;
    B_row0 = init_array(3*num_atoms);
    B_row1 = init_array(3*num_atoms);
    B_row2 = init_array(3*num_atoms);

  for(i=1;i<num_atoms;++i) {
      if(i==1) {
	  B_row0 = B_row_bond(carts, i, simple_arr[pos].get_bond()-1);
	  B[pos] = B_row0;
          ++pos;
          fprintf(outfile,"\n1"); fflush(outfile);
	}
      if(i==2) {
	  B_row0 = B_row_bond(carts, i, simple_arr[pos].get_bond()-1);
	  B_row1 = B_row_angle(carts, i, simple_arr[pos+1].get_bond()-1, simple_arr[pos+1].get_angle()-1);
	  B[pos] = B_row0;
	  B[pos+1] = B_row1;
	  pos += 2;
          fprintf(outfile,"\n2"); fflush(outfile);
	}
      if(i>2) {
	  B_row0 = B_row_bond(carts, i, simple_arr[pos].get_bond()-1);
	  B_row1 = B_row_angle(carts, i, simple_arr[pos+1].get_bond()-1, simple_arr[pos+1].get_angle()-1);
	  B_row2 = B_row_tors(carts, i, simple_arr[pos+2].get_bond()-1, simple_arr[pos+2].get_angle()-1, simple_arr[pos+2].get_tors()-1);
	  B[pos] = B_row0;
	  B[pos+1] = B_row1;
	  B[pos+2] = B_row2;
	  pos += 3;
          fprintf(outfile,"\n%d",i); fflush(outfile);
	}
    }

  /*form u*/
  for(j=0;j<num_atoms;++j) {
      fprintf(outfile,"\nMASS: %lf",masses[j]);
      u[3*j][3*j] = 1.0 / masses[j]; 
      u[3*j+1][3*j+1] = 1.0 /masses[j];
      u[3*j+2][3*j+2] = 1.0 / masses[j];
  }

  return;
}
































