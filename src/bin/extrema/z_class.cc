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
#include "opt.h"
#include "simple_internal.h"
#include "coord_base.h" 
#include "z_class.h"

extern double *B_row_bond( int atom1, int atom2 );
extern double *B_row_angle( int atom1, int atom2, int atom3 );
extern double *B_row_tors( int atom1, int atom2, int atom3, int atom4 );

/*need to declare this structure so that we can read zmat from file30
  (this is dictated by how file30 is written to, so don't change it)*/
struct z_entry {
  int bond_atom;                    /*first reference atom (bond)*/
  int angle_atom;                   /*second reference atom (angle)*/
  int tors_atom;                    /*third reference atom (torsion)*/
  int bond_opt;                     /*flag indicating to optimize variable (default=false)*/
  int angle_opt;
  int tors_opt;
  double bond_val;                  /*coordinate values*/
  double angle_val;
  double tors_val;
  };
            


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
         *B_row2;
  z_entry *z_geom;          /*array to hold z-matrix*/

  /*read z_mat from file30*/
  file30_init();
  z_geom = (z_entry *) malloc(num_atoms*sizeof(z_entry));
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

fflush(outfile);
  
  /*--------------
    form B matrix
    -------------*/

  B_mat = (double **) malloc(num_coord*sizeof(double *));
  B_row0 = init_array(3*num_atoms);
  B_row1 = init_array(3*num_atoms);
  B_row2 = init_array(3*num_atoms);
  pos = 0;

  for(i=1;i<num_atoms;++i) {
      if(i==1) {
	  B_row0 = B_row_bond(i, z_geom[i].bond_atom-1);
	  B_mat[pos] = B_row0;
	  ++pos;
	}
      if(i==2) {
	  B_row0 = B_row_bond(i, z_geom[i].bond_atom-1);
	  B_row1 = B_row_angle(i, z_geom[i].bond_atom-1, z_geom[i].angle_atom-1);
	  B_mat[pos] = B_row0;
	  B_mat[pos+1] = B_row1;
	  pos += 2;
	}
      if(i>2) {
	  B_row0 = B_row_bond(i, z_geom[i].bond_atom-1);
	  B_row1 = B_row_angle(i, z_geom[i].bond_atom-1, z_geom[i].angle_atom-1);
	  B_row2 = B_row_tors(i, z_geom[i].bond_atom-1, z_geom[i].angle_atom-1, z_geom[i].tors_atom-1);
	  B_mat[pos] = B_row0;
	  B_mat[pos+1] = B_row1;
	  B_mat[pos+2] = B_row2;
	  pos += 3;
	}
    }

  /*print B_mat*/
  print_mat(B_mat, num_coord, num_atoms*3, outfile);
  
 
  free(z_geom);
  return;
 }




