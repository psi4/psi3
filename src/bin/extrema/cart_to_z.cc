/*###################################################################
#
#  cart_to_z.cc
#
#  computes z-mat values for cartesian coordinates
#                                                     J.Kenny 10-15-00
####################################################################*/						      

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern "C" {
#include <libciomr.h>
#include <file30.h>
#include <physconst.h>
}

#define EXTERN
#include "simple_internal.h"
#include "coord_base.h"
#include "z_class.h"
#include "opt.h"

inline double *cross_pdt( double* vec1, double* vec2 );
inline double dot_pdt( double* vec1, double* vec2 );
inline double vec_normy( int atom1, int atom2 );
inline double *unit_vec( int atom1, int atom2 );

double *cart_to_z(struct z_class& z) {

    int i, j,pos=0;
    double *z_val, *temp1, *temp2, temp_num, div;

    z_val = init_array(num_coords);
    temp1 = init_array(3);
    temp2 = init_array(3);

    for(i=1;i<num_atoms;++i) {

	if(i==1) {
            fprintf(outfile,"\nn=1\n");
	    z_val[pos] = vec_normy(i,z.get_bond_atom(i)-1);
	    ++pos ;
	}

	if(i==2){
	    z_val[pos] = vec_normy(i,z.get_bond_atom(i)-1);
            temp_num = dot_pdt( unit_vec(z.get_bond_atom(i)-1,i), unit_vec(z.get_bond_atom(i)-1,z.get_angle_atom(i)-1) );
            fprintf(outfile,"\ntemp_num: %lf\n",temp_num);
	    z_val[pos+1] = acos( temp_num );
            pos += 2;
	}

	if(i>2) {
	    z_val[pos] = vec_normy(i,z.get_bond_atom(i)-1);
	    z_val[pos+1] = acos( dot_pdt(unit_vec(z.get_bond_atom(i)-1,i), unit_vec(z.get_bond_atom(i)-1,z.get_angle_atom(i)-1)));
	    temp1=cross_pdt(unit_vec(z.get_tors_atom(i)-1,z.get_angle_atom(i)-1),unit_vec(z.get_angle_atom(i)-1,z.get_bond_atom(i)-1));
	    temp2=cross_pdt(unit_vec(z.get_angle_atom(i)-1,z.get_bond_atom(i)-1),unit_vec(z.get_bond_atom(i)-1,i));
	    temp_num=dot_pdt(temp1,temp2);
            div=sin(acos(dot_pdt(unit_vec(z.get_bond_atom(i)-1,i),unit_vec(z.get_bond_atom(i)-1,z.get_angle_atom(i)-1))))*sin(acos(dot_pdt(unit_vec(z.get_angle_atom(i)-1,z.get_bond_atom(i)-1),unit_vec(z.get_angle_atom(i)-1,z.get_tors_atom(i)-1))));
	    z_val[pos+2]=acos(temp_num/div);
	    pos += 2;
	}
    }
   
    fprintf(outfile,"\nnew z-mat\n");
    for(i=0;i<num_coords;++i) {
	fprintf(outfile,"%lf\n",z_val[i]);
    }
    return z_val;
}								     



/**************************************************************
*
*  unit_vec, vec_normy, dot_pdt, cross_pdt
*
*  obvious math functions
*
***************************************************************/


  
inline double *unit_vec( int atom1, int atom2 ) {

  int i;
  double *unit_vec;

  unit_vec = init_array(3);

  for(i=0;i<3;++i) {
      fprintf(outfile,"%lf - %lf\n",cart_geom[atom2][i],cart_geom[atom1][i]);
      unit_vec[i] = (cart_geom[atom2][i]-cart_geom[atom1][i]) / vec_normy(atom1,atom2);
    }

  return unit_vec;
}



inline double vec_normy( int atom1, int atom2 ) {

  int i;
  double norm;

  norm = 0;
  for(i=0;i<3;++i) {
      norm += (cart_geom[atom2][i] - cart_geom[atom1][i]) * (cart_geom[atom2][i] - cart_geom[atom1][i]);
    }

  norm = sqrt(norm);

  fprintf(outfile,"\nnorm: %lf\n",norm);
  
  return norm;
}



inline double dot_pdt( double* vec1, double* vec2 ) {

  int i;
  double val;

  val=0;
  for(i=0;i<3;++i) {
      val += vec1[i] * vec2[i];
    }

  return val;
}
     
  
	   
inline double *cross_pdt( double* vec1, double* vec2 ) {

  double* result_vec;

  result_vec = init_array(3);

  result_vec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  result_vec[1] = -(vec1[0]*vec2[2] - vec1[2]*vec2[0]);
  result_vec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];

  return result_vec;
}













