/*##################################################################
#
#  B_mat.cc
#
#  necessary functoins to compute the B matrix 
#
#  formulas can be found in Ch4, of Molecular Vibrations by
#  Wilson, Decious and Cross
#                                                   J.Kenny 7-22-00
###################################################################*/

extern "C"{
  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <libciomr.h>
  #include <physconst.h>
}

#define EXTERN
#include "simple_internal.h"
#include "coord_base.h"
#include "z_class.h"
#include "opt.h"

double *unit_vector( int atom1, int atom2 );
double dot_pdt( double* vec1, double* vec2 );
double *cross_pdt( double* vec1, double* vec2);
double vec_norm( int atom1, int atom2 ); 
double *B_row_bond(const int atom1, const int atom2 );
double *B_row_angle(const int atom1, const int atom2, const int atom3 );
double *B_row_tors(const int atom1, const int atom2, const int atom3, const int atom4 );



/********************************************************************
*
*  form_B
*
*  calls b_row functions to form B matrix
*
*  parameters: none
*
*  returns: nothing
*********************************************************************/

void form_B(struct z_class& z) {

  int i, pos=0;
  int temp_val;                       //needed to keep z const
  double *B_row0, *B_row1, *B_row2;
  B_row0 = init_array(3*num_atoms);
  B_row1 = init_array(3*num_atoms);
  B_row2 = init_array(3*num_atoms);

  for(i=1;i<num_atoms;++i) {
      if(i==1) {
	  B_row0 = B_row_bond(i, z.get_bond_atom(i)-1);
	  B_mat[pos] = B_row0;
	  ++pos;
	}
      if(i==2) {
	  B_row0 = B_row_bond(i, z.get_bond_atom(i)-1);
	  B_row1 = B_row_angle(i, z.get_bond_atom(i)-1, z.get_angle_atom(i)-1);
	  B_mat[pos] = B_row0;
	  B_mat[pos+1] = B_row1;
	  pos += 2;
	}
      if(i>2) {
	  B_row0 = B_row_bond(i, z.get_bond_atom(i)-1);
	  B_row1 = B_row_angle(i, z.get_bond_atom(i)-1, z.get_angle_atom(i)-1);
	  B_row2 = B_row_tors(i, z.get_bond_atom(i)-1, z.get_angle_atom(i)-1, z.get_tors_atom(i)-1);
	  B_mat[pos] = B_row0;
	  B_mat[pos+1] = B_row1;
	  B_mat[pos+2] = B_row2;
	  pos += 3;
	}
    }
}



/********************************************************************
*
*  B_row_bond
*
*  calculates s_vectors for a bond length and places them in proper
*  position of row vector which eventually becomes row of B matrix
*
*  parameters: atom1, atom2
*
*  returns: *row -- pointer to the row vector
*********************************************************************/

double *B_row_bond(const int atom1, const int atom2 ) {

  int i;
  double *row, *unit12;

  row = init_array(num_atoms*3);

  unit12 = unit_vector(atom1, atom2);

  for(i=0;i<3;++i) {
      row[3*atom1+i] = -unit12[i];
      row[3*atom2+i] =  unit12[i];
    }
  return row;
}



/*******************************************************************
*
*  B_row_angle
*
*  same as above but for dihedral angles
*
********************************************************************/

double *B_row_angle(const int atom1, const int atom3, const int atom2 ) {

  int i;
  double *row, *unit31, *unit32, r31, r32, cosine, sine;
  
  row = init_array(num_atoms*3);

  unit31 = unit_vector(atom3, atom1);
  unit32 = unit_vector(atom3, atom2);
  cosine = dot_pdt(unit31,unit32);
  sine = sin(acos(cosine));
  r31 = vec_norm(atom3,atom1);
  r32 = vec_norm(atom3,atom2);

  for(i=0;i<3;++i) {
      row[3*atom1+i] = ( cosine * unit31[i] - unit32[i] ) / ( r31 * sine );
      row[3*atom2+i] = ( cosine * unit32[i] - unit31[i] ) / ( r32 * sine );
      row[3*atom3+i] = ((r31 - r32 * cosine) * unit31[i] + (r32 - r31 * cosine) * unit32[i]) / (r31 * r32 * sine); 
    }

  return row;
} 



/*******************************************************************
*
*  B_row_tors
*
*  same as above but for torsional angles
*
********************************************************************/

double *B_row_tors(const int atom1, const int atom2, const int atom3, const int atom4 ) {

  int i;
  double *row, *unit12, *unit23, *unit43, *unit32, *cross_12_23, *cross_43_32,
         cosine_an2, sine_an2, sine2_an2, cosine_an3, sine_an3, sine2_an3, r12, r23;

  row = init_array(num_atoms*3);

  unit12 = unit_vector(atom1, atom2);
  unit23 = unit_vector(atom2, atom3);
  unit43 = unit_vector(atom4, atom3);
  unit32 = unit_vector(atom3, atom2);
  cross_12_23 = cross_pdt(unit12, unit23);
  cross_43_32 = cross_pdt(unit43, unit32);
  cosine_an2 = dot_pdt( unit_vector(atom2,atom1), unit_vector(atom2,atom3));
  sine_an2 = sin(acos( cosine_an2 ));
  sine2_an2 = sine_an2 * sine_an2;
  cosine_an3 = dot_pdt( unit_vector(atom3,atom2), unit_vector(atom3,atom4));
  sine_an3 = sin(acos( cosine_an3 ));
  sine2_an3 = sine_an3 * sine_an3;
  r12 = vec_norm(atom1,atom2);
  r23 = vec_norm(atom2,atom3);

  for(i=0;i<3;++i) {
      row[3*atom1+i] = - cross_12_23[i] / ( r12 * sine2_an2 );
      row[3*atom2+i] = (r23 - r12 * cosine_an2) * cross_12_23[i] / ( r23 * r12 * sine2_an2 ) +
		       cosine_an3 * cross_43_32[i] / ( r23 * sine2_an3 );
    }

  
  /* entries for atom3 same as those for atom2 with permutation of 1 with 4 and 2 with 3,
     entries for atom4 same as those for atom1 with same permutations,
     (this may be lazy, but I write and debug less code this way) */

  /* I permute everything but the angles here */
  unit12 = unit_vector(atom4, atom3);
  unit23 = unit_vector(atom3, atom2);
  unit43 = unit_vector(atom1, atom2);
  unit32 = unit_vector(atom2, atom3);
  cross_12_23 = cross_pdt(unit12, unit23);
  cross_43_32 = cross_pdt(unit43, unit32);
  r12 = vec_norm(atom4,atom3);
  r23 = vec_norm(atom3,atom2);

  /* I permute the angles here */
  for(i=0;i<3;++i) {
      row[3*atom4+i] = - cross_12_23[i] / ( r12 * sine2_an3 );
      row[3*atom3+i] = (r23 - r12 * cosine_an3) * cross_12_23[i] / ( r23 * r12 * sine2_an3 ) +
		       cosine_an2 * cross_43_32[i] / ( r23 * sine2_an2 );
    }
  
  return row;
} 



  
/**************************************************************
*
*  unit_vector, vec_norm, dot_pdt, cross_pdt
*
*  obvious math functions
*
***************************************************************/


  
inline double *unit_vector( int atom1, int atom2 ) {

  int i;
  double *unit_vec;

  unit_vec = init_array(3);

  for(i=0;i<3;++i) {
      unit_vec[i] = (cart_geom[atom2][i]*_bohr2angstroms-cart_geom[atom1][i]*_bohr2angstroms) / vec_norm(atom1,atom2);
    }

  return unit_vec;
}



inline double vec_norm( int atom1, int atom2 ) {

  int i;
  double norm;

  norm = 0;
  for(i=0;i<3;++i) {
      norm += (cart_geom[atom2][i] - cart_geom[atom1][i]) * (cart_geom[atom2][i] - cart_geom[atom1][i])*_bohr2angstroms*_bohr2angstroms;
    }

  norm = sqrt(norm);
  
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


