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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
extern "C" {
#include <libciomr.h>
#include <ip_libv1.h>
}
#define EXTERN
#include "opt.h"

inline double *cross_pdt( double* vec1, double* vec2 );
inline double dot_pdt( double* vec1, double* vec2 );
inline double vec_norm(double* cart_arr, int atom1, int atom2 );
inline double *unit_vector(double* cart_arr, int atom1, int atom2 );

double *B_row_bond(double* c_arr, int atom1, int atom2 ) {

  int i;
  double *row, *unit12;

  row = init_array(num_atoms*3);

  unit12 = unit_vector(c_arr, atom1, atom2);

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

double *B_row_angle(double* c_arr, int atom1, int atom3, int atom2 ) {

  int i;
  double *row, *unit31, *unit32, r31, r32, cosine, sine;
  
  row = init_array(num_atoms*3);

  unit31 = unit_vector(c_arr,atom3, atom1);
  unit32 = unit_vector(c_arr,atom3, atom2);
  cosine = dot_pdt(unit31,unit32);
  sine = sin(acos(cosine));
  r31 = vec_norm(c_arr,atom3,atom1);
  r32 = vec_norm(c_arr,atom3,atom2);

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

double *B_row_tors(double* c_arr,int atom1, int atom2, int atom3, int atom4 ) {

  int i;
  double *row, *unit12, *unit23, *unit43, *unit32, *cross_12_23, *cross_43_32,
         cosine_an2, sine_an2, sine2_an2, cosine_an3, sine_an3, sine2_an3, r12, r23;

  row = init_array(num_atoms*3);

  unit12 = unit_vector(c_arr,atom1, atom2);
  unit23 = unit_vector(c_arr,atom2, atom3);
  unit43 = unit_vector(c_arr,atom4, atom3);
  unit32 = unit_vector(c_arr,atom3, atom2);
  cross_12_23 = cross_pdt(unit12, unit23);
  cross_43_32 = cross_pdt(unit43, unit32);
  cosine_an2 = dot_pdt( unit_vector(c_arr, atom2,atom1), unit_vector(c_arr, atom2,atom3));
  sine_an2 = sin(acos( cosine_an2 ));
  sine2_an2 = sine_an2 * sine_an2;
  cosine_an3 = dot_pdt( unit_vector(c_arr, atom3,atom2), unit_vector(c_arr, atom3,atom4));
  sine_an3 = sin(acos( cosine_an3 ));
  sine2_an3 = sine_an3 * sine_an3;
  r12 = vec_norm(c_arr,atom1,atom2);
  r23 = vec_norm(c_arr,atom2,atom3);

  for(i=0;i<3;++i) {
      row[3*atom1+i] = - cross_12_23[i] / ( r12 * sine2_an2 );
      row[3*atom2+i] = (r23 - r12 * cosine_an2) * cross_12_23[i] / ( r23 * r12 * sine2_an2 ) +
		       cosine_an3 * cross_43_32[i] / ( r23 * sine2_an3 );
    }

  
  /* entries for atom3 same as those for atom2 with permutation of 1 with 4 and 2 with 3,
     entries for atom4 same as those for atom1 with same permutations,
     (this may be lazy, but I write and debug less code this way) */

  /* I permute everything but the angles here */
  unit12 = unit_vector(c_arr,atom4, atom3);
  unit23 = unit_vector(c_arr,atom3, atom2);
  unit43 = unit_vector(c_arr,atom1, atom2);
  unit32 = unit_vector(c_arr,atom2, atom3);
  cross_12_23 = cross_pdt(unit12, unit23);
  cross_43_32 = cross_pdt(unit43, unit32);
  r12 = vec_norm(c_arr,atom4,atom3);
  r23 = vec_norm(c_arr,atom3,atom2);

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


  
inline double *unit_vector(double* cart_arr, int atom1, int atom2 ) {

  int i;
  double *unit_vec;

  unit_vec = init_array(3);

  for(i=0;i<3;++i) {
      unit_vec[i] = (cart_arr[3*atom2+i]-cart_arr[3*atom1+i]) / vec_norm(cart_arr,atom1,atom2);
    }

  return unit_vec;
}



inline double vec_norm(double* cart_arr, int atom1, int atom2 ) {

  int i;
  double norm, temp1, temp2;

  norm = 0;
  for(i=0;i<3;++i) {
      temp1 = cart_arr[3*atom1+i];
      temp2 = cart_arr[3*atom2+i];
      norm += (temp2 - temp1) * (temp2 - temp1);
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
