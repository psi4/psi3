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
#include <ip_libv1.h>
#include <file30.h>
#include <physconst.h>
}

#define EXTERN
#include "opt.h"

inline double norm(double* carts, int atom1, int atom2 );
inline double* unit_vec(double* carts, int atom1, int atom2 );
inline double dot_pdt( double* vec1, double* vec2 );
inline double *cross_pdt( double* vec1, double* vec2 );

double* z_class::cart_to_z() {

    int i, j,pos=0;
    double *z_val, *temp1, *temp2, temp_num, div, n1, n2;
    double tnum;

    z_val = init_array(num_coords);
    temp1 = init_array(3);
    temp2 = init_array(3);

    fprintf(outfile,"\n\n acos(1.0) = %.20lf\n\n",acos(1.0));

    for(i=1;i<num_atoms;++i) {


	if(i==1) {
	    z_val[pos] = norm(carts, i ,simple_arr[pos].get_bond()-1);
	    ++pos ;
	}


	if(i==2){

		z_val[pos] = norm( carts, i, simple_arr[pos].get_bond()-1);
		temp_num = dot_pdt( unit_vec(carts, simple_arr[pos+1].get_bond()-1, i ), 
				    unit_vec(carts, simple_arr[pos+1].get_bond()-1, simple_arr[pos+1].get_angle()-1 ) );
		z_val[pos+1] = acos( temp_num );
	    pos += 2;
	}


	if(i>2) {
	    
		z_val[pos] = norm( carts, i, simple_arr[pos].get_bond()-1);
		z_val[pos+1] = acos( dot_pdt( unit_vec( carts, simple_arr[pos+1].get_bond()-1, i ), 
					      unit_vec( carts, simple_arr[pos+1].get_bond()-1, simple_arr[pos+1].get_angle()-1 ) ) );

		temp1 = cross_pdt( unit_vec( carts, simple_arr[pos+2].get_angle()-1, simple_arr[pos+2].get_tors()-1 ), 
				   unit_vec( carts, simple_arr[pos+2].get_angle()-1, simple_arr[pos+2].get_bond()-1 ) );
		
		temp2 = cross_pdt( unit_vec( carts, simple_arr[pos+2].get_bond()-1, simple_arr[pos+2].get_angle()-1 ),
				   unit_vec( carts,simple_arr[pos+2].get_bond()-1,i ) );

		n1 = sqrt(dot_pdt(temp1,temp1));
		n2 = sqrt(dot_pdt(temp2,temp2));

		for(j=0;j<3;++j) {
		    temp1[j] /= n1;
		    temp2[j] /= n2;
		}

		temp_num = dot_pdt(temp1,temp2);

		if(temp_num>0.999999999999999)
		    z_val[pos+2] = 0.0;
		else if(temp_num<-0.999999999999999)
		    z_val[pos+2] = _pi;
		else 
		    z_val[pos+2] = acos(temp_num);
		pos += 3;
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
*  inline.h
*
*  unit_vec, norm, dot_pdt, cross_pdt
***************************************************************/


  
inline double norm(double* carts, int atom1, int atom2 ) {

  int i;
  double norm;

  norm = 0;
  for(i=0;i<3;++i) {
      norm += ((carts[3*atom2+i] - carts[3*atom1+i]) * (carts[3*atom2+i] - carts[3*atom1+i]));
    }

  norm = sqrt(norm);

//  fprintf(outfile,"\nnorm: %lf\n",norm);
  
  return norm;
}


inline double* unit_vec(double* carts, int atom1, int atom2 ) {

  int i;
  double *unit_vec;

  unit_vec = init_array(3);

  for(i=0;i<3;++i) {
      unit_vec[i] = ( carts[3*atom2+i]-carts[3*atom1+i] ) / norm(carts,atom1,atom2);
    }

  return unit_vec;
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

  
	   













