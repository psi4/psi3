/*#############################################################################
  inline.h

  inline functions
						Joseph P. Kenny 11/29/01
  ###########################################################################*/

#include<math.h>

inline double *unit_vec(double* cart_arr, int atom1, int atom2 );
inline double vec_norm(double* cart_arr, int atom1, int atom2 );
inline double dot_pdt( double* vec1, double* vec2 );
inline double *cross_pdt( double* vec1, double* vec2 );
inline double norm(double* cart_arr, int atom1, int atom2 );


inline double *unit_vec(double* cart_arr, int atom1, int atom2 ) {

  int i;
  double *u_vec;

  u_vec = init_array(3);

  for(i=0;i<3;++i) {
      u_vec[i] = (cart_arr[3*atom2+i]-cart_arr[3*atom1+i]) 
	  / vec_norm(cart_arr,atom1,atom2);
    }

  return u_vec;
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



/*  double dot_pdt(double *vec1, double *vec2, int num) { */

/*    int i; */
/*    double result=0; */

/*    for(i=0;i<num;++i) { */
/*        result += vec1[i] * vec2[i]; */
/*      } */

/*    return result; */
/*  }      */
  

   
inline double *cross_pdt( double* vec1, double* vec2 ) {

  double* result_vec;

  result_vec = init_array(3);

  result_vec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
  result_vec[1] = -(vec1[0]*vec2[2] - vec1[2]*vec2[0]);
  result_vec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];

  return result_vec;
}



inline double norm(double* cart_arr, int atom1, int atom2 ) {

  int i;
  double norm;

  norm = 0;
  for(i=0;i<3;++i) {
      norm += ((cart_arr[3*atom2+i] - cart_arr[3*atom1+i]) * (cart_arr[3*atom2+i] - cart_arr[3*atom1+i]));
    }

  norm = sqrt(norm);

   return norm;
}





























