/*#############################################################################
  internals.cc

  derived internals class;
    that which is common to both zmatrix and delocalized coodinate systems
  ###########################################################################*/

#include <string.h>

#define EXTERN
#include "extrema.h"

void internals :: back_transform() {

  int i, j, pos;
  double conv=1.0, *dq, *dx;

  //print_u(); 

  dq = init_array(num_coords);
  dx = init_array(3*num_entries);         

  int loop=0;
  int converged=0;
  double criteria = 1.0e-14;
  double dx_sum; 

  fprintf(outfile,"\n\n  Performing iterative transformation to find");
  fprintf(outfile," new cartesian coordinates");
  fprintf(outfile,"\n\n  Iter    dq (internals)          dx (cartesians)");
  fprintf(outfile,  "\n  ---- ----------------------  ----------------------");  

  while((!converged) && (loop<bt_loop) ) {
  
      /*compute A*/
      compute_G(); 
      compute_A(); 

      for(i=0;i<num_coords;++i) {
	  dq[i] = coords[i] - coord_temp[i];
      }

      /*compute dx = A dq */
      for(i=0;i<3*num_entries;++i) 
          dx[i]=0;
      for(i=0;i<3*num_entries;++i) {
	  for(j=0;j<num_coords;++j) {
	      dx[i] += A[i][j] * dq[j];
	  }
      }

      pos=0;
      dx_sum = 0.0;
      for(i=0;i<3*num_entries;++i) {
	  carts[i] += dx[i];
          dx_sum += sqrt(dx[i]*dx[i]);
      }
      dx_sum /= (3*num_entries);
      
      cart_to_internal(coord_temp);

      compute_B();
    
      if(print_lvl >= RIDICULOUS_PRINT) {
        fprintf(outfile,"\n\n  dx:");
        for(i=0;i<(3*num_entries);++i)
          fprintf(outfile,"\n   dx[%d] = %lf",i,dx[i]);
        fprintf(outfile,"\n\n  dq:");
        for(i=0;i<num_coords;++i)
          fprintf(outfile,"\n   dq[%d] = %lf - %lf = %lf",i+1,coords[i],coord_temp[i],dq[i]);
        print_B();
        print_G();
        print_A();
      }

      conv=0;
      for(i=0;i<num_coords;++i) {
	  conv += sqrt((coords[i] - coord_temp[i])*(coords[i] -coord_temp[i]));
      }
      conv /= num_coords;
      fprintf(outfile,"\n  %3d  %.20lf  %.20lf",loop+1,conv,dx_sum);
       if( (conv<BT_CONV) && (dx_sum<BT_CONV) )
	  converged = 1;
      ++loop;
  }
  
  if(!converged) 
      punt("Back transformation to cartesians has failed");
  else
      fprintf(outfile,"\n  Back transformation to cartesians completed\n");
  return;


}



inline double *cross_pdt( double* vec1, double* vec2 );
inline double dot_pdt( double* vec1, double* vec2 );
inline double vec_norm(double* cart_arr, int atom1, int atom2 );
inline double *unit_vector(double* cart_arr, int atom1, int atom2 );

double* internals :: B_row_bond(double* c_arr, int atom1, int atom2 ) {

  int i;
  double *row, *unit12;

  row = init_array(num_entries*3);

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

double* internals :: B_row_angle(double* c_arr, int atom1, int atom3, int atom2 ) {

  int i;
  double *row, *unit31, *unit32, r31, r32, cosine, sine;
 
  row = init_array(num_entries*3);

  unit31 = unit_vector(c_arr,atom3, atom1);
  unit32 = unit_vector(c_arr,atom3, atom2);
  cosine = dot_pdt(unit31,unit32);
  sine = sin(acos(cosine));
  r31 = vec_norm(c_arr,atom3,atom1);
  r32 = vec_norm(c_arr,atom3,atom2);

  for(i=0;i<3;++i) {
      row[3*atom1+i] = ( cosine * unit31[i] - unit32[i] ) / ( r31 * sine );
      row[3*atom2+i] = ( cosine * unit32[i] - unit31[i] ) / ( r32 * sine );
      row[3*atom3+i] = ((r31 - r32 * cosine) * unit31[i] + (r32 - r31 * cosine)
			* unit32[i]) / (r31 * r32 * sine); 
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

double* internals :: B_row_tors(double* c_arr,int atom1, int atom2, int atom3, int atom4 ) {

  int i;
  double *row, *unit12, *unit23, *unit43, *unit32, *cross_12_23, *cross_43_32,
         cosine_an2, sine_an2, sine2_an2, cosine_an3, sine_an3, sine2_an3, r12, r23;

  row = init_array(num_entries*3);

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
     entries for atom4 same as those for atom1 with same permutations*/

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
      unit_vec[i] = (cart_arr[3*atom2+i]-cart_arr[3*atom1+i]) 
	  / vec_norm(cart_arr,atom1,atom2);
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


void internals :: print_B() {
        fprintf(outfile,"\n  B matrix:\n");
        print_mat(B,num_coords,3*num_entries,outfile);
        return;
    }

void internals :: compute_G() {
	double **temp1;
	temp1 = init_matrix(num_coords,3*num_entries);
	mmult(B,0,u,0,temp1,0,num_coords,3*num_entries,3*num_entries,0);
	mmult(temp1,0,B,1,G,0,num_coords,3*num_entries,num_coords,0);
        free_matrix(temp1,num_coords);
	return;
    }

void internals :: print_G() {
        fprintf(outfile,"\n  G matrix:\n");
        print_mat(G,num_coords,num_coords,outfile);
        return;
    }

void internals :: compute_A() {
	double **temp1, **temp2;
	temp1 = init_matrix(num_coords,num_coords);
	temp2 = init_matrix(3*num_entries,num_coords);
	/* A = u Bt G-1 */
	temp1 = symm_matrix_invert(G, num_coords, 0, 1);
	mmult(B,1,temp1,0,temp2,0,3*num_entries,num_coords,num_coords,0);
	mmult(u,0,temp2,0,A,0,3*num_entries,3*num_entries,num_coords,0);
        free_matrix(temp1,num_coords);
        free_matrix(temp2,3*num_entries);
	return;
    }

void internals :: print_A() {
        fprintf(outfile,"\n  A matrix:\n");
        print_mat(A,3*num_entries,num_coords,outfile);
        return;
    }

void internals :: grad_trans() {
	int i,j;
	for(i=0;i<num_coords;++i) {
            grads[i] = 0;
	    for(j=0;j<3*num_entries;++j) {
		grads[i] += B[i][j] * c_grads[j];
	    }}}


void internals :: optimize_internals(internals *icrd) {
 
    (*icrd).grad_test();
    (*icrd).print_internals();
    (*icrd).compute_B();
    if(print_lvl > RIDICULOUS_PRINT)
      (*icrd).print_B();
    (*icrd).grad_trans();
    switch(iteration) {
        case 1: (*icrd).initial_H(); break;
	default: 
	    if(!strcmp(update,"BFGS"))
		(*icrd).update_bfgs();
	    else if(!strcmp(update,"MS"))
		(*icrd).update_ms();
	    break; 
}
    (*icrd).diagonalize_H();
    (*icrd).print_H();
    (*icrd).opt_step();
    (*icrd).back_transform();
    (*icrd).print_carts(1.0);
    (*icrd).print_internals();
    (*icrd).write_opt();
    (*icrd).write_file30();
    return;
}
	    

    
    

