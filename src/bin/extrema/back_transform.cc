/*###################################################################
#
#  back_transform
#
#  performs iterative transformation from internals to cartesians
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

double *cart_to_z(struct z_class& z);
double **symm_matrix_invert(double **A, int dim, int print_det, int redundant);

void z_class::back_transform() {

    fprintf(outfile,"\nNew coords:");
    

  int i, j, pos;
  double conv=1.0, *dq, *dx;

  for(i=0;i<num_coords;++i) {
      fprintf(outfile,"\ncoord: %lf",coord_arr[i]);
  }

  print_u(); 

  dq = init_array(num_coords);
  dx = init_array(3*num_atoms);         

  int loop=0;
  int converged=0;
  double criteria = 0.000000000001;
  double dx_sum; 
  while((!converged) && (loop<500) ) {
  
      fprintf(outfile,"\nbt iteration %d\n", loop+1);
      
      /*compute A*/
      //print_B();
      compute_G(); //print_G();
      compute_A(); //print_A();

      for(i=0;i<num_coords;++i) {
	  fprintf(outfile,"%lf-%lf\n",coord_arr[i],coord_temp[i]);
	  dq[i] = coord_arr[i] - coord_temp[i];
      }
      //fprintf(outfile,"\nchanges in internals\n");
      //for(i=0;i<num_coords;++i){
//	  fprintf(outfile,"%lf\n",dq[i]);
//      }

      /*compute dx = A dq */
      for(i=0;i<3*num_atoms;++i) 
          dx[i]=0;
      for(i=0;i<3*num_atoms;++i) {
	  for(j=0;j<num_coords;++j) {
	      dx[i] += A[i][j] * dq[j];
	  }
      }
      fprintf(outfile,"\nchanges in cartesians\n");
      for(i=0;i<3*num_atoms;++i){
           fprintf(outfile,"%.20lf\n",dx[i]);
           }  

      pos=0;
      dx_sum = 0.0;
      for(i=0;i<3*num_atoms;++i) {
	  carts[i] += dx[i];
          dx_sum += sqrt(dx[i]*dx[i]);
      }
      dx_sum /= (3*num_atoms);
      
      fprintf(outfile,"\nintermediate cartesian coodinates\n"); print_carts();

      coord_temp = cart_to_z();

      for(i=0;i<num_coords;++i) {
	  fprintf(outfile,"\ncoord_temp: %lf",coord_temp[i]); 
	      simple_arr[i].set_val(coord_temp[i]);
      }

      compute_B();
      //fprintf(outfile,"\nintermediate B matrix:\n"); print_B();
    
      conv=0;
      for(i=0;i<num_coords;++i) {
	  conv += sqrt((coord_arr[i] - coord_temp[i])*(coord_arr[i] -coord_temp[i]));
      }
      fprintf(outfile,"\nConvergence = %.20lf\n",conv);
      fprintf(outfile,"\ndx_conv = %.20lf\n",dx_sum);
      if((conv<criteria)&&(dx_sum<criteria))
	  converged = 1;
      ++loop;
  }
  
  if(!converged) 
      punt("Back transformation to cartesians has failed");
  else
      fprintf(outfile,"\nBack transformation to cartesians completed\n");
  return;

}

  


//  /*----------------------------------------------------------------------------

//         **SYM_MATRIX_INVERT

//         inverts a matrix by diagonalization

//         parameters:
//               **A = matrix to be inverted
//               dim = dimension of A
//               print_det = print determinant if 1, nothing if 0
//               redundant = zero eigenvalues allowed if 1

//         returns:
//               **inv_A = inverse of A

//                                                        written by Rollin King
//  ----------------------------------------------------------------------------*/

//  double **symm_matrix_invert(double **A, int dim, int print_det, int redundant) {
//    int i;
//    double **A_inv, **A_vects, *A_vals, **A_temp, det=1.0;

//    A_inv   = init_matrix(dim,dim);
//    A_temp  = init_matrix(dim,dim);
//    A_vects = init_matrix(dim,dim);
//    A_vals  = init_array(dim);

//    sq_rsp(dim,dim,A,A_vals,1,A_vects,1E-10);

//    if (redundant == 0) {
//       for (i=0;i<dim;++i) {
//          det *= A_vals[i];
//          A_inv[i][i] = 1.0/A_vals[i];
//       }
//       if (print_det)
//          fprintf(outfile,"Determinant: %10.6e\n",det);
//       if (fabs(det) < 1E-10) {
//          fprintf(outfile,"Determinant: %10.6e\n",det);
//          fprintf(outfile,"Determinant is too small...aborting.\n");
//          fclose(outfile);
//          exit(2);
//       }
//    }
//    else {
//       for (i=0;i<dim;++i) {
//          det *= A_vals[i];
//          if (fabs(A_vals[i]) > 1E-10)
//             A_inv[i][i] = 1.0/A_vals[i];
//       }
//       if (print_det)
//          fprintf(outfile,"Determinant: %10.6e\n",det);
//    }

//    mmult(A_inv,0,A_vects,1,A_temp,0,dim,dim,dim,0);
//    mmult(A_vects,0,A_temp,0,A_inv,0,dim,dim,dim,0);

//    free(A_vals);
//    free_matrix(A_vects,dim);
//    free_matrix(A_temp,dim);
//    return A_inv;
//  }
