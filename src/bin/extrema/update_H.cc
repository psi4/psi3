/*#######################################################
  UPDATE_H

  this function performs a bfgs update of the H matrix

  no parameters or returns
  #####################################################*/

#include <stdlib.h>
#include <stdio.h>

extern "C" {
  #include <libciomr.h>
}

#define EXTERN
#include "opt.h"

void update_H() {

  int i, j;
  double NUM1, NUM2, NUM12, *d_coord, *d_grad, *temp_arr1,
         **MAT1, **temp_mat1, **temp_mat2, **MAT2, **MAT3;

  /*allocate memory*/
  d_coord = init_array(num_coords);
  d_grad = init_array(num_coords);
  temp_arr1 = init_array(num_coords); 
  temp_mat1= init_matrix(num_coords,num_coords);
  temp_mat2 = init_matrix(num_coords,num_coords);
  MAT1 = init_matrix(num_coords,num_coords);
  MAT2 = init_matrix(num_coords,num_coords);
  MAT3 = init_matrix(num_coords,num_coords);

  /*the basic idea is to break up a nasty equation into pieces and put
    it the pieces together, don't even try to understand without the
    BFGS update equation in front of you, even then good luck*/

  for(i=0;i<num_coords;++i) {
      d_coord[i] = coord_vec[i] - coord_old[i];
      d_grad[i] = grad_vec[i] - grad_old[i];
    }

  fprintf(outfile,"\ncoordinate difference:\n");
  for(i=0;i<num_coords;++i) {
      fprintf(outfile,"%lf\n",d_coord[i]);
    }

  fprintf(outfile,"\ngradient difference:\n");
  for(i=0;i<num_coords;++i) {
      fprintf(outfile,"%lf\n",d_grad[i]);
    }

  for(i=0;i<num_coords;++i) {
      for(j=0;j<num_coords;++j) {
	  temp_arr1[i] += H_old[i][j] * d_grad[j];
	}
    }

  NUM1 = dot_pdt(d_grad,temp_arr1,num_coords);
  NUM2 = dot_pdt(d_coord,d_grad,num_coords);

  for(i=0;i<num_coords;++i) {
      for(j=0;j<num_coords;++j) {
        MAT1[i][j] = d_coord[i] * d_coord[j];
	}
    }

  for(i=0;i<num_coords;++i) {
      for(j=0;j<num_coords;++j) {
	temp_mat1[i][j] = d_coord[i] * d_grad[j];
	}
    }

  for(i=0;i<num_coords;++i) {
      for(j=0;j<num_coords;++j) {
	temp_mat2[i][j] = d_grad[i] * d_coord[j];
	}
    }

  mmult(temp_mat1,0,H_old,0,MAT2,0,num_coords,num_coords,num_coords,0);

  mmult(H_old,0,temp_mat2,0,MAT3,0,num_coords,num_coords,num_coords,0);

  NUM12 = 1/NUM2 + NUM1/(NUM2 * NUM2);
  
  /*finally put the pieces together*/
  for(i=0;i<num_coords;++i) {
      for(j=0;j<num_coords;++j) {

        H[i][j] = H_old[i][j] + NUM12 * MAT1[i][j]
		  - (MAT2[i][j] + MAT3[i][j]) / NUM2;
	  
	}
    }

  /*free up memory*/
  free(d_coord);
  free(d_grad);
  free(temp_arr1);
  free_matrix(temp_mat1, num_coords);
  free_matrix(temp_mat2, num_coords);
  free_matrix(MAT1, num_coords);
  free_matrix(MAT2, num_coords);
  free_matrix(MAT3, num_coords);

  return;
}















































  
  

  

  


