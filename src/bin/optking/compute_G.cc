// This function computes G via BuB^t where u is a diagonal matrix
// of inverse masses.

extern "C" {
  #include <stdio.h>
  #include <stdlib.h>
  #include <math.h>
  #include <libciomr.h>
}

#define EXTERN
#include "opt.h"
#include "cartesians.h"

double **compute_G(double **B, int num_intcos, cartesians &carts) {
  double **u, **G, **temp_mat, *masses;
  int i, dim_carts;

  dim_carts = 3*carts.get_num_atoms();
  masses = carts.get_mass();

  G = init_matrix(num_intcos,num_intcos);
  temp_mat = init_matrix(dim_carts,num_intcos);
  u = init_matrix(dim_carts,dim_carts);

  for (i=0;i<3*carts.get_num_atoms();++i)
     u[i][i] = 1.0/masses[i];

  mmult(u,0,B,1,temp_mat,0,dim_carts,dim_carts,num_intcos,0);
  mmult(B,0,temp_mat,0,G,0,num_intcos,dim_carts,num_intcos,0);

  free(masses);
  free_matrix(u,dim_carts);
  free_matrix(temp_mat,dim_carts);

  return G;
}
