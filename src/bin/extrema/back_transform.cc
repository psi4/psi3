/*###################################################################
#
#  back_transform
#
#  performs iterative transformation from internals to cartesians
#                                                     J.Kenny 10-15-00
####################################################################*/ 	          

#include <stdio.h>
#include <stdlib.h>

extern "C" {
#include <libciomr.h>
#include <file30.h>
#include <masses.h>
#include <physconst.h>
}

#define EXTERN
#include "simple_internal.h"
#include "coord_base.h"
#include "z_class.h"
#include "opt.h"

void form_B(struct z_class& z);
double *cart_to_z(struct z_class& z);
double **symm_matrix_invert(double **A, int dim, int print_det, int redundant);

void back_transform(struct z_class& z) {

  int i, j, pos;
  double conv=1.0;
  double **A, **u, **C;         
  double **temp1, **temp2, **temp3, *dq, *dx, *old_cart;

  A = init_matrix(3*num_atoms,num_coords);
  u= init_matrix(3*num_atoms,3*num_atoms);
  C= init_matrix(3*num_atoms,num_coords);
  temp1 = init_matrix(num_coords,num_coords);
  temp2 = init_matrix(num_coords,num_coords);
  temp3 = init_matrix(num_coords,3*num_atoms);
  dq = init_array(num_coords);
  dx = init_array(3*num_atoms);
  old_cart = init_array(3*num_atoms);

  for(i=0;i<3*num_atoms;++i) {
      u[i][i]= 1.0; //an2masses[8];
      //u[i+3][i+3] = 1.0/ an2masses[1];
      //u[i+6][i+6] = 1.0/ an2masses[1];
  }
 
  /*need cartesians in angstroms*/
  for(i=0;i<num_atoms;++i) {
      for(j=0;j<3;++j) {
	  cart_geom[i][j] = cart_geom[i][j] * _bohr2angstroms;
      }
  }

  while(conv>0.000000001) {
  
      fprintf(outfile,"\nbt iteration");
      /*compute G_inv*/
      fprintf(outfile,"\nBmat:\n");
      print_mat(B_mat,num_coords,3*num_atoms,outfile);
      mmult(B_mat,0,u,0,temp3,0,num_coords,3*num_atoms,3*num_atoms,0);
      mmult(temp3,0,B_mat,1,temp1,0,num_coords,3*num_atoms,num_coords,0);
      fprintf(outfile,"\nG:\n");
      print_mat(temp1,num_coords,num_coords,outfile);
      temp2 = symm_matrix_invert(temp1, num_coords, 0, 1);
      fprintf(outfile,"\nG_inv:\n");
      print_mat(temp2,num_coords,num_coords,outfile);
      mmult(B_mat,1,temp2,0,A,0,3*num_atoms,num_coords,num_coords,0);
      mmult(u,0,A,0,C,0,3*num_atoms,3*num_atoms,num_coords,0);
      fprintf(outfile,"\nA:\n");
      print_mat(C,3*num_atoms,num_coords,outfile);

      for(i=0;i<num_coords;++i) {
          fprintf(outfile,"%lf-%lf\n",coord_vec[i],old_coord_vec[i]);
	  dq[i] = coord_vec[i] - old_coord_vec[i];
      }
      fprintf(outfile,"\nchanges in internals\n");
      for(i=0;i<num_coords;++i){
	  fprintf(outfile,"%lf\n",dq[i]);
      }

      /*compute dx = A dq */
      for(i=0;i<3*num_atoms;++i) 
          dx[i]=0;
      for(i=0;i<3*num_atoms;++i) {
	  for(j=0;j<num_coords;++j) {
	      dx[i] += C[i][j] * dq[j];
	  }
      }
      fprintf(outfile,"\nchanges in cartesians\n");
      for(i=0;i<3*num_atoms;++i){
	  fprintf(outfile,"%lf\n",dx[i]);
      }  

      pos=0;
      for(i=0;i<num_atoms;++i) {
	  cart_geom[i][0] += dx[pos];
          cart_geom[i][1] += dx[pos+1];
          cart_geom[i][2] += dx[pos+2];
	  pos += 3;
      }
      
      fprintf(outfile,"\nintermediate cartesian coodinates\n");
      print_mat(cart_geom,num_atoms,3,outfile);
      fflush(outfile);

      old_coord_vec = cart_to_z(z);

      pos=0;
      for(i=1;i<num_atoms;++i) {
	  if(i==1) {
	      z.set_bond(i,old_coord_vec[pos]);
	      ++pos;
	  }
	  if(i==2) {
	      z.set_bond(i,old_coord_vec[pos]);
              z.set_angle(i,old_coord_vec[pos+1]);
	      pos += 2;
	  }
	  if(i>2) {
	      z.set_bond(i,old_coord_vec[pos]);
              z.set_angle(i,old_coord_vec[pos+1]);
              z.set_angle(i,old_coord_vec[pos+2]);
              pos += 3;
	  }
      }

      form_B(z);
      fprintf(outfile,"\nintermediate B matrix:\n");
      print_mat(B_mat,num_coords,3*num_atoms,outfile);

      conv=0;
      for(i=0;i<num_coords;++i) {
	  conv += (coord_vec[i] - old_coord_vec[i])*(coord_vec[i] - old_coord_vec[i]);
      }
  }
  fprintf(outfile,"\nBack transformation to cartesians completed\n");
  return;
}
  


/*----------------------------------------------------------------------------

       **SYM_MATRIX_INVERT

       inverts a matrix by diagonalization

       parameters:
             **A = matrix to be inverted
             dim = dimension of A
             print_det = print determinant if 1, nothing if 0
             redundant = zero eigenvalues allowed if 1

       returns:
             **inv_A = inverse of A

                                                      written by Rollin King
----------------------------------------------------------------------------*/

double **symm_matrix_invert(double **A, int dim, int print_det, int redundant) {
  int i;
  double **A_inv, **A_vects, *A_vals, **A_temp, det=1.0;

  A_inv   = init_matrix(dim,dim);
  A_temp  = init_matrix(dim,dim);
  A_vects = init_matrix(dim,dim);
  A_vals  = init_array(dim);

  sq_rsp(dim,dim,A,A_vals,1,A_vects,1E-10);

  if (redundant == 0) {
     for (i=0;i<dim;++i) {
        det *= A_vals[i];
        A_inv[i][i] = 1.0/A_vals[i];
     }
     if (print_det)
        fprintf(outfile,"Determinant: %10.6e\n",det);
     if (fabs(det) < 1E-10) {
        fprintf(outfile,"Determinant: %10.6e\n",det);
        fprintf(outfile,"Determinant is too small...aborting.\n");
        fclose(outfile);
        exit(2);
     }
  }
  else {
     for (i=0;i<dim;++i) {
        det *= A_vals[i];
        if (fabs(A_vals[i]) > 1E-10)
           A_inv[i][i] = 1.0/A_vals[i];
     }
     if (print_det)
        fprintf(outfile,"Determinant: %10.6e\n",det);
  }

  mmult(A_inv,0,A_vects,1,A_temp,0,dim,dim,dim,0);
  mmult(A_vects,0,A_temp,0,A_inv,0,dim,dim,dim,0);

  free(A_vals);
  free_matrix(A_vects,dim);
  free_matrix(A_temp,dim);
  return A_inv;
}
