// This function performs the back-transformation. It computes a
// new cartesian geometry from an old cartesian geometry and a set
// of internal coordinate displacements

extern "C" {
  #include <stdio.h>
  #include <file30.h>
  #include <stdlib.h>
  #include <string.h>
  #include <math.h>
  #include <libciomr.h>
  #include <physconst.h>
}

#define EXTERN
#include "opt.h"
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

double *compute_q(internals &simples, salc_set &symm);
double **compute_B(int num_atoms, internals &simples, salc_set &symm);
double **compute_G(double **B, int num_intcos, cartesians &carts);

void new_geom(cartesians &carts, internals &simples, salc_set &symm,
       double *dq, int print_to_geom_file, int restart_geom_file,
       char *disp_label, int disp_num, int last_disp) {
  
  int bmat_iter_done,count,i,j,dim_carts,num_atoms;
  double **A, **G, **G_inv, **B, **u, **temp_mat;
  double dx_sum, dq_sum, *dx, *new_x, *x, *new_q, *q, *masses;

  num_atoms = carts.get_num_atoms();
  dim_carts = 3*num_atoms;

  dx = init_array(dim_carts);
  new_x = init_array(dim_carts);
  q = init_array(symm.get_num());
  new_q = init_array(symm.get_num());
  masses = init_array(3*num_atoms);
  masses = carts.get_mass();

  A = init_matrix(dim_carts,symm.get_num());
  G = init_matrix(symm.get_num(),symm.get_num());
  G_inv = init_matrix(symm.get_num(),symm.get_num());
  temp_mat = init_matrix(dim_carts,symm.get_num());
  u = init_matrix(dim_carts,dim_carts);
  for (i=0;i<3*num_atoms;++i)
    u[i][i] = 1.0/masses[i];

  x = carts.get_coord();
  scalar_mult(_bohr2angstroms,x,dim_carts); // x now holds geom in Ang

  // Compute B matrix -- Isn't this slick?
  simples.compute_internals(num_atoms,carts.get_coord());
  simples.compute_s(num_atoms,carts.get_coord());
  B = compute_B(num_atoms,simples,symm);

  q = compute_q(simples,symm);
  for (i=0;i<symm.get_num();++i)
     q[i] += dq[i];

  fprintf(outfile,"\nBack-transformation to cartesian coordinates...\n");
  fprintf(outfile," Iter   RMS Delta(dx)   RMS Delta(dq)\n");

// Start back transformation iterations
  bmat_iter_done = 0;
  count = 1;
  do {
    free_matrix(G,symm.get_num());
    G = compute_G(B,symm.get_num(),carts);
    G_inv = symm_matrix_invert(G,symm.get_num(),0,optinfo.redundant);

// BMAT computes G_inv only once like the following.
// OPTKING recomputes G_inv at each iteration, which
// is slower but gives better convergence.
//   if (count == 0) {
//     G_inv = symm_matrix_invert(G,symm.get_num(),0,optinfo.redundant);
//   }

  // u B^t G_inv = A
    mmult(B,1,G_inv,0,temp_mat,0,dim_carts,symm.get_num(),symm.get_num(),0);
    mmult(u,0,temp_mat,0,A,0,dim_carts,dim_carts,symm.get_num(),0);
  // A dq = dx
    mmult(A,0,&dq,1,&dx,1,dim_carts,symm.get_num(),1,0);

   dx_sum = dq_sum = 0.0;
   if (count == 1) {
     for (i=0;i<dim_carts;++i)
       dx_sum += dx[i]*dx[i];
     dx_sum = dx_sum / dim_carts;
     for (i=0;i<symm.get_num();++i)
       dq_sum += dq[i]*dq[i];
     dq_sum = dq_sum / ((double) symm.get_num());
     dx_sum = sqrt(dx_sum);
     dq_sum = sqrt(dq_sum);
     fprintf (outfile,"%5d %15.12lf %15.12lf\n", count, dx_sum, dq_sum);
   }

  // Compute new cart coordinates, then internals and B matrix
    for (i=0;i<dim_carts;++i)
       new_x[i] = (x[i] + dx[i]) / _bohr2angstroms;
    simples.compute_internals(num_atoms,new_x);
    simples.compute_s(num_atoms,new_x);
    free_matrix(B,symm.get_num());
    B = compute_B(num_atoms,simples,symm);

    free(new_q);
    new_q = compute_q(simples,symm);

    for (i=0;i<symm.get_num();++i)
      dq[i] = q[i] - new_q[i];

  // Test for convergence of iterations
   for (i=0;i<dim_carts;++i)
     dx_sum += dx[i]*dx[i];
   dx_sum = dx_sum / dim_carts;
   for (i=0;i<symm.get_num();++i)
     dq_sum += dq[i]*dq[i];
   dq_sum = dq_sum / ((double) symm.get_num());
   dx_sum = sqrt(dx_sum);
   dq_sum = sqrt(dq_sum);
   if ((dx_sum < optinfo.bt_dx_conv) && (dq_sum < optinfo.bt_dq_conv))
      bmat_iter_done = 1;
   fprintf (outfile,"%5d %15.12lf %15.12lf\n", count+1, dx_sum, dq_sum);

    for (i=0;i<dim_carts;++i)
      x[i] = new_x[i] * _bohr2angstroms;

    ++count;
  } while( (bmat_iter_done == 0) && (count < optinfo.bt_max_iter) );

  if (count >= optinfo.bt_max_iter) {
    fprintf(outfile,"Could not converge new geometry in %d iterations.",count);
    exit(2);
  }
  else
    fprintf(outfile,
      "\nConvergence to displaced geometry took %d iterations.\n",count);

// write geometry to output.dat
  cartesians cart_temp;
  cart_temp.set_coord(x);
  cart_temp.mult(1.0/_bohr2angstroms);
  fprintf(outfile,"\n%s\n",disp_label);
  cart_temp.print(1,outfile,0,disp_label,0);


// write geometry to file30 or geom.dat
  if (print_to_geom_file == 1) {
     FILE *fp_geom;
     if (restart_geom_file) {
        ffile(&fp_geom, "geom.dat",0);
        fprintf(fp_geom, "optking: (\n");
       }
     else
        ffile(&fp_geom, "geom.dat",1);
     cart_temp.print(4,fp_geom,restart_geom_file,disp_label,disp_num);
     if(last_disp) fprintf(fp_geom,")\n");
     fclose(fp_geom);
     fflush(outfile);
  }
  else {
     cart_temp.print(30,outfile,0,disp_label,disp_num);
  }

  free(masses);
  free(q); free(new_q);
  free(x); free(dx); free(new_x);
  free_matrix(A,dim_carts);
  free_matrix(G,symm.get_num());
  free_matrix(G_inv,symm.get_num());
  free_matrix(u,dim_carts);
  free_matrix(temp_mat,dim_carts);
  return;
}

