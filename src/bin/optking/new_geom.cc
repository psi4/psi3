// This function performs the back-transformation. It computes a
// new cartesian geometry from an old cartesian geometry and a set
// of internal coordinate displacements

extern "C" {
#include <stdio.h>
#include <libchkpt/chkpt.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <physconst.h>
}

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"

double *compute_q(internals &simples, salc_set &symm);
double **compute_B(internals &simples, salc_set &symm);
double **compute_G(double **B, int num_intcos, cartesians &carts);

void new_geom(cartesians &carts, internals &simples, salc_set &symm,
    double *dq, int print_flag, int restart_geom_file,
    char *disp_label, int disp_num, int last_disp, double *return_geom) {

  int bmat_iter_done,count,i,j,dim_carts,nallatom,natom;
  double **A, **G, **G_inv, **B, **u, **temp_mat, *no_fx;
  double dx_sum, dq_sum, *dx, *new_x, *x, *new_q, *q, *masses, *fcoord;

  nallatom = optinfo.nallatom;
  natom = optinfo.natom;
  dim_carts = 3*optinfo.nallatom;

  dx = init_array(dim_carts);
  new_x = init_array(dim_carts);
  new_q = init_array(symm.get_num());

  masses = carts.get_fmass();
  u = mass_mat(masses);
  free(masses);

  A = block_matrix(dim_carts,symm.get_num());
  G = block_matrix(symm.get_num(),symm.get_num());
  G_inv = block_matrix(symm.get_num(),symm.get_num());
  temp_mat = block_matrix(dim_carts,symm.get_num());

  x = carts.get_fcoord();
  scalar_mult(_bohr2angstroms,x,dim_carts); // x now holds geom in Ang

  // Compute B matrix -- Isn't this slick?
  fcoord = carts.get_fcoord();
  simples.compute_internals(nallatom,fcoord);
  simples.compute_s(nallatom,fcoord);
  free(fcoord);

  B = compute_B(simples,symm);
  q = compute_q(simples,symm);
  for (i=0;i<symm.get_num();++i)
    q[i] += dq[i];

  fprintf(outfile,"\nBack-transformation to cartesian coordinates...\n");
  fprintf(outfile," Iter   RMS Delta(dx)   RMS Delta(dq)\n");

  // Start back transformation iterations
  bmat_iter_done = 0;
  count = 1;
  do {
    free_block(G);
    free_block(G_inv);
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

    if (count == 1) {
      // dx_sum = RMS change in cartesian coordinates 
      dx_sum = dq_sum = 0.0;
      for (i=0;i<dim_carts;++i)
        dx_sum += dx[i]*dx[i];
      dx_sum = sqrt(dx_sum / dim_carts);

      // dq_sum = RMS change in internal coordinates
      for (i=0;i<symm.get_num();++i)
        dq_sum += dq[i]*dq[i];
      dq_sum = sqrt(dq_sum / ((double) symm.get_num()));
      fprintf (outfile,"%5d %15.12lf %15.12lf\n", count, dx_sum, dq_sum);
    }

    // Compute new cart coordinates in au, then B matrix
    for (i=0;i<dim_carts;++i)
      new_x[i] = (x[i] + dx[i]) / _bohr2angstroms;
    simples.compute_internals(nallatom,new_x);
    simples.compute_s(nallatom,new_x);
    free_block(B);
    B = compute_B(simples,symm);

    // compute new internal coordinate values
    free(new_q);
    new_q = compute_q(simples,symm);

    for (i=0;i<symm.get_num();++i)
      dq[i] = q[i] - new_q[i];

    // Test for convergence of iterations
    dx_sum = dq_sum = 0.0;
    for (i=0;i<dim_carts;++i)
      dx_sum += dx[i]*dx[i];
    dx_sum = sqrt(dx_sum / ((double) dim_carts));

    for (i=0;i<symm.get_num();++i)
      dq_sum += dq[i]*dq[i];
    dq_sum = sqrt(dq_sum / ((double) symm.get_num()));

    if ((dx_sum < optinfo.bt_dx_conv) && (dq_sum < optinfo.bt_dq_conv))
      bmat_iter_done = 1;
    fprintf (outfile,"%5d %15.12lf %15.12lf\n", count+1, dx_sum, dq_sum);

    // store new x in angstroms 
    for (i=0;i<dim_carts;++i)
      x[i] = new_x[i] * _bohr2angstroms;

    ++count;
  } while( (bmat_iter_done == 0) && (count < optinfo.bt_max_iter) );

  free_block(B);

  if (count >= optinfo.bt_max_iter) {
    fprintf(outfile,"Could not converge new geometry in %d iterations.",count);
    exit(2);
  }
  else
    fprintf(outfile,
        "Convergence to displaced geometry took %d iterations.\n",count);

  // take x back to bohr
  scalar_mult(1.0/_bohr2angstroms, x, dim_carts);
 
  // set coord(x)
  no_fx = new double [3*natom];
  for (i=0;i<natom;++i) {
    no_fx[3*i+0] = x[optinfo.to_dummy[i]*3+0];
    no_fx[3*i+1] = x[optinfo.to_dummy[i]*3+1];
    no_fx[3*i+2] = x[optinfo.to_dummy[i]*3+2];
  }

  // write geometry to output.dat and return it
  cartesians cart_temp;
  cart_temp.set_fcoord(x);
  cart_temp.set_coord(no_fx);
//  cart_temp.mult(1.0/_bohr2angstroms);
  fprintf(outfile,"\n%s\n",disp_label);
  cart_temp.print(1,outfile,0,disp_label,0);

  for (i=0;i<dim_carts;++i)
    return_geom[i] = x[i];
//    return_geom[i] = x[i]/_bohr2angstroms;

  // write geometry to chkpt or geom.dat
  if (print_flag == PRINT_TO_GEOM) {
    FILE *fp_geom;
    if (restart_geom_file) {
      ffile(&fp_geom, "geom.dat",0);
      fprintf(fp_geom, "geom.dat: (\n");
    }
    else
      ffile(&fp_geom, "geom.dat",1);

    cart_temp.print(10,fp_geom,restart_geom_file,disp_label,disp_num);
    if(last_disp) fprintf(fp_geom,")\n");
    fclose(fp_geom);
    fflush(outfile);
  }
  else if (print_flag == 32) {
    cart_temp.print(32,outfile,0,disp_label,disp_num);
  }

  delete [] no_fx;
  free(q); free(new_q);
  free(x); free(dx); free(new_x);
  free_block(A);
  free_block(G);
  free_block(G_inv);
  free_block(u);
  free_block(temp_mat);
  return;
}

