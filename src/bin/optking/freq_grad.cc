/** FREQ_GRAD computes frequencies from gradients */

extern "C" {
#include <stdio.h>
#include <libchkpt/chkpt.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <physconst.h>
#include <libpsio/psio.h>
#include <psifiles.h>
}

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"
#include "salc.h"
#include "bond_lengths.h"

extern double **compute_B(internals &simples, salc_set &salcs);
extern double *compute_q(internals &simples, salc_set &symm);
extern double **compute_G(double **B, int num_intcos, cartesians &carts);

void freq_grad(cartesians &carts, internals &simples, salc_set &salcs) {

  int i,j,k,a,b, cnt, dim, dim_carts, total_num_disps;
  double **B, **G, **G_inv, *masses, **u, *geom, *forces, **force_constants;
  double energy, *energies, **displacements, cm_convert;
  double *f, **all_f_q, *f_q, *temp_arr, *temp_arr2, **all_q, tval, **geom2D;
  double **evects, *evals, **FG;
  char *disp_label;

  disp_label = new char[MAX_LINELENGTH];
  dim_carts = 3*carts.get_nallatom();

  double *micro_e, *micro_geom, *micro_grad, *grad;

  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Total num. of disp.",
      (char *) &(total_num_disps), sizeof(int));

  // needed?
  micro_e = new double[total_num_disps];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(micro_e[0]), total_num_disps*sizeof(double));

  micro_grad = new double [total_num_disps*3*carts.get_natom()];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(micro_grad[0]), total_num_disps*3*carts.get_natom()*sizeof(double));

  micro_geom = new double [total_num_disps*3*carts.get_natom()];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced geometries",
      (char *) &(micro_geom[0]), total_num_disps*3*carts.get_nallatom()*sizeof(double));

  close_PSIF();

  all_q = (double **) malloc(total_num_disps*sizeof(double *));
  all_f_q = (double **) malloc(total_num_disps*sizeof(double *));
  f = (double *) malloc(dim_carts*sizeof(double));

  // compute forces in internal coordinates, f_q = G_inv B u f
  for (i=0; i<total_num_disps; ++i) {
    simples.compute_internals(carts.get_nallatom(), &(micro_geom[i*dim_carts]));
    simples.compute_s(carts.get_nallatom(), &(micro_geom[i*dim_carts]));
    all_q[i] = compute_q(simples, salcs);
    B = compute_B(simples,salcs);
    G = compute_G(B,salcs.get_num(),carts);
    fprintf(outfile,"BuB^t ");
    G_inv = symm_matrix_invert(G,salcs.get_num(),1,optinfo.redundant);

    masses = carts.get_fmass();
    u = mass_mat(masses);

    // load up forces
    for (j=0;j<dim_carts;++j)
      f[j] = micro_grad[i*dim_carts+j] * -1.0 * _hartree2J * 1.0E18 / _bohr2angstroms;

    cnt = -1;
    /*
    fprintf(outfile,"gradient %d\n",i);
      for (j=0;j<carts.get_nallatom();++j) {
        for (k=0;k<3;++k)
          fprintf(outfile,"%15.10lf", f[++cnt]);
        fprintf(outfile,"\n");
      }
      */

    all_f_q[i] = init_array(salcs.get_num());
    temp_arr = init_array(salcs.get_num());
    temp_arr2 = init_array(dim_carts);

    mmult(u,0,&f,1,&temp_arr2,1,dim_carts,dim_carts,1,0);
    mmult(B,0,&temp_arr2,1,&temp_arr,1,salcs.get_num(),dim_carts,1,0);
    mmult(G_inv,0,&temp_arr,1,&(all_f_q[i]),1,salcs.get_num(),salcs.get_num(),1,0);

    free(temp_arr);
    free(temp_arr2);
    free_block(u);
    free_block(B);
    free_block(G);
    free_block(G_inv);
  }

  free(f);

  /*
  for (i=0;i<total_num_disps;++i) {
    fprintf(outfile,"Values of internal coordinates, displacement %d\n",i);
    for (j=0; j<salcs.get_num();++j)
      fprintf(outfile,"%15.10lf",all_q[i][j]);
    fprintf(outfile,"\n");
  }
  for (i=0;i<total_num_disps;++i) {
    fprintf(outfile,"Values of internal coordinate forces, displacement %d\n",i);
    for (j=0; j<salcs.get_num();++j)
      fprintf(outfile,"%15.10lf",all_f_q[i][j]);
    fprintf(outfile,"\n");
  }
  */

  // apply three point formula
  force_constants = block_matrix(salcs.get_num(),salcs.get_num());
  for (i=0;i<salcs.get_num();++i)
    for (j=0;j<salcs.get_num();++j)
      force_constants[i][j] = force_constants[j][i] =
        (all_f_q[2*i][j]-all_f_q[2*i+1][j]) / (2.0 * optinfo.disp_size);

  //print_mat2(all_f_q, num_disps, salcs.get_num(), outfile);
  free_block(all_f_q);

  fprintf(outfile,"\nForce Constants\n");
  print_mat(force_constants, salcs.get_num(), salcs.get_num(), outfile);
  fflush(outfile);

  // build G = BuB^t
  B = block_matrix(salcs.get_num(), 3*carts.get_natom());
  B = compute_B(simples, salcs);
  // G = block_matrix(salcs.get_num(), salcs.get_num());
  G = compute_G(B, salcs.get_num(), carts);
  free_block(B);

  // compute FG and diagonalize 
  FG = block_matrix(salcs.get_num(), salcs.get_num());
  mmult(force_constants,0,G,0,FG,0,
      salcs.get_num(),salcs.get_num(),salcs.get_num(),0);
  free_block(force_constants);
  free_block(G);

  //fprintf(outfile,"FG Matrix\n");
  //print_mat2(FG,salcs.get_num(),salcs.get_num(),outfile);
  //fflush(outfile);

  evals  = init_array(salcs.get_num());
  G = block_matrix(salcs.get_num(), salcs.get_num());
  dgeev_optking(salcs.get_num(), FG, evals, G);
  free_block(FG);
  free_block(G);

  cm_convert = 1.0/(2.0 * _pi * _c * 100.0);
  for (i=0;i<salcs.get_num();++i) {
    evals[i] = evals[i] * 1.0E-18 / ( 1.0E-20 * _amu2kg );
    evals[i] = cm_convert * sqrt( evals[i] );
  }

  fprintf(outfile,"\nHarmonic Vibrational Frequencies\n");
  for (i=salcs.get_num()-1;i>=0;--i)
    fprintf(outfile,"\t %d       %15.1lf\n",salcs.get_num()-i,evals[i]);

  free(evals);

  open_PSIF();
  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));
  close_PSIF();
}

