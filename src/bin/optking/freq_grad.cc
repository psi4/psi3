/** FREQ_GRAD computes frequencies from gradients */

#include <cmath>
extern "C" {
  #include <stdio.h>
  #include <libchkpt/chkpt.h>
  #include <stdlib.h>
  #include <string.h>
  #include <ctype.h>
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

/* FREQ_GRAD_IRREP compute frequencies from gradients for irrep block IRREP */

void freq_grad_irrep(cartesians &carts, internals &simples, salc_set &all_salcs,
    int points) {

  int i,j,ii,jj,k,a,b, cnt, dim, dim_carts, ndisps,irrep;
  int nirr_salcs, nsalcs, *irrep_salcs;
  double **B, **G, **G_inv, *masses, **u, *geom, *forces, **force_constants;
  double energy, *energies, **displacements, cm_convert;
  double *f, *f_q, *temp_arr, *temp_arr2, *q, tval, **geom2D;
  double **all_f_q; // internal coordinate forces for all unique displacements
  double **full_all_f_q; // internal coordinate forces for all displacements
  double **evects, *evals, **FG;
  double *micro_e, *micro_geom, *micro_grad, *grad, tmp;
  char *salc_lbl;

  dim_carts = 3*carts.get_natom();
  irrep_salcs = new int[all_salcs.get_num()];
  irrep = optinfo.irrep;
  nsalcs = all_salcs.get_num();

  /* count and identify IRREP salcs */
  nirr_salcs = 0;
  cnt = 0;
  for (i=0; i<all_salcs.get_num(); ++i) {
    salc_lbl = all_salcs.get_label(i);
    if ( strcmp(salc_lbl, syminfo.irrep_lbls[irrep]) == 0) {
      ++nirr_salcs;
      irrep_salcs[cnt++] = i;
    }
  }
  if (nirr_salcs == 0) {
    fprintf(outfile,"No coordinates of irrep %d\n.", irrep);
  }

fprintf(outfile,"Found %d salcs of this irrep\n",nirr_salcs);

  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Total num. of disp.",
      (char *) &(ndisps), sizeof(int));

  // needed?
  micro_e = new double[ndisps];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(micro_e[0]), ndisps*sizeof(double));

  micro_grad = new double [ndisps*3*carts.get_natom()];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(micro_grad[0]), ndisps*3*carts.get_natom()*sizeof(double));

  micro_geom = new double [ndisps*3*carts.get_natom()];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced geometries",
      (char *) &(micro_geom[0]), ndisps*3*carts.get_natom()*sizeof(double));

  close_PSIF();

  // compute forces in internal coordinates for all disps, f_q = G_inv B u f
  all_f_q = block_matrix(ndisps, nsalcs);
  f = init_array(dim_carts);
  temp_arr = init_array(nsalcs);
  temp_arr2 = init_array(dim_carts);
  masses = carts.get_fmass();
  u = mass_mat(masses);
  for (i=0; i<ndisps; ++i) {

    simples.compute_internals(carts.get_natom(), &(micro_geom[i*dim_carts]));
    simples.compute_s(carts.get_natom(), &(micro_geom[i*dim_carts]));
    q = compute_q(simples, all_salcs);

    /* fprintf(outfile,"Values of internal coordinates, displacement %d\n",i);
     for (j=0; j<salcs.get_num();++j) fprintf(outfile,"%15.10lf",all_q[i][j]);
     fprintf(outfile,"\n"); */

    B = compute_B(simples,all_salcs);
    G = compute_G(B, nsalcs, carts);
    fprintf(outfile,"BuB^t ");
    G_inv = symm_matrix_invert(G, nsalcs, 1, optinfo.redundant);

    // load up cartesian forces
    for (j=0;j<dim_carts;++j)
      f[j] = micro_grad[i*dim_carts+j] * -1.0 * _hartree2J * 1.0E18 /
        _bohr2angstroms;

    mmult(u,0,&f,1,&temp_arr2,1,dim_carts,dim_carts,1,0);
    mmult(B,0,&temp_arr2,1,&temp_arr,1, nsalcs,dim_carts,1,0);
    mmult(G_inv,0,&temp_arr,1,&(all_f_q[i]),1, nsalcs, nsalcs, 1, 0);

    free(q);
    free_block(B);
    free_block(G);
    free_block(G_inv);

  }
  free(f);
  free(temp_arr);
  free(temp_arr2);
  free(masses);
  free_block(u);

  /* expand unique displacements to redundant displacements */
  full_all_f_q = block_matrix( 2*all_salcs.get_num(), all_salcs.get_num());
  if (irrep == 0) {
    for (i=0; i<ndisps; ++i) // loop over displacements
      for (j=0; j<nsalcs; ++j)
        full_all_f_q[i][j] = all_f_q[i][j];
  }
  else {
    for (i=0; i<ndisps; ++i) { // loop over displacements
      for (j=0; j<nsalcs; ++j) {
        full_all_f_q[2*irrep_salcs[i]][j] = all_f_q[i][j]; // the - disp
        salc_lbl = all_salcs.get_label(j);
        if ( strcmp(salc_lbl, syminfo.irrep_lbls[irrep]) == 0)
          full_all_f_q[2*irrep_salcs[i]+1][j] = -1 * all_f_q[i][j];
        else
          full_all_f_q[2*irrep_salcs[i]+1][j] = all_f_q[i][j];
      }
    }
  }
  free_block(all_f_q);

  for (i=0;i<2*nsalcs;++i) {
    fprintf(outfile,
        "Redundant values of internal coordinate forces, displacement %d\n",i);
    for (j=0; j<nsalcs;++j)
      fprintf(outfile,"%15.10lf",full_all_f_q[i][j]);
    fprintf(outfile,"\n");
  }

  // apply three point formula - to generate force constants in this irrep block
  fprintf(outfile,"Applying %d-point formula\n",points);
  force_constants = block_matrix(nsalcs,nsalcs);
  for (i=0;i<nirr_salcs;++i)
    for (j=0;j<nirr_salcs;++j) {
      ii = irrep_salcs[i];
      jj = irrep_salcs[j];
      force_constants[ii][jj] = force_constants[jj][ii] =
        (full_all_f_q[2*ii][jj]-full_all_f_q[2*ii+1][jj]) / (2.0 * optinfo.disp_size);
    }

  //print_mat2(full_all_f_q, num_disps, salcs.get_num(), outfile);
  free_block(full_all_f_q);

  fprintf(outfile,"\nForce Constants\n");
  print_mat(force_constants, nsalcs, nsalcs, outfile);
  fflush(outfile);

  // build G = BuB^t
  B = compute_B(simples, all_salcs);
  G = compute_G(B, nsalcs, carts);
  free_block(B);

  // compute FG and diagonalize 
  FG = block_matrix(nsalcs, nsalcs);
  mmult(force_constants,0,G,0,FG,0,
      nsalcs,nsalcs,nsalcs,0);
  free_block(force_constants);
  free_block(G);

  //fprintf(outfile,"FG Matrix\n");
  //print_mat2(FG,salcs.get_num(),salcs.get_num(),outfile);
  //fflush(outfile);

  evals  = init_array(nsalcs);
  G = block_matrix(nsalcs, nsalcs);
  dgeev_optking(nsalcs, FG, evals, G);
  free_block(FG);
  free_block(G);

  cm_convert = 1.0/(2.0 * _pi * _c * 100.0);
  for (i=0;i<nsalcs;++i) {
    evals[i] = evals[i] * 1.0E-18 / ( 1.0E-20 * _amu2kg );
    evals[i] = cm_convert * sqrt( evals[i] );
  }

  fprintf(outfile,"\nHarmonic Vibrational Frequencies in cm^(-1) for Irrep %s\n",
      syminfo.irrep_lbls[irrep]) ;
  for (i=0; i<nirr_salcs; ++i) {
    tmp = -9999;
    for (j=0; j<nsalcs; ++j) {
      if (evals[j] > tmp) {
        tmp = evals[j];
        ii = j;
      }
    }
    fprintf(outfile,"%5d       %15.1lf\n",nirr_salcs-i,evals[ii]);
    evals[ii] = -9999;
  }
  free(evals);

  open_PSIF();
  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));
  close_PSIF();

  delete [] irrep_salcs;
}


/** FREQ_GRAD_NOSYMM compute frequencies from gradients for all
* coordinates ignoring symmetry */

void freq_grad_nosymm(cartesians &carts, internals &simples,
    salc_set &all_salcs, int points) {

  int i,j,k,a,b, ii, cnt, dim, dim_carts, ndisps;
  int nsalcs;
  double **B, **G, **G_inv, *masses, **u, *geom, *forces, **force_constants;
  double energy, *energies, **displacements, cm_convert;
  double *f, **all_f_q, *f_q, *temp_arr, *temp_arr2, **all_q, tval, **geom2D;
  double **evects, *evals, **FG, tmp;
  double *micro_e, *micro_geom, *micro_grad, *grad;

  dim_carts = 3*carts.get_natom();
  nsalcs = all_salcs.get_num();

  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Total num. of disp.",
      (char *) &(ndisps), sizeof(int));
  if (ndisps != 2*nsalcs) 
    punt("Error: number of displacements is incorrect.");

  // needed?
  micro_e = new double[ndisps];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(micro_e[0]), ndisps*sizeof(double));

  micro_grad = new double [ndisps*3*carts.get_natom()];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(micro_grad[0]), ndisps*3*carts.get_natom()*sizeof(double));

  micro_geom = new double [ndisps*3*carts.get_natom()];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced geometries",
      (char *) &(micro_geom[0]), ndisps*3*carts.get_natom()*sizeof(double));

  close_PSIF();

  all_q = (double **) malloc(ndisps*sizeof(double *));
  all_f_q = (double **) malloc(ndisps*sizeof(double *));
  f = (double *) malloc(dim_carts*sizeof(double));

  // compute forces in internal coordinates, f_q = G_inv B u f
  for (i=0; i<ndisps; ++i) {
    simples.compute_internals(carts.get_natom(), &(micro_geom[i*dim_carts]));
    simples.compute_s(carts.get_natom(), &(micro_geom[i*dim_carts]));

    all_q[i] = compute_q(simples, all_salcs);
    B = compute_B(simples,all_salcs);
    G = compute_G(B,nsalcs,carts);
    // fprintf(outfile,"BuB^t ");
    G_inv = symm_matrix_invert(G,nsalcs,0,optinfo.redundant);
    masses = carts.get_fmass();
    u = mass_mat(masses);

    // load up forces
    for (j=0;j<dim_carts;++j)
      f[j] = micro_grad[i*dim_carts+j] * -1.0 * _hartree2J * 1.0E18 /
        _bohr2angstroms;

    all_f_q[i] = init_array(nsalcs);
    temp_arr = init_array(nsalcs);
    temp_arr2 = init_array(dim_carts);

    mmult(u,0,&f,1,&temp_arr2,1,dim_carts,dim_carts,1,0);
    mmult(B,0,&temp_arr2,1,&temp_arr,1,nsalcs,dim_carts,1,0);
    mmult(G_inv,0,&temp_arr,1,&(all_f_q[i]),1,nsalcs,nsalcs,1,0);

    free(temp_arr);
    free(temp_arr2);
    free_block(u);
    free_block(B);
    free_block(G);
    free_block(G_inv);
  }

  free(f);

  /*
  fprintf(outfile,"Values of intcos for each displacement\n");
  for (i=0;i<ndisps;++i) {
    for (j=0; j<nsalcs;++j)
      fprintf(outfile,"%15.10lf",all_q[i][j]);
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"Values of intcos forces for each displacement\n");
  for (i=0;i<ndisps;++i) {
    for (j=0; j<nsalcs;++j)
      fprintf(outfile,"%15.10lf",all_f_q[i][j]);
    fprintf(outfile,"\n");
  }
  */

  // apply three point formula
  fprintf(outfile,"Applying %d-point formula\n",points);
  force_constants = block_matrix(nsalcs,nsalcs);
  for (i=0;i<nsalcs;++i)
    for (j=0;j<nsalcs;++j)
      force_constants[i][j] = force_constants[j][i] =
        (all_f_q[2*i][j]-all_f_q[2*i+1][j]) / (2.0 * optinfo.disp_size);

  //print_mat2(all_f_q, num_disps, salcs.get_num(), outfile);
  free_block(all_f_q);

  fprintf(outfile,"\nForce Constants\n");
  print_mat(force_constants, nsalcs, nsalcs, outfile);
  fflush(outfile);

  // build G = BuB^t
  B = block_matrix(nsalcs, 3*carts.get_natom());
  B = compute_B(simples, all_salcs);
  G = compute_G(B, nsalcs, carts);
  free_block(B);

  // compute FG and diagonalize 
  FG = block_matrix(nsalcs, nsalcs);
  mmult(force_constants,0,G,0,FG,0,nsalcs,nsalcs,nsalcs,0);
  free_block(force_constants);
  free_block(G);

  //fprintf(outfile,"FG Matrix\n");
  //print_mat2(FG,salcs.get_num(),salcs.get_num(),outfile);
  //fflush(outfile);

  evals  = init_array(nsalcs);
  G = block_matrix(nsalcs, nsalcs);
  dgeev_optking(nsalcs, FG, evals, G);
  free_block(FG);
  free_block(G);

  cm_convert = 1.0/(2.0 * _pi * _c * 100.0);
  for (i=0;i<nsalcs;++i) {
    evals[i] = evals[i] * 1.0E-18 / ( 1.0E-20 * _amu2kg );
    evals[i] = cm_convert * sqrt( evals[i] );
  }

  fprintf(outfile,"\nHarmonic Vibrational Frequencies\n");
  for (i=0; i<nsalcs; ++i) {
    tmp = -9999;
    for (j=0; j<nsalcs; ++j) {
      if (evals[j] > tmp) {
        tmp = evals[j];
        ii = j;
      }
    }
    fprintf(outfile,"%5d       %15.1lf\n",nsalcs-i,evals[ii]);
    evals[ii] = -9999;
  }
  free(evals);

  open_PSIF();
  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));
  close_PSIF();

}

