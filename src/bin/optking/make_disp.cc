/** displaces along all coordinates + and - by disp_size **/ 
// need to save geometry to last step for gradients by energy

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


extern void new_geom(cartesians &carts, internals &simples, salc_set &all_salcs, 
    double *dq, int print_to_geom_file, int restart_geom_file,
    char *disp_label, int disp_num, int last_disp, double *return_geom);

/* MAKE_DISP_IRREP - make displacements for modes labelled with 
* symmetry IRREP (+ and - if symmetric; otherwise, just +) */

int make_disp_irrep(cartesians &carts, internals &simples, salc_set &all_salcs, 
    int points) 
{
  int i,j,a,b, cnt,dim, dim_carts, ndisps, nsalcs, *irrep_salcs, irrep;
  double **B, **G, **G_inv, *masses, **u, *fgeom, *forces, **force_constants;
  double energy, *energies, **micro_geoms, **displacements, cm_convert;
  double *f, **all_f_q, *f_q, *temp_arr, *temp_arr2, *dq, *q, tval, **fgeom2D;
  double **evects, *evals, **FG;
  char *disp_label, **disp_irrep_lbls, *salc_lbl;

  disp_label = new char[MAX_LINELENGTH];
  dim_carts = 3*carts.get_nallatom();
  irrep_salcs = new int[all_salcs.get_num()];
  irrep = optinfo.irrep;

  /* count number of IRREP salcs */
  nsalcs = 0;
  cnt = 0;
  for (i=0; i<all_salcs.get_num(); ++i) {
    salc_lbl = all_salcs.get_label(i);
    if ( strcmp(salc_lbl, syminfo.irrep_lbls[irrep]) == 0) {
      ++nsalcs;
      irrep_salcs[cnt++] = i;
    }
  }
  fprintf(outfile,"Found %d internals of irrep %d \n", nsalcs, irrep+1);
//  for (i=0; i<nsalcs; ++i)
//    fprintf(outfile,"internal id %d\n",irrep_salcs[i]);

  if (nsalcs == 0) { 
    fprintf(outfile,"Produced 0 displacements.\n");
    return(0);
  }

  // write reference geometry (the one being replaced to PSIF)
  fgeom = carts.get_fcoord();
  energy = carts.get_energy();
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Reference geometry",
      (char *) &(fgeom[0]), dim_carts*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Reference energy",
      (char *) &(energy), sizeof(double));
  close_PSIF();

  /*** make list of internal displacements for micro_iterations ***/
  if (irrep == 0) {
    ndisps = 2*nsalcs;
    displacements = block_matrix(ndisps,all_salcs.get_num());
    cnt = 0;
    for (i=0; i<nsalcs; ++i) {
      displacements[cnt++][irrep_salcs[i]] = -1.0 * optinfo.disp_size;
      displacements[cnt++][irrep_salcs[i]] =  1.0 * optinfo.disp_size;
    }
  }
  else { // non-symmetric irrep
    ndisps = nsalcs;
    displacements = block_matrix(ndisps,all_salcs.get_num());
    cnt = 0;
    for (i=0; i<nsalcs; ++i)
      displacements[cnt++][irrep_salcs[i]] = -1.0 * optinfo.disp_size;
  }

  // count number of unique displacements; 2 disps for symm modes
  fprintf(outfile,"Generating a total of %d displacements ", ndisps);
  fprintf(outfile,"using %d-point formula for modes of irrep %d.\n", points, irrep+1);

  fprintf(outfile,"\nDisplacement Matrix\n");
  print_mat5(displacements, ndisps, all_salcs.get_num(), outfile);

  /*** generate and store Micro_iteration cartesian geometries ***/
  micro_geoms = block_matrix(ndisps, dim_carts);
  for (i=0;i<ndisps;++i)  {
    sprintf(disp_label,"Displaced geometry %d in a.u.\n",i+1);
    new_geom(carts,simples,all_salcs,displacements[i],0,
        0, disp_label, i, 0, micro_geoms[i]);
  }
  free_block(displacements);

  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced geometries",
      (char *) &(micro_geoms[0][0]), ndisps*dim_carts*sizeof(double));

  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));

  psio_write_entry(PSIF_OPTKING, "OPT: Total num. of disp.",
      (char *) &(ndisps), sizeof(int));

  close_PSIF();
  free_block(micro_geoms);

  // make room for storage of energy and gradients of displacements
  double *disp_e, *disp_grad;

  disp_e = new double[ndisps];
  disp_grad = new double [ndisps*3*carts.get_natom()*sizeof(double)];

  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(disp_grad[0]), ndisps*3*carts.get_natom()*sizeof(double));

  psio_write_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(disp_e[0]), ndisps*sizeof(double));

  close_PSIF();


  // Reset microiteration counter
  optinfo.micro_iteration = 0;
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Micro_iteration",
      (char *) &(optinfo.micro_iteration),sizeof(int));
  close_PSIF();

  fprintf(outfile,"Produced a total of %d displacements.\n",ndisps);

  delete [] disp_label;
  delete [] irrep_salcs;
  delete [] disp_e;
  delete [] disp_grad;
  return(ndisps);
}


/* MAKE_DISP_NOSYMM generate displacements - do positive and
* negative displacements along all coordinates ignorning symmetry */

int make_disp_nosymm(cartesians &carts, internals &simples,
    salc_set &all_salcs, int points) 
{
  int i,j,a,b, dim, dim_carts, ndisps, nsalcs;
  double **B, **G, **G_inv, *masses, **u, *fgeom, *forces, **force_constants;
  double energy, *energies, **micro_geoms, **displacements, cm_convert;
  double *f, **all_f_q, *f_q, *temp_arr, *temp_arr2, *dq, *q, tval, **fgeom2D;
  double **evects, *evals, **FG;
  char *disp_label, **disp_irrep_lbls;

  disp_label = new char[MAX_LINELENGTH];
  dim_carts = 3*carts.get_nallatom();

  nsalcs = all_salcs.get_num();
  /* only 3-pt formula for now */
  ndisps = 2*nsalcs;

  if (nsalcs == 0)
    punt("There are no appropriate SALCs present to displace along.\n");

  // write reference geometry (the one being replaced to PSIF)
  fgeom = carts.get_fcoord();
  energy = carts.get_energy();
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Reference geometry",
      (char *) &(fgeom[0]), dim_carts*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Reference energy",
      (char *) &(energy), sizeof(double));
  close_PSIF();

  fprintf(outfile,"Generating a total of %d displacements ", ndisps);
  fprintf(outfile,"using %d-point formula for all modes.\n", points);

  /*** make list of internal displacements for micro_iterations ***/
  displacements = block_matrix(ndisps, nsalcs);
  for (i=0;i<all_salcs.get_num();++i) {
    displacements[2*i][i] = -1.0 * optinfo.disp_size;
    displacements[2*i+1][i] = 1.0 * optinfo.disp_size;
  }

  /*
  fprintf(outfile,"\nDisplacement Matrix\n");
  print_mat2(displacements, ndisps, nsalcs, outfile);
  */

  /*** generate and store Micro_iteration cartesian geometries ***/
  micro_geoms = block_matrix(ndisps, dim_carts);
  for (i=0;i<ndisps;++i)  {
    sprintf(disp_label,"Displaced geometry %d in a.u.\n",i+1);
    new_geom(carts,simples,all_salcs,displacements[i],0,
        0, disp_label, i, 0, micro_geoms[i]);
  }
  free_block(displacements);

  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced geometries",
      (char *) &(micro_geoms[0][0]), ndisps*dim_carts*sizeof(double));

  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num),sizeof(int));

  psio_write_entry(PSIF_OPTKING, "OPT: Total num. of disp.",
      (char *) &(ndisps), sizeof(int));

  close_PSIF();
  free_block(micro_geoms);

  // write zeroes for initial energy and gradients of displacements
  double *disp_e, *disp_grad;

  disp_e = new double[ndisps];
  disp_grad = new double [ndisps*3*carts.get_natom()*sizeof(double)];

  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(disp_grad[0]), ndisps*3*carts.get_natom()*sizeof(double));

  psio_write_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(disp_e[0]), ndisps*sizeof(double));

  close_PSIF();


  // Reset microiteration counter
  optinfo.micro_iteration = 0;
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Micro_iteration",
      (char *) &(optinfo.micro_iteration),sizeof(int));
  close_PSIF();

  delete [] disp_e;
  delete [] disp_grad;
  return(ndisps);
}

