/** displaces along all coordinates + and - by disp_size **/ 
// need to save geometry to last step for gradients by energy

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


extern void new_geom(cartesians &carts, internals &simples, salc_set &symm, 
   double *dq, int print_to_geom_file, int restart_geom_file,
   char *disp_label, int disp_num, int last_disp, double *return_geom);


int disp_all(cartesians &carts, internals &simples, salc_set &salcs, 
              int points) 
{
  int i,j,a,b, dim, dim_carts, num_disps;
  double **B, **G, **G_inv, *masses, **u, *fgeom, *forces, **force_constants;
  double energy, *energies, **micro_geoms, **displacements, cm_convert;
  double *f, **all_f_q, *f_q, *temp_arr, *temp_arr2, *dq, *q, tval, **fgeom2D;
  double **evects, *evals, **FG;
  char *disp_label;

  disp_label = new char[MAX_LINELENGTH];
  dim_carts = 3*carts.get_nallatom();

  if (salcs.get_num() == 0) 
    punt("There are no SALCs present to displace along.\n");

  // write reference geometry (the one being replaced to PSIF)
  fgeom = carts.get_fcoord();
  energy = carts.get_energy();
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Reference geometry",
      (char *) &(fgeom[0]), dim_carts*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Reference energy",
      (char *) &(energy), sizeof(double));
  close_PSIF();

  num_disps = 2 * salcs.get_num();
  /*** make list of internal displacements for micro_iterations ***/
  displacements = block_matrix(num_disps,salcs.get_num());
  for (i=0;i<salcs.get_num();++i) {
    displacements[2*i][i] = -1.0 * optinfo.disp_size;
    displacements[2*i+1][i] = 1.0 * optinfo.disp_size;
  }
  fprintf(outfile,"\nDisplacement Matrix\n");
  print_mat2(displacements, num_disps, salcs.get_num(), outfile);

  /*** generate and store Micro_iteration cartesian geometries ***/
  micro_geoms = block_matrix(num_disps, dim_carts);
  for (i=0;i<num_disps;++i)  {
    sprintf(disp_label,"Displaced geometry %d in a.u.\n",i+1);
    new_geom(carts,simples,salcs,displacements[i],0,
        0, disp_label, i, 0, micro_geoms[i]);
  }
  free_block(displacements);

  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced geometries",
      (char *) &(micro_geoms[0][0]), num_disps*dim_carts*sizeof(double));

  optinfo.disp_num = 0;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
    (char *) &(optinfo.disp_num),sizeof(int));

  psio_write_entry(PSIF_OPTKING, "OPT: Total num. of disp.",
    (char *) &(num_disps), sizeof(int));

  close_PSIF();
  free_block(micro_geoms);

  // make room for storage of energy and gradients of displacements
  double *disp_e, *disp_grad;

  disp_e = new double[num_disps];
  disp_grad = new double [num_disps*3*carts.get_natom()*sizeof(double)];

  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced gradients",
    (char *) &(disp_grad[0]), num_disps*3*carts.get_natom()*sizeof(double));

  psio_write_entry(PSIF_OPTKING, "OPT: Displaced energies",
   (char *) &(disp_e[0]), num_disps*sizeof(double));

 close_PSIF();


 // Reset microiteration counter
 optinfo.micro_iteration = 0;
 open_PSIF();
 psio_write_entry(PSIF_OPTKING, "Micro_iteration",
   (char *) &(optinfo.micro_iteration),sizeof(int));
 close_PSIF();

 delete [] disp_e;
 delete [] disp_grad;
 return(num_disps);
}

