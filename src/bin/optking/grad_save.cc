/*** ENERGY_SAVE.CC Rollin King, 2002 ***/
// function executes if optinfo.mode == MODE_ENERGY_SAVE

#if HAVE_CMATH
# include <cmath>
#else
# include <math.h>
#endif

extern "C" {
#include <stdio.h>
#include <libchkpt/chkpt.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <physconst.h>
#include <psifiles.h>
}

#define EXTERN
#include "opt.h"
#undef EXTERN
#include "cartesians.h"
#include "internals.h"

void grad_save(cartesians &carts) {
  int i,j,dim_carts,total_num_disps;
  double energy, *micro_e, *micro_grad, *grad, **ggrad, **refgrad, **rref;

  dim_carts = 3*optinfo.natom;

  fprintf(outfile,"Saving gradient and energy.\n");

  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num), sizeof(int));
  psio_read_entry(PSIF_OPTKING, "OPT: Total num. of disp.",
      (char *) &(total_num_disps), sizeof(int));

  micro_e = new double[total_num_disps];
  micro_grad = new double [total_num_disps*3*carts.get_natom()*sizeof(double)];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(micro_grad[0]), total_num_disps*3*carts.get_natom()*sizeof(double));

  /*
  fprintf(outfile,"gradients\n");
  for (i=0; i<total_num_disps*3*carts.get_natom(); ++i)
    fprintf(outfile,"%15.10lf\n",micro_grad[i]);
    */

  psio_read_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(micro_e[0]), total_num_disps*sizeof(double));
  close_PSIF();

  chkpt_init(PSIO_OPEN_OLD);
  grad = chkpt_rd_grad();
  rref = chkpt_rd_rref();
  energy = chkpt_rd_etot();
  chkpt_close();

  // Rotate the gradient back to the reference frame in which all geometries were generated and stores
  int natoms = carts.get_natom();
  ggrad = block_matrix(natoms,3);
  int atomxyz=0;
  for(int atom=0; atom<natoms; atom++)
    for(int xyz=0; xyz<3; xyz++,atomxyz++)
      ggrad[atom][xyz] = grad[atomxyz];
  delete[] grad;
  refgrad = block_matrix(carts.get_natom(),3);
  mmult(ggrad,0,rref,0,refgrad,0,carts.get_natom(),3,3,0);
  free_block(rref);
  free_block(ggrad);

  micro_e[optinfo.disp_num] = energy;
  for (i=0; i<dim_carts; ++i) {
    int atom = i/3;
    int xyz = i%3;
    micro_grad[3*carts.get_natom()*(optinfo.disp_num)+i] = refgrad[atom][xyz];
  }

  open_PSIF();
  psio_write_entry(PSIF_OPTKING,"OPT: Displaced energies",
      (char *) &(micro_e[0]), total_num_disps*sizeof(double));
  psio_write_entry(PSIF_OPTKING, "OPT: Displaced gradients",
      (char *) &(micro_grad[0]), total_num_disps*3*carts.get_natom()*sizeof(double));
  close_PSIF();

  delete [] micro_e; delete [] micro_grad;

  // increment disp_num
  open_PSIF();
  optinfo.disp_num += 1;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num), sizeof(int));
  close_PSIF();

  return ;
}

