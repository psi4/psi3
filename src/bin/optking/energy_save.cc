/*** ENERGY_SAVE.CC Rollin King, 2002 ***/
// function executes if optinfo.mode == MODE_GRAD_SAVE

#include <cmath>
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

void energy_save(cartesians &carts) {
  int i,j,dim_carts,total_num_disps;
  double energy, *micro_e, *micro_grad, *grad;

  dim_carts = 3*optinfo.natom;
  fprintf(outfile,"Saving energy.\n");

  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num), sizeof(int));
  psio_read_entry(PSIF_OPTKING, "OPT: Total num. of disp.",
      (char *) &(total_num_disps), sizeof(int));

  micro_e = new double[total_num_disps];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(micro_e[0]), total_num_disps*sizeof(double));
  close_PSIF();

  chkpt_init(PSIO_OPEN_OLD);
  energy = chkpt_rd_etot();
  chkpt_close();

  micro_e[optinfo.disp_num] = energy;

  open_PSIF();
  psio_write_entry(PSIF_OPTKING,"OPT: Displaced energies",
      (char *) &(micro_e[0]), total_num_disps*sizeof(double));
  close_PSIF();
  delete [] micro_e;

  // increment disp_num
  open_PSIF();
  optinfo.disp_num += 1;
  psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
      (char *) &(optinfo.disp_num), sizeof(int));
  close_PSIF();
  return ;
}
