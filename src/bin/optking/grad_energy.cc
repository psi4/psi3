/*** GRAD_ENERGY computes a file11 entry from energies in chkpt Rollin King, 2002 ***/ 

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

double **compute_B(internals &simples, salc_set &symm);

void grad_energy(cartesians &carts, internals &simples, salc_set &symm) {

  int i,j,a,b, dim, dim_carts, num_disps, cnt;
  double **B, *geom, *forces;
  double energy, *energies, **micro_geoms, **displacements;
  double *f, *f_q, *dq, *q, tval, **geom2D;
  char *disp_label;

  disp_label = new char[MAX_LINELENGTH];
  dim_carts = 3*carts.get_natom();

  if (symm.get_num() == 0) {
    punt("No symmetric internal coordinates to optimize.\n");
  }

  num_disps = 2 * symm.get_num();

  // report energies
  open_PSIF();
  energies = new double[num_disps];
  psio_read_entry(PSIF_OPTKING, "OPT: Displaced energies",
      (char *) &(energies[0]), num_disps*sizeof(double));
  close_PSIF();
  fprintf(outfile,"Displacement energies, Check for precision!\n");
  for(i=0;i<num_disps;++i) fprintf(outfile,"%15.10lf\n",energies[i]);
  fflush(outfile);

  // Calculate forces in internal coordinates
  f_q = new double[symm.get_num()];
  for (i=0;i<symm.get_num();++i) {
    f_q[i] = (energies[2*i+1]-energies[2*i]) / (2.0 * optinfo.disp_size);
    f_q[i] = -1.0 * f_q[i] * _hartree2J * 1.0E18 ;
  }
  free(energies);

  // Print internal coordinate forces
  // fprintf(outfile,"\nInternal coordinate forces\n");
  // for (i=0;i<symm.get_num();++i)
  //   fprintf(outfile,"%13.10lf\n",f_q[i]);

  // write out approximate file11.dat
  geom = new double[dim_carts];
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Reference geometry",
      (char *) &(geom[0]), dim_carts*sizeof(double));
  psio_read_entry(PSIF_OPTKING, "OPT: Reference energy",
      (char *) &(energy), sizeof(double));
  close_PSIF();

  // Transform forces to cartesian coordinates
  simples.compute_internals(dim_carts, geom);
  simples.compute_s(dim_carts, geom);
  B = compute_B(simples, symm);
  f = new double[dim_carts];
  mmult(B,1,&f_q,1,&f,1,dim_carts,symm.get_num(),1,0);
  free_block(B);

  // change forces to gradient for writing a file11 entry
  for(i=0;i<dim_carts;++i)
    f[i] = -1.0 * f[i] / _hartree2J / 1.0E18 * _bohr2angstroms;

  ffile(&fp_11, "file11.dat", 1);
  sprintf(disp_label,"iteration: %d", optinfo.iteration+1);
  carts.set_energy(energy);
  carts.set_coord(geom);
  carts.set_grad(f);
  carts.print(11,fp_11,0,disp_label, 0);
  fclose(fp_11);

  // write out geometry, gradient and energy to chkpt file
  cnt = -1;
  geom2D = block_matrix(carts.get_natom(),3);
  for (i=0; i<carts.get_natom(); ++i)
    for (j=0; j<3; ++j)
      geom2D[i][j] = geom[++cnt];

  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_geom(geom2D);
  chkpt_wt_grad(f);
  chkpt_wt_etot(energy);
  chkpt_close();
  free_block(geom2D);
  free(f);
  free(geom);

  // recompute values of internals and s vectors -- too late!
//  simples.compute_internals(carts.get_natom(),carts.get_coord());
//  simples.compute_s(carts.get_natom(),carts.get_coord() );

  // use optking --opt_step to take a step
  // opt_step(carts, simples, symm);

  // reset microiteration value in disp_all
  /*
  optinfo.micro_iteration = 0;
  open_PSIF();
  psio_write_entry(PSIF_OPTKING, "Micro_iteration",
      (char *) &(optinfo.micro_iteration),sizeof(int));
  close_PSIF();
  */
}

