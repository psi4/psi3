/*** OPT_ENERGIES handles opt by energy points, Rollin  King, 2002 ***/ 

// this file should now be obseleted

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

double **compute_B(internals &simples, salc_set &symm);
void new_geom(cartesians &carts, internals &simples, salc_set &symm, double *dq,
    int print_to_geom_file, int restart_geom_file, 
    char *disp_label, int disp_num, int last_disp, double *return_geom);

extern void opt_step(cartesians &carts, internals &simples, salc_set &symm);

void opt_energies(cartesians &carts, internals &simples, salc_set &symm);

void opt_energies(cartesians &carts, internals &simples, salc_set &symm) {

  int i,j,a,b, dim, dim_carts, num_disps;
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

  // if first optking run at current geometry
  if (optinfo.micro_iteration == 0) {

    /** Store current step geometry and energy ***/
    geom = carts.get_coord();
    energy = carts.get_energy();
    open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Cart geom last step",
        (char *) &(geom[0]), dim_carts*sizeof(double));
    psio_write_entry(PSIF_OPTKING, "Energy last step",
        (char *) &(energy), sizeof(double));
    close_PSIF();
    free(geom);
    fprintf(outfile,"\nWriting geometry and energy of step to PSIF.\n");

    /*** make list of internal displacements for micro_iterations ***/
    displacements = block_matrix(num_disps,symm.get_num());
    for (i=0;i<symm.get_num();++i) {
      displacements[2*i][i] = -1.0 * optinfo.disp_size;
      displacements[2*i+1][i] = 1.0 * optinfo.disp_size;
    }
    //  fprintf(outfile,"\nDisplacement Matrix\n");
    //  print_mat2(displacements, num_disps, symm.get_num(), outfile);

    /*** generate and store Micro_iteration cartesian geometries ***/
    micro_geoms = block_matrix(num_disps, dim_carts);
    for (i=0;i<num_disps;++i)  {
      sprintf(disp_label,"Displaced geometry %d in a.u.\n",i+1);
      new_geom(carts,simples,symm,displacements[i],0,
          0, disp_label, i, 0, micro_geoms[i]);
    }
    free_block(displacements);

    open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Micro_iteration cart geoms",
        (char *) &(micro_geoms[0][0]), num_disps*dim_carts*sizeof(double));
    close_PSIF();
    free_block(micro_geoms);
  }

  // add energy to PSIF list of energies
  if (optinfo.micro_iteration > 0) {
    energies = new double[num_disps];
    energy = energy_chkpt();
    open_PSIF();
    if (psio_tocscan(PSIF_OPTKING, "Energy of Displacements") != NULL)
      psio_read_entry(PSIF_OPTKING, "Energy of Displacements",
          (char *) &(energies[0]), num_disps*sizeof(double));
    energies[optinfo.micro_iteration-1] = energy;
    psio_write_entry(PSIF_OPTKING, "Energy of Displacements",
        (char *) &(energies[0]), num_disps*sizeof(double));
    close_PSIF();
    free(energies);
  }

  /*** if it is NOT time to take a step, 
   ++ micro_iteration and load next geometry into chkpt  ***/
  if (optinfo.micro_iteration < num_disps) {
    micro_geoms = block_matrix(num_disps, dim_carts);

    open_PSIF();
    psio_read_entry(PSIF_OPTKING, "Micro_iteration cart geoms",
        (char *) &(micro_geoms[0][0]), num_disps*dim_carts*sizeof(double));
    close_PSIF();

    geom2D = block_matrix(carts.get_natom(),3);
    for (i=0; i<carts.get_natom(); ++i)
      for (j=0; j<3; ++j)
        geom2D[i][j] = micro_geoms[optinfo.micro_iteration][3*i+j];

    chkpt_init(PSIO_OPEN_OLD);
    chkpt_wt_geom(geom2D);
    chkpt_close();
    fprintf(outfile, "\nGeometry for displacement %d sent to chkpt.\n",optinfo.micro_iteration+1);
    free_block(micro_geoms);
    free_block(geom2D);
    open_PSIF();
    optinfo.micro_iteration += 1;
    psio_write_entry(PSIF_OPTKING, "Micro_iteration",
        (char *) &(optinfo.micro_iteration),sizeof(int));
    close_PSIF();
    psio_done();
    ip_done();
    exit(0);
  }

  /*** if time to produce a file11 and take a step ***/
  else if (optinfo.micro_iteration == num_disps) {

    // report energies
    open_PSIF();
    energies = new double[num_disps];
    psio_read_entry(PSIF_OPTKING, "Energy of Displacements",
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

    // Transform forces to cartesian coordinates
    B = compute_B(simples, symm);
    f = new double[dim_carts];
    mmult(B,1,&f_q,1,&f,1,dim_carts,symm.get_num(),1,0);
    free_block(B);

    // change forces to gradient for writing a file11 entry
    for(i=0;i<dim_carts;++i)
      f[i] = -1.0 * f[i] / _hartree2J / 1.0E18 * _bohr2angstroms;

    // write out approximate file11.dat
    geom = new double[dim_carts];
    open_PSIF();
    psio_read_entry(PSIF_OPTKING, "Cart geom last step",
        (char *) &(geom[0]), dim_carts*sizeof(double));
    psio_read_entry(PSIF_OPTKING, "Energy last step",
        (char *) &(energy), sizeof(double));
    close_PSIF();

    fp_11 = fopen("file11.dat","a");
    sprintf(disp_label,"iteration: %d", optinfo.iteration);
    carts.set_energy(energy);
    carts.set_coord(geom);
    carts.set_grad(f);
    free(f);
    free(geom);
    carts.print(11,fp_11,0,disp_label, 0);
    fclose(fp_11);

    // recompute values of internals and s vectors
    simples.compute_internals(carts.get_natom(),carts.get_coord());
    simples.compute_s(carts.get_natom(),carts.get_coord() );

    // take a step
    opt_step(carts, simples, symm);

    optinfo.micro_iteration = 0;
    open_PSIF();
    psio_write_entry(PSIF_OPTKING, "Micro_iteration",
        (char *) &(optinfo.micro_iteration),sizeof(int));
    close_PSIF();
  }

}

