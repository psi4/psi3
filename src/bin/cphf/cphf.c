/*
** CPHF: Program to solve the Coupled Perturbed Hartree-Fock equations
** for nuclear and electric field perturbations and to compute
** electric polarizabilities, harmonic vibrational frequencies, and IR
** intensities.
**
** Limitations of and future plans for this code:
**
** (1) Spin-restricted closed-shell Hartree-Fock (RHF) wave functions
** only.  Extension to ROHF and UHF cases is needed.
**
** (2) All two-electron integrals are held in core and used in the MO
** basis.  Out-of-core and AO-based algorithms are needed in mohess.c,
** cphf_X.c, and build_hessian.c in order to handle larger basis sets
** and to avoid the two-electron integral transformation.
**
** (3) Symmetry-blocking is used in most of the loops, but is not
** actually used to improve storage or computational order.  I've put
** this off because the nuclear perturbations (x, y, and z on each
** nucleus) are not yet symmetry-adapted.  Some effort in this area is
** needed.
**
** (4) Thermodynamic functions, including enthalpies, entropies, heat
** capacities, free energies, and partition functions can be computed,
** given the vibrational data computed in vibration.c.
**
** TDC, December 2001 and October 2002
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "globals.h"

void init_io(int argc, char *argv[]);
void exit_io(void);
void title(void);
void init_ioff(void);

void setup(void);
void cleanup(void);

void mohess(double **);
void cphf_X(double **, double ***);
void cphf_F(double **, double ***);
void polarize(double ***);
void vibration(double **, double **);
void build_hessian(double ***, double **);
void build_dipder(double ***, double **);

int main(int argc, char *argv[])
{
  int coord, errcod;
  double **A, ***UX, ***UF, **hessian, **dipder;

  init_io(argc, argv);
  title();
  init_ioff();

  timer_init();
  timer_on("CPHF Main");

  /* user-specified printing value */
  print_lvl = 0;
  errcod = ip_data("PRINT", "%d", &(print_lvl), 0);

  setup(); /* get useful data from chkpt and compute various lookup arrays */

  /* Grab the two-electron integrals */
  ints = init_array(ntei);
  iwl_rdtwo(PSIF_MO_TEI, ints, ioff, nmo, 0, 0, 0, outfile);

  /* Build the MO Hessian */
  A = block_matrix(num_ai,num_ai);
  mohess(A);

  /* init memory for UX matrices */
  UX = (double ***) malloc(natom*3 * sizeof(double **));
  for(coord=0; coord < natom*3; coord++) UX[coord] = block_matrix(nmo,nmo);

  /* solve CPHF for nuclear perturbations */
  cphf_X(A, UX);

  /* init memory for UF matrices */
  UF = (double ***) malloc(3 * sizeof(double **));
  for(coord=0; coord < 3; coord++) UF[coord] = block_matrix(nmo,nmo);

  /* solve CPHF for electric field perturbations */
  cphf_F(A, UF);

  /* compute the polarizability tensor */
  polarize(UF);

  /* Build the cartesian hessian */
  hessian = block_matrix(natom*3, natom*3);
  build_hessian(UX, hessian);

  /* Build the dipole moment derivatives */
  dipder = block_matrix(3, natom*3);
  build_dipder(UX, dipder);

  /* compute vibrational frequencies and ir intensities */
  vibration(hessian, dipder);

  /* Free memory */
  for(coord=0; coord < natom*3; coord++) free_block(UX[coord]);
  free(UX);
  for(coord=0; coord < 3; coord++) free_block(UF[coord]);
  free(UF);
  free_block(A); 
  free(ints);
  free_block(hessian);
  free_block(dipder);

  cleanup();  /* free memory allocated in setup(); */

  timer_off("CPHF Main");
  timer_done();

  exit_io();
  exit(0);
}

void init_io(int argc, char *argv[])
{
  extern char *gprgid(void);
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  psi_start(argc-1, argv+1,0);
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);

  psio_init();
}

void title(void)
{
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*          CPHF          *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
}

void exit_io(void)
{
  int i;
  psio_done();
  tstop(outfile);
  psi_stop();
}

char *gprgid(void)
{
   char *prgid = "CPHF";

   return(prgid);
}

void init_ioff(void)
{
  int i;
  ioff = (int *) malloc(IOFF_MAX * sizeof(int));
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) {
      ioff[i] = ioff[i-1] + i;
    }
}
