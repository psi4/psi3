/*
 * NONBONDED
 *
 * Evaluate nonbonded interactions such as van der Waals terms and
 * point-charge electrostatics
 *
 * C. David Sherrill
 * January 2008
 */

/*! 
** \file
** \ingroup NONBONDED
** \brief Evaluate empirical non-bonded interactions
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <psifiles.h>
#include <string.h>
#include <physconst.h>
#include <psifiles.h>

#include "globals.h"
#include "nonbonded.h"

namespace psi { namespace nonbonded {

/* function prototypes this module */
void start_io(int argc, char *argv[]);
void stop_io(void);
void print_intro(void);
void compute_R(int natom, double **geom, double *R);
double compute_ddisp(int natom, double *R, double *AN, double s6, double d);

}}; // close namespace decl

using namespace psi::nonbonded;

int main(int argc, char *argv[]) {

  double **geom;                  /* geometry matrix (cols are x,y,z)   */
  double *AC;                     /* atomic charges                     */
  double *AN;                     /* atomic numbers                     */
  double *R;                      /* distance matrix (lwr triangle)     */
  int natom;                      /* number of atoms                    */
  int have_partial_charges=0;     /* flag for partial charges available */
  double s6=1.0;                  /* global scaling parameter (Grimme),
                                     0.75 PBE, 1.2 BLYP, 1.05 B-P86,
                                     1.0 TPSS, 1.05 B3LYP               */
  double d=20.0;                  /* exponent for damping term (Grimme) */
  double energy_dd;               /* damped dispersion energy (J/mol)   */
  double energy_dd_hartree;       /* in hartree                         */
  int errcod;

  start_io(argc,argv);
  print_intro();
  geom = chkpt_rd_geom();
  natom = chkpt_rd_natom();
  AN = chkpt_rd_zvals();
  if (ip_exist("PARTIAL_CHARGES",0)) {
    AC = init_array(natom);
    errcod = ip_double_array("PARTIAL_CHARGES",AC,natom);
  }
  errcod = ip_data("S6","%lf",&s6,0);
  errcod = ip_data("D_DMP","%lf",&d,0);

  R = init_array((natom*(natom+1))/2);
  compute_R(natom, geom, R);
  energy_dd = compute_ddisp(natom, R, AN, s6, d);
  energy_dd_hartree = energy_dd / (_na * _hartree2J);

  fprintf(outfile, "\nDamped dispersion energy:");
  fprintf(outfile, " %12.9lf hartree (%12.9lf kcal/mol)\n\n", 
    energy_dd_hartree, (energy_dd / 4184));
  fflush(outfile);

  /* clean up */
  free(AN);
  free(R);
  if (have_partial_charges) free(AC);
  /* free geom --- how exactly is it allocated now? */
  stop_io();
}

extern "C" { char *gprgid(void) { char *prgid = "TRANSQT"; return (prgid); } }

namespace psi { namespace nonbonded {

/*!
 * start_io()
 *
 * Initiate PSI IO routines
 *
 * \ingroup NONBONDED
 */
void start_io(int argc, char *argv[])
{
  int errcod;
  
  errcod = psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0);
  if (errcod != PSI_RETURN_SUCCESS)
    abort();
  ip_cwk_add(":NONBONDED");
  tstart(outfile);
  psio_init(); psio_ipv1_config();
  chkpt_init(PSIO_OPEN_OLD);
   
  return;
}


/*!
 * stop_io()
 * 
 * Shut down PSI IO
 */
void stop_io(void)
{
  tstop(outfile);
  psio_done();
  psi_stop(infile,outfile,psi_file_prefix);
}
 


/*!
 * print_intro()
 * 
 * Print an intro for the module
 *
 * \ingroup NONBONDED
 */
void print_intro(void)
{
 fprintf(outfile,"             ------------------------------------------\n");
 fprintf(outfile,"                            NONBONDED                  \n");
 fprintf(outfile,"              Evaluate empirical non-bonded terms      \n");
 fprintf(outfile,"                                                       \n");
 fprintf(outfile,"                        C. David Sherrill              \n");
 fprintf(outfile,"             ------------------------------------------\n");
}


/*!
 * compute_R
 *
 * Compute the distances between each pair of atoms.  Convert to Angstrom
 * because those are the units Grimme's parameters are in
 *
 * \ingroup NONBONDED
 */
void compute_R(int natom, double **geom, double *R)
{
  int i, j, ij;
  double dx, dy, dz, dist;

  for (i=0,ij=0; i<natom; i++) {
    for (j=0; j<i; j++,ij++) {
      dx = geom[i][0] - geom[j][0];
      dy = geom[i][1] - geom[j][1];
      dz = geom[i][2] - geom[j][2];
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      R[ij] = dist * _bohr2angstroms;
    }
  }
 
}



/*!
 * compute_ddisp()
 *
 * Compute damped dispersion terms (ala S. Grimme, J. Comput. Chem. 27,
 * 1787, 2006)
 *
 * \ingroup NONBONDED
 */
double compute_ddisp(int natom, double *R, double *AN, double s6, double d)
{
  int i, j, ij, Zi, Zj;
  double energy=0.0, tval, r, r_i, r_j, r_vdw;
  double C6i, C6j; 

  /* loop over unique pairs of atoms and evaluate the C6 dispersion terms */
  for (i=0,ij=0; i<natom; i++) {
    Zi = (int) AN[i];
    C6i = vdw_C6_grimme[Zi];
    r_i = vdw_radii_grimme[Zi];
    for (j=0; j<i; j++,ij++) {
      Zj = (int) AN[j];
      r = R[ij];
      C6j = vdw_C6_grimme[Zj];
      r_j = vdw_radii_grimme[Zj];
      r_vdw = r_i + r_j;
      tval = sqrt(C6i * C6j) * 1000000.0; /* nm^6 -> Ang^6 */
      tval = tval/pow(r, 6.0);
      tval = tval / (1.0 + exp(-d * (r / r_vdw - 1.0))); 
      energy += tval;
    }
  } 

  energy = -s6 * energy; /* in J mol^-1 */
  return energy;
}

}} // close off psi::nonbonded
