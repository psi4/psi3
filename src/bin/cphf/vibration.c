#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <physconst.h>
#include <masses.h>
#define EXTERN
#include "globals.h"

/* vibration(): Computes the harmonic vibrational frequencies and IR
** integrated absorption coefficients (intensities) using the
** cartesian hessian from build_hessian() and dipole derivatives from
** build_dipder().  
**
** The procedure for computing the normal modes is described in:
**
**  E.B. Wilson, J.C Decius, and P.C. Cross, "Molecular Vibrations:
** The Theory of Infrared and Raman Vibrational Spectra", Dover, New
** York, 1955.
**
** The relationship between IR intensities and dipole derivatives is
** described in the same text (Ch. 7) as well as in:
**
** I. Mills, T. Cvitas, K. Homann, N. Kallay, and K. Kuchitsu,
** "Quantities, Units, and Symbols in Physical Chemistry", Oxford,
** London, 1993, pp. 33-35.
**
** B.S. Galabov and T. Dudev, "Vibrational Intensities", v.22 of
** Vibrational Spectra and Structure, J.R. Durig, ed., Elsevier, 1996,
** pp. 2-12.
**
** G. Zerbi, "Introduction to the Theory of Vibrational Frequencies
** and Vibrational Intensities," in Vibrational Intensities in
** Infrared and Raman Spectroscopy, v.20 of Studies in Physical and
** Theoretical Chemistry, W.B. Person and G. Zerbi, eds., Elsevier
** 1982, pp. 45-51.
**
** TDC, October 2002.
*/

void vibration(double **hessian, double **dipder)
{
  int i, j;
  double **M, *irint;
  double **lx, **dipder_q, **TMP;
  double *km, k_convert, cm_convert;
  double *work;
  int stat;
  double dipder_conv, ir_prefactor;

  /* mass-weight the hessian */
  M = block_matrix(natom*3, natom*3);
  for(i=0; i < natom; i++) {
    for(j=0; j < 3; j++)  {
      M[i*3+j][i*3+j] = 1/sqrt(an2masses[(int) zvals[i]]);
    }
  }

  TMP = block_matrix(natom*3,natom*3);
  C_DGEMM('n','n', natom*3, natom*3, natom*3, 1.0, &(M[0][0]), natom*3,
          &(hessian[0][0]), natom*3, 0.0, &(TMP[0][0]), natom*3);
  C_DGEMM('n','n', natom*3, natom*3, natom*3, 1.0, &(TMP[0][0]), natom*3,
	  &(M[0][0]), natom*3, 0.0, &(hessian[0][0]), natom*3);
  free_block(TMP);

  /* diagonalize mass-weighted hessian */
  km = init_array(natom*3);  /* mass-weighted force constants */
  work = init_array(natom*3*3); /* scratch array */

  if(stat = C_DSYEV('v','u',natom*3,&(hessian[0][0]),natom*3,&(km[0]),&(work[0]),natom*3*3)) {
    fprintf(outfile, "vibration(): Error in hessian diagonalization. stat = %d\n", stat);
    exit(2);
  }

  /* Construct the mass-weighted normal coordinates, lx */
  /* (note, after C_DSYEV, hessian contains the eigenvectors) */
  lx = block_matrix(natom*3,natom*3);
  C_DGEMM('n','t', natom*3,natom*3,natom*3,1,&(M[0][0]),natom*3,
	  &(hessian[0][0]),natom*3,0,&(lx[0][0]),natom*3);

  if(print_lvl & 1) {
    fprintf(outfile, "\n\tLX matrix:\n");
    print_mat(lx, natom*3, natom*3, outfile);
  }

  /* Transform dipole derivatives to normal coordinates */
  dipder_q = block_matrix(3, natom*3);
  C_DGEMM('n','n', 3, natom*3, natom*3, 1, &(dipder[0][0]),natom*3,
	  &(lx[0][0]), natom*3, 0, &(dipder_q[0][0]), natom*3);

  fprintf(outfile, "\n\tDipole Derivatives W.R.T Normal Coordinates (debye/A-amu+1/2):\n");
  print_mat(dipder_q, 3, natom*3, outfile);

  /* Compute the IR intensities in D^2/(A^2 amu) */
  irint = init_array(natom*3);
  for(i=0; i < natom*3; i++)
    for(j=0; j < 3; j++)
      irint[i] += dipder_q[j][i] * dipder_q[j][i];

  /* conversion factor from D^2/(A^2 amu) to C^2/kg */
  dipder_conv = _dipmom_debye2si*_dipmom_debye2si/(1e-20 * _amu2kg);
  for(i=0; i < natom*3; i++)
    irint[i] *= dipder_conv;

  /* IR integrated absorption coefficient prefactor */
  ir_prefactor = _na * _pi/(3.0 * _c * _c * 4.0 * _pi * _e0 * 1000.0);

/*  fprintf(outfile, "\n\tIR conversion = %20.10f\n", dipder_conv * ir_prefactor); */

  /* compute the frequencies and spit them out in a nice table */
  fprintf(outfile, "\n\t        Harmonic Frequency   Infrared Intensity\n");
  fprintf(outfile,   "\t              (cm-1)               (km/mol)    \n");
  fprintf(outfile,   "\t-----------------------------------------------\n");
  k_convert = _hartree2J/(_bohr2m * _bohr2m * _amu2kg);
  cm_convert = 1.0/(2.0 * _pi * _c * 100.0);
  for(i=natom*3-1; i >= 0; i--) {
    if(km[i] < 0.0)
      fprintf(outfile, "\t  %3d   %17.3fi       %10.4f\n", i, 
              cm_convert * sqrt(-k_convert * km[i]), irint[i]*ir_prefactor);
    else
      fprintf(outfile, "\t  %3d   %17.3f        %10.4f\n", i, 
              cm_convert * sqrt(k_convert * km[i]), irint[i]*ir_prefactor);
  }
  fprintf(outfile,   "\t-----------------------------------------------\n");

  free(work);
  free(irint);
  free_block(lx);
  free_block(dipder_q);
  free_block(M);
  free(km);
}
