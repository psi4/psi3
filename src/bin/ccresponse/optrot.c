#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <physconst.h>
#include <masses.h>
#define EXTERN
#include "globals.h"

void transmu(void);
void transp(double sign);
void sortmu(void);
void mubar(void);

void transL(double sign);
void sortL(void);
void Lbar(void);

void compute_X(char *pert, char *cart, int irrep, double omega);
double LCX(char *pert_c, char *cart_c, int irrep_c, 
	   char *pert_x, char *cart_x, int irrep_x, double omega);
double HXY(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
	   char *pert_y, char *cart_y, int irrep_y, double omega_y);
double LHX1Y1(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
	      char *pert_y, char *cart_y, int irrep_y, double omega_y);
double LHX2Y2(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
	      char *pert_y, char *cart_y, int irrep_y, double omega_y);
double LHX1Y2(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
	      char *pert_y, char *cart_y, int irrep_y, double omega_y);

void optrot(void)
{
  double ***tensor, **tensor0;
  double polar, polar_LCX, polar_HXY, polar_LHX1Y1, polar_LHX2Y2, polar_LHX1Y2;
  char **cartcomp;
  int alpha, beta, i, j, k;
  double TrG, M, nu, bohr2a4, m2a, hbar, prefactor, *rotation, *rotation_mod;

  cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");

  tensor = (double ***) malloc(params.nomega * sizeof(double **));
  for(i=0; i < params.nomega; i++) 
    tensor[i] = block_matrix(3,3);
  rotation = init_array(params.nomega);
  if(!strcmp(params.gauge,"VELOCITY")) rotation_mod = init_array(params.nomega);

  /* compute the zero-frequency Rosenfeld tensor for Koch's modified velocity optical rotation */
  if(!strcmp(params.gauge,"VELOCITY")) {

    tensor0 = block_matrix(3,3);

    fprintf(outfile, "\n\tComputing zero-frequency velocty-gauge linear response function...\n");

    /* prepare the dipole-length or dipole-velocity integrals */
    transp(+1.0);
    sortmu();
    mubar();

    /* prepare the magnetic-dipole integrals */
    transL(+1.0);
    sortL();
    Lbar();

    /* Compute the magnetic-dipole and electric-dipole CC wave functions */
    for(alpha=0; alpha < 3; alpha++) {
      compute_X("Mu", cartcomp[alpha], moinfo.l_irreps[alpha], 0);
      compute_X("L", cartcomp[alpha], moinfo.l_irreps[alpha], 0);
    }

    for(alpha=0; alpha < 3; alpha++) {
      for(beta=0; beta < 3; beta++) {

	polar_LCX = 0.0;
	polar_HXY = 0.0;
	polar_LHX1Y1 = 0.0;
	polar_LHX2Y2 = 0.0;
	polar_LHX1Y2 = 0.0;

	if(!(moinfo.mu_irreps[alpha]^moinfo.l_irreps[beta])) {

	  polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "L", cartcomp[beta],
			  moinfo.l_irreps[beta], 0.0);
	  polar_LCX += LCX("L", cartcomp[beta], moinfo.l_irreps[beta], "Mu", cartcomp[alpha],
			   moinfo.mu_irreps[alpha], 0.0);
	  polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
			  "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	  polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				"L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	  polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				"L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	  polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				"L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	  polar_LHX1Y2 += LHX1Y2("L", cartcomp[beta], moinfo.l_irreps[beta], 0.0,
				 "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0);

	  polar = polar_LCX + polar_HXY + polar_LHX1Y1 + polar_LHX2Y2 + polar_LHX1Y2;

	  tensor0[alpha][beta] = polar;
	}
      }
    }
  }

  for(i=0; i < params.nomega; i++) {

    /* prepare the dipole-length or dipole-velocity integrals */
    if(!strcmp(params.gauge,"LENGTH")) transmu();
    else if(!strcmp(params.gauge,"VELOCITY")) transp(+1.0);
    sortmu();
    mubar();

    /* prepare the magnetic-dipole integrals */
    transL(+1.0);
    sortL();
    Lbar();

    /* Compute the +omega magnetic-dipole and -omega electric-dipole CC wave functions */
    for(alpha=0; alpha < 3; alpha++) {
      compute_X("Mu", cartcomp[alpha], moinfo.l_irreps[alpha], -params.omega[i]);
      compute_X("L", cartcomp[alpha], moinfo.l_irreps[alpha], params.omega[i]);
    }

    for(alpha=0; alpha < 3; alpha++) {
      for(beta=0; beta < 3; beta++) {

	polar_LCX = 0.0;
	polar_HXY = 0.0;
	polar_LHX1Y1 = 0.0;
	polar_LHX2Y2 = 0.0;
	polar_LHX1Y2 = 0.0;

	if(!(moinfo.mu_irreps[alpha]^moinfo.l_irreps[beta])) {

	  if(params.omega[i] != 0.0) {
	    polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "L", cartcomp[beta], 
			    moinfo.l_irreps[beta], params.omega[i]);
	    polar_LCX += LCX("L", cartcomp[beta], moinfo.l_irreps[beta], "Mu", cartcomp[alpha],
			     moinfo.mu_irreps[alpha], -params.omega[i]);
	    polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i],
			    "L", cartcomp[beta], moinfo.l_irreps[beta], params.omega[i]);
	    polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i],
				  "L", cartcomp[beta], moinfo.l_irreps[beta], params.omega[i]);
	    polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i],
				  "L", cartcomp[beta], moinfo.l_irreps[beta], params.omega[i]);
	    polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i],
				  "L", cartcomp[beta], moinfo.l_irreps[beta], params.omega[i]);
	    polar_LHX1Y2 += LHX1Y2("L", cartcomp[beta], moinfo.l_irreps[beta], params.omega[i],
				   "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega[i]);
	  }
	  else {
	    polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "L", cartcomp[beta],
			    moinfo.l_irreps[beta], 0.0);
	    polar_LCX += LCX("L", cartcomp[beta], moinfo.l_irreps[beta], "Mu", cartcomp[alpha],
			     moinfo.mu_irreps[alpha], 0.0);
	    polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
			    "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	    polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				  "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	    polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				  "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	    polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				  "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	    polar_LHX1Y2 += LHX1Y2("L", cartcomp[beta], moinfo.l_irreps[beta], 0.0,
				   "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0);
	  }

	  polar = polar_LCX + polar_HXY + polar_LHX1Y1 + polar_LHX2Y2 + polar_LHX1Y2;
	  /*
	    fprintf(outfile, "polar_LCX    = %20.12f\n", polar_LCX);
	    fprintf(outfile, "polar_HXY    = %20.12f\n", polar_HXY);
	    fprintf(outfile, "polar_LHX1Y1 = %20.12f\n", polar_LHX1Y1);
	    fprintf(outfile, "polar_LHX2Y2 = %20.12f\n", polar_LHX2Y2);
	    fprintf(outfile, "polar_LHX1Y2 = %20.12f\n", polar_LHX1Y2);
	  */

	  tensor[i][alpha][beta] = 0.5 * polar;
	}
	/*      fprintf(outfile, "%1s%1s polar = %20.12f\n", cartcomp[alpha], cartcomp[beta], polar); */
      }
    }

    /* prepare the dipole-length or dipole-velocity integrals */
    if(!strcmp(params.gauge,"LENGTH")) transmu();
    else if(!strcmp(params.gauge,"VELOCITY")) transp(-1.0);
    sortmu();
    mubar();

    /* prepare the complex-conjugate of the magnetic-dipole integrals */
    transL(-1.0);
    sortL();
    Lbar();

    /* Compute the -omega magnetic-dipole and +omega electric-dipole CC wave functions */
    for(alpha=0; alpha < 3; alpha++) {
      compute_X("Mu", cartcomp[alpha], moinfo.l_irreps[alpha], params.omega[i]);
      compute_X("L", cartcomp[alpha], moinfo.l_irreps[alpha], -params.omega[i]);
    }

    for(alpha=0; alpha < 3; alpha++) {
      for(beta=0; beta < 3; beta++) {

	polar_LCX = 0.0;
	polar_HXY = 0.0;
	polar_LHX1Y1 = 0.0;
	polar_LHX2Y2 = 0.0;
	polar_LHX1Y2 = 0.0;

	if(!(moinfo.mu_irreps[alpha]^moinfo.l_irreps[beta])) {

	  if(params.omega[i] != 0.0) {
	    polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "L", cartcomp[beta], 
			    moinfo.l_irreps[beta], -params.omega[i]);
	    polar_LCX += LCX("L", cartcomp[beta], moinfo.l_irreps[beta], "Mu", cartcomp[alpha],
			     moinfo.mu_irreps[alpha], params.omega[i]);
	    polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega[i],
			    "L", cartcomp[beta], moinfo.l_irreps[beta], -params.omega[i]);
	    polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega[i],
				  "L", cartcomp[beta], moinfo.l_irreps[beta], -params.omega[i]);
	    polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega[i],
				  "L", cartcomp[beta], moinfo.l_irreps[beta], -params.omega[i]);
	    polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega[i],
				  "L", cartcomp[beta], moinfo.l_irreps[beta], -params.omega[i]);
	    polar_LHX1Y2 += LHX1Y2("L", cartcomp[beta], moinfo.l_irreps[beta], -params.omega[i],
				   "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega[i]);
	  }
	  else {
	    polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "L", cartcomp[beta],
			    moinfo.l_irreps[beta], 0.0);
	    polar_LCX += LCX("L", cartcomp[beta], moinfo.l_irreps[beta], "Mu", cartcomp[alpha],
			     moinfo.mu_irreps[alpha], 0.0);
	    polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
			    "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	    polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				  "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	    polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				  "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	    polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0,
				  "L", cartcomp[beta], moinfo.l_irreps[beta], 0.0);
	    polar_LHX1Y2 += LHX1Y2("L", cartcomp[beta], moinfo.l_irreps[beta], 0.0,
				   "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], 0.0);
	  }

	  polar = polar_LCX + polar_HXY + polar_LHX1Y1 + polar_LHX2Y2 + polar_LHX1Y2;
	  /*
	    fprintf(outfile, "polar_LCX    = %20.12f\n", polar_LCX);
	    fprintf(outfile, "polar_HXY    = %20.12f\n", polar_HXY);
	    fprintf(outfile, "polar_LHX1Y1 = %20.12f\n", polar_LHX1Y1);
	    fprintf(outfile, "polar_LHX2Y2 = %20.12f\n", polar_LHX2Y2);
	    fprintf(outfile, "polar_LHX1Y2 = %20.12f\n", polar_LHX1Y2);
	  */
	  tensor[i][alpha][beta] += 0.5 * polar;
	}
	/*      fprintf(outfile, "%1s%1s polar = %20.12f\n", cartcomp[alpha], cartcomp[beta], polar); */
      }
    }

    fprintf(outfile, "\n                      CCSD Optical Rotation Tensor:\n");
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[i],
	    (_c*_h*1e9)/(_hartree2J*params.omega[i]), _hartree2ev*params.omega[i],
	    _hartree2wavenumbers*params.omega[i]);
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    mat_print(tensor[i], 3, 3, outfile);

    /* compute the specific rotation */
    TrG = (tensor[i][0][0] + tensor[i][1][1] + tensor[i][2][2])/(3.0 * params.omega[i]);
    if(!strcmp(params.gauge,"VELOCITY")) TrG /= params.omega[i];

    for(j=0,M=0.0; j < moinfo.natom ;j++) M += an2masses[(int) moinfo.zvals[j]]; /* amu */
    nu = params.omega[i]; /* hartree */
    bohr2a4 = _bohr2angstroms * _bohr2angstroms * _bohr2angstroms * _bohr2angstroms;
    m2a = _bohr2angstroms * 1.0e-10;
    hbar = _h/(2.0 * _pi);
    prefactor = 1.0e-2 * hbar/(_c * 2.0 * _pi * _me * m2a * m2a);
    prefactor *= prefactor;
    prefactor *= 288.0e-30 * _pi * _pi * _na * bohr2a4;
    rotation[i] = prefactor * TrG * nu * nu / M;
    if(!strcmp(params.gauge,"VELOCTY"))
      fprintf(outfile, "\tSpecific rotation using simple velocity-gauge Rosenfeld tensor.\n");
    fprintf(outfile, "\n\t[alpha]_(%5.3f) = %20.12f deg/[dm (gm/cm^3)]\n", params.omega[i], rotation[i]);


    if(!strcmp(params.gauge,"VELOCITY")) {

      /* subtract the zero-frequency beta tensor */
      for(j=0; j < 3; j++)
	for(k=0; k < 3; k++)
	  tensor[i][j][k] -= tensor0[j][k];

      fprintf(outfile, "\n        CCSD Modified Velocity-Gauge Optical Rotation Tensor:\n");
      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega[i],
	      (_c*_h*1e9)/(_hartree2J*params.omega[i]), _hartree2ev*params.omega[i],
	      _hartree2wavenumbers*params.omega[i]);
      fprintf(outfile, "  -------------------------------------------------------------------------\n");
      mat_print(tensor[i], 3, 3, outfile);

      /* compute the specific rotation */
      TrG = (tensor[i][0][0] + tensor[i][1][1] + tensor[i][2][2])/(3.0 * params.omega[i]);
      if(!strcmp(params.gauge,"VELOCITY")) TrG /= params.omega[i];

      for(j=0,M=0.0; j < moinfo.natom ;j++) M += an2masses[(int) moinfo.zvals[j]]; /* amu */
      nu = params.omega[i]; /* hartree */
      bohr2a4 = _bohr2angstroms * _bohr2angstroms * _bohr2angstroms * _bohr2angstroms;
      m2a = _bohr2angstroms * 1.0e-10;
      hbar = _h/(2.0 * _pi);
      prefactor = 1.0e-2 * hbar/(_c * 2.0 * _pi * _me * m2a * m2a);
      prefactor *= prefactor;
      prefactor *= 288.0e-30 * _pi * _pi * _na * bohr2a4;
      rotation_mod[i] = prefactor * TrG * nu * nu / M;
      if(!strcmp(params.gauge,"VELOCTY"))
	fprintf(outfile, "\tSpecific rotation using zero-frequency renormalized\n\tvelocity-gauge Rosenfeld tensor.\n");
      fprintf(outfile, "\n\t[alpha]_(%5.3f) = %20.12f deg/[dm (gm/cm^3)]\n", params.omega[i], rotation_mod[i]);
    }
  }

  if(params.nomega > 1) {  /* print a summary table for multi-wavelength calcs */

    if(!strcmp(params.gauge,"VELOCITY")) {
      fprintf(outfile, "\n\t------------------------------------------------------\n");
      fprintf(outfile,   "\t         CCSD Velocity-Gauge Optical Rotation\n");
      fprintf(outfile,   "\t------------------------------------------------------\n");
      fprintf(outfile,   "\t    Omega           alpha (deg/[dm (gm/cm^3)]\n");
      fprintf(outfile, "\n\t E_h      nm   Velocity-Gauge  Modified Velocity-Gauge\n");
      fprintf(outfile,   "\t-----   ------ --------------  -----------------------\n");
      for(i=0; i < params.nomega; i++)
	fprintf(outfile, "\t%5.3f   %6.2f   %10.5f          %10.5f\n", params.omega[i], (_c*_h*1e9)/(_hartree2J*params.omega[i]), rotation[i], rotation_mod[i]);
    }
    else {
      fprintf(outfile, "\n\t------------------------------------------\n");
      fprintf(outfile,   "\t    CCSD Length-Gauge Optical Rotation\n");
      fprintf(outfile,   "\t------------------------------------------\n");
      fprintf(outfile,   "\t    Omega           alpha\n");
      fprintf(outfile,   "\t E_h      nm   deg/[dm (gm/cm^3)]\n");
      fprintf(outfile,   "\t-----   ------ ------------------\n");
       for(i=0; i < params.nomega; i++)
        fprintf(outfile, "\t%5.3f   %6.2f      %10.5f\n", params.omega[i], (_c*_h*1e9)/(_hartree2J*params.omega[i]), rotation[i]);
    }

  }

  free(rotation);
  if(!strcmp(params.gauge,"VELOCITY")) free(rotation_mod);

  for(i=0; i < params.nomega; i++)
    free_block(tensor[i]);
  free(tensor);

  if(!strcmp(params.gauge,"VELOCITY")) free_block(tensor0);

  free(cartcomp[0]);
  free(cartcomp[1]);
  free(cartcomp[2]);
  free(cartcomp);
}
