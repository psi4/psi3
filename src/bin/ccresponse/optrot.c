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
  double **tensor;
  double polar, polar_LCX, polar_HXY, polar_LHX1Y1, polar_LHX2Y2, polar_LHX1Y2;
  char **cartcomp;
  int alpha, beta, i;
  double TrG, M, nu, bohr2a4, m2a, hbar, prefactor, rotation;

  cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");

  tensor = block_matrix(3,3);

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
    compute_X("Mu", cartcomp[alpha], moinfo.l_irreps[alpha], -params.omega);
    compute_X("L", cartcomp[alpha], moinfo.l_irreps[alpha], params.omega);
  }

  for(alpha=0; alpha < 3; alpha++) {
    for(beta=0; beta < 3; beta++) {

      polar_LCX = 0.0;
      polar_HXY = 0.0;
      polar_LHX1Y1 = 0.0;
      polar_LHX2Y2 = 0.0;
      polar_LHX1Y2 = 0.0;

      if(!(moinfo.mu_irreps[alpha]^moinfo.l_irreps[beta])) {

	if(params.omega != 0.0) {
	  polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "L", cartcomp[beta], 
			  moinfo.l_irreps[beta], params.omega);
	  polar_LCX += LCX("L", cartcomp[beta], moinfo.l_irreps[beta], "Mu", cartcomp[alpha],
			   moinfo.mu_irreps[alpha], -params.omega);
	  polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
			  "L", cartcomp[beta], moinfo.l_irreps[beta], params.omega);
	  polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
				"L", cartcomp[beta], moinfo.l_irreps[beta], params.omega);
	  polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
				"L", cartcomp[beta], moinfo.l_irreps[beta], params.omega);
	  polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega,
				"L", cartcomp[beta], moinfo.l_irreps[beta], params.omega);
	  polar_LHX1Y2 += LHX1Y2("L", cartcomp[beta], moinfo.l_irreps[beta], params.omega,
				 "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], -params.omega);
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

	tensor[alpha][beta] = 0.5 * polar;
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
    compute_X("Mu", cartcomp[alpha], moinfo.l_irreps[alpha], params.omega);
    compute_X("L", cartcomp[alpha], moinfo.l_irreps[alpha], -params.omega);
  }

  for(alpha=0; alpha < 3; alpha++) {
    for(beta=0; beta < 3; beta++) {

      polar_LCX = 0.0;
      polar_HXY = 0.0;
      polar_LHX1Y1 = 0.0;
      polar_LHX2Y2 = 0.0;
      polar_LHX1Y2 = 0.0;

      if(!(moinfo.mu_irreps[alpha]^moinfo.l_irreps[beta])) {

	if(params.omega != 0.0) {
	  polar_LCX = LCX("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], "L", cartcomp[beta], 
			  moinfo.l_irreps[beta], -params.omega);
	  polar_LCX += LCX("L", cartcomp[beta], moinfo.l_irreps[beta], "Mu", cartcomp[alpha],
			   moinfo.mu_irreps[alpha], params.omega);
	  polar_HXY = HXY("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega,
			  "L", cartcomp[beta], moinfo.l_irreps[beta], -params.omega);
	  polar_LHX1Y1 = LHX1Y1("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega,
				"L", cartcomp[beta], moinfo.l_irreps[beta], -params.omega);
	  polar_LHX2Y2 = LHX2Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega,
				"L", cartcomp[beta], moinfo.l_irreps[beta], -params.omega);
	  polar_LHX1Y2 = LHX1Y2("Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega,
				"L", cartcomp[beta], moinfo.l_irreps[beta], -params.omega);
	  polar_LHX1Y2 += LHX1Y2("L", cartcomp[beta], moinfo.l_irreps[beta], -params.omega,
				 "Mu", cartcomp[alpha], moinfo.mu_irreps[alpha], params.omega);
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
	tensor[alpha][beta] += 0.5 * polar;
      }
      /*      fprintf(outfile, "%1s%1s polar = %20.12f\n", cartcomp[alpha], cartcomp[beta], polar); */
    }
  }

  fprintf(outfile, "\n                      CCSD Optical Rotation Tensor:\n");
  fprintf(outfile, "  -------------------------------------------------------------------------\n");
  fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega,
	  (_c*_h*1e9)/(_hartree2J*params.omega), _hartree2ev*params.omega,
	  _hartree2wavenumbers*params.omega);
  fprintf(outfile, "  -------------------------------------------------------------------------\n");
  mat_print(tensor, 3, 3, outfile);

  /* compute the specific rotation */
  TrG = (tensor[0][0] + tensor[1][1] + tensor[2][2])/(3.0 * params.omega);
  if(!strcmp(params.gauge,"VELOCITY")) TrG /= params.omega;

  for(i=0,M=0.0; i < moinfo.natom ;i++) M += an2masses[(int) moinfo.zvals[i]]; /* amu */
  nu = params.omega; /* hartree */
  bohr2a4 = _bohr2angstroms * _bohr2angstroms * _bohr2angstroms * _bohr2angstroms;
  m2a = _bohr2angstroms * 1.0e-10;
  hbar = _h/(2.0 * _pi);
  prefactor = 1.0e-2 * hbar/(_c * 2.0 * _pi * _me * m2a * m2a);
  prefactor *= prefactor;
  prefactor *= 288.0e-30 * _pi * _pi * _na * bohr2a4;
  rotation = prefactor * TrG * nu * nu / M;
  if(!strcmp(params.gauge,"VELOCTY"))
    fprintf(outfile, "\tSpecific rotation using simple velocity-gauge Rosenfeld tensor.\n");
  fprintf(outfile, "\n\t[alpha]_(%5.3f) = %20.12f deg/[dm (gm/cm^3)]\n", params.omega, rotation);

  /* compute the zero-frequency Rosenfeld tensor for Koch's modified velocity optical rotation */
  if(!strcmp(params.gauge,"VELOCITY")) {

    fprintf(outfile, "\n\tComputing zero-frequency velocty-gauge linear response function...\n");

    /* prepare the dipole-length or dipole-velocity integrals */
    transp(+1.0);
    sortmu();
    mubar();

    /* prepare the magnetic-dipole integrals */
    transL(+1.0);
    sortL();
    Lbar();

    /* Compute the +omega magnetic-dipole and -omega electric-dipole CC wave functions */
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

	  tensor[alpha][beta] -= polar;
	}
      }
    }

    fprintf(outfile, "\n        CCSD Modified Velocity-Gauge Optical Rotation Tensor:\n");
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    fprintf(outfile,   "   Evaluated at omega = %8.6f E_h (%6.2f nm, %5.3f eV, %8.2f cm-1)\n", params.omega,
	    (_c*_h*1e9)/(_hartree2J*params.omega), _hartree2ev*params.omega,
	    _hartree2wavenumbers*params.omega);
    fprintf(outfile, "  -------------------------------------------------------------------------\n");
    mat_print(tensor, 3, 3, outfile);

    /* compute the specific rotation */
    TrG = (tensor[0][0] + tensor[1][1] + tensor[2][2])/(3.0 * params.omega);
    if(!strcmp(params.gauge,"VELOCITY")) TrG /= params.omega;

    for(i=0,M=0.0; i < moinfo.natom ;i++) M += an2masses[(int) moinfo.zvals[i]]; /* amu */
    nu = params.omega; /* hartree */
    bohr2a4 = _bohr2angstroms * _bohr2angstroms * _bohr2angstroms * _bohr2angstroms;
    m2a = _bohr2angstroms * 1.0e-10;
    hbar = _h/(2.0 * _pi);
    prefactor = 1.0e-2 * hbar/(_c * 2.0 * _pi * _me * m2a * m2a);
    prefactor *= prefactor;
    prefactor *= 288.0e-30 * _pi * _pi * _na * bohr2a4;
    rotation = prefactor * TrG * nu * nu / M;
    if(!strcmp(params.gauge,"VELOCTY"))
      fprintf(outfile, "\tSpecific rotation using zero-frequency renormalize\n\tvelocity-gauge Rosenfeld tensor.\n");
    fprintf(outfile, "\n\t[alpha]_(%5.3f) = %20.12f deg/[dm (gm/cm^3)]\n", params.omega, rotation);
  }

  free_block(tensor);

  free(cartcomp[0]);
  free(cartcomp[1]);
  free(cartcomp[2]);
  free(cartcomp);
}
