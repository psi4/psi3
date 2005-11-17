#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#define EXTERN
#include "globals.h"

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
double cc2_LHX1Y1(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
		  char *pert_y, char *cart_y, int irrep_y, double omega_y);
double cc2_LHX1Y2(char *pert_x, char *cart_x, int irrep_x, double omega_x, 
		  char *pert_y, char *cart_y, int irrep_y, double omega_y);

void linresp(double **tensor, double A, double B,
	     char *pert_x, int *x_irreps, double omega_x,
	     char *pert_y, int *y_irreps, double omega_y)
{
  int alpha, beta;
  double polar, polar_LCX, polar_HXY, polar_LHX1Y1, polar_LHX1Y2, polar_LHX2Y2;
  char **cartcomp;

  cartcomp = (char **) malloc(3 * sizeof(char *));
  cartcomp[0] = strdup("X");
  cartcomp[1] = strdup("Y");
  cartcomp[2] = strdup("Z");

  for(alpha=0; alpha < 3; alpha++) {
    for(beta=0; beta < 3; beta++) {

      polar_LCX = 0.0;
      polar_HXY = 0.0;
      polar_LHX1Y1 = 0.0;
      polar_LHX2Y2 = 0.0;
      polar_LHX1Y2 = 0.0;

      if(!(x_irreps[alpha]^y_irreps[beta])) {

	if(omega_y != 0.0) {  /* we assume omega_x = -omega_y */
	  polar_LCX = LCX(pert_x, cartcomp[alpha], x_irreps[alpha], 
			  pert_y, cartcomp[beta], y_irreps[beta], omega_y);
	  polar_LCX += LCX(pert_y, cartcomp[beta], y_irreps[beta], 
			   pert_x, cartcomp[alpha], x_irreps[alpha], omega_x);

	  if (!strcmp(params.wfn,"CC2")) {
	    polar_HXY = HXY(pert_x, cartcomp[alpha], x_irreps[alpha], omega_x,
			    pert_y, cartcomp[beta], y_irreps[beta], omega_y);
	    polar_LHX1Y1 = cc2_LHX1Y1(pert_x, cartcomp[alpha], x_irreps[alpha], omega_x,
				      pert_y, cartcomp[beta], y_irreps[beta], omega_y);
	    polar_LHX1Y2 = cc2_LHX1Y2(pert_x, cartcomp[alpha], x_irreps[alpha], omega_x,
				      pert_y, cartcomp[beta], y_irreps[beta], omega_y);
	    polar_LHX1Y2 += cc2_LHX1Y2(pert_y, cartcomp[beta], y_irreps[beta], omega_y,
				       pert_x, cartcomp[alpha], x_irreps[alpha], omega_x);
	  }
	  else {
	    polar_LHX1Y1 = LHX1Y1(pert_x, cartcomp[alpha], x_irreps[alpha], omega_x,
				  pert_y, cartcomp[beta], y_irreps[beta], omega_y);
	    polar_LHX2Y2 = LHX2Y2(pert_x, cartcomp[alpha], x_irreps[alpha], omega_x,
				  pert_y, cartcomp[beta], y_irreps[beta], omega_y);
	    polar_LHX1Y2 = LHX1Y2(pert_x, cartcomp[alpha], x_irreps[alpha], omega_x,
				  pert_y, cartcomp[beta], y_irreps[beta], omega_y);
	    polar_LHX1Y2 += LHX1Y2(pert_y, cartcomp[beta], y_irreps[beta], omega_y,
				   pert_x, cartcomp[alpha], x_irreps[alpha], omega_x);
	  }
	}
	else {
	  polar_LCX = LCX(pert_x, cartcomp[alpha], x_irreps[alpha], pert_y, cartcomp[beta],
			  y_irreps[beta], 0.0);
	  polar_LCX += LCX(pert_y, cartcomp[beta], y_irreps[beta], pert_x, cartcomp[alpha],
			   x_irreps[alpha], 0.0);
	  if (!strcmp(params.wfn,"CC2")) {
	    polar_HXY = HXY(pert_x, cartcomp[alpha], x_irreps[alpha], 0.0,
			    pert_y, cartcomp[beta], y_irreps[beta], 0.0);
	    polar_LHX1Y1 = cc2_LHX1Y1(pert_x, cartcomp[alpha], x_irreps[alpha], 0.0,
				      pert_y, cartcomp[beta], y_irreps[beta], 0.0);
	    polar_LHX1Y2 = cc2_LHX1Y2(pert_x, cartcomp[alpha], x_irreps[alpha], 0.0,
				      pert_y, cartcomp[beta], y_irreps[beta], 0.0);
	    polar_LHX1Y2 += cc2_LHX1Y2(pert_y, cartcomp[beta], y_irreps[beta], 0.0,
				       pert_x, cartcomp[alpha], x_irreps[alpha], 0.0);
	  }
	  else {
	    polar_LHX1Y1 = LHX1Y1(pert_x, cartcomp[alpha], x_irreps[alpha], 0.0,
				  pert_y, cartcomp[beta], y_irreps[beta], 0.0);
	    polar_LHX2Y2 = LHX2Y2(pert_x, cartcomp[alpha], x_irreps[alpha], 0.0,
				  pert_y, cartcomp[beta], y_irreps[beta], 0.0);
	    polar_LHX1Y2 = LHX1Y2(pert_x, cartcomp[alpha], x_irreps[alpha], 0.0,
				  pert_y, cartcomp[beta], y_irreps[beta], 0.0);
	    polar_LHX1Y2 += LHX1Y2(pert_y, cartcomp[beta], y_irreps[beta], 0.0,
				   pert_x, cartcomp[alpha], x_irreps[alpha], 0.0);
	  }
	}

	polar = polar_LCX + polar_HXY + polar_LHX1Y1 + polar_LHX2Y2 + polar_LHX1Y2;

	if(params.print & 2) {
	  fprintf(outfile, "\tLinresp tensor[%s][%s]\n", cartcomp[alpha], cartcomp[beta]);
	  if(!strcmp(params.wfn,"CC2"))
	    fprintf(outfile, "\tpolar_LCX    = %20.12f\n", polar_LCX);
	  fprintf(outfile, "\tpolar_HXY    = %20.12f\n", polar_HXY);
	  fprintf(outfile, "\tpolar_LHX1Y1 = %20.12f\n", polar_LHX1Y1);
	  fprintf(outfile, "\tpolar_LHX1Y2 = %20.12f\n", polar_LHX1Y2);
	  fprintf(outfile, "\tpolar_LHX2Y2 = %20.12f\n\n", polar_LHX2Y2);
	  fflush(outfile);
	}

	tensor[alpha][beta] = A * polar + B * tensor[alpha][beta];
      }
    }
  }

  free(cartcomp[0]);
  free(cartcomp[1]);
  free(cartcomp[2]);
  free(cartcomp);
}