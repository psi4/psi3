#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr.h>
#include <ip_libv1.h>
#include "input.h"
#include "global.h"
#include "defines.h"

void read_zmat()
{
  int i, j, a, b, c, errcod;
  char *buffer;
  int num_entries, entry_length, atomcount;
  int A, B, C, D;
  int linearOn = 1;	/* Flag indicating the code is working on the linear fragment in the begginning of the Z-matrix */
  double rAB, rBC, rCD, thetaABC, thetaBCD, phiABCD, val, norm1, norm2;
  double cosABC, sinABC, cosBCD, sinBCD, cosABCD, sinABCD;
  double eAB[3], eBC[3], ex[3], ey[3], ez[3];
  double Z = 0.0;
  double **full_geom;   /* Full matrix of coordinates (including those of dummy atoms) */

  /* Read number of lines and count atoms in ZMAT */
  num_atoms = num_entries = 0;
  ip_count("ZMAT",&num_entries,0);
  if (num_entries == 0)
    punt("Z-matrix is empty!");
  for(i=0;i<num_entries;i++){
    errcod = ip_string("ZMAT",&buffer,2,i,0);
    if (errcod != IPE_OK)
      punt("Problem with the Z-matrix.");
    if (strcmp(buffer,"X")) {
      free(buffer);
      num_atoms++;
      if (num_atoms > MAXATOM)
	punt("There are more atoms than allowed!");
    }
  }
  if (num_atoms == 0)
    punt("Z-matrix contains no atoms!");

  /*-----------------------
    Allocate global arrays
   -----------------------*/
  geometry = init_matrix(num_atoms,3);
  element = (char **) malloc(sizeof(char *)*num_atoms);
  nuclear_charges = init_array(num_atoms);

  full_geom = init_matrix(num_entries,3);
  atomcount = 0;
  
  for (i=0;i<num_entries;++i) {
     /* Process a line of ZMAT */
     ip_count("ZMAT",&entry_length,1,i);
     if ( ((i == 0) && (entry_length != 1)) ||
          ((i == 1) && (entry_length != 3)) ||
          ((i == 2) && (entry_length != 5)) ||
          ((i  > 2) && (entry_length != 7)) ) {
       fprintf(outfile,"  Line %d of ZMAT has a wrong number of entries.\n",i+1);
       punt("Invalid ZMAT");
     }
	

     if (i == 0) {					/*	1st atom */
        full_geom[i][0] = 0.0;
        full_geom[i][1] = 0.0;
        full_geom[i][2] = 0.0;
     }

     else if (i == 1) {					/*	2nd atom */
        ip_data("ZMAT","%d",&a,2,i,1);
        if (a != 1)
          punt("Problem in line 2 in zmat.");
        ip_data("ZMAT","%lf",&rAB,2,i,2);
        if (rAB < ZERO_BOND_DISTANCE)
           punt("Invalid bond length in line 2.");
        full_geom[i][0] = 0.0;
        full_geom[i][1] = 0.0;
        full_geom[i][2] = rAB;
     }

     else if (i == 2) {					/*	3rd atom */
        ip_data("ZMAT","%d",&a,2,i,1);
        ip_data("ZMAT","%d",&b,2,i,3);
        if ( ((a == 2) && (b == 1)) ||
             ((a == 1) && (b == 2)) ) {
           ip_data("ZMAT","%lf",&rBC,2,i,2);
           if (rBC <= ZERO_BOND_DISTANCE) {
              fprintf(outfile,"  Invalid bond length in line 3.\n");
              punt("Invalid ZMAT");
           }
           
           ip_data("ZMAT","%lf",&thetaABC,2,i,4);
           if (thetaABC <= ZERO_BOND_ANGLE) {
              fprintf(outfile,"  Invalid bond angle in line 3.\n");
              punt("Invalid ZMAT");
           }
           thetaABC = thetaABC*M_PI/180.0;
           
           if (a == 2) {				/*	ABC case */
              full_geom[i][0] = full_geom[a-1][0] + rBC*sin(thetaABC);
              full_geom[i][1] = full_geom[a-1][1];
              full_geom[i][2] = full_geom[a-1][2] - rBC*cos(thetaABC);
           }
           else {					/*	BAC case */
              full_geom[i][0] = full_geom[a-1][0] + rBC*sin(thetaABC);
              full_geom[i][1] = full_geom[a-1][1];
              full_geom[i][2] = full_geom[a-1][2] + rBC*cos(thetaABC);
           }
        }
        else {
           fprintf(outfile,"  Problem in line 3 in zmat.\n");
           punt("Invalid ZMAT");;
        }
     }

     else { 
        ip_data("ZMAT","%d",&c,2,i,1);
        ip_data("ZMAT","%d",&b,2,i,3);
        ip_data("ZMAT","%d",&a,2,i,5);
        a -= 1; b -= 1; c -= 1;

        if ( (a == b) || (b == c) || (a == c) ||
             (a >= i) || (b >= i) || (c >= i) ) {
           fprintf(outfile,"  Problem in line %d of zmat.\n",i);
           punt("Invalid ZMAT");
        }

	ip_data("ZMAT","%lf",&rCD,2,i,2);
	if (rCD <= ZERO_BOND_DISTANCE) {
	   fprintf(outfile,"  Invalid bond length in line %d.\n",i+1);
	   punt("Invalid ZMAT");
	}
	
	ip_data("ZMAT","%lf",&thetaBCD,2,i,4);
	if (thetaBCD <= ZERO_BOND_ANGLE) {
	   fprintf(outfile,"  Invalid bond angle in line %d.\n",i+1);
	   punt("Invalid ZMAT");
	}
	thetaBCD = thetaBCD * M_PI/180.0;
	
	ip_data("ZMAT","%lf",&phiABCD,2,i,6);
	phiABCD = phiABCD * M_PI/180.0;


	/* If you want to have linear fragment defined in 
	   the beginning of the Z-matrix - fine, but you still have to
	   supply the third "dummy" atom and a value for the dihedral angle*/

	if ( ((full_geom[i-1][0] - LINEAR_CUTOFF) < 0.0) && linearOn ) {
	   if ((full_geom[c][2] - full_geom[b][2]) > 0.0) {
	      full_geom[i][0] = rCD*sin(thetaBCD);
	      full_geom[i][1] = 0.0;
	      full_geom[i][2] = full_geom[c][2] - rCD*cos(thetaBCD);
	   }
	   else {
	      full_geom[i][0] = rCD*sin(thetaBCD);
	      full_geom[i][1] = 0.0;
	      full_geom[i][2] = full_geom[c][2] + rCD*cos(thetaBCD);
	   }
	}

	/* Here starts a piece of code for nonlinear ABC fragment */

	else {
	   	/* Set "linearity" Flag to false so that the current "if" test 
	   	   is "false" ever since */
	   linearOn = 0;
	   
		/* Note : x,y,z - unit vectors defining a coordinate 
		   system associated with atom C; BC defines z, ABC defines 
		   the xz-plane */
	   unit_vec(full_geom[b],full_geom[a],eAB);	/* U1 in the old zmat */
	   unit_vec(full_geom[c],full_geom[b],ez);	/* U2 in the old zmat */
	   cosABC = -dot_prod(ez,eAB);	/* ez is essentially eBC */
	   sinABC = sqrt(1 - (cosABC * cosABC) );
	   
	   if ( (sinABC - LINEAR_CUTOFF) < 0.0 ) {
	     fprintf(outfile,"  Dihedral angle in line %d is defined with respect to a linear fragment.\n",i+1);
	     punt("Invalid ZMAT");
	   }
	   
	   cross_prod(eAB,ez,ey);		/* ey is U3 in the old zmat */
	   for(j=0;j<3;j++)			/* normalization of ey */
	     ey[j] /= sinABC;
	   cross_prod(ey,ez,ex);		/* ex is U4 in the old zmat */

	   /* Intermediates - to avoid calling 
	                      trig. functions multiple times */
	   cosBCD = cos(thetaBCD);
	   sinBCD = sin(thetaBCD);
	   cosABCD = cos(phiABCD);
	   sinABCD = sin(phiABCD);
	   
	   for (j=0;j<3;j++)
	     full_geom[i][j] = full_geom[c][j] + rCD * 
	                  ( - ez[j] * cosBCD + 
	                    ex[j] * sinBCD * cosABCD + 
	                    ey[j] * sinBCD * sinABCD );
	}
     }
     errcod = ip_string("ZMAT",&buffer,2,i,0);
     if (strcmp(buffer,"X")) {
       atom_num(buffer, &Z);
       free(buffer);
       nuclear_charges[atomcount] = Z;
       element[atomcount] = elem_name[(int)Z];
       geometry[atomcount][0] = full_geom[i][0]*conv_factor;
       geometry[atomcount][1] = full_geom[i][1]*conv_factor;
       geometry[atomcount][2] = full_geom[i][2]*conv_factor;
       atomcount++;
     }
  }

  free_matrix(full_geom,num_entries);
  return;
}


