/*####################################################################

  read_file11.cc

  this function reads the cartesian coordinates and gradients
  and atomic numbers into global arrays from file11.dat 

  modified code from OPTKING, written by Rollin King
                                                     J. Kenny 7-22-00
  ###################################################################*/

#include <stdio.h>
#include <stdlib.h>

extern "C" {
  #include <string.h>
  #include <libciomr.h>
  #include <ip_libv1.h>
  #include <physconst.h>
}


#define EXTERN
#include "opt.h"

void read_file11() {

  int i;
  char label[MAX_LINELENGTH], line1[MAX_LINELENGTH];
  FILE *fp_11;
  double an,x,y,z,energy;
  
  if ((fp_11 = fopen("file11.dat","r")) == NULL) {
     punt("Could not open file11.dat");
   }

  if ((fgets(label, MAX_LINELENGTH, fp_11)) == NULL) {
     punt("Touble reading first line of file11.dat");
  }

  fgets(line1, MAX_LINELENGTH, fp_11);
  if (sscanf(line1, "%d %lf", &num_atoms, &energy) != 2) {
     punt("Trouble reading natoms and energy from file11.dat");
  }

  cart_geom = init_matrix(num_atoms,3);
  cart_grad = init_matrix(num_atoms,3);
  atomic_nums = init_array(num_atoms);
 
  /*read in one chunk at a time*/ 
  for (i=0; i<num_atoms ; i++) {
      if(fscanf(fp_11, "%lf %lf %lf %lf", &an, &x, &y, &z) != 4) {
	  punt("Trouble reading cartesian coordinates from file11.dat");
	}
      atomic_nums[0]=an;
      cart_geom[i][0]=x * _bohr2angstroms;
      cart_geom[i][1]=y * _bohr2angstroms;
      cart_geom[i][2]=z * _bohr2angstroms;
      }
  
  for (i=0; i<num_atoms ; i++) {
      if(fscanf(fp_11, "%lf %lf %lf", &x, &y, &z) != 3) {
	  punt("Trouble reading gradients from file11.dat");
	}
     cart_grad[i][0] = x;
     cart_grad[i][1] = y;
     cart_grad[i][2] = z;
  } 

   fclose(fp_11);
   return;
}






