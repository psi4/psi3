/*####################################################################

  read_file11.cc

  this function reads in the cartesian coordinates and gradients
  and masses  

  modified code from OPTKING, written by Rollin King
                                                     J. Kenny 7-22-00
  ###################################################################*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern "C" {
  #include <string.h>
  #include <libciomr.h>
  #include <ip_libv1.h>
  #include <physconst.h>
  #include <masses.h>
}


#define EXTERN
#include "opt.h"

void read_file11(coord_base* coord) {

  int i, natom, count = 1, continue_flag = 1;
  char label[MAX_LINELENGTH], line1[MAX_LINELENGTH], *tmp_ptr;
  FILE *fp_11;
  double an,x,y,z,energy;
  
  if ((fp_11 = fopen("file11.dat","r")) == NULL) {
     punt("Could not open file11.dat");
   }

  tmp_ptr = fgets(label, MAX_LINELENGTH, fp_11);

  if (tmp_ptr == NULL) {
     punt("Touble reading first line of file11.dat");
  }

  fgets(line1, MAX_LINELENGTH, fp_11);
  if (sscanf(line1, "%d %lf", &natom, &energy) != 2) {
     punt("Trouble reading natoms and energy from file11.dat");
  }

  if(natom!=num_atoms)
      punt("Number of atoms differs in file11 and file30");

  rewind(fp_11);
  
  while ( fgets(label, MAX_LINELENGTH, fp_11) != NULL ) {

      fgets(line1, MAX_LINELENGTH, fp_11);
      sscanf(line1, "%d %lf", &natom, &energy);
           
      /*read in one chunk at a time*/ 
      for (i=0; i<num_atoms ; i++) {
	  if(fscanf(fp_11, "%lf %lf %lf %lf", &an, &x, &y, &z) != 4) {
	      punt("Trouble reading cartesian coordinates from file11.dat");
	    }
	  (*coord).set_mass(i,an2masses[(int) an]);
          fprintf(outfile,"\nmass: %lf",an2masses[(int) an]);
	  (*coord).set_cart(3*i,x );
	  (*coord).set_cart(3*i+1,y );
	  (*coord).set_cart(3*i+2,z );
	}
  
      for (i=0; i<num_atoms ; i++) {
	  if(fscanf(fp_11, "%lf %lf %lf", &x, &y, &z) != 3) {
	      punt("Trouble reading gradients from file11.dat");
	    }
	  (*coord).set_c_grad(3*i,x);
	  (*coord).set_c_grad(3*i+1,y);
	  (*coord).set_c_grad(3*i+2,z);
	}

       fgets(line1, MAX_LINELENGTH, fp_11);

    }

   fclose(fp_11);
   return;
}






