#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include <ip_libv1.h>
#include "input.h"
#include "global.h"
#include "defines.h"

void read_cart()
{
  int i, j, errcod;
  double Z = 0.0;
  double tmp = 0.0;
  char *atom_label;

  num_atoms = 0;
  ip_count("GEOMETRY",&num_atoms,0);
  if (num_atoms == 0)
    punt("GEOMETRY is empty!");
  else if (num_atoms > MAXATOM)
    punt("There are more atoms than allowed!");

  /*-----------------------
    Allocate global arrays
   -----------------------*/
  geometry = block_matrix(num_atoms,3);
  element = (char **) malloc(sizeof(char *)*num_atoms);
  nuclear_charges = init_array(num_atoms);

  for(i=0;i<num_atoms;i++){
    errcod = ip_string("GEOMETRY",&atom_label,2,i,0);
    if (errcod != IPE_OK)
      punt("Problem with the GEOMETRY array.");
    atom_num(atom_label, &Z);
    free(atom_label);
    nuclear_charges[i] = Z;
    element[i] = elem_name[(int)Z];
    for(j=0; j<3;j++){
      errcod = ip_data("GEOMETRY","%lf", &tmp,2,i,j+1);
      if (errcod != IPE_OK)
	punt("Problem with the GEOMETRY array.");
      else
	geometry[i][j] = tmp*conv_factor;
    }
  }

  return;
}


