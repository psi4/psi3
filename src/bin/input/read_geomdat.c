#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include "input.h"
#include "global.h"
#include "defines.h"

void read_geomdat()
{
  FILE *geomdat;
  int i, j, errcod;
  double Z = 0.0;
  double tmp = 0.0;
  char entry_name[20];

  geomdat = fopen("./geom.dat", "r");
  if (geomdat != NULL) {
    ip_append(geomdat, outfile);
    fclose(geomdat);
  }

  num_atoms = 0;
  sprintf(entry_name,"GEOMETRY%d",geomdat_entry);
  ip_count(entry_name,&num_atoms,0);
  if (num_atoms == 0)
    punt("The entry in geom.dat is empty or missing!");
  else if (num_atoms > MAXATOM)
    punt("There are more atoms than allowed!");

  /*-----------------------
    Allocate global arrays
   -----------------------*/
  geometry = block_matrix(num_atoms,3);
  element = (char **) malloc(sizeof(char *)*num_atoms);
  full_element = (char **) malloc(sizeof(char *)*num_atoms);
  elemsymb_charges = init_array(num_atoms);

  for(i=0;i<num_atoms;i++){
    errcod = ip_data(entry_name,"%lf",&Z,2,i,0);
    if (errcod != IPE_OK)
      punt("Problem with the geom.dat entry.");
    elemsymb_charges[i] = Z;
    element[i] = elem_name[(int)Z];
    full_element[i] = elem_name[(int)Z];
    for(j=0; j<3;j++){
      errcod = ip_data(entry_name,"%lf", &tmp,2,i,j+1);
      if (errcod != IPE_OK)
	punt("Problem with the geom.dat entry.");
      else
	geometry[i][j] = tmp;
    }
  }

  read_charges();

  return;
}


