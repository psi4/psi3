#define EXTERN
#include <stdio.h>
#include <stdlib.h>
#include <libciomr.h>
#include <libipv1/ip_lib.h>
#include <file30.h>
#include "input.h"
#include "global.h"
#include "defines.h"

void read_chkpt_geom()
{
  int i, j, errcod;
  double Z = 0.0;
  double tmp = 0.0;
  char *atom_label;

  file30_init();
  num_atoms = file30_rd_natom();
  if (num_atoms == 0)
    punt("GEOMETRY in the checkpoint file is empty!");
  else if (num_atoms > MAXATOM)
    punt("There are more atoms than allowed!");

  /*-----------------------
    Allocate global arrays
   -----------------------*/
  geometry = file30_rd_geom();
  disp_num = file30_rd_disp();
  nuclear_charges = file30_rd_zvals();
  element = (char **) malloc(sizeof(char *)*num_atoms);

  for(i=0;i<num_atoms;i++)
    element[i] = elem_name[(int)nuclear_charges[i]];

  file30_close();
  
  return;
}


