#include<stdio.h>
#include<libciomr/libciomr.h>
#include<libfile30/file30.h>
#include<stdlib.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"prints.h"

void read_gen_opdm()
{ 
  int natri = ioff[BasisSet.num_ao];
  PSI_FPTR next;
  double *dens, *lagr;

  dens = init_array(natri);
  lagr = init_array(natri);

  rfile(IOUnits.itapD);
  wreadw(IOUnits.itapD, (char *) dens, sizeof(double)*natri,    0, &next);
  wreadw(IOUnits.itapD, (char *) lagr, sizeof(double)*natri, next, &next);
  rclose(IOUnits.itapD,3);

  /* convert to square forms */
  Dens = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  Lagr = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
  tri_to_sq(dens,Dens,BasisSet.num_ao);
  tri_to_sq(lagr,Lagr,BasisSet.num_ao);
  free(dens);
  free(lagr);

  print_opdm();

  return;
}

