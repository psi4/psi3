#include<stdio.h>
#include<stdlib.h>
#include<memory.h>
#include<libciomr/libciomr.h>
#include<libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"


/*-------------------------------------
  Swap ket and bra of a 4-index buffer
 -------------------------------------*/
void ijkl_to_klij(double *ijkl_buf, double *klij_buf, int nij, int nkl)
{
  int ij,kl;
  int ijkl, klij;
  double temp;

  ijkl = 0;
  for(ij=0;ij<nij;ij++) {
    klij = ij;
    for(kl=0;kl<nkl;kl++,ijkl++,klij+=nij)
	klij_buf[klij] = ijkl_buf[ijkl];
  }

  return;
}
	  
