#include<stdio.h>
#include<stdlib.h>
#include<memory.h>
#include<libciomr.h>
#include<libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"


/*-------------------------------------
  Swap ket and bra of a 4-index buffer
 -------------------------------------*/
double *ijkl_to_klij(double *ijkl_buf, int nij, int nkl)
{
  int ij,kl;
  int ijkl, klij;
  double temp;
  /*--- local buffer ---*/
  static double *klij_buf;
  static int buf_size = 10000;

  if (klij_buf == NULL) {
    klij_buf = (double *) malloc(sizeof(double)*buf_size);
  }
  else if ((nij*nkl) > buf_size) {
    buf_size = (nij*nkl);
    free(klij_buf);
    klij_buf = (double *) malloc(sizeof(double)*buf_size);
  }

  ijkl = 0;
  for(ij=0;ij<nij;ij++) {
    klij = ij;
    for(kl=0;kl<nkl;kl++,ijkl++,klij+=nij)
	klij_buf[klij] = ijkl_buf[ijkl];
  }

  return klij_buf;
}
	  
