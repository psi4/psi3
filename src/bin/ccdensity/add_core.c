/*! \file add_core.c
    \ingroup (CCDENSITY)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libciomr/libciomr.h>
#include <libiwl/iwl.h>
#define EXTERN
#include "globals.h"

void add_core(struct iwlbuf *OutBuf)
{
  int p,q,m,n;
  int nmo, nfzv, nfzc;
  double value;

  nmo = moinfo.nmo;
  nfzv = moinfo.nfzv;
  nfzc = moinfo.nfzc;

  for(p=nfzc; p < (nmo - nfzv); p++) {
      for(q=nfzc; q < (nmo - nfzv); q++) {
	  value = moinfo.opdm[p][q];
	  for(m=0; m < nfzc; m++) {
	      
	      iwl_buf_wrt_val(OutBuf, p, q, m, m,value,0,outfile,0);
     	      iwl_buf_wrt_val(OutBuf, p, m, m, q,-0.5*value,0,outfile,0);
	      
	    }
	}
    }
}

