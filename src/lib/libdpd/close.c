#include <stdio.h>
#include <stdlib.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

int dpd_close(void)
{
  int h,i,j,k,cnt;

  for(i=0; i < dpd_num_subspaces; i++)  free(dpd_orboff[i]);
  free(dpd_orboff);

  for(i=0; i < dpd_num_subspaces; i++) {
      for(j=0; j < 5; j++) {
	  free_int_matrix(dpd_pairidx[5*i+j], dpd_numorbs[i]);
	  for(k=0; k < dpd_nirreps; k++) 
	      if(dpd_pairtot[5*i+j][k])
		free_int_matrix(dpd_pairorb[5*i+j][k],dpd_pairtot[5*i+j][k]);
	}
    }
  for(i=0,cnt=5*dpd_num_subspaces; i < dpd_num_subspaces; i++)
      for(j=i+1; j < dpd_num_subspaces; j++,cnt+=2) {
	  free_int_matrix(dpd_pairidx[cnt],dpd_numorbs[i]);
	  free_int_matrix(dpd_pairidx[cnt+1],dpd_numorbs[j]);
	  for(k=0; k < dpd_nirreps; k++) {
	      if(dpd_pairtot[cnt][k])
		free_int_matrix(dpd_pairorb[cnt][k],dpd_pairtot[cnt][k]);
	      if(dpd_pairtot[cnt+1][k])
		free_int_matrix(dpd_pairorb[cnt+1][k],dpd_pairtot[cnt+1][k]);
	    }
	}
  free(dpd_pairidx); free(dpd_pairorb);

  for(i=0; i < dpd_num_subspaces; i++) {
      free(dpd_oe_orbidx[i]);
      for(j=0; j < dpd_nirreps; j++) {
	  if(dpd_orbspi[i][j]) free(dpd_oe_orbs[i][j]);
	}
      free(dpd_oe_orbs[i]);
    }
  free(dpd_oe_orbidx); free(dpd_oe_orbs);

  free_int_matrix(dpd_pairtot, dpd_num_pairs);

  free(dpd_numorbs);

  for(i=0; i < dpd_num_pairs; i++)
      free(dpd_params[i]);
  free(dpd_params);
  for(i=0; i < dpd_num_subspaces; i++)
      free(dpd_oe_params[i]);
  free(dpd_oe_params);

  return 0;
}
