#include <stdio.h>
#include <stdlib.h>
extern "C" {
#include <libciomr/libciomr.h>
#include <libfile30/file30.h>
#include <libpsio/psio.h>
#include <psifiles.h>
}
#include "moinfo.h"

extern MOInfo_t MOInfo;
extern FILE *outfile;

double **eval_S_alpha()
{
  int num_mo = MOInfo.num_mo;
  int num_so = MOInfo.num_so;
  
  file30_init();
  double **rhf_evec_p = file30_rd_alpha_scf();
  file30_close();
    
  psio_open(PSIF_OLD_CHKPT, PSIO_OPEN_OLD);
  double **rhf_evec_m = block_matrix(num_so,num_mo);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:MOs alpha", (char *)&(rhf_evec_m[0][0]),
		  num_so*num_mo*sizeof(double));
  double **Spm = block_matrix(num_so,num_so);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:New-Old basis overlap", (char *)&(Spm[0][0]),
		  num_so*num_so*sizeof(double));
  psio_close(PSIF_OLD_CHKPT, 1);
    
  double **tmpmat1 = block_matrix(num_mo,num_so);
  mmult(rhf_evec_p,1,Spm,0,tmpmat1,0,
	num_mo,num_so,num_so,0);
  double **S = block_matrix(num_mo,num_mo);
  mmult(tmpmat1,0,rhf_evec_m,0,S,0,
	num_mo,num_so,num_mo,0);
  free_block(tmpmat1);
  free_block(rhf_evec_p);
  free_block(rhf_evec_m);
  free_block(Spm);

  return S;
}

double **eval_S_beta()
{
  int num_mo = MOInfo.num_mo;
  int num_so = MOInfo.num_so;
  
  file30_init();
  double **rhf_evec_p = file30_rd_beta_scf();
  file30_close();
    
  psio_open(PSIF_OLD_CHKPT, PSIO_OPEN_OLD);
  double **rhf_evec_m = block_matrix(num_so,num_mo);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:MOs beta", (char *)&(rhf_evec_m[0][0]),
		  num_so*num_mo*sizeof(double));
  double **Spm = block_matrix(num_so,num_so);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:New-Old basis overlap", (char *)&(Spm[0][0]),
		  num_so*num_so*sizeof(double));
  psio_close(PSIF_OLD_CHKPT, 1);
    
  double **tmpmat1 = block_matrix(num_mo,num_so);
  mmult(rhf_evec_p,1,Spm,0,tmpmat1,0,
	num_mo,num_so,num_so,0);
  double **S = block_matrix(num_mo,num_mo);
  mmult(tmpmat1,0,rhf_evec_m,0,S,0,
	num_mo,num_so,num_mo,0);
  free_block(tmpmat1);
  free_block(rhf_evec_p);
  free_block(rhf_evec_m);
  free_block(Spm);

  return S;
}

