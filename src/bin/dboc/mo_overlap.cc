#include <stdio.h>
#include <stdlib.h>
extern "C" {
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <psifiles.h>
}
#include "float.h"
#include "linalg.h"
#include "moinfo.h"
#include "mo_overlap.h"

extern MOInfo_t MOInfo;
extern FILE *outfile;

FLOAT **eval_S_alpha()
{
  int num_mo = MOInfo.num_mo;
  int num_so = MOInfo.num_so;

  chkpt_init(PSIO_OPEN_OLD);
  double **rhf_evec_p = chkpt_rd_scf();
  chkpt_close();
    
  psio_open(PSIF_OLD_CHKPT, PSIO_OPEN_OLD);
  double **rhf_evec_m = block_matrix(num_so,num_mo);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:MOs alpha", (char *)&(rhf_evec_m[0][0]),
		  num_so*num_mo*sizeof(double));
  FLOAT **Spm_FLOAT = create_matrix(num_so,num_so);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:New-Old basis overlap", (char *)&(Spm_FLOAT[0][0]),
		  num_so*num_so*sizeof(FLOAT));
  psio_close(PSIF_OLD_CHKPT, 1);


  //
  // Convert matrices of doubles into matrices of FLOAT's
  //
  FLOAT** rhf_evec_m_FLOAT = convert_matrix(rhf_evec_m, num_so, num_mo, 0);
  FLOAT** rhf_evec_p_FLOAT_transp = convert_matrix(rhf_evec_p, num_so, num_mo, 1);
  //  FLOAT** Spm_FLOAT = convert_matrix(Spm, num_so, num_so, 0);

  FLOAT** tmpmat1 = create_matrix(num_mo,num_so);
  if (matrix_mult(rhf_evec_p_FLOAT_transp, num_mo, num_so, Spm_FLOAT, num_so, num_so, tmpmat1))
    abort();
  FLOAT** S = create_matrix(num_mo, num_mo);
  if (matrix_mult(tmpmat1, num_mo, num_so, rhf_evec_m_FLOAT, num_so, num_mo, S))
    abort();

  delete_matrix(tmpmat1);
  delete_matrix(Spm_FLOAT);
  delete_matrix(rhf_evec_p_FLOAT_transp);
  delete_matrix(rhf_evec_m_FLOAT);
  
  /*  double **tmpmat2 = block_matrix(num_mo,num_so);
  mmult(rhf_evec_p,1,Spm,0,tmpmat2,0,
	num_mo,num_so,num_so,0);
  double **S2 = block_matrix(num_mo,num_mo);
  mmult(tmpmat2,0,rhf_evec_m,0,S2,0,
	num_mo,num_so,num_mo,0);
	free_block(tmpmat2);*/
  free_block(rhf_evec_p);
  free_block(rhf_evec_m);
  //  free_block(Spm);

  return S;
}

FLOAT **eval_S_beta()
{
  int num_mo = MOInfo.num_mo;
  int num_so = MOInfo.num_so;
  
  chkpt_init(PSIO_OPEN_OLD);
  double **rhf_evec_p = chkpt_rd_beta_scf();
  chkpt_close();
    
  psio_open(PSIF_OLD_CHKPT, PSIO_OPEN_OLD);
  double **rhf_evec_m = block_matrix(num_so,num_mo);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:MOs beta", (char *)&(rhf_evec_m[0][0]),
		  num_so*num_mo*sizeof(double));
  double **Spm = block_matrix(num_so,num_so);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:New-Old basis overlap", (char *)&(Spm[0][0]),
		  num_so*num_so*sizeof(double));
  psio_close(PSIF_OLD_CHKPT, 1);
  
  //
  // Convert matrices of doubles into matrices of FLOAT's
  //
  FLOAT** rhf_evec_m_FLOAT = convert_matrix(rhf_evec_m, num_so, num_mo, 0);
  FLOAT** rhf_evec_p_FLOAT_transp = convert_matrix(rhf_evec_p, num_so, num_mo, 1);
  FLOAT** Spm_FLOAT = convert_matrix(Spm, num_so, num_so, 0);

  FLOAT** tmpmat1 = create_matrix(num_mo,num_so);
  if (matrix_mult(rhf_evec_p_FLOAT_transp, num_mo, num_so, Spm_FLOAT, num_so, num_so, tmpmat1))
    abort();
  FLOAT** S = create_matrix(num_mo, num_mo);
  if (matrix_mult(tmpmat1, num_mo, num_so, rhf_evec_m_FLOAT, num_so, num_mo, S))
    abort();

  delete_matrix(tmpmat1);
  delete_matrix(Spm_FLOAT);
  delete_matrix(rhf_evec_p_FLOAT_transp);
  delete_matrix(rhf_evec_m_FLOAT);

  /*  double **tmpmat1 = block_matrix(num_mo,num_so);
  mmult(rhf_evec_p,1,Spm,0,tmpmat1,0,
	num_mo,num_so,num_so,0);
  double **S = block_matrix(num_mo,num_mo);
  mmult(tmpmat1,0,rhf_evec_m,0,S,0,
	num_mo,num_so,num_mo,0);
	free_block(tmpmat1);*/
  free_block(rhf_evec_p);
  free_block(rhf_evec_m);
  free_block(Spm);

  return S;
}

