#include <stdio.h>
#include <stdlib.h>
extern "C" {
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <psifiles.h>
}
#include "float.h"
#include "params.h"
#include "linalg.h"
#include "moinfo.h"
#include "mo_overlap.h"
#include <libbasis/basisset.h>
#include <libbasis/overlap.h>
#include <libbasis/rotation.h>

extern void done(const char * message);
extern MOInfo_t MOInfo;
extern FILE *outfile;
extern BasisSet *BasisDispM, *BasisDispP;
extern Params_t Params;

FLOAT **eval_S_alpha()
{
  int num_mo = MOInfo.num_mo;
  int num_so = MOInfo.num_so;

  chkpt_init(PSIO_OPEN_OLD);
  double **rhf_evec_p = chkpt_rd_scf();
  double **aotoso_p = chkpt_rd_usotao();
  double **rref_p = chkpt_rd_rref();
  int num_ao = chkpt_rd_nao();
  chkpt_close();
    
  psio_open(PSIF_OLD_CHKPT, PSIO_OPEN_OLD);

  double **rhf_evec_m = block_matrix(num_so,num_mo);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:MOs alpha", (char *)&(rhf_evec_m[0][0]),
		  num_so*num_mo*sizeof(double));

  double **aotoso_m = block_matrix(num_so, num_ao);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:AO->SO", (char *)aotoso_m[0],
		  num_so*num_ao*sizeof(double));

  double **rref_m = block_matrix(3,3);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:Transmat to reference frame", (char *)rref_m[0],
		  9*sizeof(double));

#if USE_INPUT_S
  FLOAT **Spm_FLOAT = create_matrix(num_so,num_so);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:New-Old basis overlap", (char *)&(Spm_FLOAT[0][0]),
		  num_so*num_so*sizeof(FLOAT));
#endif

  psio_close(PSIF_OLD_CHKPT, 1);


  //
  // Convert matrices of doubles into matrices of FLOAT's
  //
  FLOAT** rhf_evec_m_FLOAT = convert_matrix(rhf_evec_m, num_so, num_mo, 0);
  FLOAT** rhf_evec_p_FLOAT_transp = convert_matrix(rhf_evec_p, num_so, num_mo, 1);

  // Compute plus/minus overlap
#if !USE_INPUT_S
  OverlapEngine overlap(BasisDispP,BasisDispM);
  double** Spm_AO_AO = overlap.compute_full_matrix();

  // Rotate bases to the original coordinate systems (prior to reorientation into principal axis system)
  RotationOp Rop_m(BasisDispM);
  RotationOp Rop_p(BasisDispP);
  double** basisRref_m = Rop_m.full_rotation_mat(rref_m);
  double** basisRref_p = Rop_p.full_rotation_mat(rref_p);

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile, "  -Rotation matrix for - AO basis\n");
    print_mat(basisRref_m, num_ao, num_ao, outfile);
    
    fprintf(outfile, "  -Rotation matrix for + AO basis\n");
    print_mat(basisRref_p, num_ao, num_ao, outfile);

    double** tmp_p1 = block_matrix(num_ao, num_mo);
    double** tmp_p2 = block_matrix(num_ao, num_mo);
    double** tmp_m1 = block_matrix(num_ao, num_mo);
    double** tmp_m2 = block_matrix(num_ao, num_mo);
    mmult(aotoso_m,1,rhf_evec_m,0,tmp_m1,0,num_ao,num_so,num_mo,0);
    mmult(basisRref_m,0,tmp_m1,0,tmp_m2,0,num_ao,num_ao,num_mo,0);
    mmult(aotoso_p,1,rhf_evec_p,0,tmp_p1,0,num_ao,num_so,num_mo,0);
    mmult(basisRref_p,0,tmp_p1,0,tmp_p2,0,num_ao,num_ao,num_mo,0);
    
    fprintf(outfile, "  -Original - eigenvector\n");
    print_mat(tmp_m1, num_ao, num_mo, outfile);
    fprintf(outfile, "  -Rotated - eigenvector\n");
    print_mat(tmp_m2, num_ao, num_mo, outfile);
    
    fprintf(outfile, "  -Original + eigenvector\n");
    print_mat(tmp_p1, num_ao, num_mo, outfile);
    fprintf(outfile, "  -Rotated + eigenvector\n");
    print_mat(tmp_p2, num_ao, num_mo, outfile);
  }

  double** tmpmat = block_matrix(num_ao,num_ao);
  mmult(Spm_AO_AO,0,basisRref_m,0,tmpmat,0,num_ao,num_ao,num_ao,0);
  mmult(basisRref_p,1,tmpmat,0,Spm_AO_AO,0,num_ao,num_ao,num_ao,0);
  free_block(tmpmat);
  free_block(basisRref_p);
  free_block(basisRref_m);

  double** Spm_AO_SO = block_matrix(num_ao,num_so);
  mmult(Spm_AO_AO, 0, aotoso_m, 1, Spm_AO_SO, 0, num_ao, num_ao, num_so, 0);
  free_block(Spm_AO_AO);
  double** Spm_SO_SO = block_matrix(num_so,num_so);
  mmult(aotoso_p, 0, Spm_AO_SO, 0, Spm_SO_SO, 0, num_so, num_ao, num_so, 0);
  free_block(Spm_AO_SO);
  FLOAT** Spm_FLOAT = convert_matrix(Spm_SO_SO, num_so, num_so, 0);
  free_block(Spm_SO_SO);

#endif

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile, "  +/- overlap matrix (SO basis)\n");
    print_mat(Spm_FLOAT, num_so, num_so, outfile);
  }

  FLOAT** tmpmat1 = create_matrix(num_mo,num_so);
  if (matrix_mult(rhf_evec_p_FLOAT_transp, num_mo, num_so, Spm_FLOAT, num_so, num_so, tmpmat1))
    done("matrix_mult failed. Report the problem to the author.");
  FLOAT** S = create_matrix(num_mo, num_mo);
  if (matrix_mult(tmpmat1, num_mo, num_so, rhf_evec_m_FLOAT, num_so, num_mo, S))
    done("matrix_mult failed. Report the problem to the author.");

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile, "  +/- overlap matrix (MO basis)\n");
    print_mat(S, num_mo, num_mo, outfile);
  }

  delete_matrix(tmpmat1);
  delete_matrix(Spm_FLOAT);
  delete_matrix(rhf_evec_p_FLOAT_transp);
  delete_matrix(rhf_evec_m_FLOAT);
  
  free_block(rhf_evec_p);
  free_block(rhf_evec_m);

  return S;
}

FLOAT **eval_S_beta()
{
  int num_mo = MOInfo.num_mo;
  int num_so = MOInfo.num_so;
  
  chkpt_init(PSIO_OPEN_OLD);
  double **rhf_evec_p = chkpt_rd_beta_scf();
  double **aotoso_p = chkpt_rd_usotao();
  double **rref_p = chkpt_rd_rref();
  int num_ao = chkpt_rd_nao();
  chkpt_close();
    
  psio_open(PSIF_OLD_CHKPT, PSIO_OPEN_OLD);

  double **rhf_evec_m = block_matrix(num_so,num_mo);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:MOs beta", (char *)&(rhf_evec_m[0][0]),
		  num_so*num_mo*sizeof(double));

  double **aotoso_m = block_matrix(num_so, num_ao);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:AO->SO", (char *)aotoso_m[0],
		  num_so*num_ao*sizeof(double));

  double **rref_m = block_matrix(3,3);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:Transmat to reference frame", (char *)rref_m[0],
		  9*sizeof(double));

#if USE_INPUT_S
  FLOAT **Spm_FLOAT = create_matrix(num_so,num_so);
  psio_read_entry(PSIF_OLD_CHKPT, ":PrevCalc:New-Old basis overlap", (char *)&(Spm_FLOAT[0][0]),
		  num_so*num_so*sizeof(FLOAT));
#endif

  psio_close(PSIF_OLD_CHKPT, 1);
  
  //
  // Convert matrices of doubles into matrices of FLOAT's
  //
  FLOAT** rhf_evec_m_FLOAT = convert_matrix(rhf_evec_m, num_so, num_mo, 0);
  FLOAT** rhf_evec_p_FLOAT_transp = convert_matrix(rhf_evec_p, num_so, num_mo, 1);

  // Compute plus/minus overlap
#if !USE_INPUT_S
  OverlapEngine overlap(BasisDispP,BasisDispM);
  double** Spm_AO_AO = overlap.compute_full_matrix();

  // Rotate bases to the original coordinate systems (prior to reorientation into principal axis system)
  RotationOp Rop_m(BasisDispM);
  RotationOp Rop_p(BasisDispP);
  double** basisRref_m = Rop_m.full_rotation_mat(rref_m);
  double** basisRref_p = Rop_p.full_rotation_mat(rref_p);

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile, "  -Rotation matrix for - AO basis\n");
    print_mat(basisRref_m, num_ao, num_ao, outfile);
    
    fprintf(outfile, "  -Rotation matrix for + AO basis\n");
    print_mat(basisRref_p, num_ao, num_ao, outfile);

    double** tmp_p1 = block_matrix(num_ao, num_mo);
    double** tmp_p2 = block_matrix(num_ao, num_mo);
    double** tmp_m1 = block_matrix(num_ao, num_mo);
    double** tmp_m2 = block_matrix(num_ao, num_mo);
    mmult(aotoso_m,1,rhf_evec_m,0,tmp_m1,0,num_ao,num_so,num_mo,0);
    mmult(basisRref_m,0,tmp_m1,0,tmp_m2,0,num_ao,num_ao,num_mo,0);
    mmult(aotoso_p,1,rhf_evec_p,0,tmp_p1,0,num_ao,num_so,num_mo,0);
    mmult(basisRref_p,0,tmp_p1,0,tmp_p2,0,num_ao,num_ao,num_mo,0);
    
    fprintf(outfile, "  -Original - eigenvector\n");
    print_mat(tmp_m1, num_ao, num_mo, outfile);
    fprintf(outfile, "  -Rotated - eigenvector\n");
    print_mat(tmp_m2, num_ao, num_mo, outfile);
    
    fprintf(outfile, "  -Original + eigenvector\n");
    print_mat(tmp_p1, num_ao, num_mo, outfile);
    fprintf(outfile, "  -Rotated + eigenvector\n");
    print_mat(tmp_p2, num_ao, num_mo, outfile);
  }

  double** tmpmat = block_matrix(num_ao,num_ao);
  mmult(Spm_AO_AO,0,basisRref_m,0,tmpmat,0,num_ao,num_ao,num_ao,0);
  mmult(basisRref_p,1,tmpmat,0,Spm_AO_AO,0,num_ao,num_ao,num_ao,0);
  free_block(tmpmat);
  free_block(basisRref_p);
  free_block(basisRref_m);

  double** Spm_AO_SO = block_matrix(num_ao,num_so);
  mmult(Spm_AO_AO, 0, aotoso_m, 1, Spm_AO_SO, 0, num_ao, num_ao, num_so, 0);
  free_block(Spm_AO_AO);
  double** Spm_SO_SO = block_matrix(num_so,num_so);
  mmult(aotoso_p, 0, Spm_AO_SO, 0, Spm_SO_SO, 0, num_so, num_ao, num_so, 0);
  free_block(Spm_AO_SO);
  FLOAT** Spm_FLOAT = convert_matrix(Spm_SO_SO, num_so, num_so, 0);
  free_block(Spm_SO_SO);

#endif

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile, "  +/- overlap matrix (SO basis)\n");
    print_mat(Spm_FLOAT, num_so, num_so, outfile);
  }

  FLOAT** tmpmat1 = create_matrix(num_mo,num_so);
  if (matrix_mult(rhf_evec_p_FLOAT_transp, num_mo, num_so, Spm_FLOAT, num_so, num_so, tmpmat1))
    done("matrix_mult failed. Report the problem to the author.");
  FLOAT** S = create_matrix(num_mo, num_mo);
  if (matrix_mult(tmpmat1, num_mo, num_so, rhf_evec_m_FLOAT, num_so, num_mo, S))
    done("matrix_mult failed. Report the problem to the author.");

  if (Params.print_lvl > PrintLevels::print_contrib) {
    fprintf(outfile, "  +/- overlap matrix (MO basis)\n");
    print_mat(S, num_mo, num_mo, outfile);
  }

  delete_matrix(tmpmat1);
  delete_matrix(Spm_FLOAT);
  delete_matrix(rhf_evec_p_FLOAT_transp);
  delete_matrix(rhf_evec_m_FLOAT);

  free_block(rhf_evec_p);
  free_block(rhf_evec_m);

  return S;
}

