#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr.h>
#include <file30.h>
#include <psio.h>
#include <libint.h>

#include "defines.h"
#define EXTERN
#include "global.h"

#define OFFDIAG_DENS_FACTOR 0.5

/*-----------------------------------------------------
  This function reads in difference (diff_flag = 1) or
  total (diff_flag = 0) density matrices

  Only the HF Fock part needs the difference density
  matrix, hence when diff_flag = 0 and reftype = uhf
  I form alpha and beta densities right away.
 -----------------------------------------------------*/
void read_scf_evec(void)
{ 
  int i, j, ij;
  int nstri;
  int num_so,num_ao,num_mo;
  PSI_FPTR next;
  double **cmat, **cmato;
  double **sq_cmat, **sq_cmato;
  double **tmp_mat;
  double pfac;
  double temp, tempb;
  
  num_ao = BasisSet.num_ao;
  num_mo = MO_info.num_mo;
  num_so = Symmetry.num_so;
  nstri = num_mo*num_ao;
  
  psio_open(IOUnits.itapDSCF, PSIO_OPEN_OLD);
  cmat = block_matrix(num_mo,num_so);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf)
    cmato = block_matrix(num_mo,num_so);

  if (UserOptions.reftype == uhf) {
      psio_read_entry(IOUnits.itapDSCF, "Alpha SCF Eigenvector", (char *) &cmat[0][0], sizeof(double)*nstri);
      psio_read_entry(IOUnits.itapDSCF, "Beta SCF Eigenvector", (char *) &cmato[0][0], sizeof(double)*nstri);
  }
  else {
      psio_read_entry(IOUnits.itapDSCF, "SCF Eigenvector", (char *) &(cmat[0][0]), sizeof(double)*nstri);
  }
  psio_close(IOUnits.itapDSCF, 1);

  
  /*--------------------------
    Remove after done testing
    --------------------------*/
  fprintf(outfile,"SCF Eigenvector in SO basis :\n");
  print_mat(cmat,num_mo,num_so,outfile);
  fprintf(outfile,"\n\n");
  /*if (UserOptions.reftype == rohf) {
    fprintf(outfile,"  Total open-shell density matrix in SO basis :\n");
    print_array(denso,Symmetry.num_so,outfile);
    fprintf(outfile,"\n\n");
    }*/
  
   /*----------------------
    transform to AO basis
    ----------------------*/
  sq_cmat = init_matrix(num_ao,num_ao);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf)
      sq_cmato = init_matrix(num_ao,num_ao);

  tmp_mat = init_matrix(num_so,num_ao);
  mmult(cmat,0,Symmetry.usotao,0,tmp_mat,0,num_so,num_so,num_ao,0);
  mmult(Symmetry.usotao,1,tmp_mat,0,sq_cmat,0,num_ao,num_so,num_ao,0);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
      mmult(sq_cmato,0,Symmetry.usotao,0,tmp_mat,0,num_so,num_so,num_ao,0);
      mmult(Symmetry.usotao,1,tmp_mat,0,sq_cmato,0,BasisSet.num_ao,Symmetry.num_so,BasisSet.num_ao,0);
  }
  free_matrix(tmp_mat,Symmetry.num_so);

  /* -----------------------
     Put in the proper matrices
     ----------------------*/
  switch(UserOptions.reftype){
      
  case rhf:
      
      Cmat = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
      
      for(i=0;i<BasisSet.num_ao;i++)
	  for(j=0;j<BasisSet.num_ao;j++)
	      Cmat[i][j] = sq_cmat[i][j];
      free_matrix(sq_cmat,BasisSet.num_ao);
      break;
      
  case uhf:
      Cmata = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
      Cmatb = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
      
      for(i=0;i<BasisSet.num_ao;i++){
	  for(j=0;j<BasisSet.num_ao;i++){
	      Cmata[i][j] = sq_cmat[i][j]+sq_cmato[i][j];
	      Cmatb[i][j] = sq_cmat[i][j]-sq_cmato[i][j];
	  }
      }
      free_matrix(sq_cmat,BasisSet.num_ao);
      free_matrix(sq_cmato,BasisSet.num_ao);
      break;
  }
  
      /*--------------------------
	Remove after done testing
	--------------------------*/
  /*fprintf(outfile,"  Total density matrix in AO basis :\n");
    print_mat(Dens,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"\n\n");
    if (UserOptions.reftype == uhf) {
    fprintf(outfile,"  Total alpha-shell density matrix in AO basis :\n");
    print_mat(,BasisSet.num_ao,BasisSet.num_ao,outfile);
    fprintf(outfile,"\n\n");
    }*/
  
  return;
}


