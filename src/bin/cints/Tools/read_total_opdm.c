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


void read_total_opdm(void)
{ 
  int i, j, ij;
  int nstri = ioff[Symmetry.num_so];
  PSI_FPTR next;
  double *dens, *denso;
  double **sq_dens, **sq_denso;
  double **tmp_mat;
  double pfac;
  double temp, tempb;

  psio_open(IOUnits.itapDSCF, PSIO_OPEN_OLD);
  psio_read_entry(IOUnits.itapDSCF, "Integrals cutoff", (char *) &(UserOptions.cutoff), sizeof(double));
  dens = init_array(nstri);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf)
    denso = init_array(nstri);

  if (UserOptions.reftype == uhf) {
      denso = init_array(nstri);
      psio_read_entry(IOUnits.itapDSCF, "Total Alpha Density", (char *) dens, sizeof(double)*nstri);
      psio_read_entry(IOUnits.itapDSCF, "Total Beta Density", (char *) denso, sizeof(double)*nstri);
      /*--- Do NOT form total and open-shell (spin) density matrices ---*/
  }
  else {
      psio_read_entry(IOUnits.itapDSCF, "Total Density", (char *) dens, sizeof(double)*nstri);
      if (UserOptions.reftype == rohf) {
	  denso = init_array(nstri);
	  psio_read_entry(IOUnits.itapDSCF, "Total Open-shell Density", (char *) denso, sizeof(double)*nstri);
      }
  }
  psio_close(IOUnits.itapDSCF, 1);

  
  /*--------------------------
    Remove after done testing
   --------------------------*/
/*  fprintf(outfile,"  Total density matrix in SO basis :\n");
  print_array(dens,Symmetry.num_so,outfile);
  fprintf(outfile,"\n\n");
  if (UserOptions.reftype == rohf) {
    fprintf(outfile,"  Total open-shell density matrix in SO basis :\n");
    print_array(denso,Symmetry.num_so,outfile);
    fprintf(outfile,"\n\n");
  }*/

  
  /*------------------------
    convert to square forms
   ------------------------*/
  sq_dens = init_matrix(BasisSet.num_ao,BasisSet.num_ao);
  ij = 0;
  for(i=0;i<Symmetry.num_so;i++) {
      for(j=0;j<i;j++) {
	  sq_dens[i][j] = sq_dens[j][i] =  OFFDIAG_DENS_FACTOR*dens[ij];
	  ij++;
      }
      sq_dens[i][i] = dens[ij];
      ij++;
  }
  free(dens);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
      sq_denso = init_matrix(BasisSet.num_ao,BasisSet.num_ao);
      ij = 0;
      for(i=0;i<Symmetry.num_so;i++) {
	  for(j=0;j<i;j++) {
	      sq_denso[i][j] = sq_denso[j][i] = OFFDIAG_DENS_FACTOR*denso[ij];
	      ij++;
	  }
	  sq_denso[i][i] = denso[ij];
	  ij++;
      }
      free(denso);
  }
  
  /*----------------------
    transform to AO basis
    ----------------------*/
  tmp_mat = init_matrix(Symmetry.num_so,BasisSet.num_ao);
  mmult(sq_dens,0,Symmetry.usotao,0,tmp_mat,0,Symmetry.num_so,Symmetry.num_so,BasisSet.num_ao,0);
  mmult(Symmetry.usotao,1,tmp_mat,0,sq_dens,0,BasisSet.num_ao,Symmetry.num_so,BasisSet.num_ao,0);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
      mmult(sq_denso,0,Symmetry.usotao,0,tmp_mat,0,Symmetry.num_so,Symmetry.num_so,BasisSet.num_ao,0);
      mmult(Symmetry.usotao,1,tmp_mat,0,sq_denso,0,BasisSet.num_ao,Symmetry.num_so,BasisSet.num_ao,0);
  }
  free_matrix(tmp_mat,Symmetry.num_so);

  /* -----------------------
     Put in the proper matrices
     ----------------------*/
  switch(UserOptions.reftype){
      
  case rhf:
      
      /*Dens = init_matrix(BasisSet.num_ao,BasisSet.num_ao);*/
      Dens = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
      
      for(i=0;i<BasisSet.num_ao;i++)
	  for(j=0;j<BasisSet.num_ao;j++)
	      Dens[i][j] = sq_dens[i][j];
      free_matrix(sq_dens,BasisSet.num_ao);
      break;
      
  case uhf:
      Densa = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
      Densb = block_matrix(BasisSet.num_ao,BasisSet.num_ao);
      
      for(i=0;i<BasisSet.num_ao;i++){
	  for(j=0;j<BasisSet.num_ao;i++){
	      Densa[i][j] = sq_dens[i][j]+sq_denso[i][j];
	      Densb[i][j] = sq_dens[i][j]-sq_denso[i][j];
	  }
      }
      free_matrix(sq_dens,BasisSet.num_ao);
      free_matrix(sq_denso,BasisSet.num_ao);
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


