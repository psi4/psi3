#include<math.h>
#include<stdio.h>
#include<string.h>
#include<memory.h>
#include<stdlib.h>
#include<ip_libv1.h>
#include<libciomr.h>
#include<psio.h>
#include<libint.h>
#include<pthread.h>

#include"defines.h"
#define EXTERN
#include"global.h"

#include"read_scf_opdm.h"
#include"shell_block_matrix.h"
#include"hf_fock.h"

pthread_mutex_t fock_mutex;            /* Lock on the global AO matrix */

void fock()
{
  int dum;
  int n, num;
  int i, j, k, l;

  int nstri;
  double temp;
  double **tmpmat1;
  double *Gtri, *Gtri_o;               /* Total and open-shell G matrices and lower triagonal form
					  in SO basis */

  /*-----------------------------
    Read in the HF/DFT densities
   -----------------------------*/
  read_scf_opdm();

  /*-------------------------------------------
    Compute HF contribution to the Fock matrix
   -------------------------------------------*/
  hf_fock();

  /*--------------------------------------------------
    Compute exch+corr contribution to the Fock matrix
   --------------------------------------------------*/
/*  xc_fock();*/
  
  /*--------------------
    Print out G for now
   --------------------*/
/*  fprintf(outfile,"  Closed-shell Fock matrix in AO basis:\n");
  print_mat(G,BasisSet.num_ao,BasisSet.num_ao,outfile);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    fprintf(outfile,"  Open-shell Fock matrix in AO basis:\n");
    print_mat(Go,BasisSet.num_ao,BasisSet.num_ao,outfile);
  }*/

  /*----------------------
    Transform to SO basis
   ----------------------*/
  if (Symmetry.nirreps > 1 || BasisSet.puream) {
    tmpmat1 = block_matrix(Symmetry.num_so,BasisSet.num_ao);
    mmult(Symmetry.usotao,0,G,0,tmpmat1,0,Symmetry.num_so,BasisSet.num_ao,BasisSet.num_ao,0);
    mmult(tmpmat1,0,Symmetry.usotao,1,G,0,Symmetry.num_so,BasisSet.num_ao,Symmetry.num_so,0);
    if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
      mmult(Symmetry.usotao,0,Go,0,tmpmat1,0,Symmetry.num_so,BasisSet.num_ao,BasisSet.num_ao,0);
      mmult(tmpmat1,0,Symmetry.usotao,1,Go,0,Symmetry.num_so,BasisSet.num_ao,Symmetry.num_so,0);
    }
    free_block(tmpmat1);
  }
/*  fprintf(outfile,"  Closed-shell Fock matrix in SO basis:\n");
  print_mat(G,Symmetry.num_so,Symmetry.num_so,outfile);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf) {
    fprintf(outfile,"  Open-shell Fock matrix in SO basis:\n");
    print_mat(Go,Symmetry.num_so,Symmetry.num_so,outfile);
  }*/

  /*-------------------------
    Write G-matrices to disk
   -------------------------*/
  nstri = ioff[Symmetry.num_so];
  Gtri = init_array(nstri);
  sq_to_tri(G,Gtri,Symmetry.num_so);
  free_block(G);
  psio_open(IOUnits.itapDSCF, PSIO_OPEN_NEW);
  switch (UserOptions.reftype) {
  case rohf:
      Gtri_o = init_array(nstri);
      sq_to_tri(Go,Gtri_o,Symmetry.num_so);
      free_block(Go);
      psio_write_entry(IOUnits.itapDSCF, "Open-shell G-matrix", (char *) Gtri_o, sizeof(double)*nstri);
      free(Gtri_o);
  case rhf:
      psio_write_entry(IOUnits.itapDSCF, "Total G-matrix", (char *) Gtri, sizeof(double)*nstri);
      free(Gtri);
      break;

  case uhf:
      Gtri_o = init_array(nstri);
      sq_to_tri(Go,Gtri_o,Symmetry.num_so);
      free_block(Go);
      /*--- Form alpha and beta Fock matrices first and then write them out ---*/
      for(i=0;i<nstri;i++) {
	temp = Gtri[i] + Gtri_o[i];
	Gtri[i] = Gtri[i] - Gtri_o[i];
	Gtri_o[i] = temp;
      }
      psio_write_entry(IOUnits.itapDSCF, "Alpha G-matrix", (char *) Gtri, sizeof(double)*nstri);
      psio_write_entry(IOUnits.itapDSCF, "Beta G-matrix", (char *) Gtri_o, sizeof(double)*nstri);
      free(Gtri);
      free(Gtri_o);
      break;
  }
  psio_close(IOUnits.itapDSCF, 1);

  return;
}

