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


void read_scf_evec(void)
{ 
  int i, j, ij;
  int nstri;
  int num_so,num_ao,num_mo;
  PSI_FPTR next;
  double **cmat, **cmato;
  
  
  
  psio_open(IOUnits.itapDSCF, PSIO_OPEN_OLD);
  psio_read_entry(IOUnits.itapDSCF, "Number of MOs", (char *) &(MOInfo.num_mo), sizeof(int));
  psio_read_entry(IOUnits.itapDSCF, "Number of Doubly Occ", (char *) &(MOInfo.ndocc),sizeof(int));
  
  num_ao = BasisSet.num_ao;
  num_so = Symmetry.num_so;
  num_mo = MOInfo.num_mo;
  nstri = num_mo*num_ao;
  
  cmat = block_matrix(num_so,num_mo);
  if (UserOptions.reftype == rohf || UserOptions.reftype == uhf)
    cmato = block_matrix(num_so,num_mo);

  if (UserOptions.reftype == uhf) {
      psio_read_entry(IOUnits.itapDSCF, "Alpha SCF Eigenvector", (char *) &(cmat[0][0]), sizeof(double)*nstri);
      psio_read_entry(IOUnits.itapDSCF, "Beta SCF Eigenvector", (char *) &(cmato[0][0]), sizeof(double)*nstri);
  }
  else {
      psio_read_entry(IOUnits.itapDSCF, "SCF Eigenvector", (char *) &(cmat[0][0]), sizeof(double)*nstri);
  }
  psio_close(IOUnits.itapDSCF, 1);

  
   /*----------------------
    transform to AO basis
    ----------------------*/
  Cmat = block_matrix(num_ao,num_mo);
  mmult(cmat,0,Symmetry.usotao,0,Cmat,0,num_so,num_mo,num_ao,0);
  
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


