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


void read_scf_occ_evec(void)
{ 
    int i, j, k, l, jj, ij;
    int nstri;
    int num_so,num_ao,num_mo,ndocc;
    int bas_off;
    int shell_start,shell_end,shell_type;
    PSI_FPTR next;
    double **SO_cmat, **SO_cmato;
    double **Cocc_un;
    
    /* Read All information necessary to do DFT procedure with Eigenvector */
    
    psio_open(IOUnits.itapDSCF, PSIO_OPEN_OLD);
    
    psio_read_entry(IOUnits.itapDSCF, "Number of MOs", 
		    (char *) &(MOInfo.num_mo), sizeof(int));
    psio_read_entry(IOUnits.itapDSCF, "Number of DOCC", 
		    (char *) &(MOInfo.ndocc),sizeof(int));
    
    
/*-------------------
  Variables needed
  -------------------*/
    
    num_ao = BasisSet.num_ao;
    num_so = Symmetry.num_so;
    num_mo = MOInfo.num_mo;
    ndocc = MOInfo.ndocc;
    nstri = num_mo*ndocc;
    
    SO_cmat = block_matrix(num_so,ndocc);
      if (UserOptions.reftype == rohf || UserOptions.reftype == uhf)
	  SO_cmato = block_matrix(num_so,ndocc);
      
      if (UserOptions.reftype == uhf) {
	  psio_read_entry(IOUnits.itapDSCF, "Alpha SCF Eigenvector", 
		      (char *) &(SO_cmat[0][0]), sizeof(double)*nstri);
	  psio_read_entry(IOUnits.itapDSCF, "Beta SCF Eigenvector", 
			  (char *) &(SO_cmato[0][0]), sizeof(double)*nstri);
      }
      else {
	  psio_read_entry(IOUnits.itapDSCF, "Occupied SCF Eigenvector", 
			  (char *) &(SO_cmat[0][0]), sizeof(double)*nstri);
      }
      psio_close(IOUnits.itapDSCF, 1);
      
  
   /*----------------------
    transform to AO basis
    ----------------------*/
      /*fprintf(outfile,"\nUSOTAO matrix");
	print_mat(Symmetry.usotao,num_so,num_ao,outfile);
	fprintf(outfile,"\nSO Cmat");
	print_mat(SO_cmat,num_so,num_mo,outfile);
      */
      /*Cocc = block_matrix(num_ao,ndocc);
      mmult(Symmetry.usotao,1,SO_cmat,0,Cocc,0,num_ao,num_so,ndocc,0);
      free_block(SO_cmat);
      */
      /* ---------------------
	 Order according to 
	 angular momentum   
	 --------------------*/
      Cocc_un = block_matrix(num_ao,ndocc);
      mmult(Symmetry.usotao,1,SO_cmat,0,Cocc_un,0,num_ao,num_so,ndocc,0);
      free_block(SO_cmat);
      
      Cocc = (double **)malloc(sizeof(double *)*num_ao);
      
      k=0;
      for(i=0;i<BasisSet.num_shells;i++){
	  bas_off = BasisSet.am2shell[i];
	  shell_type = BasisSet.shells[bas_off].am;
	  shell_start = BasisSet.shells[bas_off].fao-1;
	  shell_end = shell_start+ioff[shell_type]; 
	  for(j=shell_start;j<shell_end;j++){
	      Cocc[k]=Cocc_un[j];
	      k++;
	  }
      }
      
      /*--------------------------
	Remove after done testing
	--------------------------*/   
      
      /*fprintf(outfile,"\nCocc");
	print_mat(Cocc,num_ao,MOInfo.ndocc,outfile);*/
      
      return;
}


