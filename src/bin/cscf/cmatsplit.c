/* ------------------------------------

   cmatsplit - takes cmat and splits
   the scf_info struct into two spin
   struct cmats.  Also intend to add
   a mixing here to force convergence
   to UHF solutions.

   -------------------------------------*/

#define EXTERN
#include "includes.h"
#include "common.h"

void cmatsplit(void){
    int i,j,k,l;
    int nn;
    double temp;
    
    for(i=0;i<num_ir;i++){
	nn=scf_info[i].num_so;
	for(j=0;j<nn;j++){
	    for(k=0;k<nn;k++){
		temp = scf_info[i].cmat[j][k];
		spin_info[0].scf_spin[i].cmat[j][k] = temp;
		spin_info[1].scf_spin[i].cmat[j][k] = temp;
	    }
	}
	/*fprintf(outfile,"\n\nCmat for spin = %d and irrep = %s\n"
		,0,scf_info[i].irrep_label);
	print_mat(spin_info[0].scf_spin[i].cmat,nn,nn,outfile);
	fprintf(outfile,"\n\nCmat for spin = %d and irrep = %s\n"
		,1,scf_info[i].irrep_label);
		print_mat(spin_info[1].scf_spin[i].cmat,nn,nn,outfile);*/
    }
}
		
    
