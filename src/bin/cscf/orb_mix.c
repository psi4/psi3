/*-----------------------------------------------------------

  orb_mix(); - A function to take the alpha and beta orbitals
  and destroy their alpha and beta symmetry to find unique 
  UHF solutions.

  The manner in which this will be done is to take the LUMO
  and make them mix them in the following way

  LUMO (alpha) = LUMO + HOMO
  LUMO (beta) = LUMO - HOMO

  Shawn Brown - (11/4/99)
  
  ---------------------------------------------------------*/

#define EXTERN
#include "includes.h"
#include "common.h"

void orb_mix(void){
    
    int i,j,k,l,m,n;
    int *lumo;
    int *homo;
    double lumo1,lumo2;
    double homo1,homo2;
    double homom,lumom;
    double mixed;
    int num_mo;
    
    fprintf(outfile,"\n***** Orbital Mixing being used to find unique UHF solution *****\n");

    /*First we must find out which orbital is actually the HOMO and LUMO */
    
    lumo = init_int_array(2);
    homo = init_int_array(2);
    
    for(n = 1; n < 2; n++){
	for(m = 1; m < num_ir; m++){
	    if(spin_info[n].scf_spin[m].noccup){
		homo1 = spin_info[n].scf_spin[m].fock_evals[spin_info[n].scf_spin[m].noccup-1];
		homo2 = spin_info[n].scf_spin[homo[n]].fock_evals[spin_info[n].scf_spin[m].noccup-1];
		lumo1 = spin_info[n].scf_spin[m].fock_evals[spin_info[n].scf_spin[m].noccup];
		lumo2 = spin_info[n].scf_spin[lumo[n]].fock_evals[spin_info[n].scf_spin[m].noccup];
		
		if(lumo1 < lumo2) lumo[n] == m;
		if(homo1 > homo2) homo[n] == m;
	    }
	}
	num_mo = scf_info[lumo[n]].num_mo;
	/*fprintf(outfile,"\n CMAT before mixing\n");
	  print_mat(spin_info[n].scf_spin[lumo[n]].cmat,num_mo,num_mo,outfile);*/
    }
    
    /* Now that we know which irrep is the lumo lets alter its orbital */
    
    
    for(n = 1; n < 2; n++){
	num_mo = scf_info[lumo[n]].num_mo;
	for(i = 0; i < num_mo;i++){
	    homom = spin_info[n].scf_spin[homo[n]].cmat[i][spin_info[n].scf_spin[homo[n]].noccup-1];
	    lumom = spin_info[n].scf_spin[lumo[n]].cmat[i][spin_info[n].scf_spin[lumo[n]].noccup];
	    if(n == 0)
		mixed = lumom + homom;
	    else
		mixed = lumom - homom;
	    
	    spin_info[n].scf_spin[lumo[n]].cmat[i][spin_info[n].scf_spin[lumo[n]].noccup] = 
		mixed;
	}
	/*fprintf(outfile,"\n CMAT after mixing\n");
	  print_mat(spin_info[n].scf_spin[lumo[n]].cmat,num_mo,num_mo,outfile);*/
    }
}
	    
	
	    
	    
