/*****************************************************/
/*                                                   */
/*      occ_init - This function will take care of   */
/*                   all of the input parameters     */
/*                   concerning occupations          */
/*      occ_det - Determines the occupations from    */
/*                multiplicity and charge            */
/*                                                   */
/*      By Shawn Brown - 6/23/99                     */
/*****************************************************/

#define EXTERN
#include "includes.h"
#include "common.h"
#include <ip_libv1.h>

/* Read the Opentype */

void occ_init(void){
    int i;
    int errcod; 
    int mguess=0;
    int open;

    iopen=0;
    mflag=0;
    uhf=special=twocon=hsos=singlet=0;
    
    num_ir = file30_rd_nirreps();
    
    multp=1;
    errcod = ip_data("MULTP","%d",&multp,0);
    if(errcod != IPE_OK){
	/*open = (int *) init_array(num_ir);*/
        if(ip_exist("SOCC",0)){
	    for(i=0;i<num_ir;i++){
		errcod = ip_data("SOCC","%d",&open,1,i);
		mguess += open;
	    }
	    
	    if(mguess == 1) multp = 2;
	    else if(mguess == 0) multp = 1;
	    else if(mguess == 2){
		fprintf(outfile,"\nYou must specify MULTP");
		fprintf(outfile,"\nI have no way of discerning between");
		fprintf(outfile,"\na triplet and open shell singlet\n");
		exit(1);
	    }
	    else{
		fprintf(outfile,"\nThe multiplicity will be highspin\n");
		multp = mguess + 1;
	    }
	}
	fprintf(outfile,"\nAnd I think the multiplicity is %d\n\n",multp);
	fprintf(outfile,"\nIf this is wrong, the specify MULTP keyword\n\n");
    }
	    
    
    reference = "RHF";
    refnum = ref_rhf; /* RHF for file30 flag */
    ksdft = 0;  /* do Kohn-Sham DFT? Default - no */
    errcod = ip_string("REFERENCE",&reference,0);
    if(!strcmp(reference,"ROHF")){
	refnum = ref_rohf; /* ROHF for file30 flag */
	iopen = 1;
	if(multp == 1)
	    singlet = 1;
	else if(multp > 1)
	    hsos = 1;
	
    }
    else if(!strcmp(reference,"TWOCON")){
	refnum = ref_tcscf; /*TCSCF for file30 flag */
	iopen = 2;
	twocon = 1;
    }
    else if(!strcmp(reference,"SPECIAL")){
	/* NO FILE30 flag for Special */
	iopen = 1;
	special = 1;
    }
    else if(!strcmp(reference,"UHF")){
	refnum = ref_uhf; /* UHF for file30 flag */
	uhf = 1;
    }
    else if(!strcmp(reference,"RKS")){
	refnum = ref_rks; /* flag for spin-restricted Kohn-Sham DFT */
	ksdft = 1;
    }
    else if(!strcmp(reference,"UKS")){
	refnum = ref_uks; /* flag for spin-unrestricted Kohn-Sham DFT */
	uhf = 1;
	ksdft = 1;
    }
    else{
	if(multp != 1){
	    fprintf(outfile,
		    "\n Please specify an open shell reference\n");
	    fprintf(outfile," with multpicity > 1\n");
	    fprintf(outfile," multiplicity = %d\n",multp);
	    fprintf(outfile," meference    = %s\n",reference);
	    exit(1);
	}
    }
	
/* Read Charge same as above */
    
    charge=0;
    errcod = ip_data("CHARGE","%d",&charge,0);

/* read in number of atoms and nuclear charges and total number of MO*/
    natom = file30_rd_natom();
    zvals = file30_rd_zvals();
    nbfso = file30_rd_nso();
    
/* Let's make sure that the molecule
   can have the specified multiplicity */
    
    nelec = 0;
    for(i=0; i < natom;i++)
	{
	    nelec = nelec + zvals[i];
	}
    nelec = nelec - charge;
    
    if(multp == 1 && iopen == 0) {
	if(nelec%2==1) {
	    fprintf(outfile,"\n Impossible multiplicity with charge");
	    fprintf(outfile," and # of electrons specified\n\n");
	    fprintf(outfile,"\nMultiplicity = %d\nCharge = %d",multp,charge);
	    fprintf(outfile,"\nNumber of Electrons = %d\n",nelec);
	    exit(1); }
    }
    else if(multp == 3 || (multp == 1 && singlet == 1)) {
	if(nelec%2==1) {
	    fprintf(outfile,"\n Impossible multiplicity with charge");
	    fprintf(outfile," and # of electrons specified\n\n");
	    fprintf(outfile,"\nMultiplicity = %d\nCharge = %d",multp,charge);
	    fprintf(outfile,"\nNumber of Electrons = %d\n",nelec);
	    exit(1); }
    }
    else if(multp%2 == 1 && nelec%2 ==1){
	fprintf(outfile,"\nImpossible multiplicity with charge");
	fprintf(outfile," and # of electrons specified\n");
	fprintf(outfile,"\nMultiplicity = %d\nCharge = %d",multp,charge);
	fprintf(outfile,"\nNumber of Electrons = %d\n",nelec);
	exit(1); 
    }
    else if(multp%2 == 0 && nelec%2 ==0){
	fprintf(outfile,"\nImpossible multiplicity with charge");
	fprintf(outfile," and # of electrons specified\n");
	fprintf(outfile,"\nMultiplicity = %d\nCharge = %d",multp,charge);
	fprintf(outfile,"\nNumber of Electrons = %d\n",nelec);
	exit(1); 
    }
    else {
	fprintf(outfile,
		"\nCannot check consistency of the multiplicity\n");
	fprintf(outfile,"\nand number of electrons, double check\n");
	fprintf(outfile,"your occupations\n\n");
    }
}

void occ_calc(void){
    int i,j;
    int a,b;
    
        
    if(multp == 1 && iopen == 0) {
	for(i = 0;i < (nelec/2);i++) 
		  scf_info[symm_tot[i]].nclosed++; 
    }
    
    else if(multp == 3 || (multp == 1 && singlet == 1)) {
	for(i = 0;i < ((nelec/2)-1);i++)
	    scf_info[symm_tot[i]].nclosed++;
	for(j = (nelec/2)-1;j < (nelec/2)+1;j++)
	    scf_info[symm_tot[j]].nopen++; 
    }
    
    else if(multp == 2) {
	for(i = 0;i < (nelec-1)/2;i++)
	    scf_info[symm_tot[i]].nclosed++;
	scf_info[symm_tot[((nelec-1)/2)]].nopen++;
    }
    
    else if(multp == 4) {
	if(nelec < 3){
	    fprintf(outfile,"\nNot enough electrons for a quartet\n\n");
	    exit(1);
	}
	fprintf(outfile,
		"\n*****CSCF3.0 can only guess at highspin quartets*****\n");
	for(i=0;i < (nelec/2)-1;i++)
	    scf_info[symm_tot[i]].nclosed++;
	for(i=(nelec/2)-1;i < (nelec/2)+2;i++)
	    scf_info[symm_tot[i]].nopen++;
	    }
    /* STB (10/29/99) Can define all other highspin multiplicities by
       whether it is odd or even */
    
    else if(multp%2 == 0){
	if (nelec < multp - 1){
	    fprintf(outfile,
		    "\n  Not enough electrons for a Multiplicity %d \n\n"
		    ,multp);
	    exit(1);
	}
	fprintf(outfile,
		"\n*****CSCF3.0 can only guess at highspin \n");
	fprintf(outfile,"     for Multplicity %d  *****\n\n",multp);
	
	a = (nelec/2)-((multp)-3);
	b = a+(multp-1);
		
	for(i=0;i<a;i++)
	    scf_info[symm_tot[i]].nclosed++;
	for(i=a;i<b;i++)
	    scf_info[symm_tot[i]].nopen++;
	    }
    
    else if(multp%2 == 1){
	if (nelec < multp - 1){
	    fprintf(outfile,
		    "\n  Not enough electrons for a Multiplicity %d \n\n"
		    ,multp);
	    exit(1);
	}
	
	a = (nelec/2)-((multp)-3);
	b = a+(multp-1);
			
	for(i=0;i<a;i++)
	    scf_info[symm_tot[i]].nclosed++;
	for(i=a;i<b;i++)
	    scf_info[symm_tot[i]].nopen++;
    }
}

void occ_read(int *cldpi, int *openpi){
    int i;	
    
    fprintf(outfile,"\n  Reading Occupations from file30\n");
    cldpi = file30_rd_clsdpi();
    openpi = file30_rd_openpi();
    
    for(i = 0;i < num_ir;i++){
	scf_info[i].nclosed=cldpi[i];
	scf_info[i].nopen=openpi[i];
    }
}

void occ_out(void){
    
    int i;
    
    fprintf(outfile,"\n  Symmetry block:  ");
    for(i=0;i<num_ir;i++) {
	fprintf(outfile,"%4s  ",scf_info[i].irrep_label);
    }
    fprintf(outfile,"\n  DOCC:            ");
    for(i=0;i<num_ir;i++) {
	fprintf(outfile,"%3d   ",scf_info[i].nclosed);
    }
    fprintf(outfile,"\n  SOCC:            ");
    for(i=0;i<num_ir;i++) {
	fprintf(outfile,"%3d   ",scf_info[i].nopen);
    }
    fprintf(outfile,"\n\n");
    fflush(outfile);
}
    
    
