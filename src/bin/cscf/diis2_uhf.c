/* $Log$
 * Revision 1.1  2000/02/04 22:52:33  evaleev
 * Initial revision
 *
/* Revision 1.1  1999/11/02 23:55:56  localpsi
/* Shawn Brown - (11/2/99) Modified to the code in a few major ways.
/*
/* 1.  Added the capability to do UHF.  All of the features available with the
/* other refrences have been added for UHF.
/*
/* 2.  For UHF, I had to alter the structure of file30. (See cleanup.c for a
/* map)  This entailed adding a pointer array right after the header in the SCF
/* section of file30 that pointed to all of the data for the SCF caclulation.
/* Functions were added to libfile30 to account for this and they are
/* incorporated in this code.
/*
/* 3.  Updated and fixed all of the problems associated with my previous
/* guessing code.  The code no longer uses OPENTYPE to specify the type of
/* occupation.  The keword REFERENCE and MULTP can now be used to indicate any
/* type of calculation.  (e.g. ROHF with MULTP of 1 is an open shell singlet
/* ROHF calculation)  This code was moved to occ_fun.c.  The code can also
/* guess at any multplicity in a highspin case, provided enough electrons.
/*
/* Revision 1.1.1.1  1999/04/12 16:59:25  evaleev
/* Added a version of CSCF that can work with CINTS.
/* -Ed
/*
 * Revision 1.1  1991/06/15  20:22:20  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

extern double delta;

static double *btemp, **bold, **bmat;

static struct diis_mats {
    double ****fock_c;
    double ****error;
    double ***error_tot;
} *diism,dtemp; 

void diis_uhf(scr1,scr2,scr3)
    double **scr1, **scr2, **scr3;
{
    int i,j,k,ij;
    int errcod;
    int a,b,c,e;
    int m,n,nn,mm,num_mo;
    int try = 0;
    int last = iter-1;
    int col = iter+1;
    double etemp, dotp, norm, determ;
    double scale;
    struct symm *s;
    struct diis_mats *d;
    int diis_print=0;
    
    if(diism == NULL){
	/*bold = (double ***) malloc(sizeof(double ***)*2);*/
	bmat = (double **) init_matrix(ndiis+1,ndiis+1);
	bold = (double **) init_matrix(ndiis,ndiis);
	/*bold[1] = (double **) init_matrix(ndiis,ndiis);*/
	btemp = (double *) init_array(ndiis+1);
	
	diism = (struct diis_mats *) malloc(sizeof(struct diis_mats)*ndiis);
	for(m=0; m < ndiis ; m++) {
	    d = &diism[m];
	    d->fock_c = 
		(double ****) malloc(sizeof(double ***)*2);
	    d->error = 
		(double ****) malloc(sizeof(double ***)*2);
	    d->error_tot =
		(double ***) malloc(sizeof(double ***)*num_ir);
            for (j=0;j<num_ir;j++){
		if(nn=scf_info[j].num_so) {
		d->error_tot[j] = (double **) init_matrix(nn,nn);
		}
	    }
	    for(n = 0; n < 2;n++){
	       d->fock_c[n] = 
		   (double ***) malloc(sizeof(double **)*num_ir);
	       d->error[n] =
		   (double ***) malloc(sizeof(double **)*num_ir);
	       for(j=0; j < num_ir ; j++) {
		   if(nn=scf_info[j].num_so) {
		       d->fock_c[n][j] = 
			   (double **) init_matrix(nn,nn);
		       d->error[n][j] = 
			   (double **) init_matrix(nn,nn);
		   }
	       }
	    }
	}
    }	
    scale = 1.0 + dampsv;
    
    for(n = 0; n < 2; n++){
	if (iter > ndiis) {
	    last = ndiis-1;
	    col = ndiis+1;
	    if(n == 0){
		dtemp = diism[0];
		for (i=0; i < last ; i++) {
		    diism[i] = diism[i+1];
		}
		diism[last] = dtemp;
	    }
	} 
	
	
	errcod = ip_boolean("DIIS_PRINT",&diis_print,0);
	if(iter == 1 && diis_print)
	    ffile(&diis_out,"diis_out.dat",0);
	
/* save ao fock matrices in fock_save */
	
	d = &diism[last];
	for (m=0; m < num_ir ; m++) {
	    s = &scf_info[m];
	    if(nn=s->num_so) {
		num_mo = s->num_mo;
		tri_to_sq(spin_info[n].scf_spin[m].fock_pac
			  ,d->fock_c[n][m],nn);
		
		/*fprintf(diis_out,"\nAO FOCK");
		  print_mat(d->fock_c[n][m],nn,nn,diis_out);*/
		
		/*mxmb(spin_info[n].scf_spin[m].cmat,nn,1
		  ,d->fock_c[n][m],1,nn,scr1,1,nn,nn,nn,nn);
		  mxmb(scr1,1,nn,spin_info[n].scf_spin[m].cmat,1
		  ,nn,scr2,1,nn,nn,nn,nn);*/

		mmult(spin_info[n].scf_spin[m].cmat,1
		      ,d->fock_c[n][m],0,scr1,0,num_mo,nn,nn,0);
		mmult(scr1,0,spin_info[n].scf_spin[m].cmat,0
		,scr2,0,num_mo,nn,num_mo,0);
		
		  /*  fprintf(diis_out,"\nMO FOCK");
		      print_mat(scr2,nn,nn,diis_out);*/

		  mm = spin_info[n].scf_spin[m].noccup;
		  zero_mat(scr1,nn,nn);
		  for (i=0; i < mm; i++) {
		      for (j=mm; j < num_mo ; j++ ) {
			  scr1[i][j]= scr2[i][j];
			  scr1[j][i]= scr2[i][j];
		      }
		  }
	    }
	    /*    fprintf(diis_out,"\nMO ERROR");
		  print_mat(scr1,nn,nn,diis_out);*/
	    
	    /* transform error matrix into an orthonormal basis */
	    
	    /*mmult(scf_info[m].sahalf,1,scr1,0,scr2,0,nn,nn,nn,0);
	      mmult(scr2,0,scf_info[m].sahalf,0,d->error[n][m],0,nn,nn,nn,0);*/
	    
	    /* transform error matrix into ao basis */
	     mxmb(spin_info[n].scf_spin[m].cmat,1,nn
		 ,scr1,1,nn,scr2,1,nn,nn,nn,nn);
	    mxmb(scr2,1,nn,spin_info[n].scf_spin[m].cmat,nn,1
	    ,d->error[n][m],1,nn,nn,nn,nn);
	    
	    /*mmult(spin_info[n].scf_spin[m].cmat,0
		  ,scr1,0,scr2,0,nn,num_mo,num_mo,0);
	    mmult(scr2,0,spin_info[n].scf_spin[m].cmat,1
	    ,d->error[n][m],0,nn,num_mo,nn,0);*/
	    
	    for(i=0; i < nn ; i++) {
		for(j=0; j <= i ; j++) {
		    etemp=fabs(scr1[i][j]);
		    diiser = MAX0(diiser,etemp);
		}
	    }
	}
    }
    
    for(m = 0;m < num_ir;m++){
	if(nn=scf_info[m].num_so){
	    for(i=0;i<nn;i++){
		for(j=0;j<nn;j++){
		    diism[last].error_tot[m][i][j] = 
			(diism[last].error[0][m][i][j]
			 +diism[last].error[1][m][i][j])/2;
		}
	    }
	}
    }
	
    /* then set up b matrix */
    if (iter > ndiis) {
	for (i=0; i < last ; i++) {
	    for (j=0; j <= i ; j++) {
		bold[i][j]=bold[j][i]=bold[i+1][j+1];
	    }
	}
    }
    for (i=0; i <= last ; i++) {
	etemp=0.0;
	for (m=0; m < num_ir ; m++) {
	    s = &scf_info[m];
	    if(nn=s->num_so) {
		sdot(diism[i].error_tot[m]
		     ,diism[last].error_tot[m],nn,&dotp);
		etemp += dotp;
	    }
	}
	bold[i][last]=bold[last][i] = etemp;
    }
    
    bmat[0][0] = 0.0;
    btemp[0] = -1.0;
    norm = 1.0/bold[0][0];
    /*norm = 1.0;*/
    for (i=1; i <= last+1 ; i++) {
	bmat[i][0]=bmat[0][i] = -1.0;
	btemp[i] = 0.0;
	for (j=1; j <= i ; j++) {
	    bmat[i][j]=bmat[j][i] = bold[i-1][j-1]*norm;
	    if(i==j) bmat[i][j] *= scale;
	}
    }
    
    if(diis_print){
	fprintf(diis_out,"\nBMAT for iter %d",iter);
	print_mat(bmat,ndiis+1,ndiis+1,diis_out);
    }
    
    if (iter-1) {
	flin(bmat,btemp,col,1,&determ);
	
	
	/* test for poorly conditioned equations */
	try = 0;
	while (fabs(determ) < 1.0e-19 && try < last) {
	    
	    try++;
	    col--;
	    
	    bmat[0][0] = 0.0;
	    btemp[0] = -1.0;
	    norm=1.0/bold[try][try];
	    /*norm = 1.0;*/
	    for (i=1; i <= ndiis-try ; i++) {
		bmat[i][0]=bmat[0][i] = -1.0;
		for (j=1; j <= i ; j++) {
		    bmat[i][j]=bmat[j][i]=bold[i+try-1][j+try-1]*norm;
		    if(i==j) bmat[i][j] *= scale;
		}
		btemp[i] = 0.0;
	    }
	    if(diis_print){
		fprintf(diis_out,"\nCorrected BMAT for iter %d",iter);
		print_mat(bmat,col,col,diis_out);
	    }
	    
	    flin(bmat,btemp,col,1,&determ);
	}
	
	if(fabs(determ) < 10.0e-20) {
	    printf(" try %d no good\n",try);
	    return;
	}
	
	if((iter >= it_diis)) {
	    for(n=0;n<2;n++){
		for (m=0; m < num_ir ; m++) {
		    s = &scf_info[m];
		    if(nn=s->num_so) {
			for (i=ij=0; i < nn ; i++) {
			    for (j=0; j <= i ; j++,ij++) {
				int kk=1;
				etemp=0.0;
				for (k=try; k < last+1 ; k++) {
				    etemp += btemp[kk]
					*diism[k].fock_c[n][m][i][j];
				    kk++;
				}
				spin_info[n].scf_spin[m].fock_pac[ij] = etemp;
			    }
			}
		    }
		}
	    }
	}
    }
}
