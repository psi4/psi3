/*----------------------------------------------------
  
  get_density.c
  
  By Shawn Brown
  
  The code contains functions that will retrieve 
  and process the density at a given x,y,z coord
  in space
  
  --------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"

struct den_info_s calc_density(struct coordinates geom){
    
    int i,j,k;
    int si,sj;
    int starti,startj;
    int shell_type;
    int si_n_ao,sj_n_ao;
    double x,y,z;
    double xa,ya,za;
    double rr;
    double bastmp;
    double den_sum=0.0;
    double den_sum_o = 0.0;
    double dena,denb;
    double bas1,bas2;
    int  shell_start;
    int  shell_end;
    double shell_den;
    double coeff;
    double expon;
    int n_shells;
    double **dens;
    
    struct den_info_s den_info;
    struct shell_pair *sp;
    
    x = geom.x;
    y = geom.y;
    z = geom.z;
    
    n_shells = BasisSet.num_shells;

    for(i=k=0;i<n_shells;i++){
	shell_type = BasisSet.shells[i].am;
	xa = x-Molecule.centers[BasisSet.shells[i].center-1].x;
	ya = y-Molecule.centers[BasisSet.shells[i].center-1].y;
	za = z-Molecule.centers[BasisSet.shells[i].center-1].z;
	
	rr = xa*xa+ya*ya+za*za;
	
	bastmp = 0.0;
	
	shell_start = BasisSet.shells[i].fprim-1;

	shell_end = shell_start+BasisSet.shells[i].n_prims;
	
	for(j=shell_start;j<shell_end;j++){
	    
	    expon = -BasisSet.cgtos[j].exp;
	    coeff = BasisSet.cgtos[j].ccoeff[shell_type-1];
	    bastmp += coeff*exp(expon*rr);
	}
	
	switch (shell_type) {
	    
	case 1:
	    DFT_options.basis[k] = bastmp;
	    k++;
	    break;
	    
	case 2:
	    DFT_options.basis[k] = bastmp*xa;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*ya;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*za;
	    k++;
	    break;
	case 3:
	    DFT_options.basis[k] = bastmp*xa*xa;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*ya*ya;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*za*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*ya;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*ya*za;
	    k++;
	    break;
	case 4:
	    DFT_options.basis[k] = bastmp*xa*xa*xa;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*ya*ya*ya;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*za*za*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*xa*ya;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*xa*za;
	    k++;

	    DFT_options.basis[k] = bastmp*xa*ya*ya;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*ya*ya*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*za*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*ya*za*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*ya*za;
	    k++;
	    break;
	case 5:
	    DFT_options.basis[k] = bastmp*xa*xa*xa*xa;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*ya*ya*ya*ya;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*za*za*za*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*xa*xa*ya;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*xa*xa*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*ya*ya*ya;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*ya*ya*ya*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*za*za*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*ya*za*za*za;
	    k++;

	    DFT_options.basis[k] = bastmp*xa*xa*ya*ya;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*xa*za*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*ya*ya*za*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*xa*ya*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*ya*ya*za;
	    k++;
	    
	    DFT_options.basis[k] = bastmp*xa*ya*za*za;
	    k++;
	    break;
	default:
	    fprintf(outfile,"\nBasis Functions of Angular Momentum %d not implemented in get_density function\n",shell_type);
	      exit(1);
	}
    }
    
    /* Now contract the basis functions with the AO density matrix elements */
    
    if(UserOptions.reftype == rhf){
	den_sum = 0.0;
	for(si=0;si<n_shells;si++){
	    si_n_ao = ioff[BasisSet.shells[si].am];
	    starti = BasisSet.shells[si].fao-1;
	    for(sj=0;sj<n_shells;sj++){
		sj_n_ao = ioff[BasisSet.shells[sj].am];
		startj = BasisSet.shells[sj].fao-1;		
		sp = &(BasisSet.shell_pairs[si][sj]);
		    for(i=0;i<si_n_ao;i++){
			bas1 = DFT_options.basis[starti+i];
			for(j=0;j<sj_n_ao;j++){
			    bas2 = DFT_options.basis[startj+j];
			    den_sum += sp->dmat[i][j]*bas1*bas2;
			}
		    }
	    }
	    
	}
	den_info.den = 0.5*den_sum;
    }
    
    /*if(UserOptions.reftype == uhf){
	for(si=0;si<n_shells;si++){
	    si_n_ao = ioff[BasisSet.shells[si].am];
	    starti = BasisSet.shells[si].fao-1;
	    for(sj=si;sj<n_shells;sj++){
		sj_n_ao = ioff[BasisSet.shells[sj].am];
		startj = BasisSet.shells[sj].fao-1;
		sp = &(BasisSet.shell_pairs[si][sj]);
		    for(i=starti;i<starti+si_n_ao;i++){
			bas1 = den_info.basis_arr[starti+i];
			for(j=startj;j<startj+sj_n_ao;j++){
			    dena = sp->dmata[i][j];
			    denb = sp->dmatb[i][j];
			    bas2 = den_info.basis_arr[startj+j];
			    den_sum += dena*bas1*bas2;
			    den_sum_o += denb*bas1*bas2; 
			}
		    }
	    }
	}
	den_info.den_a = den_sum;
	den_info.den_b = den_sum_o;
	}*/
    return den_info;
}

	
	
