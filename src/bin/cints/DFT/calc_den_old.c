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
#include <qt.h>
#include <libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"

struct den_info_s calc_density_old(struct coordinates geom){

    
    int i,j,k;
    int si,sj;
    int starti,startj;
    int shell_type;
    int si_n_ao,sj_n_ao;
    double x,y,z;
    double xa,ya,za;
    double rr;
    double rrtmp;
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
    int shell_center;
    double **dens;
    double *norm_ptr;
    double *dist_atom;
    struct coordinates *dist_coord;
    
    struct den_info_s den_info;
    struct shell_pair *sp;
    
    x = geom.x;
    y = geom.y;
    z = geom.z;
    
    dist_atom = init_array(Molecule.num_atoms);
    
    dist_coord = (struct coordinates *)malloc(sizeof(struct coordinates)*Molecule.num_atoms);
    timer_on("distance");
    for(i=0;i<Molecule.num_atoms;i++){
        dist_coord[i].x = x-Molecule.centers[i].x;
        dist_coord[i].y = y-Molecule.centers[i].y;
	dist_coord[i].z = z-Molecule.centers[i].z;
	dist_atom[i] = dist_coord[i].x*dist_coord[i].x
	    +dist_coord[i].y*dist_coord[i].y
	    +dist_coord[i].z*dist_coord[i].z;
    }
    n_shells = BasisSet.num_shells;
timer_off("distance");
timer_on("basis");
    for(i=k=0;i<n_shells;i++){
        
	shell_type = BasisSet.shells[i].am;
	shell_center = BasisSet.shells[i].center-1;
	xa = dist_coord[shell_center].x;
        ya = dist_coord[shell_center].y;
	za = dist_coord[shell_center].z;
	rr = dist_atom[shell_center];

	
	shell_start = BasisSet.shells[i].fprim-1;
	shell_end = shell_start+BasisSet.shells[i].n_prims;
	
	norm_ptr = GTOs.bf_norm[shell_type-1];
	timer_on("exponent");
	bastmp = 0.0;
	for(j=shell_start;j<shell_end;j++){
	    expon = -BasisSet.cgtos[j].exp;
	    coeff = BasisSet.cgtos[j].ccoeff[shell_type-1];
	    bastmp += coeff*exp(expon*rr);
	}
timer_off("exponent");	
	/*----------------------------------
	  Compute values of basis functions

	  NOTE: using Psi 3 ordering of
	        functions within shells
	 ----------------------------------*/
	switch (shell_type) {
	    
	case 1:
	    DFT_options.basis[k] = norm_ptr[0]*bastmp;
	    k++;
	    break;
	    
	case 2:
	    DFT_options.basis[k] = norm_ptr[0]*bastmp*xa;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[1]*bastmp*ya;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[2]*bastmp*za;
	    k++;
	    break;
	case 3:
	    DFT_options.basis[k] = norm_ptr[0]*bastmp*xa*xa;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[1]*bastmp*xa*ya;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[2]*bastmp*xa*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[3]*bastmp*ya*ya;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[4]*bastmp*ya*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[5]*bastmp*za*za;
	    k++;
	    break;
	case 4:
	    DFT_options.basis[k] = norm_ptr[0]*bastmp*xa*xa*xa;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[1]*bastmp*xa*xa*ya;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[2]*bastmp*xa*xa*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[3]*bastmp*xa*ya*ya;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[4]*bastmp*xa*ya*za;
	    k++;

	    DFT_options.basis[k] = norm_ptr[5]*bastmp*xa*za*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[6]*bastmp*ya*ya*ya;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[7]*bastmp*ya*ya*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[8]*bastmp*ya*za*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[9]*bastmp*za*za*za;
	    k++;
	    break;
	case 5:
	    DFT_options.basis[k] = norm_ptr[0]*bastmp*xa*xa*xa*xa;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[1]*bastmp*xa*xa*xa*ya;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[2]*bastmp*xa*xa*xa*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[3]*bastmp*xa*xa*ya*ya;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[4]*bastmp*xa*xa*ya*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[5]*bastmp*xa*xa*za*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[6]*bastmp*xa*ya*ya*ya;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[7]*bastmp*xa*ya*ya*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[8]*bastmp*xa*ya*za*za;
	    k++;

	    DFT_options.basis[k] = norm_ptr[9]*bastmp*xa*za*za*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[10]*bastmp*ya*ya*ya*ya;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[11]*bastmp*ya*ya*ya*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[12]*bastmp*ya*ya*za*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[13]*bastmp*ya*za*za*za;
	    k++;
	    
	    DFT_options.basis[k] = norm_ptr[14]*bastmp*za*za*za*za;
	    k++;
	    break;
	default:
	    fprintf(stderr,"Basis Functions of Angular Momentum %d not implemented in get_density function",shell_type);
	    fprintf(outfile,"Basis Functions of Angular Momentum %d not implemented in get_density function",shell_type);
	    punt("");
	}
    }
	timer_off("basis"); 
    /* Now contract the basis functions with the AO density matrix elements */
   timer_on("density"); 
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
   timer_off("density"); 
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
   free(dist_coord);
   free(dist_atom);
   return den_info;
}

	
	
