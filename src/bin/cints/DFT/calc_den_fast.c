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

struct den_info_s calc_density_fast(struct coordinates geom){
    
    int i,j,k,l,m;
    int am2shell;
    int shell_type;
    int shell_start;
    int shell_end;
    int n_shells;
    int num_ao,ndocc;
    int shell_center;
    double x,y,z;
    double xa,ya,za;
    double rr;
    double rrtmp;
    double bastmp;
    double den_sum;
    double coeff;
    double expon;
    double *norm_ptr;
    double *dist_atom;
    double *temp_arr;
    double norm_ptr0,norm_ptr1,norm_ptr2;
    double norm_ptr3,norm_ptr4,norm_ptr5;
    double norm_ptr6,norm_ptr7,norm_ptr8;
    double norm_ptr9;
    
    struct coordinates *dist_coord;
    struct den_info_s den_info;

    x = geom.x;
    y = geom.y;
    z = geom.z;
    
    num_ao = BasisSet.num_ao;
    ndocc = MOInfo.ndocc;
    temp_arr = init_array(ndocc);
    
    /* ---------------------------------
       Compute distances from atom that 
       the basis function is centered on
       to the grid point
       --------------------------------*/
    
    dist_atom = init_array(Molecule.num_atoms);
    dist_coord = (struct coordinates *)
	malloc(sizeof(struct coordinates)*Molecule.num_atoms);
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
    
    l=0;
    m=0;
    for(i=0;i<BasisSet.max_am;i++){
	
	norm_ptr = GTOs.bf_norm[i];
	shell_type = i;
	switch(i){
	    
	    /* ------------------------------------------------------*/
	    /* S-type functions */
	    
	case 0:
	    norm_ptr0 = norm_ptr[0];
	    for(j=0;j<BasisSet.shells_per_am[i];j++){
		
		am2shell = BasisSet.am2shell[m];
		shell_center = BasisSet.shells[am2shell].center-1;
		xa = dist_coord[shell_center].x;
		ya = dist_coord[shell_center].y;
		za = dist_coord[shell_center].z;
		rr = dist_atom[shell_center];
		
		shell_start = BasisSet.shells[am2shell].fprim-1;
		shell_end = shell_start
		    +BasisSet.shells[am2shell].n_prims;
		
		
		timer_on("exponent");
		bastmp = 0.0;
		for(k=shell_start;k<shell_end;k++){
		    expon = -BasisSet.cgtos[k].exp;
		    coeff = BasisSet.cgtos[k].ccoeff[shell_type];
		    bastmp += coeff*exp(expon*rr);
		}
		
		timer_off("exponent");
		
		DFT_options.basis[l] = norm_ptr0*bastmp;
		l++;
		
		m++;
	    }
 
	    break;
	    /* ------------------------------------------------------*/
	    
	    
	    /* ------------------------------------------------------*/
	    /* P-type functions */
	    
	case 1:
	    norm_ptr0 = norm_ptr[0];
	    norm_ptr1 = norm_ptr[1];
	    norm_ptr2 = norm_ptr[2];
	    for(j=0;j<BasisSet.shells_per_am[i];j++){
		
		am2shell = BasisSet.am2shell[m];
		shell_center = BasisSet.shells[am2shell].center-1;
		xa = dist_coord[shell_center].x;
		ya = dist_coord[shell_center].y;
		za = dist_coord[shell_center].z;
		rr = dist_atom[shell_center];
		
		shell_start = BasisSet.shells[am2shell].fprim-1;
		shell_end = shell_start
		    +BasisSet.shells[am2shell].n_prims;
		
		
		timer_on("exponent");
		bastmp = 0.0;
		for(k=shell_start;k<shell_end;k++){
		    expon = -BasisSet.cgtos[k].exp;
		    coeff = BasisSet.cgtos[k].ccoeff[shell_type];
		    bastmp += coeff*exp(expon*rr);
		}
		
		timer_off("exponent");
		
		DFT_options.basis[l] = norm_ptr0*bastmp*xa;
		l++;
		
		DFT_options.basis[l] = norm_ptr1*bastmp*ya;
		l++;
		
		DFT_options.basis[l] = norm_ptr2*bastmp*za;
		l++;
		m++;
	    }
	    break;
	    /* ------------------------------------------------------*/
	    
	    /* ------------------------------------------------------*/
	    /* D-type functions */ 
	
	case 2:
	    norm_ptr0 = norm_ptr[0];
	    norm_ptr1 = norm_ptr[1];
	    norm_ptr2 = norm_ptr[2];
	    norm_ptr3 = norm_ptr[3];
	    norm_ptr4 = norm_ptr[4];
	    norm_ptr5 = norm_ptr[5];
	    for(j=0;j<BasisSet.shells_per_am[i];j++){
		am2shell = BasisSet.am2shell[m];
		shell_center = BasisSet.shells[am2shell].center-1;
		xa = dist_coord[shell_center].x;
		ya = dist_coord[shell_center].y;
		za = dist_coord[shell_center].z;
		rr = dist_atom[shell_center];
		
		shell_start = BasisSet.shells[am2shell].fprim-1;
		shell_end = shell_start
		    +BasisSet.shells[am2shell].n_prims;
		
		
		timer_on("exponent");
		bastmp = 0.0;
		for(k=shell_start;k<shell_end;k++){
		    expon = -BasisSet.cgtos[k].exp;
		    coeff = BasisSet.cgtos[k].ccoeff[shell_type];
		    bastmp += coeff*exp(expon*rr);
		}
		timer_off("exponent");
		
		DFT_options.basis[l] = norm_ptr0*bastmp*xa*xa;
		l++;
		
		DFT_options.basis[l] = norm_ptr1*bastmp*xa*ya;
		l++;
		
		DFT_options.basis[l] = norm_ptr2*bastmp*xa*za;
		l++;
		
		DFT_options.basis[l] = norm_ptr3*bastmp*ya*ya;
		l++;
		
		DFT_options.basis[l] = norm_ptr4*bastmp*ya*za;
		l++;
		
		DFT_options.basis[l] = norm_ptr5*bastmp*za*za;
		l++;
		m++;
	    }
	    break;
	    /* ------------------------------------------------------*/
	    
	    /* ------------------------------------------------------*/
	    /* F-type functions */
	    
	case 3:
	    norm_ptr0 = norm_ptr[0];
	    norm_ptr1 = norm_ptr[1];
	    norm_ptr2 = norm_ptr[2];
	    norm_ptr3 = norm_ptr[3];
	    norm_ptr4 = norm_ptr[4];
	    norm_ptr5 = norm_ptr[5];
	    norm_ptr6 = norm_ptr[6];
	    norm_ptr7 = norm_ptr[7];
	    norm_ptr8 = norm_ptr[8];
	    norm_ptr9 = norm_ptr[9];
	    for(j=0;j<BasisSet.shells_per_am[i];j++){
		am2shell = BasisSet.am2shell[m];
		shell_center = BasisSet.shells[am2shell].center-1;
		xa = dist_coord[shell_center].x;
		ya = dist_coord[shell_center].y;
		za = dist_coord[shell_center].z;
		rr = dist_atom[shell_center];
		
		shell_start = BasisSet.shells[am2shell].fprim-1;
		shell_end = shell_start
		    +BasisSet.shells[am2shell].n_prims;
		
		
		timer_on("exponent");
		bastmp = 0.0;
		for(k=shell_start;k<shell_end;k++){
		    expon = -BasisSet.cgtos[k].exp;
		    coeff = BasisSet.cgtos[k].ccoeff[shell_type];
		    bastmp += coeff*exp(expon*rr);
		}
		timer_off("exponent");
		
		DFT_options.basis[l] = norm_ptr0*bastmp*xa*xa*xa;
		l++;
		
		DFT_options.basis[l] = norm_ptr1*bastmp*xa*xa*ya;
		l++;
		
		DFT_options.basis[l] = norm_ptr2*bastmp*xa*xa*za;
		l++;
		
		DFT_options.basis[l] = norm_ptr3*bastmp*xa*ya*ya;
		l++;
		
		DFT_options.basis[l] = norm_ptr4*bastmp*xa*ya*za;
		l++;
		
		DFT_options.basis[l] = norm_ptr5*bastmp*xa*za*za;
		l++;
		
		DFT_options.basis[l] = norm_ptr6*bastmp*ya*ya*ya;
		l++;
		
		DFT_options.basis[l] = norm_ptr7*bastmp*ya*ya*za;
		l++;
		
		DFT_options.basis[l] = norm_ptr8*bastmp*ya*za*za;
		l++;
		
		DFT_options.basis[l] = norm_ptr9*bastmp*za*za*za;
		l++;
		m++;
	    }
	    
	    break;
	    
	    /* ------------------------------------------------------*/
	    
	default:
	    fprintf(stderr,"Basis Functions of Angular Momentum %d not implemented in get_density function",shell_type);
	    fprintf(outfile,"Basis Functions of Angular Momentum %d not implemented in get_density function",shell_type);
	    punt("");
	}
	
    }
   
    timer_off("basis"); 
/* Now contract the basis functions with the AO density matrix elements */
    timer_on("density"); 
    /*for(i=0;i<num_ao;i++)
	fprintf(outfile,"\nBasis[%d] = %10.10lf",i,DFT_options.basis[i]);
    fprintf(outfile,"\n");*/
    
    if(UserOptions.reftype == rhf){
	den_sum = 0.0;    
	for(i=0;i<ndocc;i++){
	    for(j=0;j<num_ao;j++){
		temp_arr[i] += Cocc[j][i]*DFT_options.basis[j];
	    }
	}
	dot_arr(temp_arr,temp_arr,MOInfo.ndocc,&den_sum);
	den_info.den = den_sum;
    }
    free(temp_arr);
    timer_off("density");
    free(dist_coord);
    free(dist_atom);
    return den_info;
}

	
	
