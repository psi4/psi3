#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <qt.h>

#include"defines.h"
#define EXTERN
#include"global.h"

double calc_exp_basis(int shell_num, double rr){
    
    int i;
    int shell_type;
    int shell_start;
    int shell_end;
    double expon,coeff;
    double bastmp;
    
    shell_type = BasisSet.shells[shell_num].am;
    shell_start = BasisSet.shells[shell_num].fprim-1;
    shell_end = shell_start+BasisSet.shells[shell_num].n_prims;
    
    bastmp = 0.0;
    for(i=shell_start;i<shell_end;i++){
	expon = -BasisSet.cgtos[i].exp;
	coeff = BasisSet.cgtos[i].ccoeff[shell_type-1];
	bastmp += coeff*exp(expon*rr);
    }
    
    return bastmp;
}

double calc_radial_bas(int shell_num, double rr, double r){
    
    int i;
    int shell_type;
    int end;
    double bastmp;
   
    shell_type = BasisSet.shells[BasisSet.am2shell[shell_num]].am;
    end = ioff[shell_type-1];
    bastmp = calc_exp_basis(shell_num,rr);
    
    for(i=0;i<end;i++)
	bastmp *= r;
    
    return bastmp;
}
