/* ------------------------------------------
   
   functionals go here

   ----------------------------------------- */ 

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <stdlib.h>
#include<ip_libv1.h>
#include<libciomr.h>
#include<psio.h>
#include<libint.h>
#include<pthread.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"physconst.h"

double slater(struct den_info_s den_info){
    double Cx;
    double a=(2.0/3.0);
    
    Cx = -(9.0/4.0)*a*pow(4.0/(3.0*_pi),(-1.0/3.0));
    
    return Cx*pow(den_info.den,4.0/3.0);
}


double d_slater(struct den_info_s den_info){
    double dCx;
    double a=(2.0/3.0);
 
    dCx = -3.0*a*pow(3.0/(4.0*_pi),(1.0/3.0));
    
    return dCx*pow(den_info.den,1.0/3.0);
}

double density(struct den_info_s den_info){
    return den_info.den;
}

double no_funct(struct den_info_s den_info){
    return 0.0;
}
