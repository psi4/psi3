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
    
    Cx = -(9.0/4.0)*a*2.0*pow(3.0/(4.0*_pi),(1.0/3.0));
    
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

/* Implementation of the VWN functionals -- Can J. Phys., 58, 1980
   pp. 1200 */

double Pade_int(double p, double x0, double b, double c, double A){

    double x,Q,X,X0,invQ;
    double term1,term2,mult1,term3,term4;
    double temp1,temp2,temp3,temp4,temp5,temp6,temp7;
    double ec;

    temp7 = 3.0/(4.0*_pi*2.0*p);
    x = pow(temp7,0.166666666667);
    Q = sqrt(4.0*c-b*b);
    invQ = 1.0/Q;
    X = 1.0/(x*x+b*x+c);
    X0 = 1.0/(x0*x0+b*x0+c);

    temp1 = x*x*X;
    term1 = log(temp1);

    temp2 = Q/(2.0*x+b);
    temp3 = atan(temp2);
    term2 = 2.0*b*invQ*temp3;

    mult1 = b*x0*X0;

    temp4 = (x-x0)*(x-x0)*X;
    term3 = log(temp4);

    temp5 = 2.0*(2.0*x0+b)*invQ;
    temp6 = Q/(2.0*x+b);
    term4 = temp5*atan(temp6);

    ec = A*(term1+term2-mult1*(term3+term4));

    return ec;

}
double d_Pade_int(double p, double x0, double b, double c, double A){
    
    double x,Q,X,X0,invQ;
    double term1,term2,term3,term4,term5,term6;
    double mult1;
    double temp1,temp2,temp3,temp4,temp5;
    double dec;

    temp5 = 3.0/(4.0*_pi*2.0*p);
    x = pow(temp5,0.166666666667);
    Q = sqrt(4*c-b*b);
    invQ = 1.0/Q;
    X = 1.0/(x*x+b*x+c);
    X0 = 1.0/(x0*x0+b*x0+c);

    term1 = 2.0/x;

    temp1 = -2.0*x+b;
    term2 = temp1*X;

    temp2 = (2.0*x+b)*(2.0*x+b)+Q*Q;
    term3 = -4.0*b/temp2;

    mult1 = -b*x0*X0;

    temp3 = x-x0;
    term4 = 2.0/temp3;

    temp4 = 2.0*x0+b;
    term5 = -4.0*temp4/temp2;

    dec = A*(term1+term2+term3+mult1*(term4+term2+term5));

    return dec;
}                                            

double VWN5_r(struct den_info_s den_info){
    double A = 0.06218414;
    double x0 = -0.10498;
    double b = 3.72744;
    double c = 12.9352;

    double ec;

    ec = Pade_int(den_info.den,x0,b,c,A);

    return den_info.den*ec;
}

double d_VWN5_r(struct den_info_s den_info){

    double A = 0.06218414;
    double x0 = -0.10498;
    double b = 3.72744;
    double c = 12.9352;
    double x, temp1;
    double ec,dec;
    double p;

    p = den_info.den;

    temp1 = 3.0/(4.0*_pi*2.0*p);
    x = pow(temp1,0.166666666667);

    ec = Pade_int(den_info.den,x0,b,c,A);
    dec = Pade_int(den_info.den,x0,b,c,A);

    temp1 = -x/6.0;

    return ec+temp1*dec;
}   

 double VWN4_r(struct den_info_s den_info){
    double A = 0.06218414;
    double x0 = -0.409286;
    double b = 13.0720;
    double c = 42.7198;

    double ec;

    ec = Pade_int(den_info.den,x0,b,c,A);

    return den_info.den*ec;
}

double d_VWN4_r(struct den_info_s den_info){

    double A = 0.06218414;
    double x0 = -0.409286;
    double b = 13.0720;
    double c = 42.7198;
    double x, temp1;
    double ec,dec;
    double p;

    p = den_info.den;

    temp1 = 3.0/(4.0*_pi*p);
    x = pow(temp1,0.166666666667);

    ec = Pade_int(den_info.den,x0,b,c,A);
    dec = Pade_int(den_info.den,x0,b,c,A);

    temp1 = -x/6.0;

    return ec+temp1*dec;
}                              
