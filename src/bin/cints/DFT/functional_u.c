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
#include<qt.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"physconst.h"

double slater_u(struct den_info_s den_info){
    double Cx;
    
    /*Cx = -(9.0/4.0)*(2.0/3.0)*pow(3.0/(4.0*_pi),(1.0/3.0));*/    
    Cx = -0.930525736349100;
    
    return Cx*(pow(den_info.dena,4.0/3.0)+pow(den_info.denb,4.0/3.0));
}


double d_slater_a(struct den_info_s den_info){
    double dCx;
  
 
    /*dCx = -3.0*(2.0/3.0)*pow(3.0/(4.0*_pi),(1.0/3.0));*/

    dCx = -1.240700981798800;
    
    return dCx*pow(den_info.dena,1.0/3.0);
}

double d_slater_b(struct den_info_s den_info){
    double dCx;
 
    /*dCx = -3.0*(2.0/3.0)*pow(3.0/(4.0*_pi),(1.0/3.0));*/

    dCx = -1.240700981798800;
    
    return dCx*pow(den_info.denb,1.0/3.0);
}


double density_a(struct den_info_s den_info){
    return den_info.dena;
}

double density_b(struct den_info_s den_info){
    return den_info.denb;
}

/* Implementation of the VWN functionals -- Can J. Phys., 58, 1980
   pp. 1200 */

/*double Pade_int(double p, double x0, double b, double c, double A, double Q){

    double x,X,X0,invQ;
    double txpb,Qd2xpb;
    double term1,term2,mult1,term3,term4;
    double temp1,temp2,temp3,temp4,temp5,temp6,temp7;
    double ec;
  
    x = sqrt(p);
    txpb = 2.0*x+b;
    Qd2xpb = Q/txpb;
    invQ = 1.0/Q;
    X = 1.0/(x*x+b*x+c);
    X0 = 1.0/(x0*x0+b*x0+c);

    temp1 = x*x*X;
    term1 = log(temp1);

    temp3 = atan(Qd2xpb);
    term2 = 2.0*b*invQ*temp3;

    mult1 = b*x0*X0;

    temp4 = (x-x0)*(x-x0)*X;
    term3 = log(temp4);

    temp5 = 2.0*(2.0*x0+b)*invQ;
    term4 = temp5*atan(Qd2xpb);

    ec = A*(term1+term2-mult1*(term3+term4));
    return ec;

}

double d_Pade_int(double p, double x0, double b, double c, double A, double Q){
    
    double x,Q2,X,X0,invQ;
    double txpb;
    double Qsptxpbs;
    double term1,term2,term3,term4,term5,term6;
    double mult1;
    double temp1,temp2,temp3,temp4,temp5;
    double dec;

    x=sqrt(p);
    A=A/(2*x);
    Q2 = Q*Q;
    txpb = 2.0*x+b;
    Qsptxpbs = txpb*txpb+Q2;
    X = 1.0/(x*x+b*x+c);
    X0 = 1.0/(x0*x0+b*x0+c);
    
    term1 = 2.0/x;
   
    term2 = -txpb*X;

    term3 = -4.0*b/Qsptxpbs;

    mult1 = -b*x0*X0;

    temp3 = x-x0;
    term4 = 2.0/temp3;

    temp4 = 2.0*x0+b;
    term5 = -4.0*temp4/Qsptxpbs;

    dec = A*(term1+term2+term3+mult1*(term4+term2+term5));
    return dec;
}                                            

double VWN5_r(struct den_info_s den_info){
    double A = 0.0621814/2.0;
    double x0 = -0.10498;
    double b = 3.72744;
    double c = 12.9352;
    double Q = 6.1519908198;
    double p;
    double ec;
    double x;
    double temp;

    p = 2.0*den_info.den;
    temp = 3.0/(4.0*_pi*p);
    x = pow(temp,0.3333333333333);
    ec = Pade_int(x,x0,b,c,A,Q);
    return p*ec;
}


double d_VWN5_r(struct den_info_s den_info){

    double A = 0.0621814/2.0;
    double x0 = -0.10498;
    double b = 3.72744;
    double c = 12.9352;
    double Q = 6.1519908198;
    double x, temp1;
    double dxdp;
    double ec,dec;
    double ec_sum,dec_sum;
    double p;

    p = 2.0*den_info.den;

    temp1 = 3.0/(4.0*_pi*p);
    x = pow(temp1,0.3333333333);
    dxdp = -x/(3.0*p);
    ec = Pade_int(x,x0,b,c,A,Q);
    dec = d_Pade_int(x,x0,b,c,A,Q);
    return ec+p*dxdp*dec;
}   
*/
/* This is the functional in which Gaussian uses */

/* double VWN4_r(struct den_info_s den_info){
    double A = 0.0621814/2.0;
    double x0 = -0.409286;
    double b = 13.0720;
    double c = 42.7198;
    double Q = 0.04489988864;
    double p;
    double ec;
    double x;
    double temp;

    p = 2.0*den_info.den;
    temp = 3.0/(4.0*_pi*p);
    x = pow(temp,0.3333333333333);
    ec = Pade_int(x,x0,b,c,A,Q);
    return p*ec;
}

double d_VWN4_r(struct den_info_s den_info){

    double A = 0.0621814/2.0;
    double x0 = -0.409286;
    double b = 13.0720;
    double c = 42.7198;
    double Q = 0.04489988864;
    double x, temp1;
    double ec,dec;
    double p;
    double dxdp;

    p = 2.0*den_info.den;
    
    temp1 = 3.0/(4.0*_pi*p);
    x = pow(temp1,0.3333333333);
    dxdp = -x/(3.0*p);
    ec = Pade_int(x,x0,b,c,A,Q);
    dec = d_Pade_int(x,x0,b,c,A,Q);
    
    return ec+p*dxdp*dec;
}     */                         


