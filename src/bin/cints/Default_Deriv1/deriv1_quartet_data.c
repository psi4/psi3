#include<math.h>
#include<stdio.h>
#include<string.h>
#include<memory.h>
#include<stdlib.h>
#include<libint.h>
#include<libderiv.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"int_fjt.h"


/*--------------------------------------------------------------------------------
  This function computes constants used in OSRR for a given quartet of primitives
 --------------------------------------------------------------------------------*/
void deriv1_quartet_data(prim_data *Data, double_array_t *fjt_table, double AB2, double CD2,
			 struct shell_pair *sp1, struct shell_pair *sp2, 
			 int am, int pi, int pj, int pk, int pl, double scale)
{

  /*------------------------------------------------------
    External data necessary to compute auxiliary function
   ------------------------------------------------------*/
  static double F0[40] = {1.0,  1.0/3.0,  1.0/5.0,  1.0/7.0,  1.0/9.0,
                  1.0/11.0, 1.0/13.0, 1.0/15.0, 1.0/17.0, 1.0/19.0,
                  1.0/21.0, 1.0/23.0, 1.0/25.0, 1.0/27.0, 1.0/29.0,
                  1.0/31.0, 1.0/33.0, 1.0/35.0, 1.0/37.0, 1.0/39.0,
                  1.0/41.0, 1.0/43.0, 1.0/45.0, 1.0/47.0, 1.0/49.0,
                  1.0/51.0, 1.0/53.0, 1.0/55.0, 1.0/57.0, 1.0/59.0,
                  1.0/61.0, 1.0/63.0, 1.0/65.0, 1.0/67.0, 1.0/69.0,
                  1.0/71.0, 1.0/73.0, 1.0/75.0, 1.0/77.0, 1.0/79.0};

  /*----------------
    Local variables
   ----------------*/
  struct coordinates PQ, W;
  int i;
  double coef1;
  double PQ2;
  double oozn;
  double zeta, eta, rho;

  zeta = sp1->gamma[pi][pj];
  eta = sp2->gamma[pk][pl];
  oozn = 1.0/(zeta+eta);
  Data->twozeta_a = 2.0*sp1->a1[pi];
  Data->twozeta_b = 2.0*sp1->a2[pj];
  Data->twozeta_c = 2.0*sp2->a1[pk];
  Data->twozeta_d = 2.0*sp2->a2[pl];
  Data->poz = eta*oozn;
  rho = zeta*Data->poz;
  coef1 = 2.0*sqrt(rho*M_1_PI)*scale*sp1->Sovlp[pi][pj]*sp2->Sovlp[pk][pl];
  PQ.x = sp1->P[pi][pj][0] - sp2->P[pk][pl][0];
  PQ.y = sp1->P[pi][pj][1] - sp2->P[pk][pl][1];
  PQ.z = sp1->P[pi][pj][2] - sp2->P[pk][pl][2];
  PQ2 = PQ.x*PQ.x;
  PQ2 += PQ.y*PQ.y;
  PQ2 += PQ.z*PQ.z;

  Data->oo2zn = 0.5*oozn;
  Data->pon = zeta*oozn;
  Data->oo2z = 0.5/zeta;
  Data->oo2n = 0.5/eta;
  W.x = (sp1->P[pi][pj][0]*zeta+sp2->P[pk][pl][0]*eta)*oozn;
  W.y = (sp1->P[pi][pj][1]*zeta+sp2->P[pk][pl][1]*eta)*oozn;
  W.z = (sp1->P[pi][pj][2]*zeta+sp2->P[pk][pl][2]*eta)*oozn;

  if(fabs(PQ2)<ZERO){ 
    for(i=0; i<=am+DERIV_LVL; i++) 
      Data->F[i] = F0[i]*coef1;
    }
  else {
    int_fjt(fjt_table,am+DERIV_LVL,rho*PQ2);
    for(i=0;i<=am+DERIV_LVL;i++)
      Data->F[i] = fjt_table->d[i]*coef1;
    }

  /* PA */
  Data->U[0][0] = sp1->PA[pi][pj][0];
  Data->U[0][1] = sp1->PA[pi][pj][1];
  Data->U[0][2] = sp1->PA[pi][pj][2];
  /* PB */
  Data->U[1][0] = sp1->PB[pi][pj][0];
  Data->U[1][1] = sp1->PB[pi][pj][1];
  Data->U[1][2] = sp1->PB[pi][pj][2];
  /* QC */
  Data->U[2][0] = sp2->PA[pk][pl][0];
  Data->U[2][1] = sp2->PA[pk][pl][1];
  Data->U[2][2] = sp2->PA[pk][pl][2];
  /* QD */
  Data->U[3][0] = sp2->PB[pk][pl][0];
  Data->U[3][1] = sp2->PB[pk][pl][1];
  Data->U[3][2] = sp2->PB[pk][pl][2];
  /* WP */
  Data->U[4][0] = W.x - sp1->P[pi][pj][0];
  Data->U[4][1] = W.y - sp1->P[pi][pj][1];
  Data->U[4][2] = W.z - sp1->P[pi][pj][2];
  /* WQ */
  Data->U[5][0] = W.x - sp2->P[pk][pl][0];
  Data->U[5][1] = W.y - sp2->P[pk][pl][1];
  Data->U[5][2] = W.z - sp2->P[pk][pl][2];

  return;
}


