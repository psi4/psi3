#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libciomr.h>
#include <libint.h>
#include <libderiv.h>

#include "defines.h"
#define EXTERN
#include "global.h"
#include "oe_osrr.h"
#include "oe_deriv1_osrr.h"
#include "small_fns.h"

/*--- These frequently used numbers are to avoid costs of passing parameters ---*/
static double oo2g, oog, gam;

/*-------------------------------
  Explicit function declarations
 -------------------------------*/
static double overlap_int(double a1, int l1, int m1, int n1, double norm1,
			  double a2, int l2, int m2, int n2, double norm2,
			  struct coordinates AB,
			  struct coordinates PA,
			  struct coordinates PB);
static double ke_int(double a1, int l1, int m1, int n1, double norm1,
		     double a2, int l2, int m2, int n2, double norm2,
		     struct coordinates AB,
		     struct coordinates PA,
		     struct coordinates PB);
static double f_n(int k, int l1, int l2, double A, double B);
static double int_pow(double a, int p);
static double ***init_box(int a, int b, int c);
static void free_box(double ***box, int a, int b);

/*-------------------------------------------------------------
  This function computes derivatives of one-electron integrals
 -------------------------------------------------------------*/
void oe_deriv1()
{
   struct coordinates PA, PB, AB, PC;
   struct shell_pair *sp;
   int i, j, k, l, ii, jj, kk, ll;
   int count;
   int si, sj;
   int np_i, np_j;
   int si_fao, sj_fao;
   int sz;
   int l1, l2, m1, m2, n1, n2;
   int ioffset, joffset ;
   int ij;
   int ijpack;
   int h1;
   int am;
   int dimension ;
   int ni,li,nj,lj,ai,aj;
   int am_i, am_j;
   int ixm, iym, izm, jxm, jym, jzm;
   int ixm1, iym1, izm1, jxm1, jym1, jzm1;
   int indmax, iind,jind;
   int atom, atom1, atom2;
   int bf;
   int iimax, jjmax;
   double a1, a2;
   double ab2;
   double tmp;
   double inorm, jnorm, over_pf;
   double dens_pf, wdens_pf;
   double *ptr1, *ptr2, norm1, norm12;
   double ***AI0;
   double ***AIX, ***AIY, ***AIZ;
   double **grad_oe, **grad_ov;

  indmax = (BasisSet.max_am+DERIV_LVL-1)*(BasisSet.max_am+DERIV_LVL)*(BasisSet.max_am+DERIV_LVL)+1;
  AI0 = init_box(indmax,indmax,2*(BasisSet.max_am+DERIV_LVL)+1);
  AIX = init_box(indmax,indmax,2*(BasisSet.max_am+DERIV_LVL)+1);
  AIY = init_box(indmax,indmax,2*(BasisSet.max_am+DERIV_LVL)+1);
  AIZ = init_box(indmax,indmax,2*(BasisSet.max_am+DERIV_LVL)+1);
  grad_oe = block_matrix(Molecule.num_atoms,3);
  grad_ov = block_matrix(Molecule.num_atoms,3);
  

  for (si=0; si<BasisSet.num_shells; si++){
    am_i = BasisSet.shells[si].am-1;
    izm = 1;
    iym = am_i+1;
    ixm = iym*iym;
    izm1 = 1;
    iym1 = am_i+DERIV_LVL+1;
    ixm1 = iym1*iym1;
    atom1 = BasisSet.shells[si].center-1;
    si_fao = BasisSet.shells[si].fao-1;
    for (sj=0; sj<=si; sj++){
      ni = ioff[BasisSet.shells[si].am];
      nj = ioff[BasisSet.shells[sj].am];
      am_j = BasisSet.shells[sj].am-1;
      jzm = 1;
      jym = am_j+1;
      jxm = jym*jym;
      jzm1 = 1;
      jym1 = am_j+DERIV_LVL+1;
      jxm1 = jym1*jym1;
      atom2 = BasisSet.shells[sj].center-1;
      sj_fao = BasisSet.shells[sj].fao - 1;
	
      sp = &(BasisSet.shell_pairs[si][sj]);
      AB.x = sp->AB[0];
      AB.y = sp->AB[1];
      AB.z = sp->AB[2];
      ab2 = AB.x * AB.x;
      ab2 += AB.y * AB.y;
      ab2 += AB.z * AB.z;
	
	/*--- contract by primitives here ---*/
	for (i = 0; i < BasisSet.shells[si].n_prims; i++) {
	  a1 = sp->a1[i];
	  inorm = sp->inorm[i];
	  for (j = 0; j < BasisSet.shells[sj].n_prims; j++) {
	    a2 = sp->a2[j];
	    gam = sp->gamma[i][j];
	    jnorm = sp->jnorm[j];
	    PA.x = sp->PA[i][j][0];
	    PA.y = sp->PA[i][j][1];
	    PA.z = sp->PA[i][j][2];
	    PB.x = sp->PB[i][j][0];
	    PB.y = sp->PB[i][j][1];
	    PB.z = sp->PB[i][j][2];
	    oog = 1.0/gam;
	    over_pf = exp(-a1*a2*ab2*oog)*sqrt(M_PI*oog)*M_PI*oog*inorm*jnorm;

	    /*--- create all am components of si ---*/
	    ai = 0;
	    for(ii = 0; ii <= am_i; ii++){
	      l1 = am_i - ii;
	      for(jj = 0; jj <= ii; jj++){
		m1 = ii - jj;
		n1 = jj ;
		/*--- create all am components of sj ---*/
		aj = 0;
		for(kk = 0; kk <= am_j; kk++){
		  l2 = am_j - kk;
		  for(ll = 0; ll <= kk; ll++){
		    m2 = kk - ll;
		    n2 = ll;

		    if (si == sj && ai < aj)
		      break;

		    dens_pf = Dens[si_fao+ai][sj_fao+aj];
		    dens_pf *= (GTOs.bf_norm[am_i][ai] * GTOs.bf_norm[am_j][aj]);
		    wdens_pf = (-1.0)*Lagr[si_fao+ai][sj_fao+aj];
		    wdens_pf *= (GTOs.bf_norm[am_i][ai] * GTOs.bf_norm[am_j][aj]);
		    if (si_fao+ai != sj_fao+aj) {
		      dens_pf *= 2.0;
		      wdens_pf *= 2.0;
		    }
		    /*--- A factor of 1/2 for correlated Lagrangians ---*/
		    if (strcmp(UserOptions.wfn,"SCF"))
		      wdens_pf *= 0.5;

		    tmp = 2.0*a1*overlap_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, 
					     n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= l1*overlap_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, 
					    n2, jnorm, AB, PA, PB);
		    grad_ov[atom1][0] += tmp*wdens_pf;

		    tmp = 2.0*a1*overlap_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, 
					     n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= m1*overlap_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, 
					    n2, jnorm, AB, PA, PB);
		    grad_ov[atom1][1] += tmp*wdens_pf;

		    tmp = 2.0*a1*overlap_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, 
					     n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= n1*overlap_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, 
					    n2, jnorm, AB, PA, PB);
		    grad_ov[atom1][2] += tmp*wdens_pf;


		    tmp = 2.0*a2*overlap_int(a1, l1, m1, n1, inorm, a2, l2+1, m2, 
					     n2, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= l2*overlap_int(a1, l1, m1, n1, inorm, a2, l2-1, m2, 
					    n2, jnorm, AB, PA, PB);
		    grad_ov[atom2][0] += tmp*wdens_pf;

		    tmp = 2.0*a2*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2+1, 
					     n2, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= m2*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2-1, 
					    n2, jnorm, AB, PA, PB);
		    grad_ov[atom2][1] += tmp*wdens_pf;

		    tmp = 2.0*a2*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
					     n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= n2*overlap_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
					    n2-1, jnorm, AB, PA, PB);
		    grad_ov[atom2][2] += tmp*wdens_pf;


		    tmp = 2.0*a1*ke_int(a1, l1+1, m1, n1, inorm, a2, l2, m2, 
					n2, jnorm, AB, PA, PB);
		    if (l1)
		      tmp -= l1*ke_int(a1, l1-1, m1, n1, inorm, a2, l2, m2, 
				       n2, jnorm, AB, PA, PB);
		    grad_oe[atom1][0] += tmp*dens_pf;

		    tmp = 2.0*a1*ke_int(a1, l1, m1+1, n1, inorm, a2, l2, m2, 
					n2, jnorm, AB, PA, PB);
		    if (m1)
		      tmp -= m1*ke_int(a1, l1, m1-1, n1, inorm, a2, l2, m2, 
				       n2, jnorm, AB, PA, PB);
		    grad_oe[atom1][1] += tmp*dens_pf;

		    tmp = 2.0*a1*ke_int(a1, l1, m1, n1+1, inorm, a2, l2, m2, 
					n2, jnorm, AB, PA, PB);
		    if (n1)
		      tmp -= n1*ke_int(a1, l1, m1, n1-1, inorm, a2, l2, m2, 
				       n2, jnorm, AB, PA, PB);
		    grad_oe[atom1][2] += tmp*dens_pf;


		    tmp = 2.0*a2*ke_int(a1, l1, m1, n1, inorm, a2, l2+1, m2, 
					n2, jnorm, AB, PA, PB);
		    if (l2)
		      tmp -= l2*ke_int(a1, l1, m1, n1, inorm, a2, l2-1, m2, 
				       n2, jnorm, AB, PA, PB);
		    grad_oe[atom2][0] += tmp*dens_pf;

		    tmp = 2.0*a2*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2+1, 
					n2, jnorm, AB, PA, PB);
		    if (m2)
		      tmp -= m2*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2-1, 
				       n2, jnorm, AB, PA, PB);
		    grad_oe[atom2][1] += tmp*dens_pf;

		    tmp = 2.0*a2*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
					n2+1, jnorm, AB, PA, PB);
		    if (n2)
		      tmp -= n2*ke_int(a1, l1, m1, n1, inorm, a2, l2, m2, 
				       n2-1, jnorm, AB, PA, PB);
		    grad_oe[atom2][2] += tmp*dens_pf;
		    
		    aj++;
		  }
		}  
		ai++;
	      }
	    } /*--- end cartesian components for (si,sj) with primitives (i,j) ---*/

      /*--- create all am components of si ---*/
	    for(atom=0;atom<Molecule.num_atoms;atom++) {
	      PC.x = sp->P[i][j][0] - Molecule.centers[atom].x;
	      PC.y = sp->P[i][j][1] - Molecule.centers[atom].y;
	      PC.z = sp->P[i][j][2] - Molecule.centers[atom].z;
	      AI_Deriv1_OSrecurs(AI0,AIX,AIY,AIZ,PA,PB,PC,gam,am_i+DERIV_LVL,am_j+DERIV_LVL);
	      ai = 0;
	      for(ii = 0; ii <= am_i; ii++){
		l1 = am_i - ii;
		for(jj = 0; jj <= ii; jj++){
		  m1 = ii - jj;
		  n1 = jj;
		  iind = n1*izm1 + m1*iym1 + l1*ixm1;
		  /*--- create all am components of sj ---*/
		  aj = 0;
		  for(kk = 0; kk <= am_j; kk++){
		    l2 = am_j - kk;
		    for(ll = 0; ll <= kk; ll++){
		      m2 = kk - ll;
		      n2 = ll ;

		      if (si == sj && ai < aj)
		      break;

		      dens_pf = Dens[si_fao+ai][sj_fao+aj];
		      dens_pf *= (GTOs.bf_norm[am_i][ai] * GTOs.bf_norm[am_j][aj]);
		      if (si_fao+ai != sj_fao+aj)
			dens_pf *= 2.0;
		      
		      jind = n2*jzm1 + m2*jym1 + l2*jxm1;

		      tmp = 2.0*a1*AI0[iind+ixm1][jind][0];
		      if (l1)
			tmp -= l1*AI0[iind-ixm1][jind][0];
		      grad_oe[atom1][0] -= tmp * Molecule.centers[atom].Z_nuc * (over_pf * dens_pf);

		      tmp = 2.0*a2*AI0[iind][jind+jxm1][0];
		      if (l2)
			tmp -= l2*AI0[iind][jind-jxm1][0];
		      grad_oe[atom2][0] -= tmp * Molecule.centers[atom].Z_nuc * (over_pf * dens_pf);

		      grad_oe[atom][0] -= AIX[iind][jind][0] * Molecule.centers[atom].Z_nuc * (over_pf * dens_pf);

		      tmp = 2.0*a1*AI0[iind+iym1][jind][0];
		      if (m1)
			tmp -= m1*AI0[iind-iym1][jind][0];
		      grad_oe[atom1][1] -= tmp * Molecule.centers[atom].Z_nuc * (over_pf * dens_pf);

		      tmp = 2.0*a2*AI0[iind][jind+jym1][0];
		      if (m2)
			tmp -= m2*AI0[iind][jind-jym1][0];
		      grad_oe[atom2][1] -= tmp * Molecule.centers[atom].Z_nuc * (over_pf * dens_pf);

		      grad_oe[atom][1] -= AIY[iind][jind][0] * Molecule.centers[atom].Z_nuc * (over_pf * dens_pf);

		      tmp = 2.0*a1*AI0[iind+izm1][jind][0];
		      if (n1)
			tmp -= n1*AI0[iind-izm1][jind][0];
		      grad_oe[atom1][2] -= tmp * Molecule.centers[atom].Z_nuc * (over_pf * dens_pf);

		      tmp = 2.0*a2*AI0[iind][jind+jzm1][0];
		      if (n2)
			tmp -= n2*AI0[iind][jind-jzm1][0];
		      grad_oe[atom2][2] -= tmp * Molecule.centers[atom].Z_nuc * (over_pf * dens_pf);

		      grad_oe[atom][2] -= AIZ[iind][jind][0] * Molecule.centers[atom].Z_nuc * (over_pf * dens_pf);
		      

		      aj++;
		    }
		  }  
		  ai++;
		}
	      } /*--- end cartesian components for (si,sj) with primitives (i,j) ---*/
	    }
	  }
	} /*--- end primitive contraction ---*/
    

    }
  }

  if (UserOptions.print_lvl >= PRINT_OEDERIV) {
    print_atomvec("One-electron contribution to the forces (a.u.)",grad_oe);
    print_atomvec("Overlap contribution to the forces (a.u.)",grad_ov);
  }

  for(i=0;i<Molecule.num_atoms;i++) {
    Grad[i][0] += grad_oe[i][0] + grad_ov[i][0];
    Grad[i][1] += grad_oe[i][1] + grad_ov[i][1];
    Grad[i][2] += grad_oe[i][2] + grad_ov[i][2];
  }
  
  /*--- flush it all away ---*/
  fflush(outfile);
  free_block(grad_oe);
  free_block(grad_ov);
  free_box(AI0,indmax,indmax);
  free_box(AIX,indmax,indmax);
  free_box(AIY,indmax,indmax);
  free_box(AIZ,indmax,indmax);
}   


/*-----------------------------------------------
  This computes overlap of 2 primitive gaussians
 -----------------------------------------------*/
double overlap_int(double a1, int l1, int m1, int n1, double norm1,
		   double a2, int l2, int m2, int n2, double norm2,
		   struct coordinates AB,
		   struct coordinates PA,
		   struct coordinates PB)
{
  int i, j, k, l;
  int imax, jmax, kmax;
  int print = 0;
  double Ix, Iy, Iz;
  double I;
  double gam;
  double AB2;
  double tval, tval1, tval2 ;
  double norm_fact ;

  
  AB2 = AB.x*AB.x+AB.y*AB.y+AB.z*AB.z;
  gam = a1+a2;
  norm_fact = norm1*norm2; 

  tval1 = 2*gam;
  imax = (l1+l2)/2;
  Ix = 0.0;
  for (i = 0; i<= imax; i++){
    tval = f_n(i*2, l1, l2, PA.x, PB.x);
    tval2 = int_pow(tval1, i);
    Ix += tval*(num_ser[i])/(tval2);
    }
  jmax = (m1+m2)/2;
  Iy = 0.0;
  for (j = 0; j<= jmax; j++){
    tval = f_n(j*2, m1, m2, PA.y, PB.y);
    tval2 = int_pow(tval1, j);
    Iy += tval*num_ser[j]/(tval2);
    }
  kmax = (n1+n2)/2;
  Iz = 0.0;
  for (k = 0; k<= kmax; k++){
    tval = f_n(k*2, n1, n2, PA.z, PB.z);
    tval2 = int_pow(tval1, k);
    Iz += tval*num_ser[k]/(tval2);
    }
 
  I=exp(-1*a1*a2*AB2/gam)*Ix*Iy*Iz*sqrt(M_PI/gam)*(M_PI/gam);

  return I*norm_fact;
}


/*-----------------------------------------------------------------------------
  This computes matrix element of kinetic energy between 2 primitive gaussians
 -----------------------------------------------------------------------------*/
double ke_int(double a1, int l1, int m1, int n1, double norm1,
	      double a2, int l2, int m2, int n2, double norm2,
	      struct coordinates AB,
	      struct coordinates PA,
	      struct coordinates PB)
{
  int localprint = 0 ;
  double Ix, Iy, Iz;
  double I1 = 0.0, I2 = 0.0, I3 = 0.0, I4 = 0.0 ;
  double gam;
  double norm_fact;

  I2 = overlap_int(a1, l1+1, m1, n1, norm1, a2, l2+1, m2, n2, norm2, 
                    AB, PA, PB);
  I1 = overlap_int(a1, l1-1, m1, n1, norm1, a2, l2-1, m2, n2, norm2, 
                    AB, PA, PB);
  I3 = overlap_int(a1, l1+1, m1, n1, norm1, a2, l2-1, m2, n2, norm2, 
                    AB, PA, PB);
  I4 = overlap_int(a1, l1-1, m1, n1, norm1, a2, l2+1, m2, n2, norm2, 
                    AB, PA, PB);
 
  Ix = (l1*l2*I1/2.0 + 2*a1*a2*I2 - a1*l2*I3 - a2*l1*I4);

  I1 = 0.0;
  I2 = 0.0;
  I3 = 0.0;
  I4 = 0.0;

  I2 = overlap_int(a1, l1, m1+1, n1, norm1, a2, l2, m2+1, n2, norm2, 
                    AB, PA, PB);
  I1 = overlap_int(a1, l1, m1-1, n1, norm1, a2, l2, m2-1, n2, norm2, 
                    AB, PA, PB);
  I3 = overlap_int(a1, l1, m1+1, n1, norm1, a2, l2, m2-1, n2, norm2, 
                    AB, PA, PB);
  I4 = overlap_int(a1, l1, m1-1, n1, norm1, a2, l2, m2+1, n2, norm2, 
                    AB, PA, PB);
  
  Iy = (m1*m2*I1/2.0 + 2*a1*a2*I2 - a1*m2*I3 - a2*m1*I4);

  I1 = 0.0;
  I2 = 0.0;
  I3 = 0.0;
  I4 = 0.0;

  I2 = overlap_int(a1, l1, m1, n1+1, norm1, a2, l2, m2, n2+1, norm2, 
                    AB, PA, PB);
  I1 = overlap_int(a1, l1, m1, n1-1, norm1, a2, l2, m2, n2-1, norm2, 
                    AB, PA, PB);
  I3 = overlap_int(a1, l1, m1, n1+1, norm1, a2, l2, m2, n2-1, norm2, 
                    AB, PA, PB);
  I4 = overlap_int(a1, l1, m1, n1-1, norm1, a2, l2, m2, n2+1, norm2, 
                    AB, PA, PB);
  
  Iz = (n1*n2*I1/2.0 + 2*a1*a2*I2 - a1*n2*I3 - a2*n1*I4);

  return((Ix+Iy+Iz));

}

double f_n(int k, int l1, int l2, double A, double B)
{
  int i, j;
  double sum = 0.0;
  double tval;
  double tmp1, tmp2, tmp3, tmp4 ;

  for(i=0; i<=l1; i++){
    for(j=0; j<=l2; j++){
      tmp1 = tmp2 = tmp3 = tmp4 = 0.0 ;
      if((i+j) == k){
        tmp1 = int_pow(A, (l1-i));
        tmp2 = int_pow(B, (l2-j));
        tmp3 = bc[l1][i];
        tmp4 = bc[l2][j];
        tval = tmp1 * tmp2 * tmp3 * tmp4 ;
        sum += tval ;
        }
      }
    }
  return sum;
}

double int_pow(double a, int p)
{
  register int i;
  double b = 1.0;
  
  for(i=0; i<p; i++) b = b*a;
  return b;
}


double ***init_box(int a, int b, int c)
{
  int i,j,k;
  double ***box;

  box = (double ***) malloc(sizeof(double **)*a);
  for(i=0;i<a;i++)
    box[i] = (double **) malloc(sizeof(double *)*b);
  for(i=0;i<a;i++)
    for(j=0;j<b;j++) {
	box[i][j] = (double *) malloc(sizeof(double)*c);
	bzero((char *) box[i][j],sizeof(double)*c);
    }

  return box;

}


void free_box(double ***box, int a, int b)
{
  int i,j;

  for(i=0;i<a;i++)
    for(j=0;j<b;j++)
      free(box[i][j]);

  for(i=0;i<a;i++)
    free(box[i]); 

  free(box);

  return;
}

