#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<libipv1/ip_lib.h>
#include<iwl.h>
#include<libciomr.h>
#include<libint.h>
#include<psifiles.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#include"oe_osrr.h"

/*-----------------------------------------------------------------------------------
  This function computes AO overlap and moment integrals and writes them out to disk
 -----------------------------------------------------------------------------------*/
void moment_ints()
{
   struct coordinates PA, PB, AB, A, B;
   struct shell_pair *sp;
   struct unique_shell_pair *usp;
   register int i, j, k, l, ii, jj, kk, ll;
   int count;
   int si, sj;
   int np_i, np_j;
   int sz;
   int l1, l2, m1, m2, n1, n2;
   int ioffset, joffset ;
   int ij;
   int h1;
   int am;
   int dimension ;
   int ni,nj,ai,aj;
   int am_i, am_j;
   double a1, a2;
   double ab2, oog, gam;
   double x0, y0, z0;
   double x1, y1, z1;
   double **stemp, **mxtemp, **mytemp, **mztemp, **temp, **tmp_dptr;
   double *S, *MX, *MY, *MZ;
   double inorm, jnorm, over_pf;
   double *ptr1, *ptr2, norm1, norm12;
   double **OIX, **OIY, **OIZ;


  /*--- allocate room for the one-e matrices ---*/
  dimension = ioff[BasisSet.num_ao];
  S = init_array(dimension);
  MX = init_array(dimension);
  MY = init_array(dimension);
  MZ = init_array(dimension);

  /*--- allocate storage for shell blocks of one electron integrals ---*/
  dimension = ioff[BasisSet.max_am];
  stemp = block_matrix(dimension,dimension);
  mxtemp = block_matrix(dimension,dimension);
  mytemp = block_matrix(dimension,dimension);
  mztemp = block_matrix(dimension,dimension);

  OIX = block_matrix(BasisSet.max_am+1,BasisSet.max_am+1);
  OIY = block_matrix(BasisSet.max_am+1,BasisSet.max_am+1);
  OIZ = block_matrix(BasisSet.max_am+1,BasisSet.max_am+1);
  
  for (si=0; si<BasisSet.num_shells; si++){
    am_i = BasisSet.shells[si].am-1;
    ni = ioff[BasisSet.shells[si].am];
    A = Molecule.centers[BasisSet.shells[si].center-1];
    ioffset = BasisSet.shells[si].fao - 1;
    for (sj=0; sj<=si; sj++){
      nj = ioff[BasisSet.shells[sj].am];
      am_j = BasisSet.shells[sj].am-1;
      B = Molecule.centers[BasisSet.shells[sj].center-1];
      joffset = BasisSet.shells[sj].fao - 1;

      sp = &(BasisSet.shell_pairs[si][sj]);
      AB.x = sp->AB[0];
      AB.y = sp->AB[1];
      AB.z = sp->AB[2];
      ab2 = AB.x * AB.x;
      ab2 += AB.y * AB.y;
      ab2 += AB.z * AB.z;
	
      /*--- zero the temporary storage for accumulating contractions ---*/
      memset(stemp[0],0,sizeof(double)*dimension*dimension);
      memset(mxtemp[0],0,sizeof(double)*dimension*dimension);
      memset(mytemp[0],0,sizeof(double)*dimension*dimension);
      memset(mztemp[0],0,sizeof(double)*dimension*dimension);
      
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

	  OI_OSrecurs(OIX,OIY,OIZ,PA,PB,gam,am_i+1,am_j+1);

	  /*--- create all am components of si ---*/
	  ai = 0;
	  for(ii = 0; ii <= am_i; ii++){
	    l1 = am_i - ii;
	    for(jj = 0; jj <= ii; jj++){
	      m1 = ii - jj;
	      n1 = jj;
	      /*--- create all am components of sj ---*/
	      aj = 0;
	      for(kk = 0; kk <= am_j; kk++){
		l2 = am_j - kk;
		for(ll = 0; ll <= kk; ll++){
		  m2 = kk - ll;
		  n2 = ll;

		  x0 = OIX[l1][l2];   y0 = OIY[m1][m2];    z0 = OIZ[n1][n2];
		  x1 = OIX[l1][l2+1]; y1 = OIY[m1][m2+1];  z1 = OIZ[n1][n2+1];
		  stemp[ai][aj] += over_pf*x0*y0*z0;
		  /*--- electrons have negative charge ---*/
		  mxtemp[ai][aj] -= over_pf*(x1+x0*B.x)*y0*z0;
		  mytemp[ai][aj] -= over_pf*x0*(y1+y0*B.y)*z0;
		  mztemp[ai][aj] -= over_pf*x0*y0*(z1+z0*B.z);
		  
		  aj++;
		}
	      }
	      ai++;
	    }
	  } /*--- end cartesian components for (si,sj) with primitives (i,j) ---*/
	}
      } /*--- end primitive contraction ---*/
    
      /*--- Normalize the contracted integrals ---*/
      ptr1 = GTOs.bf_norm[am_i];
      ptr2 = GTOs.bf_norm[am_j];
      for(i=0; i<ni; i++) {
	norm1 = ptr1[i];
	for(j=0; j<nj; j++) {
	  norm12 = norm1*ptr2[j];
	  stemp[i][j] *= norm12;
	  mxtemp[i][j] *= norm12;
	  mytemp[i][j] *= norm12;
	  mztemp[i][j] *= norm12;
	}
      }

      for(i=0;i<ni;i++)
	for(j=0;j<nj;j++) {
	  ij = INDEX(ioffset + i, joffset + j);
	  S[ij] = stemp[i][j];
	  MX[ij] = mxtemp[i][j];
	  MY[ij] = mytemp[i][j];
	  MZ[ij] = mztemp[i][j];
	}
    }
  }  /*--- This shell pair is done ---*/

  /*--- flush it all away ---*/
  fflush(outfile);
  free_block(OIX);
  free_block(OIY);
  free_block(OIZ);
  dimension=ioff[BasisSet.num_ao];
  iwl_wrtone(IOUnits.itapS_AO, PSIF_AO_S, dimension,S);
  iwl_wrtone(IOUnits.itapMX_AO,PSIF_AO_MX,dimension,MX);
  iwl_wrtone(IOUnits.itapMY_AO,PSIF_AO_MY,dimension,MY);
  iwl_wrtone(IOUnits.itapMX_AO,PSIF_AO_MZ,dimension,MZ);
  if (UserOptions.print_lvl >= PRINT_OEI) {
    fprintf(outfile,"  -Overlap AO integrals:\n\n");
    print_array(S,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -mu(x) AO integrals:\n\n");
    print_array(MX,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -mu(y) AO integrals:\n\n");
    print_array(MY,BasisSet.num_ao,outfile);
    fprintf(outfile,"  -mu(z) AO integrals:\n\n");
    print_array(MZ,BasisSet.num_ao,outfile);
    fprintf(outfile,"\n");
  }

  free(S);
  free(MX);
  free(MY);
  free(MZ);

  return;
}   

