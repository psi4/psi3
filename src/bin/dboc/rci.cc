#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
extern "C" {
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libqt/slaterd.h>
#include <psifiles.h>
}
#include "moinfo.h"
#include "float.h"
#include "linalg.h"
#include "mo_overlap.h"

using namespace std;

// Wrap a,b indices into one composite index assuming S2 symmetry
#define INDEX2(a,b) ((a) > (b)) ? ( (((a)*(a+1)) >> 1) + (b) ) : ( (((b)*(b+1)) >> 1) + (a) )
// Wrap a>=b indices into one composite index assuming S2 symmetry
#define INDEX2_ORD(a,b) ( (((a)*(a+1))/2) + (b) )
// Wrap a>=b>=c indices into one composite index assuming S3 symmetry
#define INDEX3_ORD(a,b,c) ( ((a)*(((a)+4)*((a)-1)+6)/6) + (((b)*(b+1))/2) + (c) )

extern MOInfo_t MOInfo;
extern FILE *outfile;
extern void done(const char *);

double eval_rci_derwfn_overlap()
{
  int ndocc = MOInfo.ndocc;
  FLOAT **CSC_full = eval_S_alpha();
  FLOAT **CSC = create_matrix(ndocc,ndocc);
  int *tmpintvec = new int[ndocc];

  // Read in CI vectors
  SlaterDetVector *vecm, *vecp;
  slaterdetvector_read(PSIF_CIVECT,"Old CI vector",&vecm);
  slaterdetvector_read(PSIF_CIVECT,"CI vector",&vecp);

  // Compute overlap between strings for alpha spin case (beta is the same)
  StringSet *ssetm;
  ssetm = vecm->sdset->alphastrings;
  int nstr_a = ssetm->size;
  int nalpha = ndocc;
  FLOAT **S_a = create_matrix(nstr_a,nstr_a);
  // Assume the order of strings is the same for - and + displacements
  for(int jp=0; jp<nstr_a; jp++) {
    String *str_j = &ssetm->strings[jp];
    for(int im=0; im<nstr_a; im++) {
      String *str_i = &ssetm->strings[im];

      for(int j=0;j<ndocc;j++)
	for(int i=0;i<ndocc;i++)
	  CSC[j][i] = CSC_full[str_j->occ[j]][str_i->occ[i]];

      // Compute the determinant
      //      C_DGETRF(ndocc,ndocc,&(CSC[0][0]),ndocc,tmpintvec);
      FLOAT sign;
      lu_decom(CSC, ndocc, tmpintvec, &sign);
      FLOAT deter1 = 1.0;
      for(int i=0;i<ndocc;i++)
	deter1 *= CSC[i][i];

      S_a[jp][im] = sign*deter1;
    }
  }

  // Evaluate total overlap in the highest available precision
  int ndets = vecm->size;
  FLOAT S_tot = 0.0;
  for(int I=0; I<ndets; I++) {
    SlaterDet *detI = vecm->sdset->dets + I;
    int Istra = detI->alphastring;
    int Istrb = detI->betastring;
    FLOAT cI = vecm->coeffs[I];

    for(int J=0; J<ndets; J++) {
      SlaterDet *detJ = vecp->sdset->dets + J;
      int Jstra = detJ->alphastring;
      int Jstrb = detJ->betastring;
      FLOAT cJ = vecp->coeffs[J];
      
      FLOAT S = S_a[Jstra][Istra] * S_a[Jstrb][Istrb];

      FLOAT contrib = cI * S * cJ;
      S_tot += cI * S * cJ;
      /* fprintf(outfile,"  %3d %3d %+15.10Le", I, J, cI);
      fprintf(outfile," %+15.10Le", cJ);
      fprintf(outfile," %+25.15Le", S);
      fprintf(outfile," %+25.15Le", contrib);
      fprintf(outfile," %+25.15Le\n", S_tot); */

    }
  }
  

  // Cleanup
  slaterdetvector_delete_full(vecm);
  slaterdetvector_delete_full(vecp);
  delete[] tmpintvec;
  delete_matrix(CSC);
  delete_matrix(CSC_full);
  delete_matrix(S_a);
  double S_tot_double = (double) S_tot;
  return fabs(S_tot_double);
}

