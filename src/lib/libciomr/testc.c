#include "includes.h"
#include "common.h"
#include "iomrparam.h"

extern double *init_array(int);
extern double **init_matrix(int,int);

#define RSP 0
#define MXM 1

int main()
   {
      int i,j;
      double k=1.0/100000;
#if RSP
      double *arr, **vecs;
      double *vals;
#else
      double **arr, **vecs, **res;
#endif
      double tol = 10.0e-14;
      int n = 500;
      int m = 5;
      int nn;

      nn = n*(n+1)/2;

      ffile(&outfile,"outc",0);
      tstart(outfile);

#if RSP
      arr = (double *) init_array(nn);
      vals = (double *) init_array(nn);
#else
      arr = (double **) init_matrix(n,n);
      res = (double **) init_matrix(n,n);
#endif
      vecs = (double **) init_matrix(n,n);

#if RSP
      for (i=0; i < nn; i++) {
         arr[i] = (double) i + 1.0;
         }
#endif

#if MXM
      for (i=0; i<n ; i++ ) {
         for (j=0; j <= i; j++) {
            arr[i][j] = k;
            arr[j][i] = k;
            k = (k*100000.0 + 1.0)/100000.0;
            vecs[i][j] = vecs[j][i] = k;
            }
         }
#endif
         
#if RSP
      resource_command();
      rsp(n,n,nn,arr,vals,1,vecs,tol);
      resource_command();
      eivout(vecs,vals,n,5,outfile);
#endif

#if MXM
      resource_command();
      mmult(arr,0,vecs,1,res,0,n,n,n,0);
      resource_command();
      print_mat(res,m,m,outfile); 
#endif

      tstop(outfile);
      }
