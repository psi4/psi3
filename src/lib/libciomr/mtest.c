#include "includes.h"
#include "common.h"
#include "iomrparam.h"

extern double **init_matrix(int, int);
extern void resource_command(void);
extern void mmult(double **AF, int ta, double **BF, int tb, double **CF, int tc,
		  int nr, int nl, int nc, int add);
extern void print_mat(double **, int, int, FILE*);

main()
   {
      int i,j,l;
      double k=1.0/100000;
      double t,t1;
      double **arr, **vecs, **res;
      double *ai,*bj;
      int n = 400;
      int m = 5;
      char name[MAX_STRING];

      strcat(name,"outc");
      ffile(&outfile,name,0);


      arr = (double **) init_matrix(n,n);
      res = (double **) init_matrix(n,n);
      vecs = (double **) init_matrix(n,n);


      for (i=0; i<n ; i++ ) {
         for (j=0; j <= i; j++) {
            arr[i][j] = k;
            arr[j][i] = k;
            k = (k*100000.0 + 1.0)/100000.0;
            vecs[i][j] = vecs[j][i] = k;
            }
         }

      resource_command();
      for(j=0; j < n ; j++) {
         for(i=0; i < n ; i++) {
            for(l=0,t=0.0; l < n ; l++) 
               t += arr[i][l]*vecs[l][j];
            res[i][j] = t;
            }
         }
         
      resource_command();
      mmult(arr,0,vecs,0,res,0,n,n,n,0);
      resource_command();
      print_mat(res,m,m,outfile); 
      }
