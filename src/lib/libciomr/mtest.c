
/* $Log$
 * Revision 1.1  2000/02/04 22:53:21  evaleev
 * Initial revision
 *
/* Revision 2.3  1999/11/01 20:10:57  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.2  1997/06/23 12:25:52  crawdad
/* Multiple changes to libciomr: Moved "param.h" to "iomrparam.h" to avoid
/*     conflicts with similarly named system file under linux.  Corrected type
/*    casting in rread(), rwrit(), sread(), and swrit() functions.  Corrected
/*     unclosed tmpdisks.dat file in sequential.c.  Corrected block_matrix() to
/*    avoid malloc'ing zero-length arrays.
/*
/* -Daniel
/*
 * Revision 2.1  1991/06/15  18:29:30  seidl
 * *** empty log message ***
 * */

static char *rcsid = "$Id$";

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
