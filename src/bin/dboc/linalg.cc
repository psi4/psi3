//
// Local definitions for matrix operations: (de)allocation, product, LU decomposition
//

#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linalg.h"

extern void done(const char * message);

FLOAT** create_matrix(int a, int b)
{
  FLOAT** M;

  if (a>=0 && b>=0) {
    M = new FLOAT*[a];
    if (M == NULL) {
      done("create_matrix failed -- probably not enough memory.");
    }
    M[0] = new FLOAT[a*b];
    if (M[0] == NULL) {
      done("create_matrix failed -- probably not enough memory.");
    }
    for(int i=1; i<a; i++)
      M[i] = M[i-1] + b;
  }
  else
    M = NULL;

  return M;
}

void delete_matrix(FLOAT** M)
{
  if (M) {
    delete[] M[0];
    delete[] M;
    M = NULL;
  }
}

void print_mat(FLOAT** a, int m, int n, FILE* out)
{
  int ii,jj,kk,nn,ll;
  int i,j,k;

  ii=0;jj=0;
L200:
  ii++;
  jj++;
  kk=10*jj;
  nn=n;
  if (nn > kk) nn=kk;
  ll = 2*(nn-ii+1)+1;
  fprintf (out,"\n");
  for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
  fprintf (out,"\n");
  for (i=0; i < m; i++) {
    fprintf (out,"\n%5d",i+1);
    for (j=ii-1; j < nn; j++) {
#if LONG_DOUBLE
      fprintf (out,"%12.7Lf",a[i][j]);
#else
      fprintf (out,"%12.7lf",a[i][j]);
#endif
    }
  }
  fprintf (out,"\n");
  if (n <= kk) {
    fflush(out);
    return;
  }
  ii=kk; goto L200;
}

FLOAT** convert_matrix(double **m, int a, int b, int transpose)
{
  int nrow = transpose ? b : a;
  int ncol = transpose ? a : b;

  FLOAT** M = create_matrix(nrow, ncol);
  if (!transpose) {
    int nelem = nrow*ncol;
    FLOAT* Melem = M[0];
    double* melem = m[0];
    for(int elem=0; elem<nelem; ++elem) {
      (*Melem) = (FLOAT) (*melem);
      ++Melem; ++melem;
    }
  }
  else {
    FLOAT* Melem = M[0];
    for(int row=0; row<nrow; row++)
      for(int col=0; col<ncol; col++) {
	(*Melem) = (FLOAT) m[col][row];
	++Melem;
      }
  }

  return M;
}

//
// C = A*B, return 1(fail) or 0(success)
//
int matrix_mult(FLOAT** A, int arow, int acol, FLOAT** B, int brow, int bcol, FLOAT** C)
{
  int nlink = acol;
  if (acol != brow)
    return 1;

  for(int row=0; row<arow; row++) {
    for(int col=0; col<bcol; col++) {
      FLOAT* aelem = A[row];
      FLOAT* belem = B[0] + col;
      FLOAT tmp = 0.0;
      for(int link=0; link<nlink; link++) {
	tmp += (*aelem) * (*belem);
	++aelem;
	belem += bcol;
      }
      C[row][col] = tmp;
    }
  }

  return 0;
}


#define TINY 1.0E-20

void lu_decom(FLOAT** a, int n, int* indx, FLOAT* d)
{
  int i,imax,j,k;
  FLOAT big,dum,sum,temp;
  FLOAT* vv = new FLOAT[n];

  *d = 1.0;

  for (i=0; i < n ; i++) {
    big=0.0;
    for (j=0; j < n; j++) {
      if ((temp=FABS(a[i][j])) > big) big=temp;
    }
    if (big == 0.0) {
      *d = 0.0;
      return;
    }
    vv[i] = 1.0/big;
  }
  for (j=0; j < n ; j++) {
    for (i=0; i < j ; i++) {
      sum = a[i][j];
      for (k=0; k < i ; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i=j ; i < n ; i++) {
      sum=a[i][j];
      for (k=0; k < j ; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if ((dum=vv[i]*fabs(sum)) >= big) {
	big = dum;
	imax = i;
      }
    }
    if (j != imax) {
      for (k=0; k < n; k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j] = TINY;
    if (j != n-1) {
      dum = 1.0/a[j][j];
      for (i=j+1; i < n ; i++) a[i][j] *= dum;
    }
  }
  delete[] vv;
}
