/*
** Interface to the BLAS routines
**
** C. David Sherrill
** Anna I. Krylov
**
** May 1998
**
** Additions by TD Crawford and EF Valeev, June 1999.
*/

#include <stdio.h>

#if FCLINK==1
#define F_DAXPY daxpy_
#define F_DCOPY dcopy_
#define F_DGEMM dgemm_
#define F_DROT drot_
#define F_DSCAL dscal_
#define F_DGEMV dgemv_
#define F_DDOT  ddot_
#elif FCLINK==2
#define F_DAXPY daxpy
#define F_DCOPY dcopy
#define F_DGEMM dgemm
#define F_DROT drot
#define F_DSCAL dscal
#define F_DGEMV dgemv
#define F_DDOT  ddot
#else
#define F_DAXPY DAXPY
#define F_DCOPY DCOPY
#define F_DGEMM DGEMM
#define F_DROT DROT
#define F_DSCAL DSCAL
#define F_DGEMV DGEMV
#define F_DDOT  DDOT
#endif

extern void F_DAXPY(int *length, double *a, double *x, int *inc_x, 
                    double *y, int *inc_y);
extern void F_DCOPY(int *length, double *x, int *inc_x, 
                    double *y, int *inc_y);
extern void F_DGEMM(char *transa, char *transb, int *m, int *n, int *k, 
                    double *alpha, double *A, int *lda, double *B, int *ldb, 
                    double *beta, double *C, int *ldc);
extern void F_DROT(int *ntot,double *x, int *incx,double *y, int *incy,
                  double *cotheta,double *sintheta);
extern void F_DSCAL(int *n, double *alpha, double *vec, int *inc);
extern void F_DGEMV(char *transa, int *m, int *n, double *alpha, double *A, 
                    int *lda, double *X, int *inc_x, double *beta, 
                    double *Y, int *inc_y);
extern double F_DDOT(int *n, double *x, int *incx, double *y, int *incy);

/*
** C_DAXPY()
**
** This function performs y = a * x + y
** Steps every inc_x in x and every inc_y in y (normally both 1)
**
*/
void C_DAXPY(int length, double a, double *x, int inc_x, 
             double *y, int inc_y)
{
  F_DAXPY(&length, &a, x, &inc_x, y, &inc_y);
}

/*
** C_DCOPY()
**
** This function performs copies x into y
** Steps every inc_x in x and every inc_y in y (normally both 1)
**
*/
void C_DCOPY(int length, double *x, int inc_x, 
             double *y, int inc_y)
{
  F_DCOPY(&length, x, &inc_x, y, &inc_y);
}


/*
** C_DSCAL()
** 
** This function scales a vector by a real scalar
*/
void C_DSCAL(int n, double alpha, double *vec, int inc)
{
  F_DSCAL(&n, &alpha, vec, &inc);
}


/*
**
** void C_DROT(int ntot, double *x, int incx, double *y, int incy,
**             double costheta, double sintheta);
** This function calculates plane Givens rotation for vectors
** x,y and angle theta:  x=x*cos+y*sin, y=-x*sin+y*cos;
** int ntot: length of x,y
** int incx,incy: increments for x,y
*/
void C_DROT(int ntot, double *x, int incx, double *y, int incy,
            double costheta, double sintheta)
{

  F_DROT(&ntot,x,&incx,y,&incy,&costheta,&sintheta);
}


/*
**
** void C_DGEMM(char transa, char transb, int m, int n, int k, double alpha,
**            double *A, int nca, double *B, int ncb, double beta, double *C,
**            int ncc);
** This function calculates C(m,n)=alpha*(opT)A(m,k)*(opT)B(k,n)+ beta*C(m,n)
**
** Below is the Fortran explanation of arguments; I have reversed 
** everything under nca, ncb, ncc, A, B, C,  to make it correct for C.
**
** char transa:       On entry, specifies the form of (op)A used in the
**                    matrix multiplication:
**                    If transa = 'N' or 'n', (op)A = A
**                    If transa = 'T' or 't', (op)A = transp(A)
**                    If transa = 'R' or 'r', (op)A = conjugate(A)
**                    If transa = 'C' or 'c', (op)A = conjug_transp(A)
**                    On exit, transa is unchanged.
** char transb:       On entry, specifies the form of (op)B used in the
**                    matrix multiplication:
**                    If transb = 'N' or 'n', (op)B = B
**                    If transb = 'T' or 't', (op)B = transp(B)
**                    If transb = 'R' or 'r', (op)B = conjugate(B)
** int m:             On entry, the number of rows of the matrix (op)A and of
**                    the matrix C; m >= 0. On exit, m is unchanged.
** int n:             On entry, the number of columns of the matrix (op)B and
**                    of the matrix C; n >= 0. On exit, n is unchanged.
** int k:             On entry, the number of columns of the matrix (op)A and
**                    the number of rows of the matrix (op)B; k >= 0. On exit,
**                    k is unchanged.
** double alpha:      On entry, specifies the scalar alpha. On exit, alpha is
**                    unchanged.
** double *A:         On entry, a two-dimensional array A with dimensions ka
**                    by nca. For (op)A = A  or  conjugate(A), nca >= k and the
**                    leading m by k portion of the array A contains the matrix
**                    A. For (op)A = transp(A) or conjug_transp(A), nca >= m
**                    and the leading k by m part of the array A contains the
**                    matrix A. On exit, a is unchanged.
** int nca:           On entry, the second dimension of array A.
**                    For (op)A = A  or conjugate(A), nca >= MAX(1,k).
**                    For (op)A=transp(A) or conjug_transp(A), nca >= MAX(1,m).
**                    On exit, nca is unchanged.
** double *B:         On entry, a two-dimensional array B with dimensions kb
**                    by ncb. For (op)B = B or conjugate(B), kb >= k and the
**                    leading k by n portion of the array contains the matrix
**                    B. For (op)B = transp(B) or conjug_transp(B), ncb >= k and
**                    the leading n by k part of the array contains the matrix
**                    B. On exit, B is unchanged.
** int ncb:           On entry, the second dimension of array B.
**                    For (op)B = B or <conjugate(B), ncb >= MAX(1,n).
**                    For (op)B = transp(B) or conjug_transp(B), ncb >=
**                    MAX(1,k). On exit, ncb is unchanged.
** double beta:       On entry, specifies the scalar beta. On exit, beta is
**                    unchanged.
** double *C:         On entry, a two-dimensional array with the dimension
**                    at least m by ncc. On exit,  the leading  m by n part of
**                    array C is overwritten by the matrix alpha*(op)A*(op)B +
**                    beta*C.
** int ncc:           On entry, the second dimension  of array C; ncc >=MAX(1,n)
**                    On exit, ncc is unchanged.
**
*/
void C_DGEMM(char transa, char transb, int m, int n, int k, double alpha,
           double *A, int nca, double *B, int ncb, double beta, double *C,
           int ncc)
{

  /* the only strange thing we need to do is reverse everything
     since the stride runs differently in C vs. Fortran
   */
  
  /* also, do nothing if a dimension is 0 */
  if (m == 0 || n == 0 || k == 0) return;

  F_DGEMM(&transb,&transa,&n,&m,&k,&alpha,B,&ncb,A,&nca,&beta,C,&ncc);

}

/*
** void C_DGEMV(char transa, int m, int n, double alpha, double *A, 
**              int nca, double *X, int inc_x, double beta, double *Y,
**              int inc_y)
**
** This function calculates the matrix-vector product:
**
** Y = alpha * A * X + beta * Y
**
** where X and Y are vectors, A is a matrix, and alpha and beta are
** constants. 
**
*****
**
** char transa:       Indicates whether the matrix A should be
**                    transposed ('t') or left alone ('n').
** int m:             The row dimension of A (regardless of transa).
** int n:             The column dimension of A (regardless of transa).
** double alpha:      The scalar alpha.
** double *A:         A pointer to the beginning of the data in A.
** int nca:           The number of columns *actually* in A.  This is
**                    useful if one only wishes to multiply the first
**                    n columns of A times X even though A
**                    contains nca columns.
** double *X:         A pointer to the beginning of the data in X.
** int inc_x:         The desired stride for X.  Useful for skipping
**                    sections of data to treat only one column of a
**                    complete matrix.  Usually 1, though.
** double beta:       The scalar beta.
** double *Y:         A pointer to the beginning of the data in Y.
** int inc_y:         The desired stride for Y.
**
** Interface written by TD Crawford and EF Valeev.
** June 1999.
**
*/

void C_DGEMV(char transa, int m, int n, double alpha, double *A, 
             int nca, double *X, int inc_x, double beta, double *Y,
             int inc_y)
{
  if (m == 0 || n == 0) return;

  if(transa == 'n') transa = 't';
  else transa = 'n';

  F_DGEMV(&transa,&n,&m,&alpha,A,&nca,X,&inc_x,&beta,Y,&inc_y);

}
/*
** void C_DDOT(int n, double *X, int inc_x, double *Y, int inc_y) 
**            
**
** This function returns the dot product of two vectors, X and Y.
**
*****
**
** int n:             Number of elements in X and Y.
** 
** double *X:         A pointer to the beginning of the data in X.
**                    Must be of at least length (1+(N-1)*abs(inc_x).
**
** int inc_x:         The desired stride of X. Useful for skipping
**                    around to certain data stored in X.
**
** double *Y:         A pointer to the beginning of the data in Y.
** 
** int inc_y:         The desired stride for Y.
**
** Interface written by ST Brown.
** July 2000
*/

double C_DDOT(int n, double *X, int inc_x, double *Y, int inc_y)
{
   if(n == 0) return 0.0;

   return F_DDOT(&n,X,&inc_x,Y,&inc_y);
}

