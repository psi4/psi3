/*!
  \file dirprd_block.c
  \ingroup (QT)
*/

/*!

  dirprd_block()

  This function takes two block matrices A and B and multiplies
  each element of B by the corresponding element of A

  \param **A: block matrix A
  \param **B: block matrix B 
  \param nrows: number of rows of A and B
  \param ncols: number of columns of A and B

  \ingroup (QT)
*/
void dirprd_block(double **A, double **B, int rows, int cols)
{
  register int i;
  double *a, *b;

  if(!(rows*cols)) return;

  a = A[0]; b= B[0];

  for(i=0; i < rows*cols; i++, a++, b++) (*b) = (*a) * (*b);
}
