/*!
  \file dot_block.c
  \ingroup (QT)
*/

/*!
** dot_block()
** This function takes two block matrices A and B and finds
** the dot product.
**
**  \param double **A: block matrix A
**  \param double **B: block matrix B 
**  \param int nrows: number of rows of A and B
**  \param int ncols: number of columns of A and B
**  \param double alpha: scale factor by which the dot product is multiplied
** \ingroup (QT)
*/
double dot_block(double **A, double **B, int rows, int cols, double alpha)
{
  register int i;
  double *a, *b;
  double value;

  if(!(rows*cols)) return 0.0;

  a = A[0]; b = B[0];

  value = 0.0;
  for(i=0; i < rows*cols; i++,a++,b++) value += (*a) * (*b);

  return alpha*value;
}

