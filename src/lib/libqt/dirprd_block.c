void dirprd_block(double **A, double **B, int rows, int cols)
{
  register int i;
  double *a, *b;

  if(!(rows*cols)) return;

  a = A[0]; b= B[0];

  for(i=0; i < rows*cols; i++, a++, b++) (*b) = (*a) * (*b);
}
