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
