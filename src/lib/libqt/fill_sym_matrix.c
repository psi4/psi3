/*
** FILL_SYM_MATRIX(): This function fills a symmetric matrix by
**    placing the elements of the lower triangle into the upper
**    triangle.
**
** Arguments:
**    A    = matrix to symmetrize
**    size = number of rows or columns in matrix (assume square)
**
** Returns: 
**    none
*/
void fill_sym_matrix(double **A, int size)
{
   double **row, *col; 
   int rc, cc;

   row = A ;
   for (rc = 0; rc < (size-1); rc++) {
      col = *row;
      for (cc = 0; cc < size; cc++) {
         if (cc > rc) {
            *col = A[cc][rc] ;
            }
         col++ ;
         }
      row++ ;
      }
}

