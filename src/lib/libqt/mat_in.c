#include <stdio.h>

/*
** MAT_IN(): Function to read in a matrix.  Simple version for now.
**
** Parameters:
**    fp         =  file pointer to input stream
**    array      =  matrix to hold data
**    width      =  number of columns to read
**    max_length = maximum number of rows to read
**    stat       = pointer to int to hold status flag (0=read ok, 1=error)
**
** Returns: 
**    number of rows read
**    Also modifies stat to = error code (0 = ok, 1 = error)
*/

int mat_in(FILE *fp, double **array, int width, int max_length, int *stat) 
{
   int i=0, j, errcod=0 ;
   int nr ;
   double data ;

   while ( (i < max_length) && (!errcod) ) {
      for (j=0; j<width; j++) {
         nr = fscanf(fp, "%lf", &data) ;
         if (feof(fp)) break ;
         if (nr != 1) {
            errcod = 1 ;
            break ;
            }
         else {
            array[i][j] = data ;
            }
         }
      if (feof(fp)) break ;
      if (!errcod) i++ ;
      }

   *stat = errcod ;
   return(i) ;
}

