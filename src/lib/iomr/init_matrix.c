#include <stdio.h>
#if 0
# include <malloc.h>
#endif

/* allocates memory for an n x m matrix */
/* returns pointer to pointer to 1st element */

double ** init_matrix(n,m)
   int n,m;

   {
    double **array=NULL;
    int i;
    char *malloc();

    if ((array = (double **) malloc(sizeof(double *)*n))==NULL) {
         fprintf(stderr,"init_matrix: trouble allocating memory \n");
         fprintf(stderr,"want %d double words of core\n",n*m);
         exit(2);
         }

    for (i = 0; i < n; i++) {
        if ((array[i] = (double *) malloc(sizeof(double)*m))==NULL) {
           fprintf(stderr,"init_matrix: trouble allocating memory \n");
           fprintf(stderr,"want %d double words of core\n",n*m);
           exit(3);
           }
        }
    return(array);
    }

void free_matrix(array,size)
   double **array;
   int size;

   {
      int i;

      for (i=0; i < size ; i++) {
         free(array[i]);
         }

      free(array);
      }
