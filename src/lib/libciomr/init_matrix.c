/*!
  \file init_matrix.c
  \ingroup (CIOMR)
*/ 

/* $Log$
 * Revision 1.5  2002/06/01 18:23:54  sherrill
 * Upgrade doxygen documentation
 *
/* Revision 1.4  2002/04/19 21:48:06  sherrill
/* Remove some unused functions and do doxygen markup of libciomr.
/*
/* Revision 1.3  2002/04/04 22:24:38  evaleev
/* Converted allocation functions (init_array, etc.) to take unsigned long ints
/* to be able to allocate properly 2GB+ chunks). Some declarations cleaned up.
/*
/* Revision 1.2  2002/03/25 02:43:45  sherrill
/* Update documentation
/*
/* Revision 1.1.1.1  2000/02/04 22:53:19  evaleev
/* Started PSI 3 repository
/*
/* Revision 2.2  1999/11/01 20:10:56  evaleev
/* Added explicit extern declarations of functions within the library.
/*
/* Revision 2.1  1991/06/15 18:29:20  seidl
/* *** empty log message ***
/* */

static char *rcsid = "$Id$";

#include "includes.h"

extern void resource_command(void);

/* allocates memory for an n x m matrix */
/* returns pointer to pointer to 1st element */

/*!
** init_matrix(): Initialize a matrix of doubles and return a pointer to
** the first row.  Note that this does not form a matrix which is 
** necessarily contiguous in memory.  Use block_matrix() for that.
** 
** \param n = number of rows
** \param m = number of columns
** \ingroup (CIOMR)
*/

double ** init_matrix(unsigned long int n, unsigned long int m)
   {
    double **array=NULL;
    unsigned long int i;

    if ((array = (double **) malloc(n*(unsigned long int)sizeof(double *)))
        ==NULL) {
         fprintf(stderr,"init_matrix: trouble allocating memory \n");
         fprintf(stderr,"n = %ld\n",n);
         resource_command();
         exit(2);
         }

    for (i = 0; i < n; i++) {
        if ((array[i] = (double *) malloc(m*(unsigned long int)sizeof(double)))
            ==NULL) {
           fprintf(stderr,"init_matrix: trouble allocating memory \n");
           fprintf(stderr,"i = %ld m = %ld\n",i,m);
           resource_command();
           exit(3);
           }
        bzero(array[i],m*(unsigned long int)sizeof(double));
        }
    return(array);
    }

/*!
** free_matrix(): Free a 2D matrix allocated with init_matrix().
**
** \param array = matrix to free
** \param size = number of rows
** \ingroup (CIOMR)
*/
void free_matrix(double **array, unsigned long int size)
   {
      unsigned long int i;

      for (i=0; i < size ; i++) {
         free(array[i]);
         }

      free(array);
      }
