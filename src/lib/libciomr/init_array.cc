/*
** \file
** \brief Initialize an array of doubles
** \ingroup CIOMR
*/

#include <psifiles.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <strings.h>

extern "C" {

/*!
** init_array(): This function initializes an array of doubles of
** length 'size' and returns a pointer to the first element
**
** \param size = length of array (unsigned long to allow large arrays)
**
** Returns: pointer to new array
**
** \ingroup CIOMR
*/
double * init_array(unsigned long int size)
{
  double *array;

  if ((array = (double *) malloc(size*(unsigned long int)sizeof(double)))
    == NULL) {
    fprintf(stderr,"init_array: trouble allocating memory \n");
    fprintf(stderr,"size = %ld\n",size);
    exit(PSI_RETURN_FAILURE);
  }
  memset(array,'\0',size*(unsigned long int)sizeof(double));
  //bzero(array,size*(unsigned long int)sizeof(double));
  return(array);
}

} /* extern "C" */
