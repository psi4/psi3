/*!
** \file ffile.c
** \ingroup (CIOMR)
*/ 

#include <stdio.h>
#include <psifiles.h>
extern char *psi_file_prefix;

/*!
** ffile(): Open a PSI3 ASCII file for reading/writing.  Returns a
** pointer to the new file.
**
** \param suffix = name of the file, not including automatic prefix
** \param code = 0 (write), 1 (write/append), 2 (read)
** \ingroup (CIOMR)
*/
FILE *ffile(char *suffix, int code)
{
  FILE *fptr;
  char name[100];

  /* build the standard file name */
  sprintf(name, "%s.%s", psi_file_prefix, suffix);

  switch (code) {
  case 0:
    fptr = fopen(name,"w+");
    break;
  case 1:
    fptr = fopen(name,"a+");
    break;
  case 2:
    fptr = fopen(name,"r+");
    break;
  default:
    fprintf(stderr,"error in ffile: invalid code %d\n",code);
  }
  if (fptr == NULL) { 
    fprintf(stderr,"error in ffile: cannot open file %s\n", suffix);
    exit(PSI_RETURN_FAILURE);
  }

  return(fptr);
}

