/*!
  \file init_in_out
  \ingroup (QT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include "qt.h"

extern FILE *infile, *outfile;

/*!
** init_in_out()
**
** This function initializes the input and output files
**
** David Sherrill and Micah Abrams, December 2002
**
** Arguments: 
**  \param nstrings   = number of strings passed
**  \param fnames     = file name strings
**
** Returns: none
** \ingroup (QT)
*/
void init_in_out(int nstrings, char *fnames[])
{
  char *ifname, *ofname;

  if (nstrings < 0 || nstrings > 2) {
    fprintf(stderr, "Error: improper number (%d) of filename arguments\n",
      nstrings);
    exit(PSI_RETURN_FAILURE);
  }
  if (nstrings == 1) {
    fprintf(stderr, "Usage: (module) [options] input output\n");
    exit(PSI_RETURN_FAILURE);
  }
  if (nstrings == 2) {
    ffile(&infile,fnames[0],2);
    ffile(&outfile,fnames[1],1);
  }
  else {
    ifname = getenv("PSI_INPUT");
    if (ifname == NULL)
      ffile(&infile,"input.dat",2);
    else
      ffile(&infile,ifname,2);
    ofname = getenv("PSI_OUTPUT");
    if (ofname == NULL)
      ffile(&outfile,"output.dat",1);
    else
      ffile(&outfile,ofname,1);
  }
}

