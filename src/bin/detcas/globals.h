/*! \file 
    \ingroup (DETCAS)
    \brief Enter brief description of file here 
*/
/*
** GLOBALS.H
**
** List of all the global data used by the program
**
** Note that these are given as "extern", so they must be defined
** in the main program!
**
** C. David Sherrill
** University of California, Berkeley
*/

#include "calcinfo.h"
#include "params.h"

extern struct calcinfo CalcInfo;
extern struct params Params;
extern int *ioff;
extern FILE *infile, *outfile;
extern char *psi_file_prefix;
