/*! \file 
    \ingroup (CCLAMBDA)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#define EXTERN
#include "globals.h"

void update(void)
{
  fprintf(outfile,"\t%4d      %20.15f    %4.3e\n",moinfo.iter,moinfo.lcc,
	  moinfo.conv);
  fflush(outfile);
}
