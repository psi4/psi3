#include <stdio.h>
#define EXTERN
#include "globals.h"

void update(void)
{
  if(!moinfo.iopen)
      fprintf(outfile,"\t%4d      %20.15f    %4.3e    %7.6f    %7.6f\n",
	      moinfo.iter,moinfo.ecc,moinfo.conv,moinfo.t1diag,moinfo.d1diag);
  else
      fprintf(outfile,"\t%4d      %20.15f    %4.3e    %7.6f\n",
	      moinfo.iter,moinfo.ecc,moinfo.conv,moinfo.t1diag);
  fflush(outfile);
}
