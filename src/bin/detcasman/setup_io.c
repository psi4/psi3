/*
** SETUP_IO.C
**
** Contains some setup/shutdown I/O stuff
**
** C. David Sherrill
** University of California, Berkeley
** April 1998
*/

#include <stdio.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "globals.h"

/*
** init_io(): Function opens input and output files
*/
void init_io(int argc, char *argv[])
{
  int i;
  int parsed=1;
  
  for (i=0; i<ncasiter && !converged; i++) {
    ci_conv = calc_ci_conv(scale_conv);
    if (ci_conv > 1.0E-7) {
      sprintf(detci_string, "detci --quiet -c %12.9lf\n", ci_conv);
      parsed+=2;
    }
    else {
      sprintf(detci_string, "detci --quiet\n");
      parsed++;
    }
  }
  
  init_in_out(argc-parsed, argv+parsed);

  ip_set_uppercase(1);
  ip_initialize(infile, outfile);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":DETCASMAN");
  ip_cwk_add(":DETCAS"); 
  tstart(outfile);
}



/*
** close_io(): Function closes down I/O and exits
*/
void close_io(void)
{
   fclose(infile);
   tstop(outfile);
   fclose(outfile);
}



/*
** check() acts as an abort function if the condition 'a' is not true;
**   it shuts down i/o and returns an error
*/
void check(int a, char *errmsg)
{
  if (!a) {
    fprintf(outfile, "%s\n", errmsg);
    close_io();
    exit(1);
  }
}

/*
** The gprgid() function required by the PSI I/O libraries
*/
char *gprgid()
{
   char *prgid = "DETCASMAN";

   return(prgid);
}


