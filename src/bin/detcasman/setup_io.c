/*! \file 
    \ingroup (DETCASMAN)
    \brief Enter brief description of file here 
*/
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
#include <stdlib.h>
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
  int num_extra_args=0;
  char **extra_args;

  extra_args = (char **) malloc(argc*sizeof(char *));
  
  for (i=1; i<argc; i++) {
    extra_args[num_extra_args++] = argv[i];
  }
  
  psi_start(num_extra_args, extra_args, 0);
  ip_cwk_add(":DETCASMAN");
  ip_cwk_add(":DETCAS"); 
  tstart(outfile);
  free(extra_args);
}



/*
** close_io(): Function closes down I/O and exits
*/
void close_io(void)
{
  tstop(outfile);
  psi_stop();
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


