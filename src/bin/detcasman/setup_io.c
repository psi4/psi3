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
#include <libciomr.h>
#include "globals.h"

/*
** init_io(): Function opens input and output files
*/
void init_io(void)
{
   ffile(&infile,"input.dat",2) ;
   ffile(&outfile,"output.dat",1);
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


