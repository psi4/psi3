/*
** SETUP_IO.C
**
** Contains some setup/shutdown I/O stuff
**
** C. David Sherrill
** University of California, Berkeley
** April 1998
*/

#include <stdlib.h>
#include <stdio.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include "globals.h"

/*
** init_io(): Function opens input and output files
*/
void init_io(void)
{
   ffile(&infile,"input.dat",2) ;
   ffile(&outfile,Params.ofname,1);
   if (Params.print_lvl) tstart(outfile);
   ip_set_uppercase(1);
   ip_initialize(infile, outfile);
   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":DETCAS");
   psio_init();
}



/*
** close_io(): Function closes down I/O and exits
*/
void close_io(void)
{
   psio_done();
   fclose(infile);
   if (Params.print_lvl) tstop(outfile);
   ip_done();
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
   char *prgid = "DETCAS";

   return(prgid);
}


