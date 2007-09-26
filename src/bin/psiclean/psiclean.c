/*! \file 
    \ingroup (PSICLEAN)
    \brief Enter brief description of file here 
*/

/*! \defgroup PSICLEAN Add a description of the group PSICLEAN */

/*
** PSICLEAN
**
** Utility program to delete scratch files.  Generalization of earlier
** PSI2.0 shell script which was limited to scratch files being put
** in /tmp[0-9]/$user/$name.* .  Here we will search the default path
** instead.
**
** C. David Sherrill
** 
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>


FILE *infile, *outfile;
char *psi_file_prefix;
void exit_bad(void);

int main(int argc, char *argv[])
{
  ULI i, nvol;
  int errcod;
  char *vpath;
  char *basename;
  char fileslist[MAX_STRING];
  char cmdstring[MAX_STRING];

  psi_start(argc-1,argv+1,0);
  
  /* Initialize the I/O system */
  psio_init(); psio_ipv1_config();

  /* Get the number of volumes */
  nvol = psio_get_numvols_default();

  errcod = psio_get_filename_default(&basename);

  for (i=0; i<nvol; i++) {
    errcod = psio_get_volpath_default(i, &vpath);

    /* errcod == 1 now means that the default of /tmp/ is used for volpath 
    **  -TDC, 8/03 */
    /*
    if (errcod) {
      fprintf(outfile, "psiclean: Trouble reading volume path %d\n", nvol);
      exit_bad();
    }
    */

    sprintf(fileslist,"%s%s.*",vpath,basename);
    sprintf(cmdstring,"echo Removing files %s%s",vpath,basename);
    system(cmdstring);
    sprintf(cmdstring,"ls -l %s",fileslist);
    system(cmdstring);
    sprintf(cmdstring,"/bin/rm %s",fileslist);
    system(cmdstring);
    free(vpath);
  }
  free(basename);

  /* we're done, clean up */
  psio_done();
  psi_stop();
  exit(0);
}


void exit_bad(void)
{
  psio_done();
  exit(1);
}


char *gprgid()
{
   char *prgid = "PSICLEAN";

   return(prgid);
}

  
