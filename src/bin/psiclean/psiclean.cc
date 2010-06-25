/*! \defgroup PSICLEAN psiclean: Delete scratch files */

/*! 
** \file
** \ingroup PSICLEAN
** \brief Delete scratch files
**
** Utility program to delete scratch files.  Generalization of earlier
** PSI2.0 shell script which was limited to scratch files being put
** in /tmp[0-9]/$user/$name.* .  Here we will search the default path
** instead.
**
** C. David Sherrill
** 
*/ 

#include <cstdio>
#include <cstdlib>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

#define MAX_STRING 300

extern "C" {
  FILE *infile, *outfile;
  char *psi_file_prefix;
}

namespace psi { namespace psiclean {
  void exit_bad(void)
  {
    psio_done();
    exit(1);
  }
}}

int main(int argc, char *argv[])
{
  ULI i, nvol;
  int errcod;
  char *vpath;
  char *basename;
  char fileslist[MAX_STRING];
  char cmdstring[MAX_STRING];

  psi_start(&infile,&outfile,&psi_file_prefix,argc-1,argv+1,0);
  
  /* Initialize the I/O system */
  psio_init(); psio_ipv1_config();

  // clean out all psio files
  psio_purge();

  /* we're done, clean up */
  psio_done();
  psi_stop(infile,outfile,psi_file_prefix);
  exit(0);
}


extern "C" const char *gprgid()
{
   const char *prgid = "PSICLEAN";

   return(prgid);
}

  
