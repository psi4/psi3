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
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>


FILE *infile, *outfile;
void exit_bad(void);


int main()
{
  ULI i, nvol;
  int errcod;
  char vpath[MAX_STRING];
  char basename[MAX_STRING];
  char fileslist[MAX_STRING];
  char cmdstring[MAX_STRING];

  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_add(":DEFAULT");


  /* Initialize the I/O system */
  psio_init();

  /* Get the number of volumes */
  nvol = psio_get_numvols_default();

  errcod = psio_get_filename_default(basename);
  if (errcod) {
    fprintf(outfile, "psiclean: Trouble reading default filename\n");
    exit_bad();
  }

  for (i=0; i<nvol; i++) {
      errcod = psio_get_volpath_default(i, vpath);

      if (errcod) {
        fprintf(outfile, "psiclean: Trouble reading volume path %d\n", nvol);
        exit_bad();
      }

      sprintf(fileslist,"%s%s.*",vpath,basename);
      sprintf(cmdstring,"echo Removing files %s%s",vpath,basename);
      system(cmdstring);
      sprintf(cmdstring,"ls -l %s",fileslist);
      system(cmdstring);
      sprintf(cmdstring,"/bin/rm %s",fileslist);
      system(cmdstring);
  }

  /* we're done, clean up */
  psio_done();
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

  
