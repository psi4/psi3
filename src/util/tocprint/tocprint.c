#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/* Function prototypes */
void init_io(void);
void exit_io(void);
int cc2unit(char *);

FILE *infile, *outfile;

int main(int argc, char *argv[])
{
  int unit=0, i=0;
  init_io();

  while(--argc > 0) {
       i++;
       if(!strcmp(argv[i], "-unit") || !strcmp(argv[i], "-u")) {
           unit = atoi(argv[i+1]);
           argc--;
         }
       else if(!strcmp(argv[i], "-cc")) {
           unit = cc2unit(argv[i+1]);
           argc--;
         }
    }

  if(!unit) { printf("Bad unit number.\n"); exit(1); }

  psio_open(unit,PSIO_OPEN_OLD);
  psio_tocprint(unit,stdout);
  psio_close(unit,1);

  exit_io();
  exit(0);
}

void init_io(void)
{
  char *gprgid();
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  ffile(&infile,"input.dat",2);
  ffile(&outfile,"output.dat",1);
  ip_set_uppercase(1);
  ip_initialize(infile,outfile);
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(progid);

  free(progid);

  psio_init();
}

void exit_io(void)
{
  psio_done();
  ip_done();
  fclose(infile);
  fclose(outfile);
}

char *gprgid()
{
   char *prgid = "TOCPRINT";

   return(prgid);
}
