#include <stdio.h>
#include <stdlib.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>

/* Function prototypes */
void init_io(int argc, char *argv[]);
void exit_io(void);
int cc2unit(char *);

FILE *infile, *outfile;
char *psi_file_prefix;

int main(int argc, char *argv[])
{
  int unit=0, i=0;

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

  init_io(argc-1, argv+i);

  psio_open(unit,PSIO_OPEN_OLD);
  psio_tocprint(unit,stdout);
  psio_close(unit,1);

  exit_io();
  exit(0);
}

void init_io(int argc, char *argv[])
{
  extern char *gprgid();
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  psi_start(argc-1, argv+1, 0);
  ip_cwk_add(progid);
  free(progid);

  psio_init();
}

void exit_io(void)
{
  psio_done();
  psi_stop();
}

char *gprgid()
{
   char *prgid = "TOCPRINT";

   return(prgid);
}
