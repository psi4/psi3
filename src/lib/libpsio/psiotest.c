#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include "psio.h"

/* Function prototypes */
void init_ip(void);
void exit_ip(void);

FILE *infile, *outfile;

int main(int argc, char *argv[])
{
  double junk1, junk4;
  
  init_ip();
  psio_init();

  junk1 = 10*acos(-1);

  psio_open(30,PSIO_OPEN_NEW);
  fprintf(outfile, "ENUC = %20.10f\n", junk1);
  psio_write_entry(30,"Nuclear Repulsion Energy",(char *) &junk1, 
                   sizeof(double));
  psio_read_entry(30, "Nuclear Repulsion Energy", (char *) &junk4,
                  sizeof(double));
  fprintf(outfile, "ENUC = %20.10f\n", junk4);
  psio_close(30,1);
  
  psio_done();
  exit_ip();
  exit(0);
}

void init_ip(void)
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
}

void exit_ip(void)
{
  ip_done();
}

char *gprgid()
{
   char *prgid = "IOTEST";

   return(prgid);
}

