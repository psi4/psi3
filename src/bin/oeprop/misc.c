/*! \file misc.c
    \ingroup (OEPROP)
    \brief Enter brief description of file here 
*/
#define EXTERN
#include "includes.h"
#include "prototypes.h"
#include "globals.h"

void start_io(int argc, char *argv[])
{
  int errcod;

  errcod = psi_start(argc-1,argv+1,0);
  if (errcod != PSI_RETURN_SUCCESS)
    abort();
  ip_cwk_add(":OEPROP");
  tstart(outfile);
  psio_init();
  chkpt_init(PSIO_OPEN_OLD);

  return;
}

void stop_io()
{
  tstop(outfile);
  psio_done();
  psi_stop();

  return;
}


void punt(char *errmsg)
{
  fprintf(stderr, "Error: %s\n", errmsg);
  stop_io();
  abort();
}
