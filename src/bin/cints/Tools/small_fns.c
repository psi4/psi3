#include<stdio.h>
#include<stdlib.h>
#include<libipv1/ip_lib.h>
#include<math.h>
#include<libciomr/libciomr.h>
#include<libchkpt/chkpt.h>
#include <psifiles.h>

#include"defines.h"
#define EXTERN
#include"global.h"

void setup()
{
   register int i, j;

   ioff[0] = 0;
   for (i=1; i<IOFFMAX; i++){
      ioff[i] = ioff[i-1] + i;
      }

/* Factorials */
   fac[0] = 1;
   fac[1] = 1;
   for(i=2;i<CINTS_MAX_AM*2;i++)
     fac[i] = fac[i-1]*i;

/* Binomial coefficients */
   for(i=0;i<=CINTS_MAX_AM;i++)
     for(j=0;j<=i;j++)
       bc[i][j] = (int)(fac[i]/(fac[i-j]*fac[j]));
 
/* df[i] gives (i-1)!!, so that (-1)!! is defined... */
/* we shouldn't need both this and lci with the range needed on df[] */
  df[0] = 1.0;
  df[1] = 1.0;
  df[2] = 1.0;
  for(i=3; i<MAXFACT*2; i++){
    df[i] = (i-1)*df[i-2];
    }

  /* (2i-1)!! a useful thing, num_ser in the integral expansion */
   num_ser[0] = 1;
   for (i=1; i<CINTS_MAX_AM+1; i++)
      num_ser[i] = (2*i-1)*num_ser[i-1];

}

void start_io(int argc, char *argv[])
{
  int i, errcod;
  int num_extra_args = 0;
  char **extra_args;
  extra_args = (char **) malloc(argc*sizeof(char *));

  /* Filter out known options */
  for (i=1; i<argc; i++) {
    if ( strcmp(argv[i], "--fock") &&
	 strcmp(argv[i], "--oeints") &&
	 strcmp(argv[i], "--teints") &&
	 strcmp(argv[i], "--deriv1") &&
	 strcmp(argv[i], "--deriv2") &&
	 strcmp(argv[i], "--oeprop") &&
	 strcmp(argv[i], "--mp2") &&
	 strcmp(argv[i], "--r12ints") &&
	 strcmp(argv[i], "--mp2r12") )
      extra_args[num_extra_args++] = argv[i];
  }
  
  errcod = psi_start(num_extra_args, extra_args, 0);
  if (errcod != PSI_RETURN_SUCCESS)
    exit(PSI_RETURN_FAILURE);
  ip_cwk_add(":CINTS");
  psio_init();
  chkpt_init(PSIO_OPEN_OLD);

  free(extra_args);
  return;
}

void stop_io()
{
  chkpt_close();
  if(UserOptions.print_lvl)
    tstop(outfile);
  psio_done();
  psi_stop();
}

void punt(char *mess)
{
  fprintf(outfile, "  error: %s\n", mess);
  fprintf(stderr, "  CINTS error: %s\n", mess);
  stop_io();
  exit(PSI_RETURN_FAILURE);
}

double distance_calc(struct coordinates g1, struct coordinates g2)
{
  return sqrt((g1.x-g2.x)*(g1.x-g2.x) +
	      (g1.y-g2.y)*(g1.y-g2.y) +
	      (g1.z-g2.z)*(g1.z-g2.z));
}
