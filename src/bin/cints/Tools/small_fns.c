#include<stdio.h>
#include<stdlib.h>
#include<ip_libv1.h>
#include<math.h>
#include<libciomr.h>
#include<malloc.h>
#include<libint.h>

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
   for(i=2;i<MAX_AM*2;i++)
     fac[i] = fac[i-1]*i;

/* Binomial coefficients */
   for(i=0;i<=MAX_AM;i++)
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
   for (i=1; i<MAX_AM+1; i++)
      num_ser[i] = (2*i-1)*num_ser[i-1];

}

void start_io()
{
  ffile(&infile, "input.dat", 2);
  ffile(&outfile, "output.dat", 1);
  ip_set_uppercase(1);
  ip_initialize(infile, outfile);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":CINTS");
  /*--- Initialize new IO system ---*/
  psio_init();
  /*--- Initialize libfile30 library ---*/
  file30_init();

  return;
}

void stop_io()
{
  file30_close();
  if(UserOptions.print_lvl)
    tstop(outfile);
  psio_done();
  fclose(outfile);
  fclose(infile);

  return;
}

void punt(char *mess)
{
  fprintf(outfile, "  error: %s\n", mess);
  fprintf(stderr, "  CINTS error: %s\n", mess);
  stop_io();
  exit(1);
}

