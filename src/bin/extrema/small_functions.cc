/*##################################################################
#
#  small_functions.cc
#
#  various small functions for extrema
#                                                  J.Kenny 7-22-00
###################################################################*/						  

#include <stdio.h>

#define EXTERN
#include "opt.h"


/*this needs to be in C*/
extern "C" {
char *gprgid()
{
   char *prgid = "EXTREMA";
 
   return(prgid);
   }
}                      


void punt(char *mess)
{
  fprintf(outfile, "  error: %s\n", mess);
  fprintf(stderr, "  EXTREMA error: %s\n", mess);
  // stop_io();
  exit(1);
}                    



