/*##################################################################
#
#  small_functions.cc
#
#  various small functions for extrema
#                                                  J.Kenny 7-22-00
###################################################################*/						  

#include <stdio.h>
#include <stdlib.h>

extern "C" {
  #include <libciomr.h>
}

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

void global_allocate() {
 
  coord_vec = init_array(num_coords);
  grad_vec = init_array(num_coords);
  H = init_matrix(num_coords,num_coords);  
  coord_old = init_array(num_coords);
  grad_old = init_array(num_coords);
  H_old = init_matrix(num_coords,num_coords);
 
  return;
}
 
void global_free() {
 
  free(coord_vec);
  free(grad_vec);
  free_matrix(H,num_coords);
  free(coord_old);
  free(grad_old);
  free_matrix(H_old,num_coords);
 
  return;
}           
                    

double dot_pdt(double *vec1, double *vec2, int num) {

  int i;
  double result=0;

  for(i=0;i<num;++i) {
      result += vec1[i] * vec2[i];
    }

  return result;
}




