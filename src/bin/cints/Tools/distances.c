/* -------------------------------------------------
   
   distances.c
   
   by Shawn Brown
   
   Calculates the distance matrix for the 
   molecule.  Kid stuff.
   
   Also has a function that will return the distance
   given two geometry structs
   
   ------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include<memory.h>
#include<stdlib.h>
#include<ip_libv1.h>
#include<libciomr.h>
#include<libint.h>

#include"data_structs.h"

double distance_calc(struct coordinates g1, struct coordinates g2){
    
    return sqrt(pow(g1.x-g2.x,2)
		+pow(g1.y-g2.y,2)
		+pow(g1.z-g2.z,2));
}
