/*##########################################################################*/
/*! \file extrema.h
  \brief Included header files, function declarations, and variables. */
/*##########################################################################*/

#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else 
# define EXTERN
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern "C" {
#include <libciomr.h>
#include <ip_libv1.h>
#include <file30.h>
#include <physconst.h>
} 

#include "defines.h"

void punt(char*);
double **symm_matrix_invert(double**, int, int, int);

EXTERN FILE *infile, *outfile;

EXTERN int coord_type;
EXTERN int errcod;

/*this needs to be in C*/
extern "C" {
    char *gprgid();
}

#include "math_tools.h"
#include "coord_base.h"
#include "simple.h"
#include "internals.h"
#include "zmat.h"


