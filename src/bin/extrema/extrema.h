/*##########################################################################*/
/*! \file 
    \ingroup (EXTREMA)
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
#include <ctype.h>
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.h>
#include <physconst.h>
} 

#include "defines.h"

void punt(char*);
double **symm_matrix_invert(double**, int, int, int);

EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;

EXTERN int coord_type;
EXTERN int errcod, error;

/*this needs to be in C*/
extern "C" {
    char *gprgid();
}

#include "math_tools.h"
#include "coord_base.h"
#include "simple.h"
#include "internals.h"
#include "zmat.h"
#include "deloc.h"

/*inline function declarations*/
inline double *unit_vec(double* cart_arr, int atom1, int atom2 );
inline double vec_norm(double* cart_arr, int atom1, int atom2 );
inline double dot_pdt( double* vec1, double* vec2 );
inline double *cross_pdt( double* vec1, double* vec2 );
inline double norm(double* cart_arr, int atom1, int atom2 );
inline double compute_bond(double *car, int atm, int bnd);
inline double compute_angle(double *car, int atm, int bnd, int ang);
inline double compute_torsion(double *car, int atm, int bnd, int ang, int tor);


