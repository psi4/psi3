#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
#else 
# define EXTERN
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern "C" {
#include <libciomr.h>
#include <ip_libv1.h>
#include <file30.h>
#include <physconst.h>
} 

#include "defines.h"

void input(void);
void punt(char*);
void parsing(void);



EXTERN FILE *infile, *outfile;

EXTERN int coord_type;
EXTERN int errcod;
EXTERN int num_entries, num_atoms, num_coords, iteration;
EXTERN double bond_lim, angle_lim;
EXTERN int converged;
EXTERN int grad_max;
EXTERN int print_lvl;
EXTERN int dummy;
EXTERN int bt_loop;
EXTERN char *update;

/*this needs to be in C*/
extern "C" {
    char *gprgid(); 
}



#include "coord_base.h"
#include "carts.h"
#include "simple.h"
#include "internals.h"
#include "zmat.h"




