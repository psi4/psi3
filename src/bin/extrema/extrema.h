#ifdef EXTERN
# undef EXTERN
# define EXTERN extern
extern int MAX_LINELENGTH;
#else 
# define EXTERN
int MAX_LINELENGTH = 133;
#endif

#include<new>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern "C" {
#include <libciomr.h>
#include <ip_libv1.h>
#include <physconst.h>
} 

void input(void);
void punt(char*);
void parsing(void);



EXTERN FILE *infile, *outfile;

EXTERN int coord_type;
EXTERN int errcod;
EXTERN int num_entries, num_atoms, num_coords, iteration;

/*this needs to be in C*/
extern "C" {
    char *gprgid(); 
}



#include "coord_base.h"
#include "carts.h"
#include "simple.h"
#include "internals.h"
#include "zmat.h"




