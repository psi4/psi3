#include <ccfiles.h>
#include "MOInfo.h"
#include "Params.h"

#define IOFF_MAX 32641

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;

EXTERN int *ioff;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;

#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))