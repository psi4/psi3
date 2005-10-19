#include <ccfiles.h>
#include "moinfo.h"
#include "params.h"

#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

#define MAXIOFF 32641
#define INDEX(i,j) ((i>j) ? (ioff[(i)]+(j)) : (ioff[(j)]+(i)))

EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
EXTERN struct moinfo mo;
EXTERN struct params params;
EXTERN int* ioff;
