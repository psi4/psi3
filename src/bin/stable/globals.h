#include <ccfiles.h>
#include "MOInfo.h"
#include "Params.h"

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

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))

