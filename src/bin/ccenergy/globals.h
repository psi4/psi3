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

#ifdef DMALLOC
#include <dmalloc.h>
#endif

EXTERN FILE *infile, *outfile;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
