#include <ccfiles.h>
#include <dpd.h>
#include "MOInfo.h"
#include "Params.h"

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

// #define EOM_DEBUG (1)

EXTERN char *progid;
EXTERN FILE *infile, *outfile;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN struct Eom_params eom_params;
EXTERN int ***dpd_dp;


