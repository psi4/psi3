#include <ccfiles.h>
#include <libdpd/dpd.h>
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
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN int L_irr;
EXTERN int R_irr;
