#include <ccfiles.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

/* #define EOM_DEBUG (1) */
#define TIME_CCEOM (1)

#define H_IRR (0)
#define MAX(I,J) ((I>J) ? I : J)
#define MIN(I,J) ((I<J) ? I : J)

EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN struct Eom_params eom_params;
EXTERN struct Local local;
EXTERN int ***dpd_dp;


