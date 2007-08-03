/*! \file 
    \ingroup (CCLAMBDA)
    \brief Enter brief description of file here 
*/
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

EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN struct L_Params *pL_params;
EXTERN struct Local local;
void check_sum(char *lbl, int L_irr);
