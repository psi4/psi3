#include <ccfiles.h>

/* Global variables */
#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

#define CC_TEI_HALFT 91

EXTERN FILE *infile, *outfile;
EXTERN char *psi_file_prefix;
EXTERN int *ioff;
EXTERN struct MOInfo moinfo;
EXTERN struct Params params;
EXTERN struct Local local;
