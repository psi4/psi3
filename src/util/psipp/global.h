
/* #define DEBUG */

#include "parse.h"

#ifndef EXTERN
#  ifdef ALLOCATE_GLOBALS
#    define EXTERN
#    define INITIALIZE(x,y) x=y
#  else
#    define EXTERN extern
#    define INITIALIZE(x,y) x
#  endif
#endif

#define BUFSIZE 512

EXTERN char *INITIALIZE(dir,NULL), *INITIALIZE(args,NULL);
EXTERN FILE *INITIALIZE(input,NULL), *INITIALIZE(output,NULL);
EXTERN char *progname;
EXTERN int lineno;
extern YYSTYPE yylval;

char *name_to_val();
char *eval_name_to_str();

