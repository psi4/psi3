
#include <stdio.h>
#include "global.h"

yyerror(s)
char *s;
{
  syntax_error(s);
  }

syntax_error(s)
char *s;
{
  include_traceback();
  scan_dump();
  fprintf(stderr,"%s: %s\n",progname,s);
  exit(1);
  }
