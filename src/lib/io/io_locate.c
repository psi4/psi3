
#include <stdio.h>
#include "param.h"

int
io_locate(input,loc_token)
FILE *input;
char loc_token[MAX_STRING];
{
  char line[MAX_STRING];

  fseek(input,0L,0);

  for (;;) {
    if (io_getline(input,line) != 0) return(-1);
    if (!strncmp(loc_token,line,10)) {
      return(0);
      }
    }
  }

