
#include <stdio.h>
#include "param.h"

int
io_getline(input,line)
FILE *input;
char line[MAX_STRING];
{
  int i;
  int ch;

  for (i=0; i<MAX_STRING; i++) {
    if ((ch=fgetc(input)) == EOF) return(-1);
    line[i] = ch;
    if (line[i] == '\n') {
      line[i] = '\0';
      return(0);
      }
    }

  fprintf(stdout,"io_getline: buffer size exceeded\n");
  fprintf(stderr,"io_getline: buffer size exceeded\n");
  io_q_abort();
  return -1;
  }
