
/* This inverts the transformation done by psipp.  A source file for a *
 * particular machine will be converted into one or more master source *
 * files.                                                              */

/* If an argument is given it is taken to be the input file name.  If  *
 * any of the master source files came from the standard input when    *
 * the given source file was built, then they will go to the standard  *
 * output.  If any master source files already exist with write access *
 * they will be overwritten.  If a master source file did not change,  *
 * then it will not be touched.                                        */

#include <stdio.h>
#include <string.h>

#define BUFSIZE 512

static char *progname;
static int lineno=0;

char *getline();

main(argc,argv)
int argc;
char **argv;
{
  int included=0;
  FILE *input;
  FILE *output;
  char *line;
  char *filename;

  progname = argv[0];

  if (argc>2) death("too many arguments");
  else if (argc==2) {
    input = fopen(argv[1],"r");
    if (!input) death("could not open the input file");
    }

  while(line=getline(input)) {

    /* Check for psipp directives. */
    if ((!strncmp(line,"c_*",3))||(!strncmp(line,"C_*",3))) {
      if (!strncmp(&line[3],"BEGIN_FILE",10)) {
        filename = strtok(line," \t\n");
        filename = strtok(NULL," \t\n");
        if (!strcmp(filename,"stdin")) output = stdout;
        else {
          output = fopen(filename,"w");
          if (!output) death("could open requested file");
          }
        }
      else if (!strncmp(&line[3],"BEGIN_INCLUDE",13)) {
        included++;
        }
      else if (!strncmp(&line[3],"BEGIN_PSI_INCLUDE",17)) {
        included++;
        }
      else if (!strncmp(&line[3],"END_INCLUDE",11)) {
        included--;
        if (included<0) death("Too many END_INCLUDEs");
        }
      else {
        death("bad psipp directive");
        }
      }
    else if ((!strncmp(line,"c_#",3))||(!strncmp(line,"C_#",3))) {
      if (!included) fprintf(output,"%s",&line[3]);
      }
    else {
      if (!included) fprintf(output,"%s",line);
      }
    }
  }

death(s)
char *s;
{
  fprintf(stderr,"%s (at %d): %s\n",progname,lineno,s);
  exit(1);
  }

/* Get a single line from the input, newline included.     *
 * When this routine is called, the old line is wiped out. */
char *
getline(input)
FILE *input;
{
  int i, ch;
  static char buf[BUFSIZE];

  for (i=0; ; i++) {
    if (i>=BUFSIZE) {
      fprintf(stderr,"%s: input buffer size has been exceeded\n",progname);
      exit(1);
      }
    ch = getc(input);
    if (ch == EOF) {
      if (i==0) return(NULL);
      buf[i] = '\0';
      break;
      }
    else if (ch == '\n') {
      buf[i++] = ch;
      buf[i] = '\0';
      break;
      }
    buf[i] = ch;
    }
  /* Increment the global line number counter. */
  lineno++;
  return(buf);
  }

