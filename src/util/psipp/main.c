
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define ALLOCATE_GLOBALS
#include "global.h"

main(argc,argv)
int argc;
char **argv;
{
  int i;

  progname = argv[0];
  lineno = 0;
  input = stdin;
  output = stdout;

  /* Command line processing. */
  for (i=1; i<argc; i++) {
    if (!strncmp(argv[i],"-D",2)) command_line_define(argv[i]);
    else if (!strncmp(argv[i],"-U",2)) command_line_undef(argv[i]);
    else if (!strncmp(argv[i],"-I",2)) include_directory(&(argv[i][2]));
    else if (!strcmp(argv[i],"-o")) {
      FILE *fp;
      i++; if (i>=argc) syntax_error("-o requires an argument");
      fp = fopen(argv[i],"w");
      if (!fp) syntax_error("trouble opening output file");
      if (output!=stdout) fclose(output);
      output = fp;
      }
    else if (!strncmp(argv[i],"-",1)) syntax_error("unknown option");
    else {
      /* Anything else is taken to be a filename. */
      /* If the file has a path, then make the path an include dir. */
      char *last_slash = strrchr(argv[i],'/');
      if (last_slash) {
          char *path = (char*)malloc(last_slash - argv[i] + 2);
          strncpy(path,argv[i],last_slash - argv[i] + 2);
          path[last_slash - argv[i]] = '\0';
          include_source_path(path);
        }
      else {
          include_source_path(0);
        }
      initial_file(argv[i]);
      yyparse();
      end_initial_file();
      }
    }

  /* If the input is still stdin, then I have yet to call yyparse(). */
  if (input==stdin) {
    initial_file("stdin");
    yyparse();
    end_initial_file();
    }

  return(0);
  }

yywrap(){return(1);}
