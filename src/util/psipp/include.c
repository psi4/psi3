
#include <stdio.h>
#include <stdlib.h>
#include "global.h"

  /* The maximum include depth. */
#define MAX_INCLUDE 16

static int n_include=0;
static char *names[MAX_INCLUDE];
static int lines[MAX_INCLUDE];
static FILE *filep[MAX_INCLUDE];

#define MAX_INC_DIR 16
static int n_inc_dir=1;
static char* inc_dir[MAX_INC_DIR];

int
include_source_path(char * p)
{
  if (inc_dir[0]) free(inc_dir[0]);
  if (p) {
      inc_dir[0] = p;
    }
  else {
      inc_dir[0] = (char*)malloc(2);
      strcpy(inc_dir[0],".");
    }
  return 0;
}

int
include_directory(dir)
    char* dir;
{
  if (n_inc_dir >= MAX_INC_DIR - 1) {
      fprintf(stderr,"psipp: increase MAX_INC_DIR\n");
      exit(1);
    }
  inc_dir[n_inc_dir] = dir;
  n_inc_dir++;
  return 0;
}

include_traceback()
{
  int i;

  /* The lines array does not contain an entry for the top level file. */
  lines[n_include-1] = lineno;
  fprintf(stderr,"Traceback of included files:\n");
  fprintf(stderr,"<bottom of stack>\n");
  for (i=0; i<n_include; i++) {
    fprintf(stderr,"In file %s at line number %d (roughly).\n",names[i],lines[i]);
    }
  }

/* Sets up a source file.  The fname "stdin" is recognized and is used *
 * to make stdin as the input file.                                    */
initial_file(fname)
char *fname;
{
  fprintf(output,"c_*BEGIN_FILE %s\n",fname);
  if (n_include!=0) syntax_error("initial_file called with n_include nonzero");
  names[n_include] = fname;
  if (input!=stdin) fclose(input);
  if (strcmp(fname,"stdin")) {
    input = fopen(fname,"r");
    filep[n_include] = input;
    if (!filep[n_include]) syntax_error("problem opening intial file");
    }
  else {
    filep[n_include] = stdin;
    }
  lines[n_include] = 0;
  n_include++;

 /* if the file 'fname' is not in the current dir, we'll need to add the
  * directory which contains fname to the include path
  */
  if (strchr(fname,'/')) {
    int i;
    char *path = (char *) malloc(strlen(fname)+1);
    strcpy(path,fname);
    for (i=strlen(path); i>-1; i--) if(path[i]== '/' ) {
      path[i]='\0';
      break;
      }
    include_directory(path);
    }
  }

/* Resets n_include to the appropiate number. */
end_initial_file()
{
  n_include--;
  }

begin_include()
{
  int i;
  int len;
  int quoted=0, psi_include=0;
  char *name = args;
  char *tmp;
  FILE *fp;

  if (n_include<=0) syntax_error("weird... begin_include has n_include <= 0");

  if (n_include>=MAX_INCLUDE) {
    syntax_error("include stack overflow... recursive include?");
    }

  /* The name of the include file is in arg. */
  if (name[0]=='<') {
    psi_include = 1;
    name++;
    }
  if (name[0]=='"') {
    quoted = 1;
    name++;
    }

  len = strcspn(name," \t\n");

  if (quoted&&name[len-1]=='"') name[len-1]='\0';
  else if (quoted) syntax_error("closing \" in include name not found");
  else if (psi_include&&name[len-1]=='>') name[len-1]='\0';
  else if (psi_include) syntax_error("closing > in include name not found");
  else name[len]='\0';

  /* If psi_include is false then look in "." first, so start the loop at 0. */
  for (i=(psi_include?1:0); i<n_inc_dir; i++) {
      tmp = (char *) malloc(strlen(inc_dir[i])+1+strlen(name)+1);
      strcpy(tmp,inc_dir[i]);
      strcat(tmp,"/");
      strcat(tmp,name);

      fp = fopen(tmp,"r");

      if (!fp)
        free(tmp);
      else
        break;
    }

  name = tmp;
  names[n_include] = name;
  lines[n_include] = lineno;

  fprintf(output,"c_*BEGIN_INCLUDE %s\n",name);

  if (!fp) {
      fprintf(stderr,"looking for \"%s\"\n",name);
      syntax_error("problem opening an include file");
    }
  filep[n_include] = fp;

  input = fp;

  n_include++;

  }

end_include()
{
  /* This will tell the scanner when the include stack is empty. */
  if (n_include == 1) return(0);
  n_include--;
  fprintf(output,"c_*END_INCLUDE %s\n",names[n_include]);
  fclose(filep[n_include]);
  input = filep[n_include-1];
  lineno = lines[n_include-1];
  return(1);
  }
