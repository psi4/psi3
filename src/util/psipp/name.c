
/* This file contains the name table manipulation routines. */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "global.h"

/* This linked list will hold the name table */
struct name_struct {
  char *name;
  char *value;
  struct name_struct *p;
  };
typedef struct name_struct name_t;

static name_t *names = NULL;

/* Evaluates name as an integer (abends if it is not an integer). */
int
eval_name_to_int(name)
char *name;
{
  int i, integer;
  char *value;

  /* If the name consists of digits 0-9 and +-, then convert it to a integer. */
  integer = 1;
  for (i=0; i<strlen(name); i++) {
    if (!(  isspace(name[i])
          ||isdigit(name[i])
          ||(i==0&&(name[i]=='+'||name[i]=='-')))) {
      integer = 0;
      }
    }
  if (integer) return(atoi(name));

  value = eval_name_to_str(name);

  for (i=0; i<strlen(value); i++) {
    if (!(  isspace(value[i])
          ||isdigit(value[i])
          ||(i==0&&(value[i]=='+'||value[i]=='-')))) {
      syntax_error("Ahrg!  Names in expressions must evaluate to integers.");
      }
    }

  return(atoi(value));
  }

/* Evaluates name as a character (abends if it has no value). */
char *
eval_name_to_str(name)
char *name;
{
  name_t *l;

  for (l=names; l!=NULL; l=l->p) {
    if (!strcmp(l->name,name)) {
      if (l->value==NULL) {
        syntax_error("Arhg! Tried to evaluate a name with no value.");
        }
      return(l->value);
      }
    }

  syntax_error("Argh! Tried to evaluate a name which wasn't in the name list.");
  return 0;
  }

/* Returns true if the name is defined. */
int
is_name(name)
char *name;
{
  name_t *l;

  for (l=names; l!=NULL; l=l->p) {
#   ifdef DEBUG
    fprintf(stderr,"l->name=%s, name=%s\n",l->name,name);
#   endif
    if (!strcmp(l->name,name)) return(1);
    }
  return(0);
  }

/* Remove the current name/value pair (in the global variable args) from the
 * names linked list.  Junk following the name is ignored.
 */
undef_name()
{
  char *name,*tmp;

  /* Remove the blanks preceeding name. */
  name = args;
  while((*name==' '||*name=='\t')&&*name!='\0') name++;
  if (*name == '\0') {
    syntax_error("Ouch! Undef directive is missing something.");
    }

  /* Remove any junk after name. */
  tmp=name;
  while (*tmp!=' '&&*tmp!='\t'&&*tmp!='\n'&&*tmp!='\0') tmp++;
  *tmp = '\0';

  undef_name_(name);
  }

undef_name_(name)
char *name;
{
  name_t *l,*back_ptr;

  /* Look thru the names list for name. */
  back_ptr = NULL;
  for (l=names; l!=NULL; l=l->p) {
    /* If name is found delete it. */
    if (!strcmp(l->name,name)) {
      free(l->name);
      free(l->value);
      if (back_ptr) {
        back_ptr->p = l->p;
        }
      else names = l->p;
      free(l);
      }
    else back_ptr = l;
    }

  }

/* Put the current name/value pair (in the global variable args) into the
 * names linked list.
 */
define_name()
{
  char *name;
  char *value;
  char *tmp;

  /* Break args into a name and a value. */
  name = args;
  for (value=name;
       *value!='\n'&&*value!=' '&&*value!='\0'&&*value!='\t';
       value++);

  /* This has the side effect of terminating name. */
  if (*value != '\0') *(value++) = '\0';

  /* Proceed by skipping whitespaces. */
  while((*value==' '||*value=='\t')&&*value!='\0') value++;

  /* If value ends with a newline, get rid of it. */
  if (*value != '\0') {
    if (value[strlen(value)-1]=='\n') value[strlen(value)-1] = '\0';
    }

  /* See if name is legal. */
  if (!(isupper(name[0])||islower(name[0])||(name[0]=='_'))) {
    syntax_error("Fed define_name an illegal name.");
    }

  tmp = (char *) malloc(strlen(name)+1);
  strcpy(tmp,name);
  name = tmp;
  tmp = (char *) malloc(strlen(value)+1);
  strcpy(tmp,value);
  value = tmp;

  define_name_(name,value);
  }

define_name_(name,value)
char *name,*value;
{
  name_t *l;
  name_t *tmp;

  /* See if the name has already been defined. */
  for (l=names; l!=NULL; l=l->p) {
    if (!strcmp(names->name,name)) {
      syntax_error("Hey!  The name in a define is already defined.");
      }
    }

  /* Add the new name to the linked list of names. */
  tmp = (name_t *) malloc(sizeof(name_t));
  tmp->p = names;
  tmp->name = name;
  tmp->value = value;
  names = tmp;

  }

/* defines a name given in the form -Dname=value or -Dname */
command_line_define(com)
char *com;
{
  char *name;
  char *value;
  int index_eq, length;

  if (strncmp(com,"-D",2)) syntax_error("command_line_define given bad string");

  length = strlen(com);
  index_eq = strcspn(com,"=");
  if (index_eq<=2) syntax_error("No name found in command line define");

  name = (char *) malloc(index_eq-1);
  strncpy(name,&com[2],index_eq-2);
  name[index_eq-2] = '\0';

  if (index_eq>=length-1) {
    value = (char *) malloc(1);
    value[0] = '\0';
    }
  else {
    value = (char *) malloc(length-index_eq+1);
    strncpy(value,&com[index_eq],length-index_eq);
    value[length-index_eq] = '\0';
    }

  define_name_(name,value);
  }

/* undefs a name given in the form -Uname */
command_line_undef(com)
char *com;
{
  char *name;

  if (strncmp(com,"-U",2)) syntax_error("command_line_undef given bad string");

  name = &com[2];

  if (strlen(name)==0) {
    syntax_error("No name found in a command line undef");
    }

  undef_name_(name);
  }
