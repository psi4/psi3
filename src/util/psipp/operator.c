
/* Here are the special operators that psipp recognizes which are
 * not handled in parse.y.
 */

#include <stdio.h>
#include <stdlib.h>
#include "global.h"

/* This returns true if the value of the two names are the same.  If
 * one (or both) of the names is quoted, then the string with the
 * quotes removed is used in the comparison.
 */
int
streq(name1,name2)
char *name1,*name2;
{
  int result;
  char *val1,*val2;

  val1 = name_to_val(name1);
  val2 = name_to_val(name2);

  result = !strcmp(val1,val2);
  free(val1);
  free(val2);
  return(result);
  }

char *
name_to_val(name)
char *name;
{
  char *result,*tmp;

  if (name[0]=='"') {
    if (name[strlen(name)-1]!='"') {
      syntax_error("Yikes! Couldn't find the closing quote on a string.");
      }
    result = (char *) malloc(strlen(name));
    strcpy(result,&name[1]);
    result[strlen(result)-1] = '\0';
    }
  else {
    tmp = eval_name_to_str(name);
    result = (char *) malloc(sizeof(tmp)+1);
    strcpy(result,tmp);
    }

  return(result);
  }
