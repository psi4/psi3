
#include <stdio.h>
#include "global.h"

print_token(token)
int token;
{
  char *name;
  char tmp[2];

  switch (token) {
    case T_INACTIVE_IF:
      name = "inactive_if";
      break;
    case T_FORTRAN:
      name = "fortran";
      break;
    case T_INCLUDE:
      name = "include";
      break;
    case T_IF:
      name = "if";
      break;
    case T_ELIF:
      name = "elif";
      break;
    case T_IFDEF:
      name = "ifdef";
      break;
    case T_IFNDEF:
      name = "ifndef";
      break;
    case T_ELSE:
      name = "else";
      break;
    case T_ENDIF:
      name = "endif";
      break;
    case T_DEFINE:
      name = "define";
      break;
    case T_UNDEF:
      name = "undef";
      break;
    case T_OR:
      name = "||";
      break;
    case T_AND:
      name = "&&";
      break;
    case T_EQ:
      name = "==";
      break;
    case T_NOT_EQ:
      name = "!=";
      break;
    case T_GREA_EQ:
      name = ">=";
      break;
    case T_LESS_EQ:
      name = "<=";
      break;
    case T_NAME:
      name = "name";
      break;
    case T_DEFINED:
      name = "defined";
      break;
    case T_STREQ:
      name = "streq";
      break;
    default:
      if (token >= 256) name = "unknown";
      else {
        name = tmp;
        tmp[0] = token;
        tmp[1] = '\0';
        }
    }

  fprintf(stderr,"token: %d; name: %s\n",token,name);

  }
