/*!
  \file keyword.c
  \ingrpu (CHKPT)
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "chkpt.h"
#include <libpsio/psio.h>

char *chkpt_build_keyword(char *key)
{
  char *keyword;
  int keylen;

  keylen = strlen(key) + strlen(chkpt_prefix) + 2;
  if(keylen > PSIO_KEYLEN) {
    printf("LIBCHKPT: requested key exceeds allowed LIBPSIO length: :%s:%s\n", 
	   chkpt_prefix, key);
    exit(2);
  }

  keyword = (char *) malloc((keylen+1)*sizeof(char));
  sprintf(keyword, ":%s:%s", chkpt_prefix, key);
  keyword[keylen] = '\0';

  return keyword;
}

