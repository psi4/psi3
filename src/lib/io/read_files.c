
#include <stdio.h>
#include "param.h"

static char *format1 = ":%s:FILES:FILE%d:%s";
static char *format2 = ":%s:FILES:DEFAULT:%s";
static char *format3 = ":DEFAULT:FILES:FILE%d:%s";
static char *format4 = ":DEFAULT:FILES:DEFAULT:%s";

int
oldstyleinput()
{
  int input,ierr;

  /* Use fortran routines to see if 'FILES' is in unit 5. */
  input=5;
  C_LOCATE (&input,"# FILES ##",&ierr);
  if (ierr == 0) return(1);
  return 0;
  }

int
read_files_string(unit,shortkey,value)
int unit;
char *shortkey;
char *value;
{
  char keyword[MAX_STRING];
  char progid[80];
  char *tmp;

  bzero(keyword,MAX_STRING);
  bzero(progid,80);

  /* Get the program id (from a F77 routine). */
  C_CPRGID(progid);
  f77unpad(progid);

  /* Go thru the search path trying to read it in. */
  sprintf(keyword,format1,progid,unit,shortkey);
  if (read_string(keyword,value)) return 1;

  sprintf(keyword,format2,progid,shortkey);
  if (read_string(keyword,value)) return 1;

  sprintf(keyword,format3,unit,shortkey);
  if (read_string(keyword,value)) return 1;

  sprintf(keyword,format4,shortkey);
  if (read_string(keyword,value)) return 1;

  return 0;
  }

int
read_files_boolean(unit,shortkey,value)
int unit;
char *shortkey;
int *value;
{
  char keyword[MAX_STRING];
  char progid[MAX_STRING];
  char *tmp;

  /* Get the program id (from a F77 routine). */
  C_CPRGID(progid);
  f77unpad(progid);

  /* Go thru the search path trying to read it in. */
  sprintf(keyword,format1,progid,unit,shortkey);
  if (read_boolean(keyword,value)) return 1;

  sprintf(keyword,format2,progid,shortkey);
  if (read_boolean(keyword,value)) return 1;

  sprintf(keyword,format3,unit,shortkey);
  if (read_boolean(keyword,value)) return 1;

  sprintf(keyword,format4,shortkey);
  if (read_boolean(keyword,value)) return 1;

  return 0;
  }

int
read_files_integer(unit,shortkey,value)
int unit;
char *shortkey;
int *value;
{
  char keyword[MAX_STRING];
  char progid[MAX_STRING];
  char *tmp;

  /* Get the program id (from a F77 routine). */
  C_CPRGID(progid);
  f77unpad(progid);

  /* Go thru the search path trying to read it in. */
  sprintf(keyword,format1,progid,unit,shortkey);
  if (read_integer(keyword,value)) return 1;

  sprintf(keyword,format2,progid,shortkey);
  if (read_integer(keyword,value)) return 1;

  sprintf(keyword,format3,unit,shortkey);
  if (read_integer(keyword,value)) return 1;

  sprintf(keyword,format4,shortkey);
  if (read_integer(keyword,value)) return 1;

  return 0;
  }

int
read_string(keyword,value)
char *keyword;
char *value;
{
  int errcod;

  /* Read the string using a F77 parsing routine. */
  f77pad(keyword,80);
  errcod = C_FFRDC(keyword,value);
  if (errcod != 0) return 0;
  f77unpad(value);
  return 1;
  }

int
read_boolean(keyword,value)
char *keyword;
int *value;
{
  int errcod;

  /* Read the string using a F77 parsing routine. */
  f77pad(keyword,80);
  errcod = C_FFRDBO(keyword,value);
  if (errcod != 0) return 0;
  return 1;
  }

int
read_integer(keyword,value)
char *keyword;
int *value;
{
  int errcod;

  /* Read the string using a F77 parsing routine. */
  f77pad(keyword,80);
  errcod = C_FFRDI(keyword,value);
  if (errcod != 0) return 0;
  return 1;
  }

int
f77pad(string,length)
char *string;
int length;
{
  int i;
  for (i=strlen(string); i<length; i++) string[i] = ' ';
  return 0;
  }

int
f77unpad(string)
char *string;
{
  char *tmp;

  for (tmp=string; *tmp!=' ' && *tmp!='\0'; tmp++);
  *tmp = '\0';
  return 0;
  }
