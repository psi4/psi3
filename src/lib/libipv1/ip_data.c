/* $Log$
 * Revision 1.2  2002/01/04 18:50:26  evaleev
 * Added ip_double_array to read in an array of floating point numbers.
 *
/* Revision 1.1.1.1  2000/02/04 22:53:26  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.11  1995/11/27 15:32:53  sherrill
/* Make error handling of ip_int_array completely standardized.
/*
 * Revision 1.10  1995/11/09  14:31:08  sherrill
 * Patch up ip_int_array().
 *
 * Revision 1.9  1995/11/09  14:27:24  sherrill
 * Added routine ip_int_array().
 *
 * Revision 1.8  1995/01/16  23:03:49  cljanss
 * Minor changes so the SGI compiler won't complain.
 *
 * Revision 1.7  1994/08/09  22:33:57  crawdad
 * Added check for AIX before including globals.
 *
 * Revision 1.6  1994/08/04  17:33:24  dcrawfrd
 * Added vars for $(LEX) = lex or flex and $(LEXLIB) = -ll or -fl.  These are
 * necessary for portability.  This will have to be thoroughly tested for IBM'
 * AIX.
 *
 * Revision 1.5  1994/06/02  02:22:24  seidl
 * using new tmpl now...change .global to .gbl and .local to .lcl
 *
 * Revision 1.1.1.1  1994/05/02  17:05:52  cljanss
 * The May 1, 1994 version of psi as on the CCQC machines.
 *
 * Revision 1.4  1991/09/18  20:49:48  seidl
 * changes made for DEC
 *
 * Revision 1.3  1991/07/30  03:28:45  seidl
 * add rcs log and id
 * */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <tmpl.h>

#ifdef DEC
#ifdef VOID
#undef VOID
#define VOID char
#endif
#endif

#include "ip_types.h"
#include "ip_global.h"
#include "ip_error.h"


static char rcs_id[] = "$Id$";

/* The way the AIX xlc compiler handles vararg decs isn't compatible with
 * the way the tmpl file is formed, so ip_data.global cannot be used.
 * Everything global should be defined before it is used here. */
/* NOTE: ip_data.global will work fine (and be correct) for other source
 * files.  It's just this one that is a problem because xlc thinks my
 * declarations below are different than the ones in the global file
 * (and they are-sort of). */
#ifndef AIX 
#include "ip_data.gbl"
#endif

#include "ip_data.lcl"

#include "ip_cwk.gbl"

GLOBAL_VA_FUNCTION int
ip_count(keyword,count,n)
char *keyword;
int *count;
int n;
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_count_v(keyword,count,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return IPE_MALLOC;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_count_v(keyword,count,n,v);
    free(v);
    return r;
    }
  }

GLOBAL_FUNCTION int
ip_count_v(keyword,count,n,v)
char *keyword;
int *count;
int n;
int *v;
{
  ip_value_t *val;
  int errcod;

  if (errcod = ip_value_v(keyword,&val,n,v)) return errcod;

  if (val->type != IP_ARRAY) return IPE_NOT_AN_ARRAY;

  *count = val->v.array->n;
  return IPE_OK;
  }

GLOBAL_VA_FUNCTION int
ip_boolean(keyword,boolean,n)
char *keyword;
int *boolean;
int n;
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_boolean_v(keyword,boolean,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return IPE_MALLOC;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_boolean_v(keyword,boolean,n,v);
    free(v);
    return r;
    }
  }

GLOBAL_FUNCTION int
ip_boolean_v(keyword,boolean,n,v)
char *keyword;
int *boolean;
int n;
int *v;
{
  ip_value_t *val;
  int errcod;
  char copy[10],*s;

  if (errcod = ip_value_v(keyword,&val,n,v)) return errcod;

  if (val->type != IP_SCALAR) return IPE_NOT_A_SCALAR;

  strncpy(copy,val->v.scalar,10);
  copy[9] = '\0';

  /* Convert the string to uppercase. */
  for (s=copy; *s!='\0'; s++) {
    if (*s>='a' && *s <='z') *s = *s + 'A' - 'a';
    }
  
  if (!strcmp(copy,"YES")) *boolean = 1;
  else if (!strcmp(copy,"NO")) *boolean = 0;
  else if (!strcmp(copy,"1")) *boolean = 1;
  else if (!strcmp(copy,"0")) *boolean = 0;
  else if (!strcmp(copy,"TRUE")) *boolean = 1;
  else if (!strcmp(copy,"FALSE")) *boolean = 0;
  else return IPE_TYPE;

  return IPE_OK;
  }

/* n should always be zero in this version of libip. */
GLOBAL_VA_FUNCTION int
ip_exist(keyword,n)
char *keyword;
int n;
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_exist_v(keyword,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) {
      ip_warn("ip_exist: problem mallocing %d integers",n);
      return 0;
      }
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_exist_v(keyword,n,v);
    free(v);
    return r;
    }
  }

/* n should always be zero in this version of libip. */
GLOBAL_FUNCTION int
ip_exist_v(keyword,n,v)
char *keyword;
int n;
int *v;
{
  if (ip_cwk_descend_tree(keyword)) return 1;

  return 0;
  }

GLOBAL_VA_FUNCTION int
ip_data(keyword,conv,value,n)
char *keyword;
char *conv;
VOID *value;
int n;
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_data_v(keyword,conv,value,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return IPE_MALLOC;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_data_v(keyword,conv,value,n,v);
    free(v);
    return r;
    }
  }

GLOBAL_FUNCTION int
ip_data_v(keyword,conv,value,n,v)
char *keyword;
char *conv;
VOID *value;
int n;
int *v;
{
  ip_value_t *val;
  int errcod;

  if (errcod = ip_value_v(keyword,&val,n,v)) return errcod;

  if (val->type != IP_SCALAR) return IPE_NOT_A_SCALAR;

  if (sscanf(val->v.scalar,conv,value) != 1) return IPE_TYPE;

  return IPE_OK;
  }

GLOBAL_VA_FUNCTION int
ip_string(keyword,value,n)
char *keyword;
char **value;
int n;
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_string_v(keyword,value,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return IPE_MALLOC;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_string_v(keyword,value,n,v);
    free(v);
    return r;
    }
  }

GLOBAL_FUNCTION int
ip_string_v(keyword,value,n,v)
char *keyword;
char **value;
int n;
int *v;
{
  ip_value_t *val;
  int errcod;

  if (errcod = ip_value_v(keyword,&val,n,v)) return errcod;

  if (val->type != IP_SCALAR) return IPE_NOT_A_SCALAR;

  *value = (char *) malloc(sizeof(char)*(strlen(val->v.scalar)+1));
  if (! *value) return IPE_MALLOC;
  strcpy(*value,val->v.scalar);
  return IPE_OK;
  }

GLOBAL_VA_FUNCTION int
ip_value(keyword,value,n)
char *keyword;
ip_value_t **value;
int n;
{
  va_list args;
  int i;
  int *v;
  int r;

  if (n==0) {
    return ip_value_v(keyword,value,n,NULL);
    }
  else {
    v = (int *) malloc(sizeof(int)*n);
    if (!v) return IPE_MALLOC;
    va_start(args, n);
    for (i=0; i<n; i++) {
      v[i] = va_arg(args,int);
      }
    va_end(args);
    r = ip_value_v(keyword,value,n,v);
    free(v);
    return r;
    }
  }

GLOBAL_FUNCTION int
ip_value_v(keyword,value,n,v)
char *keyword;
ip_value_t **value;
int n;
int *v;
{
  int i;
  ip_value_t *val;

  /* Use the cwk list to obtain the value associated with the keyword. */
  val = ip_key_value(keyword);
  if (!val) return IPE_KEY_NOT_FOUND;

  /* Descend thru val to find the subarray that were are interested in. */
  for (i=0; i<n; i++) {
    if (val->type != IP_ARRAY) return IPE_NOT_AN_ARRAY;
    if (v[i] < 0) return IPE_OUT_OF_BOUNDS;
    if (v[i] >= val->v.array->n) return IPE_OUT_OF_BOUNDS;
    val = val->v.array->values[v[i]];
    }

  *value = val;
  return IPE_OK;
  }


/*
** ip_int_array()
**
** Function reads in an integer array using the PSI input parser.
** It checks for errors at all stages and makes sure that the array
** has the proper length.
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** August 1995
**
** Parameters:
**    keyword = string containing the keyword for the input parser
**    arr     = array to hold results
**    len     = length of array
**
** Returns: IP Error code
**
** Note: keyword should ordinarily be an uppercase string.
*/

GLOBAL_FUNCTION int
ip_int_array(keyword, arr, len)
char *keyword;
int *arr;
int len;
{
  int i, errcod, cnt;

  errcod = ip_count(keyword,&cnt,0);
  if (errcod != IPE_OK) return(errcod);
  if (cnt != len) {
    fprintf(ip_out," (ip_int_array): Trouble parsing %s array.\n", keyword);
    fprintf(ip_out," Length is %d should be %d\n",cnt,len);
    return(IPE_OUT_OF_BOUNDS);
    }
  for (i=0; i<len; i++) {
    errcod = ip_data(keyword,"%d",arr+i,1,i);
    if (errcod != IPE_OK) {
      fprintf(ip_out," (ip_int_array): Trouble parsing %s array element %d\n",
        keyword,i+1);
      ip_warn(ip_error_message(errcod));
      return(errcod);
      }
   }

  return(IPE_OK);
}


/*
** ip_double_array()
**
** Function reads in an array of doubles using the PSI input parser.
** It checks for errors at all stages and makes sure that the array
** has the proper length.
**
** Based on ip_int_array by C. David Sherrill
**
** Parameters:
**    keyword = string containing the keyword for the input parser
**    arr     = array to hold results
**    len     = length of array
**
** Returns: IP Error code
*/

GLOBAL_FUNCTION int ip_double_array(keyword, arr, len)
char *keyword;
double *arr;
int len;
{
  int i, errcod, cnt;

  errcod = ip_count(keyword,&cnt,0);
  if (errcod != IPE_OK) return(errcod);
  if (cnt != len) {
    fprintf(ip_out," (ip_array): Trouble parsing %s array.\n", keyword);
    fprintf(ip_out," Length is %d should be %d\n",cnt,len);
    return(IPE_OUT_OF_BOUNDS);
    }
  for (i=0; i<len; i++) {
    errcod = ip_data(keyword,"%lf",arr+i,1,i);
    if (errcod != IPE_OK) {
      fprintf(ip_out," (ip_array): Trouble parsing %s array element %d\n",
        keyword,i+1);
      ip_warn(ip_error_message(errcod));
      return(errcod);
      }
   }

  return(IPE_OK);
}
