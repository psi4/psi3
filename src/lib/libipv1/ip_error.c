/* $Log$
 * Revision 1.1  2000/02/04 22:53:26  evaleev
 * Initial revision
 *
/* Revision 1.5  1994/08/04 17:33:33  dcrawfrd
/* Added vars for $(LEX) = lex or flex and $(LEXLIB) = -ll or -fl.  These are
/* necessary for portability.  This will have to be thoroughly tested for IBM'
/* AIX.
/*
 * Revision 1.4  1994/06/02  02:22:25  seidl
 * using new tmpl now...change .global to .gbl and .local to .lcl
 *
 * Revision 1.1.1.1  1994/05/02  17:05:52  cljanss
 * The May 1, 1994 version of psi as on the CCQC machines.
 *
 * Revision 1.3  1991/07/30  03:28:45  seidl
 * add rcs log and id
 * */

#include <stdio.h>
#include <stdarg.h>
#include <tmpl.h>
#include "ip_types.h"
#include "ip_global.h"

static char rcs_id[] = "$Id$";

/* Cannot include ip_error.global due to xlc's handling of varargs. */
#include "ip_error.gbl"
#include "ip_error.lcl"
#include "ip_error.h"

#include "scan.gbl"

/* Returns some text for an errcod. */
GLOBAL_FUNCTION char *
ip_error_message(errcod)
int errcod;
{
  static char *ipe_ok = "No problem has been detected.";
  static char *ipe_key_not_found = "No match was found for the given keyword.";
  static char *ipe_out_of_bounds = "An array index is out of bounds.";
  static char *ipe_malloc = "Memory allocation failed.";
  static char *ipe_not_an_array = "An index was given for a scalar quantity.";
  static char *ipe_not_a_scalar = "Expected a scalar, but found an array.";
  static char *ipe_type = "The datum is not of the appropiate type.";
  static char *huh = "The nature of the problem is unknown.";

  if (errcod == IPE_OK) return ipe_ok;
  if (errcod == IPE_KEY_NOT_FOUND) return ipe_key_not_found;
  if (errcod == IPE_OUT_OF_BOUNDS) return ipe_out_of_bounds;
  if (errcod == IPE_MALLOC) return ipe_malloc;
  if (errcod == IPE_NOT_AN_ARRAY) return ipe_not_an_array;
  if (errcod == IPE_NOT_A_SCALAR) return ipe_not_a_scalar;
  if (errcod == IPE_TYPE) return ipe_type;
  return huh;
  }

GLOBAL_VA_FUNCTION VOID
ip_error(msg)
char *msg;
{
  va_list args;
  va_start(args,msg);
  fprintf(ip_out,"IP_ERROR: ");
  vfprintf(ip_out,msg,args);
  fprintf(ip_out,"\n");
  va_end(args);
  showpos();
  exit(1);
  }

GLOBAL_VA_FUNCTION VOID
ip_warn(msg)
char *msg;
{
  va_list args;
  va_start(args,msg);
  fprintf(ip_out,"IP_WARN: ");
  vfprintf(ip_out,msg,args);
  fprintf(ip_out,"\n");
  va_end(args);
  }
