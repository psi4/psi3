
/* $Log$
 * Revision 1.1  2000/02/04 22:53:26  evaleev
 * Initial revision
 *
/* Revision 1.5  1995/01/16 22:58:56  cljanss
/* Minor changes so the SGI compiler won't complain.
/*
 * Revision 1.4  1994/06/02  02:22:21  seidl
 * using new tmpl now...change .global to .gbl and .local to .lcl
 *
 * Revision 1.1.1.1  1994/05/02  17:05:52  cljanss
 * The May 1, 1994 version of psi as on the CCQC machines.
 *
 * Revision 1.3  1991/07/30  03:28:45  seidl
 * add rcs log and id
 * */

#include <stdio.h>
#include <stdlib.h>
#include <tmpl.h>
#include "ip_types.h"
#include "ip_global.h"

#include "ip_alloc.gbl"
#include "ip_alloc.lcl"

#include "ip_error.gbl"

static char rcs_id[] = "$Id$";

GLOBAL_FUNCTION ip_keyword_tree_t *
ip_alloc_keyword_tree()
{
  ip_keyword_tree_t *result;

  result = (ip_keyword_tree_t *) malloc(sizeof(ip_keyword_tree_t));
  if (!result) {
    perror("ip_alloc_keyword_tree: malloc failed");
    ip_error(NULL);
    }

  result->up = NULL;
  result->down = NULL;
  result->across = NULL;
  result->keyword = NULL;
  result->value = NULL;

  return result;
  }

GLOBAL_FUNCTION VOID
ip_free_keyword_tree(tree)
ip_keyword_tree_t *tree;
{
  ip_keyword_tree_t *I,*nextI;

  if (!tree) return;

  /* Free all the keyword_trees in the across list. */
  I=tree;
  do {
    /* Free the sub trees first. */
    ip_free_keyword_tree(I->down);
    free(I->keyword);
    nextI = I->across;

    /* Zero out I (so if I accidently use it again I'll get SEGV). */
    I->down = NULL;
    I->keyword = NULL;
    I->across = NULL;
    I->up = NULL;

    free(I);
    } while ((I = nextI) != tree);

  }

GLOBAL_FUNCTION ip_value_t *
ip_alloc_value()
{
  ip_value_t *result;

  result = (ip_value_t *) malloc(sizeof(ip_value_t));
  if (!result) {
    perror("ip_alloc_value: malloc failed");
    ip_error(NULL);
    }
  result->type = IP_UNDEFINED;
  result->v.scalar = NULL;
  return result;
  }

GLOBAL_FUNCTION VOID
ip_free_value(value)
ip_value_t *value;
{
  if (!value) return;

  if (value->type == IP_ARRAY) {
    ip_free_array(value->v.array);
    value->v.array = NULL;
    }
  else if (value->type == IP_SCALAR) {
    free(value->v.scalar);
    value->v.scalar = NULL;
    }

  free(value);
  }

GLOBAL_FUNCTION VOID
ip_free_array(array)
ip_array_t *array;
{
  int i;

  for (i=0; i<array->n; i++) {
    ip_free_value(array->values[i]);
    }

  free(array->values);
  free(array);
  }

