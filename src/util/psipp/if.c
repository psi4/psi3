
/* This handles conditionally executed code.  Routine names imply use of
 * a stack; however, it is not necessary to really use a stack anywhere.
 */

#include <stdio.h>
#include "global.h"

static int active_level = 0; /* the last level with active code */
static int if_level = 0;     /* the number of nexted ifs */

static int if_satisfied;

/*************************************************************************/

/* This is called for if, ifdef, and ifndef. */
if_action(test)
int test;
{
  if (is_previous_level_active() && test) satisfy_if();
  }

/* This is called for elif and else. */
else_action(test)
int test;
{

  /* Perhaps we are in the middle of an active if. */
  if ( is_level_active() ) {
    active_level--;
    if (active_level<0)
      syntax_error("#endif out of place (parser should have caught this)");
    if (!if_satisfied) syntax_error("else_action: weird error 1");
    }

  if (   is_previous_level_active()
      && if_not_satisfied()
      && test ) {
    satisfy_if();
    }

  }


/*************************************************************************/

push_if()
{
# ifdef DEBUG
  fprintf(stderr,"Push if\n");
# endif
  if_satisfied = 0;
  if_level++;
  }

pop_if()
{
# ifdef DEBUG
  fprintf(stderr,"Pop if\n");
# endif
  if (active_level==if_level) {
    active_level--;
    if_satisfied = 1;
    }
  else if (active_level==if_level-1) {
    if_satisfied = 1;
    }
  if_level--;
  }

int
if_not_satisfied()
{
  return (!if_satisfied);
  }

int
is_level_active()
{
# ifdef DEBUG
  fprintf(stderr,"is_level_active: active_level=%d, if_level=%d\n",
          active_level,if_level);
# endif
  return (active_level==if_level);
  }

int
is_previous_level_active()
{
# ifdef DEBUG
  fprintf(stderr,"is_previous_level_active: active_level=%d, if_level=%d\n",
          active_level,if_level);
# endif
  return ((active_level==if_level)||(active_level==(if_level-1)));
  }

satisfy_if()
{
# ifdef DEBUG
  fprintf(stderr,"satisfy_if: entered\n");
# endif
  active_level = if_level;
  if_satisfied=1;
  }

