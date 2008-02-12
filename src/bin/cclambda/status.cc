/*! \file status.cc
    \ingroup CCLAMBDA
    \brief Enter brief description of file here 
*/
#include <stdio.h>

namespace psi { namespace cclambda {

void status(char *s, FILE *out)
{
  fprintf(out, "     %-15s...complete\n", s);
  fflush(out);
}

}} // namespace psi::cclambda
