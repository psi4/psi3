#define EXTERN
#include "globals.h"

/* DEANTI(): Convert the two-particle density from Dirac to Mulliken
** ordering.  The original, Fock-adjusted density (see the comments in
** fock.c) corresponds to a two-electron energy (or energy derivative)
** expression of the form:
**
** E(TWO) = 1/4 sum_pqrs Gpqrs <pq||rs>
**
** However, the derivative two-electron integrals are produced in
** Mulliken-ordered, symmetric form rather than Dirac-ordered
** antisymmetric form.  This code alters the two-particle density
** matrix ordering for the energy expression of the form:
**
** E(TWO) = 1/2 sum_pqrs Gpqrs <pq|rs>
**
** The final conversion to Mulliken ordering is taken care of in
** dump.c
**
** The second equation above may be derived via
**
** E(TWO) = 1/4 sum_pqrs Gpqrs (<pq|rs> - <pq|sr>)
**        = 1/4 sum_pqrs Gpqrs <pq|rs> - 1/4 sum_pqrs Gpqrs <pq|sr>
**        = 1/4 sum_pqrs Gpqrs <pq|rs> - 1/4 sum_pqrs Gpqsr <pq|rs>
**        = 1/4 sum_pqrs (Gpqrs - Gpqsr) <pq|rs>
**        = 1/2 sum_pqrs Gpqrs <pq|rs>
**
** Equations for the individual components are given explicitly in
** comments below.
** */

void deanti_ROHF(void);
void deanti_UHF(void);

void deanti(void)
{
  if(params.ref == 0 || params.ref == 1) deanti_ROHF();
  else if(params.ref == 2) deanti_UHF();
}
