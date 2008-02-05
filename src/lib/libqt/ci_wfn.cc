/*!
  \file ci_wfn.cc
  \brief Check if wavefunction is CI-type
  \ingroup QT
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <psifiles.h>

extern "C" {
	
/*!
** ci_wfn(): Examine the wavefunction type and return 1 if a CI/MCSCF-type,
** otherwise 0
**
** \param *wfn = wavefunction string
**
** Returns: 1 if a CI/MCSCF-type wavefunction, otherwise 0
**
** \ingroup QT
*/
int ci_wfn(char *wfn)
{

  if (strcmp(wfn, "CI")==0     || strcmp(wfn, "DETCAS")==0 || 
      strcmp(wfn, "CASSCF")==0 || strcmp(wfn, "RASSCF")==0 ||
      strcmp(wfn, "DETCI")==0 || strcmp(wfn, "MCSCF")==0 ||
      strcmp(wfn, "OOCCD")==0 || strcmp(wfn,"ZAPTN")==0) 
  {
    return(1);
  }
  else  {
    return 0;
  }
}

} /* extern "C" */
