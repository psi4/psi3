#define EXTERN
#include <stdio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include "input.h"
#include "global.h"
#include "defines.h"


void freeze_core()
{
  int large = 0;
  int atom;

  if (!strcmp(frozen_core,"FALSE") ||
      !strcmp(frozen_core,"NO") ||
      !strcmp(frozen_core,"0")) {
    nfzc = 0;
    return;
  }
  else if (!strcmp(frozen_core,"TRUE") ||
	   !strcmp(frozen_core,"YES") ||
	   !strcmp(frozen_core,"1") ||
	   !strcmp(frozen_core,"SMALL") ||
	   !strcmp(frozen_core,"LARGE")) {

      if (!strcmp(frozen_core,"LARGE"))
	large = 1;

      nfzc = 0;
      for(atom=0; atom<num_atoms; atom++) {
	/* H - Be */
	if (nuclear_charges[atom] < 4.1)
	  continue;
	/* B - Ne */
	else if (nuclear_charges[atom] > 4.9 && nuclear_charges[atom] < 10.1)
	  nfzc++;
	/* Na - Ar */
	else if (nuclear_charges[atom] > 10.9 && nuclear_charges[atom] < 18.1)
	  if (large)
	    nfzc += 5;
	  else
	    nfzc++;
	else
	  punt("Cannot freeze core automatically for fourth and higher row elements yet");
      }
  }
  else
    punt("Invalid value for FROZEN_CORE");

}
