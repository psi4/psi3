#define EXTERN
#include <stdio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"

void enmo_setup()
{
   int i;
   int errcod;
   char *nucleus_type;
   int num_isotopes;

   num_isotopes = 0;
   errcod = ip_count("ENMO_ISOTOPES",&num_isotopes,0);
   if (num_isotopes == 0) {
     fprintf(outfile,"\nENMO_ISOTOPES not specified -- using default of most abundant isotope for each nucleus\n");
     /* isotopes were not given - default is most abundant isotope */
   } 
   else {
     errcod = ip_count("ENMO_ISOTOPES",&num_isotopes,1,0);
     if (num_isotopes != num_atoms || errcod != IPE_NOT_AN_ARRAY) {
       fprintf(outfile,"\nKeyword ENMO_ISOTOPES should be a vector of num_atoms elements\n");
       punt("Keyword ENMO_ISOTOPES should be a vector of num_atoms elements");
     }
     /*----------
      Nuclei sets for each atom is specified, e.g.
      ENM0_ISOTOPES = (h d)
      ---------*/
     for(i=0;i<num_atoms;i++) {
       errcod = ip_string("ENMO_ISOTOPES",&nucleus_type,1,i);
       if (errcod != IPE_OK) 
         punt("There is a problem with the BASIS array!");
         isotope[i] = nucleus_type;
       }
       fprintf(outfile,"\n   ENMO_ISOTOPES = ( ");
       for(i=0; i<num_atoms; i++) 
         fprintf(outfile,"%s ",isotope[i]);
       fprintf(outfile,")\n");
     } 
}
