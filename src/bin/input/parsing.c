#define EXTERN
#include <stdio.h>
#include <ip_libv1.h>
#include <libciomr.h>
#include <string.h>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"


void parsing()
{
   int errcod;

     /*Setting some defaults*/
     puream = 0;
     shownorm = 0;
     print_lvl = 1;
     units = strdup("BOHR");
     subgroup = NULL;
     unique_axis = NULL;
     no_reorient = 0;

     errcod = ip_string("LABEL", &label,0);
     if(errcod != IPE_OK)
       punt("\nERROR: Where is the label?\n");


     /*------------------------------------------
       Parse Some boolean information from input
      ------------------------------------------*/

     /*Default = false*/
     /*Flag if the user wants to see that his/her basis set is normalized*/
     errcod = ip_boolean("SHOWNORM",&shownorm,0);

     /*Get the subgroup label */
     errcod = ip_string("SUBGROUP", &subgroup,0);
     if (subgroup != NULL && strlen(subgroup) != 0)
       if (strcmp(subgroup,"C1") && strcmp(subgroup,"C2") && strcmp(subgroup,"CS") && strcmp(subgroup,"CI") &&
	   strcmp(subgroup,"C2V") && strcmp(subgroup,"C2H") && strcmp(subgroup,"D2"))
	 subgroup = NULL;

     /*Get the unique axis*/
     errcod = ip_string("UNIQUE_AXIS", &unique_axis, 0);
     if (unique_axis != NULL && strlen(unique_axis) != 0)
       if (strcmp(unique_axis,"X") && strcmp(unique_axis,"Y") && strcmp(unique_axis,"Z"))
	 unique_axis = NULL;
     
     /*Default = false*/
     errcod = ip_boolean("PUREAM",&puream,0);

     /*Default = 1*/
     errcod = ip_data("PRINT","%d",&print_lvl,0);

     /*No default :-) */
     if (ip_exist("ZMAT",0) == 1)
       cartOn = 0;
     else if (ip_exist("GEOMETRY",0) == 1)
       cartOn = 1;
     else
       punt("\nERROR: Both ZMAT and GEOMETRY are missing!\n\n");

     /*Default = BOHR*/
     errcod = ip_string("UNITS",&units,0);
     if (!strcmp(units,"BOHR") || !strcmp(units,"AU"))
       conv_factor = 1.0;
     else if (!strcmp(units,"ANGSTROMS") || !strcmp(units,"ANGSTROM"))
       conv_factor = 1.0 / _bohr2angstroms;
     else
       punt("\nERROR: unrecognized UNITS\n\n");

     /*No reorientation?*/
     errcod = ip_boolean("NO_REORIENT",&no_reorient,0);
}

     
