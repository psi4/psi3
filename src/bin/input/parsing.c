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
   char *guess;

     errcod = ip_string("LABEL", &label,0);
     if(errcod != IPE_OK)
       punt("Where is the label?");

     /*Default = false*/
     /*Flag if the user wants to see that his/her basis set is normalized*/
     shownorm = 0;
     errcod = ip_boolean("SHOWNORM",&shownorm,0);
	 
     /*Default = false*/
     puream = 0;
     errcod = ip_boolean("PUREAM",&puream,0);
	 
     /*------------------------------------------
       Parse Some boolean information from input
      ------------------------------------------*/

     /* If not reading geometry from chkpt file, we will
	need some options specified in input.dat */
     if (chkpt_geom == 0) {
	 /*Default = 1*/
	 print_lvl = 1;
	 errcod = ip_data("PRINT","%d",&print_lvl,0);

	 /*Get the subgroup label */
	 subgroup = NULL;
	 errcod = ip_string("SUBGROUP", &subgroup,0);
	 if (subgroup != NULL && strlen(subgroup) != 0)
	     if (strcmp(subgroup,"C1") && strcmp(subgroup,"C2") && strcmp(subgroup,"CS") && strcmp(subgroup,"CI") &&
		 strcmp(subgroup,"C2V") && strcmp(subgroup,"C2H") && strcmp(subgroup,"D2"))
		 subgroup = NULL;
	 
	 /*Get the unique axis*/
	 unique_axis = NULL;
	 errcod = ip_string("UNIQUE_AXIS", &unique_axis, 0);
	 if (unique_axis != NULL && strlen(unique_axis) != 0)
	     if (strcmp(unique_axis,"X") && strcmp(unique_axis,"Y") && strcmp(unique_axis,"Z"))
		 unique_axis = NULL;

	 /*No default for these two unless running a findif procedure*/
	 if (ip_exist("ZMAT",0) == 1)
	     cartOn = 0;
	 else if (ip_exist("GEOMETRY",0) == 1)
	     cartOn = 1;
	 else
	     punt("Both ZMAT and GEOMETRY are missing!");

	 /*Default = BOHR*/
	 units = strdup("BOHR");
	 errcod = ip_string("UNITS",&units,0);
	 if (!strcmp(units,"BOHR") || !strcmp(units,"AU"))
	     conv_factor = 1.0;
	 else if (!strcmp(units,"ANGSTROMS") || !strcmp(units,"ANGSTROM"))
	     conv_factor = 1.0 / _bohr2angstroms;
	 else
	     punt("Unrecognized UNITS");
	 
	 /*No reorientation?*/
	 no_reorient = 0;
	 errcod = ip_boolean("NO_REORIENT",&no_reorient,0);
     }

     return;
}

     
void parsing_cmdline(int argc, char *argv[])
{
   int i;

   chkpt_mos = 0;
   chkpt_geom = 0;
   overwrite_output = 1;
   
   for (i=1; i<argc; i++) {
       
       /*--- read MOs and project onto new basis? ---*/
       if (strcmp(argv[i], "--getmos") == 0) {
	   chkpt_mos = 1;
       }

       /*--- read geometry from file30? (in findif calculations) ---*/
       if (strcmp(argv[i], "--getgeom") == 0) {
	   chkpt_geom = 1;
	   print_lvl = 0;
	   cartOn = 1;
	   overwrite_output = 0;
       }

   }

   return;
}
