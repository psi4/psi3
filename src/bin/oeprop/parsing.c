#define EXTERN
#include "includes.h"
#include "globals.h"
#include "prototypes.h"



void parsing()
{
  int i,errcod;
  double xmin,xmax,ymin,ymax,zmin,zmax;

        /* Set some defaults for certain wavefunctions */

  errcod = ip_string("WFN", &wfn, 0);

  if (errcod == IPE_OK) {
          /* wfn = one of CC types */
    /* this seems to not work, at least for NO's ---CDS 10/30/02 */
    /*
    if ( (strcmp(wfn, "CCSD")==0) ) {
      read_opdm = 1;
      opdm_file = PSIF_AO_OPDM;
      corr = 0;
      opdm_basis = (char *) malloc(3*sizeof(char));
      strcpy(opdm_basis,"AO");
      opdm_format = (char *) malloc(7*sizeof(char));
      strcpy(opdm_format,"TRIANG");
    }
    */
          /* wfn = ooccd, detci, detcas */
    /*
    if (strcmp(wfn, "CI")==0 || strcmp(wfn, "OOCCD")==0 ||
	strcmp(wfn, "DETCI")==0 || strcmp(wfn, "DETCAS")==0) {
      read_opdm = 1;
      opdm_file = PSIF_AO_OPDM;
      corr = 0;
      opdm_basis = (char *) malloc(3*sizeof(char));
      strcpy(opdm_basis,"SO");
      opdm_format = (char *) malloc(7*sizeof(char));
      strcpy(opdm_format,"TRIANG");
    }
    */
    if (!strcmp(wfn, "CI") || !strcmp(wfn, "DETCI") ||
        !strcmp(wfn, "CCSD") || !strcmp(wfn, "DETCAS") ||
	!strcmp(wfn, "MP2"))  {
      read_opdm = 1;
      opdm_file = PSIF_MO_OPDM;
      corr = 0;
      opdm_basis = (char *) malloc(3*sizeof(char));
      strcpy(opdm_basis,"MO");
      opdm_format = (char *) malloc(7*sizeof(char));
      strcpy(opdm_format,"SQUARE");
    }

  }

  
	/* Parsing section */

  errcod = ip_boolean("READ_OPDM",&read_opdm,0);
  if (read_opdm) {
    errcod = ip_data("OPDM_FILE","%d",&opdm_file,0);
    if ((opdm_file >= MAX_UNIT) || (opdm_file <= 0))
      punt("OPDM_FILE out of range");
    if ((opdm_file != 40) && (opdm_file != 79) && 
        (opdm_file != PSIF_MO_OPDM)) {
      errcod = ip_string("OPDM_BASIS",&opdm_basis,0);
      if ((errcod == IPE_KEY_NOT_FOUND) || (errcod != IPE_OK)) {
        opdm_basis = (char *) malloc(3*sizeof(char));
        strcpy(opdm_basis,"AO");
      }
      errcod = ip_string("OPDM_FORMAT",&opdm_format,0);
      if ((errcod == IPE_KEY_NOT_FOUND) || (errcod != IPE_OK)) { 
        opdm_format = (char *) malloc(7*sizeof(char));
        strcpy(opdm_format,"TRIANG");
      }
    }

    errcod = ip_boolean("WRTNOS",&wrtnos,0);
    errcod = ip_boolean("ASYMM_OPDM",&asymm_opdm,0);
  }

  errcod = ip_boolean("SPIN_PROP",&spin_prop,0);
  if (spin_prop && read_opdm)
    spin_prop = 0;
  if (iopen == 0)
    spin_prop = 0;

  errcod = ip_data("PRINT","%d",&print_lvl,0);
  if (print_lvl < 0)
    print_lvl = 1;
  errcod = ip_boolean("PRINT_NOS",&print_nos,0);

  errcod = ip_boolean("CORREL_CORR",&corr,0);
  /*--- corr should be zero since we are not using Psi2 any longer ---*/
  corr = 0;
  if (corr) {
    errcod = ip_data("ZVEC_FILE","%d",&zvec_file,0);
    if ((zvec_file >= MAX_UNIT) || (zvec_file <= 0))
      punt("ZVEC_FILE out of range");
    errcod = ip_boolean("DELETE_ZVEC",&delete_zvec,0);
  }

  errcod = ip_data("MPMAX","%d",&mpmax,0);
  if (mpmax < 1)
    mpmax = 1;
  else if (mpmax > 3)
         mpmax = 3;


  if (ip_exist("MP_REF_XYZ",0)) {
    ip_count("MP_REF_XYZ",&i,0);
    if (i != 3)
      punt("MP_REF_XYZ must have 3 components");
    for (i=0;i<3;i++) {
      errcod = ip_data("MP_REF_XYZ","%lf",&mp_ref_xyz[i],1,i);
      if (errcod != IPE_OK)
        punt("Error in the definition of MP_REF_XYZ");
    }
    mp_ref = -1;         /* mp_ref = -1 means that mp_ref_xyz specified by user */
  }
  else {
    errcod = ip_data("MP_REF","%d",&mp_ref,0);
    if (mp_ref <= 0)             /* Default is COM */
      mp_ref = 1;
  }

  if (ip_exist("LM_REF_XYZ",0)) {
    ip_count("LM_REF_XYZ",&i,0);
    if (i != 3)
      punt("LM_REF_XYZ must have 3 components");
    for (i=0;i<3;i++) {
      errcod = ip_data("LM_REF_XYZ","%lf",&Lm_ref_xyz[i],1,i);
      if (errcod != IPE_OK)
        punt("Error in the definition of LM_REF_XYZ");
    }
  }
  else {
    Lm_ref_xyz[0] = Lm_ref_xyz[1] = Lm_ref_xyz[2] = 0.0;
  }

  errcod = ip_boolean("NUC_ESP",&nuc_esp,0);
  if (spin_prop)
    nuc_esp = 1;

  errcod = ip_data("GRID","%d",&grid,0);
  if (grid < 0)
    punt("GRID type must be positive");
  if (grid == 5 || grid == 6)
    grid3d = 1;
  if (grid == 5) {
    mo_to_plot = 1;
    read_opdm = 0;
    errcod = ip_data("MO_TO_PLOT","%d",&mo_to_plot,0);
    if (mo_to_plot <= 0 || mo_to_plot > nmo)
      punt("MO_TO_PLOT out of range");
    mo_to_plot--;
  }
  
  if (grid !=0) {
    if (ip_exist("GRID_ORIGIN",0)) {
      ip_count("GRID_ORIGIN",&i,0);
      if (i != 3)
        punt("GRID_ORIGIN must have 3 components");
      for (i=0;i<3;i++) {
        errcod = ip_data("GRID_ORIGIN","%lf",&grid_origin[i],1,i);
        if (errcod != IPE_OK)
          punt("Error in the definition of GRID_ORIGIN");
      }
    }
    else
      punt("GRID_ORIGIN is not defined");
    if (ip_exist("GRID_UNIT_X",0)) {
      ip_count("GRID_UNIT_X",&i,0);
      if (i != 3)
        punt("GRID_UNIT_X must have 3 components");
      for (i=0;i<3;i++) {
        errcod = ip_data("GRID_UNIT_X","%lf",&grid_unit_x[i],1,i);
        if (errcod != IPE_OK)
          punt("Error in the definition of GRID_UNIT_X");
      }
    }
    else
      punt("GRID_UNIT_X is not defined");
    if (ip_exist("GRID_UNIT_Y",0)) {
      ip_count("GRID_UNIT_Y",&i,0);
      if (i != 3)
        punt("GRID_UNIT_Y must have 3 components");
      for (i=0;i<3;i++) {
        errcod = ip_data("GRID_UNIT_Y","%lf",&grid_unit_y[i],1,i);
        if (errcod != IPE_OK)
          punt("Error in the definition of GRID_UNIT_Y");
      }
    }
    else
      punt("GRID_UNIT_Y is not defined");
    
    if (grid3d == 0) {
      if (ip_exist("GRID_XY0",0)) {
	ip_count("GRID_XY0",&i,0);
	if (i != 2)
	  punt("GRID_XY0 must have 2 components");
	for (i=0;i<2;i++) {
	  errcod = ip_data("GRID_XY0","%lf",&grid_xy0[i],1,i);
	  if (errcod != IPE_OK)
	    punt("Error in the definition of GRID_XY0");
	}
      }
      else
	punt("GRID_XY0 is not defined");
    
      if (ip_exist("GRID_XY1",0)) {
	ip_count("GRID_XY1",&i,0);
	if (i != 2)
	  punt("GRID_XY1 must have 2 components");
	for (i=0;i<2;i++) {
	  errcod = ip_data("GRID_XY1","%lf",&grid_xy1[i],1,i);
	  if (errcod != IPE_OK)
	    punt("Error in the definition of GRID_XY1");
	  else if (grid_xy1[i] <= grid_xy0[i])
	    punt("GRID_XY1 must point to the upper right corner of the grid");
	}
      }
      else
	punt("GRID_XY1 is not defined");
    }
    else {
      if (ip_exist("GRID_XYZ0",0)) {
	ip_count("GRID_XYZ0",&i,0);
	if (i != 3)
	  punt("GRID_XYZ0 must have 3 components");
	for (i=0;i<3;i++) {
	  errcod = ip_data("GRID_XYZ0","%lf",&grid_xyz0[i],1,i);
	  if (errcod != IPE_OK)
	    punt("Error in the definition of GRID_XYZ0");
	}
      }
      else
	punt("GRID_XYZ0 is not defined");
    
      if (ip_exist("GRID_XYZ1",0)) {
	ip_count("GRID_XYZ1",&i,0);
	if (i != 3)
	  punt("GRID_XYZ1 must have 3 components");
	for (i=0;i<3;i++) {
	  errcod = ip_data("GRID_XYZ1","%lf",&grid_xyz1[i],1,i);
	  if (errcod != IPE_OK)
	    punt("Error in the definition of GRID_XYZ1");
	  else if (grid_xyz1[i] <= grid_xyz0[i])
	    punt("GRID_XYZ1 must point to the upper right corner of the grid parallelepiped");
	}
      }
      else
	punt("GRID_XYZ1 is not defined");
    }

    errcod = ip_data("NIX","%d",&nix,0);
    if (nix <= 0)
      punt("NIX must be positive integer");
    errcod = ip_data("NIY","%d",&niy,0);
    if (niy <= 0)
      punt("NIY must be positive integer");
    if (grid3d) {
      errcod = ip_data("NIZ","%d",&niz,0);
      if (niz <= 0)
	punt("NIZ must be positive integer");
    }

    if (grid3d == 0) {
	errcod = ip_data("GRID_ZMIN","%lf",&grid_zmin,0);
	errcod = ip_data("GRID_ZMAX","%lf",&grid_zmax,0);
	if (grid_zmin >= grid_zmax)
	  punt("GRID_ZMIN must be less than GRID_ZMAX");
	errcod = ip_data("EDGRAD_LOGSCALE","%d",&edgrad_logscale,0);
    }
  }

}
