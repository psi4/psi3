#define EXTERN
#include <stdio.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <stdlib.h>
#include <psifiles.h>
#include "input.h"
#include <physconst.h>
#include "global.h"
#include "defines.h"

void start_io(int argc, char *argv[])
{
  int i;
  int parsed = 1; /* account for program name */
  char *ifname, *ofname;
  
  read_chkpt = 0;
  chkpt_mos = 0;
  chkpt_geom = 0; 
  dont_project_mos = 0; 
  geomdat_geom = 0;
  save_oldcalc = 0;
  overwrite_output = 1;
  no_comshift = 0;
  no_reorient = 0;
			  
  for (i=1; i<argc; i++) {
    /*--- read MOs from checkpoint file and project onto new basis ---*/
    if (strcmp(argv[i], "--chkptmos") == 0) {
      read_chkpt = 1;
      chkpt_mos = 1;
      parsed++;
    }
    /*--- read MOs from checkpoint file and save to a separate file ---*/
    else if (strcmp(argv[i], "--savemos") == 0) {
      read_chkpt = 1;
      save_oldcalc = 1;
      parsed++;
    }
    /*--- don't project MOs but simply keep them ---*/
    else if (strcmp(argv[i], "--noproject") == 0) {
      dont_project_mos = 1;
      parsed++;
    }
    /*--- read geometry from checkpoint file (in findif calculations) ---*/
    else if (strcmp(argv[i], "--chkptgeom") == 0) {
      read_chkpt = 0;
      chkpt_geom = 1;
      /* preserve the information about the original reference frame
       * so that properties can be rotated back later by other codes */
      keep_ref_frame = 1;     
      print_lvl = 0;
      cartOn = 1;
      overwrite_output = 0;
      parsed++;
    }
    /*--- read geometry from geom.dat file (in findif calculations) ---*/
    else if (strcmp(argv[i], "--geomdat") == 0) {
      geomdat_geom = 1;
      geomdat_entry = atoi(argv[i+1]);  i++;
      keep_ref_frame = 1;
      cartOn = 1;
      overwrite_output = 0;
      parsed++;
    }
    else if (strcmp(argv[i], "--nocomshift") == 0) {
      no_comshift = 1;
      parsed++;
    }
    else if (strcmp(argv[i], "--noreorient") == 0) {
      no_reorient = 1;
      parsed++;
    }
  }
  
  /* initialize input and output files.  We want to pass the list
   * of arguments starting after the last one we have parsed, so
   * we do pointer arithmetic on argv.  Since this module opens
   * the output file using different flags determined by the 
   * overwrite option, the init_in_out function will not be used. */
  
  if (argc-parsed < 0 || argc-parsed > 2 ) {
    fprintf(stderr, "Error: improper number (%d) of filename arguments\n",
            argc-parsed);
    exit(PSI_RETURN_FAILURE);
  }
  if (argc-parsed == 1) {
    fprintf(stderr, "Usage: (module) [options] input output\n");
    exit(PSI_RETURN_FAILURE);
  }
  if (argc-parsed == 2) {
    if (overwrite_output) { 
      ffile(&infile, argv[parsed+0], 2);
      outfile = fopen(argv[parsed+1], "w+"); 
    }  
    else {
      ffile(&infile, argv[parsed+0], 2);
      outfile = fopen(argv[parsed+1], "a+");
    }      
  }
  else {
    ifname = getenv("PSI_INPUT");
    if (ifname == NULL)
      ffile(&infile,"input.dat",2);
    else
      ffile(&infile,ifname,2);
    ofname = getenv("PSI_OUTPUT");
    if (ofname == NULL) {
      if (overwrite_output) 
	outfile = fopen("./output.dat", "w+");
      else
        outfile = fopen("./output.dat", "a+");
    }
    else {
      if (overwrite_output)
	outfile = fopen("./ofname", "w+");
      else
        outfile = fopen("./ofname", "a+");
    }
  }
  
  ip_set_uppercase(1);
  ip_initialize(infile, outfile);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");
  ip_cwk_add(":PSI");
  ip_cwk_add(":INPUT");
  tstart(outfile);
  /*--- Initialize new IO system ---*/
  psio_init();

  return;
}

void stop_io()
{
  tstop(outfile);
  ip_done();
  psio_done();
  fclose(outfile);
  fclose(infile);
}


