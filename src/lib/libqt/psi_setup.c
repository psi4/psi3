/*!
  \file psi_setup
  \ingroup (QT)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include "qt.h"

extern FILE *infile, *outfile;

/*!
** psi_setup()
**
** This function initializes the input, output files, file prefix, etc.,
** by checking command line arguments and environmental variables
**
** Arguments: 
**  \param argc       = number of command-line arguments passed
**  \param argv       = command-line arguments
**  \param overwrite  = whether to overwrite output file (1) or append to it (0).
**                      Most PSI modules will want to append.
**
** Returns: one of standard PSI error codes
** \ingroup (QT)
*/

void psi_setup(int args, char *argv[], int overwrite_output)
{
  int i;
                                /* state flags */
  int found_if_np = 0;          /* found input file name without -i */
  int found_of_np = 0;          /* found output file name without -o */
  int found_fp_np = 0;          /* found file prefix name without -p */
  int found_if_p = 0;           /* found input file name with -i */
  int found_of_p = 0;           /* found output file name with -o */
  int found_fp_p = 0;           /* found file prefix name with -p */
  char *ifname, *ofname, *fprefix, *arg;

  /* process command-line arguments in sequence */
  for(i=0; i<argc; i++) {
    arg = argv[i];
    if (!strcmp(arg,"-f") && !found_if_p) {
      ifname = argv[++i];
      found_if_p = 1;
    }
    else if (!strcmp(arg,"-o") && !found_of_p) {
      ofname = argv[++i];
      found_of_p = 1;
    }
    else if (!strcmp(arg,"-p") && !found_fp_p) {
      fprefix = argv[++i];
      found_fp_p = 1;
    }
    else if (arg[0] == '-') {
      fprintf(stderr, "Error: unrecognized command-line argument %s\n", arg);
      return(PSI_RETURN_FAILURE);
    }
    else if (!found_if_np) {
      ifname = arg;
      found_if_np = 1;
    }
    else if (!found_of_np) {
      ofname = arg;
      found_of_np = 1;
    }
    else if (!found_fp_np) {
      fprefix = arg;
      found_fp_np = 1;
    }
    else {
      fprintf(stderr, "Error: too many command-line arguments given\n");
      return(PSI_RETURN_FAILURE);
    }
  }
  

  /* check if some arguments were specified in both prefixed and nonprefixed form */
  if (found_if_p && found_if_np) {
    fprintf(stderr, "Error: input file name specified both with and without -f\n");
    fprintf(stderr, "Usage: (module) [options] -f input -o output [-p prefix]  OR\n");
    fprintf(stderr, "       (module) [options] input output [prefix]\n");
    return(PSI_RETURN_FAILURE);
  }
  if (found_of_p && found_of_np) {
    fprintf(stderr, "Error: output file name specified both with and without -o\n");
    fprintf(stderr, "Usage: (module) [options] -f input -o output [-p prefix]  OR\n");
    fprintf(stderr, "       (module) [options] input output [prefix]\n");
    return(PSI_RETURN_FAILURE);
  }
  if (found_fp_p && found_fp_np) {
    fprintf(stderr, "Error: file prefix specified both with and without -p\n");
    fprintf(stderr, "Usage: (module) [options] -f input -o output -p prefix  OR\n");
    fprintf(stderr, "       (module) [options] input output prefix\n");
    return(PSI_RETURN_FAILURE);
  }

  
  /* if some arguments were not specified on command-line - check the environment */
  if (ifname == NULL)
    ifname = getenv("PSI_INPUT");
  if (ofname == NULL)
    ofname = getenv("PSI_OUTPUT");
  if (fpname == NULL)
    fpname = getenv("PSI_PREFIX");

  /* if some arguments still not defined - assign default values */
  if (ifname == NULL)
    ifname = strdup("input.dat");
  if (ofname == NULL)
    ofname = strdup("output.dat");
  /* prefix is not assigned here yet because someone else should check
     input file for it??? */

  /* copy over file prefix, etc. into their appropriate variables */
  /* NOTHING TO DO YET */

  /* finally, open input and output files */
  infile = fopen(ifname, "r");
  if (overwrite_output)
    outfile = fopen(ofname, "w+");
  else
    outfile = fopen(ofname, "a+");

  return(PSI_RETURN_SUCCESS);
}
