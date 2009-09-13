/*! \file
    \ingroup OPTKING
    \brief DISP_USER only performs input-specified displacements
*/

#define EXTERN
#include "globals.h"
#undef EXTERN
#include "cartesians.h"
#include "simples.h"
#include "salc.h"
#include "opt.h"

#include <libipv1/ip_lib.h>

namespace psi { namespace optking {

void disp_user(const cartesians &carts, simples_class & simples, const salc_set &all_salcs) {
  int i,j,a,b,success;
  int  num_disps = 0, disp_length = 0, restart_geom_file, line_length_count;
  double *geom, *djunk, *dq, **displacements, disp = 0;
  char *disp_label, *ch_temp;

  disp_label = new char[MAX_LINELENGTH];
  djunk = new double[3*carts.get_natom()];

  dq = init_array(all_salcs.get_num());

  ip_count("DISPLACEMENTS",&num_disps,0);

  if (num_disps < 1) {
    punt("No DISPLACEMENTS vector found in input.");
  }

  displacements = init_matrix(num_disps,all_salcs.get_num());
  for (i=0;i<num_disps;++i) {
    disp_length = 0;
    ip_count("DISPLACEMENTS",&disp_length,1,i);
    if (div_int(disp_length,2) && disp_length != 0) {
      for (j=0;j<disp_length;j+=2) {
        ip_data("DISPLACEMENTS","%d",&a,2,i,j);
        ip_data("DISPLACEMENTS","%lf",&disp,2,i,j+1);
        if ((a>0) && (a<=all_salcs.get_num()))
          displacements[i][a-1] = disp; 
        else {
          punt("internal to be displaced does not exist");
        }
      }
    }
    else {
      punt("displacement vector has wrong number of elements");
      exit(1);
    }
  }

  fprintf(outfile,"Displacement Matrix\n");
  print_mat5(displacements, num_disps, all_salcs.get_num(), outfile);

  fprintf(outfile,"\nDisplaced geometries in a.u.\n");
  for (i=0;i<num_disps;++i)  {
    sprintf(disp_label, "Disp:");
    ch_temp = disp_label + 5 ;
    line_length_count = 5;
    for (j=0; j<all_salcs.get_num(); ++j) {
      if (fabs(displacements[i][j]) > 1.0E-8) {
        if (line_length_count < (MAX_LINELENGTH - 10)) {
          sprintf(ch_temp, " (%d %5.3lf)", j+1, displacements[i][j]);
          ch_temp += 10;
          line_length_count += 10;
        }
      }
    }
    restart_geom_file = 0;
    if (i == 0) restart_geom_file = 1;
    success = new_geom(carts,simples,all_salcs,displacements[i],PRINT_TO_GEOM,
        restart_geom_file,disp_label,i, (i==(num_disps-1)?1:0), djunk);
    if (!success) {
      fprintf(outfile,"Unable to generate displaced geometry.\n");
      fclose(outfile);
      exit(PSI_RETURN_FAILURE);
    }
  }
  free_matrix(displacements);
  free(disp_label);
  free(djunk);
}

}} /* namespace psi::optking */

