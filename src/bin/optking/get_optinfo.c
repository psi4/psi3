/* GET_OPTINFO   reads optimization parameters from input.dat */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libpsio/psio.h>
#include <psifiles.h>

#define EXTERN
#define C_CODE
#define C_EXTERN
#include "opt.h"
#undef C_EXTERN
#undef C_CODE
#undef EXTERN

double power(double x, int y);

void get_optinfo() {
  int a, i, cnt, cnt2, natom, nallatom;

  optinfo.iteration = 0;
  optinfo.micro_iteration = 0;
  open_PSIF();
  if (psio_tocscan(PSIF_OPTKING, "Iteration") != NULL)
    psio_read_entry(PSIF_OPTKING, "Iteration",
        (char *) &(optinfo.iteration),sizeof(int));
  if (psio_tocscan(PSIF_OPTKING, "Micro_iteration") != NULL)
    psio_read_entry(PSIF_OPTKING, "Micro_iteration",
        (char *) &(optinfo.micro_iteration),sizeof(int));
  close_PSIF();

  optinfo.dertype = 0;
  optinfo.numerical_dertype = 0;

  /* print options */
  optinfo.print_simples = 0;
  ip_boolean("PRINT_SIMPLES", &(optinfo.print_simples),0);
  optinfo.print_params = 0;
  ip_boolean("PRINT_PARAMS", &(optinfo.print_params),0);
  optinfo.print_delocalize = 0;
  ip_boolean("PRINT_DELOCALIZE", &(optinfo.print_delocalize),0);
  optinfo.print_symmetry = 0;
  ip_boolean("PRINT_SYMMETRY", &(optinfo.print_symmetry),0);
  optinfo.print_hessian = 0;
  ip_boolean("PRINT_HESSIAN", &(optinfo.print_hessian),0);

  /* optimization parameters */
  optinfo.optimize = 1;
  if (ip_exist("DISPLACEMENTS",0)) {
    optinfo.do_only_displacements = 1;
    optinfo.mode = MODE_DISP_USER;
    optinfo.optimize = 0;
  }

  optinfo.bfgs = 1;
  ip_boolean("BFGS",&(optinfo.bfgs),0);
  optinfo.bfgs_use_last = 6;
  ip_data("BFGS_USE_LAST","%d",&(optinfo.bfgs_use_last),0);
  optinfo.mix_types = 1;
  ip_boolean("MIX_TYPES",&(optinfo.mix_types),0);

  optinfo.redundant = 1;
  optinfo.delocalize = 0;
  ip_boolean("DELOCALIZE", &(optinfo.delocalize),0);
  if (optinfo.delocalize)
    optinfo.redundant = 0;
  if ((optinfo.mode == MODE_DISP_IRREP) || (optinfo.mode == MODE_DISP_NOSYMM) ) 
    {  optinfo.redundant = 0; optinfo.delocalize =1; }

  /* takes values of 1,2,3 for x,y,z for location of first dummy of linear bend*/
  optinfo.dummy_axis_1 = 1;
  ip_data("DUMMY_AXIS_1","%d",&(optinfo.dummy_axis_1),0);
  optinfo.dummy_axis_1 -= 1;
  optinfo.dummy_axis_2 = 2;
  ip_data("DUMMY_AXIS_2","%d",&(optinfo.dummy_axis_2),0);
  optinfo.dummy_axis_2 -= 1;

  optinfo.zmat = 0;
  if (ip_exist("ZMAT",0)) optinfo.zmat = 1;

  optinfo.zmat_simples = 0;
  ip_boolean("ZMAT_SIMPLES",&(optinfo.zmat_simples),0);

  a = 5;
  ip_data("CONV","%d",&a,0);
  optinfo.conv = power(10.0, -1*a);

  a= 5;
  ip_data("EV_TOL","%d",&a,0);
  optinfo.ev_tol = power(10.0, -1*a);

  optinfo.scale_connectivity = 1.3;
  ip_data("SCALE_CONNECTIVITY","%lf",&(optinfo.scale_connectivity),0);

  optinfo.disp_size = 0.0010;
  ip_data("EDISP","%lf",&(optinfo.disp_size),0);

  /* back-transformation parameters */
  optinfo.bt_max_iter = 60;
  ip_data("BT_MAX_ITER","%d",&(optinfo.bt_max_iter),0);
  a = 2; /* some minimal level of this seems necessary for butane */
  ip_data("BT_DQ_CONV","%d",&a,0);
  optinfo.bt_dq_conv  = power(10.0, -1*a);
  a = 7;
  ip_data("BT_DX_CONV","%d",&a,0);
  optinfo.bt_dx_conv  = power(10.0, -1*a);


  /* Obscure limits in intco evaluation */
  optinfo.cos_tors_near_1_tol = 0.999999999;
  ip_data("COS_TORS_NEAR_1_TOL","%lf",&(optinfo.cos_tors_near_1_tol),0);
  optinfo.cos_tors_near_neg1_tol = -0.999999999;
  ip_data("COS_TORS_NEAR_NEG1_TOL","%lf",&(optinfo.cos_tors_near_neg1_tol),0);
  optinfo.sin_phi_denominator_tol = 0.000000001;
  ip_data("SIN_PHI_DENOMINATOR_TOL","%lf",&(optinfo.sin_phi_denominator_tol),0);


  /* Compute dummy atom lookup arrays */

  chkpt_init(PSIO_OPEN_OLD);
  optinfo.natom = natom = chkpt_rd_natom();
  optinfo.nallatom = nallatom = chkpt_rd_nallatom();
  optinfo.atom_dummy = chkpt_rd_atom_dummy();
  chkpt_close();

  optinfo.to_dummy = init_int_array(natom);
  optinfo.to_nodummy = init_int_array(nallatom);

  cnt=0; cnt2=0;
  for (i=0; i<nallatom; ++i) {
    if (optinfo.atom_dummy[i]) continue;
    else optinfo.to_dummy[cnt++] = i;
    
    if (optinfo.atom_dummy[i]) continue;
    else optinfo.to_nodummy[i] = cnt2++;
  }

  if (optinfo.print_params) {
    for (i=0;i<nallatom;++i)
      fprintf(outfile,"atom_dummy[%d]: %d\n",i,optinfo.atom_dummy[i]);
    for (i=0;i<nallatom;++i)
      fprintf(outfile,"to_nodummy[%d]: %d\n",i,optinfo.to_nodummy[i]);
    for (i=0;i<natom;++i)
      fprintf(outfile,"to_dummy[%d]: %d\n",i,optinfo.to_dummy[i]);
    fflush(outfile);
  }

  /*
  fprintf(outfile,"\nIteration: %d     ",optinfo.iteration+1);
  if (optinfo.numerical_dertype > 0)
    fprintf(outfile,"Micro_iteration: %d",optinfo.micro_iteration+1);
  fprintf(outfile,"\n");
  */
  if (optinfo.print_params) {
    fprintf(outfile,"\n+++ Optinfo Parameters +++\n");
    fprintf(outfile,"print_params:  %d\n",optinfo.print_params);
    fprintf(outfile,"print_simples: %d\n",optinfo.print_simples);
    fprintf(outfile,"print_delocalize %d\n",optinfo.print_delocalize);
    fprintf(outfile,"print_symmetry %d\n",optinfo.print_symmetry);
    fprintf(outfile,"optimize:      %d\n",optinfo.optimize);
    fprintf(outfile,"zmat:          %d\n",optinfo.zmat);
    fprintf(outfile,"dummy_axis_1:    %d\n",optinfo.dummy_axis_1);
    fprintf(outfile,"dummy_axis_2:    %d\n",optinfo.dummy_axis_2);
    fprintf(outfile,"points:        %d\n",optinfo.points);
    fprintf(outfile,"zmat_simples:  %d\n",optinfo.zmat_simples);
    fprintf(outfile,"redundant:     %d\n",optinfo.redundant);
    fprintf(outfile,"bfgs:          %d\n",optinfo.bfgs);
    fprintf(outfile,"bfgs_use_last: %d\n",optinfo.bfgs_use_last);
    fprintf(outfile,"mix_types:     %d\n",optinfo.mix_types);
    fprintf(outfile,"delocalize:    %d\n",optinfo.delocalize);
    fprintf(outfile,"conv:          %.1e\n",optinfo.conv);
    fprintf(outfile,"dertype:       %d\n",optinfo.dertype);
    fprintf(outfile,"numerical dertype: %d\n",optinfo.numerical_dertype);
    fprintf(outfile,"iteration:       %d\n",optinfo.iteration);
    fprintf(outfile,"micro_iteration: %d\n",optinfo.micro_iteration);
    fprintf(outfile,"scale_connectivity: %.3lf\n",optinfo.scale_connectivity);
    fprintf(outfile,"disp_size: %.4lf\n",optinfo.disp_size);
    fprintf(outfile,"bt_max_iter:   %d\n",optinfo.bt_max_iter);
    fprintf(outfile,"bt_max_dq_conv:    %.1e\n",optinfo.bt_dq_conv);
    fprintf(outfile,"bt_max_dx_conv:    %.1e\n",optinfo.bt_dx_conv);
    fprintf(outfile,"cos_tors_near_1_tol:     %10.6e\n",
        optinfo.cos_tors_near_1_tol);
    fprintf(outfile,"cos_tors_near_neg1_tol: %10.6e\n",
        optinfo.cos_tors_near_neg1_tol);
    fprintf(outfile,"sin_phi_denominator_tol: %10.6e\n",
        optinfo.sin_phi_denominator_tol);
    fflush(outfile);
  }
}

/* POWER raises number to a power */
double power(double x, int y) {
  double tval = 1.0;
  int invert = 0;

  if (y < 0) {
    invert = 1;
    y = abs(y);
  }
  for ( ; y>0; --y) 
    tval *= x;
  if (invert) tval = 1.0/tval;
  return tval;
}

