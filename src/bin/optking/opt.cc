/*** OPT.CC Rollin King, begun 1999 ***/

/*
command-line      internal specifier   what it does
--opt_step         MODE_OPT_STEP        take an ordinary geometry optimization step
--disp_nosymm      MODE_DISP_NOSYMM     displaces along all internals, assumes no sym
--disp_irrep       MODE_DISP_IRREP      displaces along all internals labeled as IRREP
--disp_load        MODE_DISP_LOAD       load a previous displacement from PSIF to chkpt
--disp_user        MODE_DISP_USER       perform displacements given in input
--load_ref         MODE_LOAD_REF        load undisplaced reference geometry into chkpt
--freq_energy      MODE_FREQ_ENERGY     not yet supported
--grad_energy      MODE_GRAD_ENERGY     use energies in chkpt to compute a gradient
--freq_grad_nosymm MODE_FREQ_GRAD_NOSYMM use grads in chkpt to compute freqs, assumes no sym
--freq_grad_irrep  MODE_FREQ_GRAD_IRREP use grads in chkpt to compute IRREP freqs
--grad_save        MODE_GRAD_SAVE       save the gradient in chkpt to PSIF list
--energy_save      MODE_ENERGY_SAVE     save the energy in chkpt to PSIF list
--disp_num         optinfo.disp_num     displacement index (by default read from PSIF)
--points           optinfo.points      3 or 5 pt formula (5pt not yet supported)
--irrep            optinfo.irrep       the irrep (1,2,...) being displaced or computed
*/

#include <cmath>
extern "C" {
  #include <stdio.h>
  #include <libchkpt/chkpt.h>
  #include <stdlib.h>
  #include <string.h>
  #include <ctype.h>
  #include <libciomr/libciomr.h>
  #include <libipv1/ip_lib.h>
  #include <physconst.h>
  #include <libpsio/psio.h>
  #include <libqt/qt.h>	  
  #include <psifiles.h>
}

#include "opt.h"
#include "cartesians.h"
#include "internals.h"
#include "salc.h"
#include "bond_lengths.h"

void intro(int argc, char **argv);
void free_info(int nsimples);
int *count_internals(cartesians &cart_geom, int intco_given);
extern "C" void zmat_to_intco(void);
extern "C" void get_optinfo();
extern void get_syminfo(internals &simples);
extern void delocalize(internals &simples, cartesians &carts);
extern void disp_user(cartesians &carts, internals &simples, 
                      salc_set &all_salcs);
extern int make_disp_irrep(cartesians &carts, internals &simples, 
                    salc_set &all_salcs, int points);
extern int make_disp_nosymm(cartesians &carts, internals &simples, 
                    salc_set &all_salcs, int points);
extern void load_ref(cartesians &carts);
extern void freq_grad_irrep(cartesians &carts, internals &simples, 
    salc_set &all_salcs, int points);
extern void freq_grad_nosymm(cartesians &carts, internals &simples, 
    salc_set &all_salcs, int points);
extern void grad_energy(cartesians &carts, internals &simples, 
                        salc_set &all_salcs);
extern void grad_save(cartesians &carts);
extern void energy_save(cartesians &carts);
extern int opt_step(cartesians &carts, internals &simples, salc_set &symm_salcs);

int main(int argc, char **argv) {

    int i,j,a,b,dim,count,dim_carts,user_intcos;
    int parsed=1, num_disps, disp_length;
    double *f, *coord; 
    char *disp_label, *buffer, *err;

    // Set defaults & read command-line arguments
    optinfo.mode = MODE_OPT_STEP;
    optinfo.disp_num = 0;
    optinfo.points = 3;
    for (i=1; i<argc; ++i) {
      if (!strcmp(argv[i],"--disp_nosymm")) {
        optinfo.mode = MODE_DISP_NOSYMM;
        parsed++;
      }
      else if (!strcmp(argv[i],"--disp_irrep")) {
        optinfo.mode = MODE_DISP_IRREP;
        parsed++;
      }
      else if (!strcmp(argv[i],"--disp_load")) {
        optinfo.mode = MODE_DISP_LOAD;
        parsed++;
      }
      else if (!strcmp(argv[i],"--disp_user")) {
        optinfo.mode = MODE_DISP_USER;
        parsed++;
      }
      else if (!strcmp(argv[i],"--load_ref")) {
        optinfo.mode = MODE_LOAD_REF;
        parsed++;
      }
      else if (!strcmp(argv[i],"--opt_step")) {
        optinfo.mode = MODE_OPT_STEP;
        parsed++;
      }
      else if (!strcmp(argv[i],"--freq_energy")) {
        optinfo.mode = MODE_FREQ_ENERGY;
        parsed++;
      }
      else if (!strcmp(argv[i],"--grad_energy")) {
        optinfo.mode = MODE_GRAD_ENERGY;
        parsed++;
      }
      else if (!strcmp(argv[i],"--freq_grad_nosymm")) {
        optinfo.mode = MODE_FREQ_GRAD_NOSYMM;
        parsed++;
      }
      else if (!strcmp(argv[i],"--freq_grad_irrep")) {
        optinfo.mode = MODE_FREQ_GRAD_IRREP;
        parsed++;
      }
      else if (!strcmp(argv[i],"--grad_save")) {
        optinfo.mode = MODE_GRAD_SAVE;
        parsed++;
      }
      else if (!strcmp(argv[i],"--energy_save")) {
        optinfo.mode = MODE_ENERGY_SAVE;
        parsed++;
      }
      else if (!strcmp(argv[i],"--disp_num")) {
        sscanf(argv[++i], "%d", &optinfo.disp_num);
        parsed+=2;
      }
      else if (!strcmp(argv[i],"--points")) {
        sscanf(argv[++i], "%d", &optinfo.points);
        parsed+=2;
      }
      else if (!strcmp(argv[i],"--irrep")) {
        sscanf(argv[++i], "%d", &optinfo.irrep);
        optinfo.irrep -= 1;
        parsed+=2;
      }
      else {
        printf("command line argument not understood.\n");
        exit(1);
      }
    }

    psi_start(argc-parsed,argv+parsed,0);
    /* init_in_out() sets the value of "infile", so we need to save it */
    fp_input = infile;
    
    intro(argc, argv);

    ip_cwk_add(":OPTKING");

    // determine if simples and salcs are present in intco.dat 
    optinfo.simples_present = 0;
    optinfo.salcs_present = 0;
    ffile_noexit(&fp_intco, "intco.dat", 2);
    if (fp_intco != NULL) {
      ip_append(fp_intco, outfile) ;
      if ( ip_exist(":INTCO",0) ) {
        ip_cwk_add(":INTCO");
        optinfo.simples_present = 1;
      }
      if ( ip_exist("SYMM",0) || ip_exist("ASYMM",0) )
        optinfo.salcs_present = 1;
      fclose(fp_intco);
    }
    // printf("simples_present already? %d\n",optinfo.simples_present);
    // printf("salcs_present already? %d\n",optinfo.salcs_present);

    psio_init();
    get_optinfo();

    // generate simples from zmat
    if ( optinfo.zmat_simples && !(optinfo.simples_present)) {
      zmat_to_intco();
      optinfo.simples_present = 1;
      ffile(&fp_intco, "intco.dat", 2);
      ip_append(fp_intco, outfile) ;
      ip_cwk_add(":INTCO");
    }

    disp_label = new char[MAX_LINELENGTH];
    /* loads up geometry, (possibly gradient and energy) from chkpt */
    cartesians carts;
    dim_carts = carts.get_natom()*3;
    if ((optinfo.mode == MODE_OPT_STEP) || (optinfo.mode == MODE_GRAD_SAVE)) {
      fprintf(outfile,
      "\nCartesian geometry and possibly gradient in a.u. with masses\n");
      carts.print(3,outfile,0,disp_label,0);
    }
    else if ((optinfo.mode != MODE_DISP_LOAD)
        && (optinfo.mode != MODE_LOAD_REF)) {
      fprintf(outfile,
      "\nCartesian geometry and possibly gradient in a.u. with masses\n");
      carts.print(2,outfile,0,disp_label,0);
    }
    delete [] disp_label;
    fflush(outfile);

    internals simples(carts, optinfo.simples_present);
    coord = carts.get_coord();
    simples.compute_internals(carts.get_natom(), coord);
    simples.compute_s(carts.get_natom(), coord);
    // simples.print_s();
    free(coord);
    if ( (optinfo.mode != MODE_DISP_LOAD) && (optinfo.mode != MODE_LOAD_REF)) {
      fprintf(outfile,"\nSimple Internal Coordinates and Values\n");
      simples.print(outfile,1);
    }
    fflush(outfile);

    /* obtain symmetry info, including simple transformation matrix */
    get_syminfo(simples);

    /*** If SYMM is not user given, produce SYMM containing delocalized \
     *** internal coordinates or else use redundant simples          ***/
    buffer = new char[MAX_LINELENGTH];
    if (!(optinfo.salcs_present)) {
      if (optinfo.delocalize) {
        fprintf(outfile,"\nForming delocalized internal coordinates.\n");
        // delocalize(carts.get_natom(),simples);
        delocalize(simples, carts);
      }
      else {
        fprintf(outfile,"\nUsing simple, possibly redundant, internal coordinates.\n");
        ffile(&fp_intco, "intco.dat", 2);
        count = 0;
        for( ; ; ) {
          err = fgets(buffer, MAX_LINELENGTH, fp_intco);
          if (err == NULL) break;
          ++count;
        }
        rewind(fp_intco);
        for(j=0;j<(count-1);++j)
          err = fgets(buffer, MAX_LINELENGTH, fp_intco);
        fflush(fp_intco);
        fprintf(fp_intco,"  symm = (\n");
            for (j=0;j<simples.get_num();++j)
            fprintf(fp_intco,"    (\" \" (%d))\n",simples.index_to_id(j));
            fprintf(fp_intco,"  )\n");
        fprintf(fp_intco,")\n");
        fclose(fp_intco);
      }
      /* Add the new intco information to the keyword tree */
      ffile(&fp_intco, "intco.dat", 2);
      if (fp_intco != NULL) {
        ip_append(fp_intco, outfile) ;
        if (ip_exist(":INTCO",0)) {
          user_intcos = 1;
          ip_cwk_add(":INTCO");
        }
        fclose(fp_intco);
      }
    }
    delete [] buffer;
    fflush(outfile);

    // do optimization step by gradients
    if (optinfo.mode == MODE_OPT_STEP) {
      fprintf(outfile," \n ** Taking normal optimization step. **\n");
      salc_set symm_salcs("SYMM");
      if (!optinfo.redundant)
        symm_salcs.print();
      a = opt_step(carts, simples, symm_salcs);
      free_info(simples.get_num());
      exit_io();
      // fprintf(stderr,"Optking returning value %d\n", a);
      return a;
    }

    // only execute a user-given displacements vector
    if (optinfo.mode == MODE_DISP_USER) {
      fprintf(outfile," \n ** Executing user-given displacements vector. **\n");
      salc_set all_salcs;
      all_salcs.print();
      disp_user(carts, simples, all_salcs);
      free_info(simples.get_num());
      exit_io();
      return 0;
    }

    if ((optinfo.mode == MODE_DISP_NOSYMM) || (optinfo.mode == MODE_DISP_IRREP)) {
      // generate unique displacements
      salc_set all_salcs;
      all_salcs.print();
      if (optinfo.mode == MODE_DISP_IRREP) {
        num_disps = make_disp_irrep(carts, simples, all_salcs, optinfo.points);
      }
      else {
        num_disps = make_disp_nosymm(carts, simples, all_salcs,
            optinfo.points);
      }
      free_info(simples.get_num());
      exit_io();
      return(num_disps);
    }

    if (optinfo.mode == MODE_LOAD_REF) {
      fprintf(outfile,"Loading undisplaced reference geometry to chkpt.\n");
      load_ref(carts);
      free_info(simples.get_num());
      exit_io();
      return 0;
    }

    // save the gradient and increment disp_num
    if (optinfo.mode == MODE_GRAD_SAVE) {
      grad_save(carts);
      free_info(simples.get_num());
      exit_io();
      return 0;
    }

    // save the energy and increment disp_num
    if (optinfo.mode == MODE_ENERGY_SAVE) {
      energy_save(carts);
      free_info(simples.get_num());
      exit_io();
      return 0;
    }

    int total_num_disps;

    if (optinfo.mode == MODE_DISP_LOAD) {
      double **micro_geoms, **geom2D;

      open_PSIF();
      psio_read_entry(PSIF_OPTKING, "OPT: Current disp_num",
          (char *) &(optinfo.disp_num), sizeof(int));
      psio_read_entry(PSIF_OPTKING, "OPT: Total num. of disp.",
          (char *) &(total_num_disps), sizeof(int));

      micro_geoms = block_matrix(total_num_disps, dim_carts);
      psio_read_entry(PSIF_OPTKING, "OPT: Displaced geometries",
          (char *) &(micro_geoms[0][0]), total_num_disps*dim_carts*
          sizeof(double));
      close_PSIF();

      geom2D = block_matrix(carts.get_natom(),3);
      for (i=0; i<carts.get_natom(); ++i)
        for (j=0; j<3; ++j)
          geom2D[i][j] = micro_geoms[optinfo.disp_num][3*i+j];

      chkpt_init(PSIO_OPEN_OLD);
      chkpt_wt_geom(geom2D);
      chkpt_close();
      fprintf(outfile,"\n ** Geometry for displacement %d sent to chkpt. **\n", 
              optinfo.disp_num+1);
      free_block(micro_geoms);
      free_block(geom2D);

      free_info(simples.get_num());
      exit_io();
      return(0);
    }

    if (optinfo.mode == MODE_GRAD_ENERGY) {
      fprintf(outfile,"\n ** Calculating file11.dat from energies. **\n");
      salc_set symm_salcs("SYMM");
      symm_salcs.print();
      grad_energy(carts, simples, symm_salcs);
      free_info(simples.get_num());
      exit_io();
      return(0);
    }

    if (optinfo.mode == MODE_FREQ_ENERGY) {
      fprintf(outfile,"\n ** frequencies by energies not yet implemented. **\n");
      free_info(simples.get_num());
      exit_io();
      return(0);
    }

    if ((optinfo.mode==MODE_FREQ_GRAD_IRREP) || (optinfo.mode==MODE_FREQ_GRAD_NOSYMM)) {
      fprintf(outfile,"\n ** Calculating frequencies from gradients. **\n");
      salc_set all_salcs;
      all_salcs.print();
      if (optinfo.mode == MODE_FREQ_GRAD_IRREP)
        freq_grad_irrep(carts, simples, all_salcs, optinfo.points);
      else
        freq_grad_nosymm(carts, simples, all_salcs, optinfo.points);
      free_info(simples.get_num());
      exit_io();
      return(0);
    }

    return(0);
  }

/***  INTRO   prints into ***/
void intro(int argc, char **argv) {
  int i;
  if (optinfo.mode == MODE_OPT_STEP) {
    fprintf(outfile,
            "\n\t------------------------------------------------------\n");
    fprintf(outfile,
              "\t    OPTKING: for internal coordinate optimizations    \n");
    fprintf(outfile,
              "\t------------------------------------------------------\n");
  }
  else {
    fprintf(outfile,"\n******* OPTKING: ");
    for (i=1; i<argc; ++i)
      fprintf(outfile,"%s ", argv[i]);
    fprintf(outfile,"\n");
  }
}

void free_info(int nsimples) {
  int i,j,natom;
  natom = optinfo.natom;

 // free syminfo
  free(syminfo.symmetry);
  for (i=0; i<syminfo.nirreps; ++i) {
    free(syminfo.ct[i]);
    free(syminfo.ict[i]);
    free(syminfo.fict[i]);
    free(syminfo.irrep_lbls[i]);
    free(syminfo.clean_irrep_lbls[i]);
  }
  free(syminfo.ct);
  free(syminfo.ict);
  free(syminfo.fict);
  free(syminfo.irrep_lbls);
  free(syminfo.clean_irrep_lbls);

  for (i=0; i<nsimples; ++i) {
    free(syminfo.ict_ops[i]);
    free(syminfo.ict_ops_sign[i]);
  }
  free(syminfo.ict_ops);
  free(syminfo.ict_ops_sign);

  // free optinfo
  free(optinfo.to_dummy);
  free(optinfo.to_nodummy);
  return;
}


void load_ref(cartesians &carts) {
  int i,j,cnt;
  double *geom, **geom2D,*grad, energy;

  geom = init_array(carts.get_natom() * 3);
  open_PSIF();
  psio_read_entry(PSIF_OPTKING, "OPT: Reference geometry",
           (char *) &(geom[0]), (int) 3*carts.get_natom()*sizeof(double));
  psio_read_entry(PSIF_OPTKING, "OPT: Reference energy",
                (char *) &(energy), sizeof(double));
  // reference gradient assumed 0
  grad = init_array(3*carts.get_natom());
  close_PSIF();

  cnt = 0;
  geom2D = block_matrix(carts.get_natom(), 3);
  for (i=0;i<carts.get_natom();++i)
    for (j=0;j<3;++j)
      geom2D[i][j] = geom[cnt++];

  /*
  for (i=0;i<carts.get_natom();++i)
    for (j=0;j<3;++j)
      fprintf(outfile,"%15.10lf\n", geom2D[i][j]);
      */

  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_geom(geom2D);
  chkpt_wt_grad(grad);
  chkpt_wt_etot(energy);
//  fprintf(outfile,"%15.10lf\n",energy);
  chkpt_close();
}
