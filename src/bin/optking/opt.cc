/*** OPT.CC Rollin King, 1999, 2002 ***/

/*
command-line      internal specifier   what it does
--opt_step        MODE_OPT_STEP        take an ordinary geometry optimization step
--disp_symm       MODE_DISP_SYMM       displaces along all symmetric internals
--disp_all        MODE_DISP_ALL        displaces along all internals
--disp_load       MODE_DISP_LOAD       load a previous displacement from PSIF to chkpt
--freq_energy     MODE_FREQ_ENERGY     not yet supported
--grad_energy     MODE_GRAD_ENERGY     use energies in chkpt to compute a gradient
--freq_grad       MODE_FREQ_GRAD       use gradients in chkpt to compute frequencies
--freq_grad_symm  MODE_FREQ_GRAD_SYMM  use gradients in chkpt to compute symm freqs
--grad_save       MODE_GRAD_SAVE       save the gradient in chkpt to PSIF list
--energy_save     MODE_ENERGY_SAVE     save the energy in chkpt to PSIF list
--disp_num                             displacement index (by default read from PSIF)
--points                               3 or 5 pt formula (5pt not yet supported)
*/

  extern "C" {
#include <stdio.h>
#include <libchkpt/chkpt.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

  void intro();
  void free_info(int nsimples);
  int *count_internals(cartesians &cart_geom, int intco_given);
  extern "C" void zmat_to_intco(void);
  extern "C" void get_optinfo();
  extern void get_syminfo(internals &simples);
  extern void delocalize(internals &simples, cartesians &carts);
  extern void disp_user(cartesians &carts, internals &simples, 
                        salc_set &all_salcs);
  extern int disp_all(cartesians &carts, internals &simples, 
                      salc_set &all_salcs, int points);
  extern void freq_grad(cartesians &carts, internals &simples, 
                        salc_set &all_salcs);
  extern void grad_energy(cartesians &carts, internals &simples, 
                          salc_set &all_salcs);
  extern void grad_save(cartesians &carts);
  extern void energy_save(cartesians &carts);
  extern int opt_step(cartesians &carts, internals &simples, salc_set &symm);

  int main(int argc, char *argv[]) {

    int i,j,a,b,dim,count,dim_carts,user_intcos;
    int parsed=1;
    int num_disps, disp_length, *number_internals;
    double *fcoord, *f;
    char *disp_label, *buffer, *err;

    optinfo.mode = MODE_OPT_STEP;
    optinfo.disp_num = 0;
    optinfo.points = 3;
    for (i=1; i<argc; ++i) {
      if (!strcmp(argv[i],"--disp_symm")) {
        optinfo.mode = MODE_DISP_SYMM;
	parsed++;
      }
      else if (!strcmp(argv[i],"--disp_all")) {
        optinfo.mode = MODE_DISP_ALL;
	parsed++;
      }
      else if (!strcmp(argv[i],"--disp_load")) {
        optinfo.mode = MODE_DISP_LOAD;
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
      else if (!strcmp(argv[i],"--freq_grad")) {
        optinfo.mode = MODE_FREQ_GRAD;
	parsed++;
      }
      else if (!strcmp(argv[i],"--freq_grad_symm")) {
        optinfo.mode = MODE_FREQ_GRAD_SYMM;
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
      else {
        printf("command line argument not understood.\n");
        exit(1);
      }
    }
    // printf("optinfo.mode %d\n",optinfo.mode);
    // printf("optinfo.points %d\n",optinfo.points);

    psi_start(argc-parsed,argv+parsed,0);
    /* init_in_out() sets the value of "infile", so we need to save it */
    fp_input = infile;
    
    intro();

    ip_cwk_add(":OPTKING");

    // determine if simples and salcs are present in intco.dat 
    optinfo.simples_present = 0;
    optinfo.salcs_present = 0;
    if ((fp_intco = fopen("intco.dat","r")) != NULL) {
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

    // generate simples from zmat if optinfo.z_matrix == 1
    if (optinfo.zmat) zmat_to_intco();

    disp_label = new char[MAX_LINELENGTH];
    /* loads up geometry, (possibly gradient and energy) from chkpt */
    cartesians carts;
    dim_carts = carts.get_nallatom()*3;
    fprintf(outfile,
      "\nCartesian geometry and gradient (if available) in a.u. with masses\n");
    if (optinfo.mode == MODE_GRAD_ENERGY || optinfo.mode == MODE_DISP_SYMM ||
        optinfo.mode == MODE_DISP_ALL || optinfo.mode == MODE_DISP_LOAD ||
        optinfo.mode == MODE_DISP_USER || optinfo.mode == MODE_ENERGY_SAVE ) 
      carts.print(2,outfile,0,disp_label,0);
    else
      carts.print(3,outfile,0,disp_label,0);

    delete [] disp_label;
    fflush(outfile);

    /* the following function and constructor have a lot of identical code */
    number_internals = count_internals(carts, optinfo.simples_present); 
    //  fprintf(outfile,"stre: %d bend: %d tors: %d oop: %d\n",
    //          number_internals[0], number_internals[1], number_internals[2], 
    //          number_internals[3]);
    internals simples(carts, optinfo.simples_present,number_internals);
    // simples.compute_internals(carts.get_natom(),carts.get_coord());
    fcoord = carts.get_fcoord();
    simples.compute_internals(carts.get_nallatom(), fcoord);
    //simples.compute_s(carts.get_natom(),carts.get_coord() );
    simples.compute_s(carts.get_nallatom(), fcoord);
    free(fcoord);
    fprintf(outfile,"\nSimple Internal Coordinates and Values\n");
    simples.print(outfile,1);
    fflush(outfile);

    /* obtain symmetry info, including simple transformation matrix */
    get_syminfo(simples);

    /*** If SYMM is not user given, produce SYMM containing delocalized \
     *** internal coordinates or else use redundant simples          ***/
    buffer = new char[MAX_LINELENGTH];
    i = 1;
    if (ip_exist("SYMM",0)) i = 0;
    if (i) {
      if (optinfo.delocalize) {
        fprintf(outfile,"\nForming delocalized internal coordinates.\n");
        // delocalize(carts.get_natom(),simples);
        delocalize(simples, carts);
      }
      else {
        fprintf(outfile,"\nUsing simple, possibly redundant, internal \
            coordinates.\n");
        fp_intco = fopen("intco.dat", "r+");
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
      fp_intco = fopen("intco.dat","r");
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
      salc_set symm("SYMM");
      symm.print();
      a = opt_step(carts, simples, symm);
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

    if (optinfo.mode == MODE_DISP_SYMM) {
      fprintf(outfile,
              " \n ** Displacing symmetric modes for 3-point formula\n");
      salc_set symm("SYMM");
      symm.print();
      // disp_all is therefore misnamed
      num_disps = disp_all(carts, simples, symm, optinfo.points);
      free_info(simples.get_num());
      exit_io();
      return(num_disps);
    }

    if (optinfo.mode == MODE_DISP_ALL) {
      fprintf(outfile,"\n ** Displacing all modes for 3-point formula\n");
      salc_set all_salcs;
      all_salcs.print();
      num_disps = disp_all(carts, simples, all_salcs, optinfo.points);
      free_info(simples.get_num());
      exit_io();
      return(num_disps);
    }

    int total_num_disps;

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

      geom2D = block_matrix(carts.get_nallatom(),3);
      for (i=0; i<carts.get_nallatom(); ++i)
        for (j=0; j<3; ++j)
          geom2D[i][j] = micro_geoms[optinfo.disp_num][3*i+j];

      chkpt_init(PSIO_OPEN_OLD);
      chkpt_wt_fgeom(geom2D);
      chkpt_close();
      fprintf(outfile,"\n ** Geometry for displacement %d sent to chkpt. **\n", 
              optinfo.disp_num);
      free_block(micro_geoms);
      free_block(geom2D);

      free_info(simples.get_num());
      exit_io();
      return(0);
    }

    if (optinfo.mode == MODE_GRAD_ENERGY) {
      fprintf(outfile,"\n ** Calculating file11.dat from energies. **\n");
      salc_set symm("SYMM");
      symm.print();
      grad_energy(carts, simples, symm);
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

    if (optinfo.mode == MODE_FREQ_GRAD) {
      fprintf(outfile,"\n ** Calculating frequencies from gradients. **\n");
      salc_set all_salcs;
      all_salcs.print();
      freq_grad(carts, simples, all_salcs);
      free_info(simples.get_num());
      exit_io();
      return(0);
    }

    if (optinfo.mode == MODE_FREQ_GRAD_SYMM) {
      fprintf(outfile,"\n ** Calculating frequencies from gradients. **\n");
      salc_set symm("SYMM");
      symm.print();
      freq_grad(carts, simples, symm);
      free_info(simples.get_num());
      exit_io();
      return 0;
    }

    free(number_internals);
    return(0);

    free(number_internals);
    return(0);
  }

  /***  INTRO   prints into ***/
  void intro() {
    fprintf(outfile,
            "\t________________________________________________________\n");
    fprintf(outfile,
            "\t|                                                      |\n");
    fprintf(outfile,
            "\t| OPTKING: For internal coordinate optimizations and   |\n");
    fprintf(outfile,
            "\t| and vibrational frequencies, Rollin King, 1999, 2002 |\n");
    fprintf(outfile,
            "\t|______________________________________________________|\n");
  }


  /*** COUNT_INTERNALS : counts the simple internals present in intco or
   * that will be generated - code reproduces lots of stuff done in the
   * internals constructor ***/

  int *count_internals(cartesians &cart_geom, int intco_given) {

    int i,j,k,l,a,b, natom, num, count, Zmax,Zmin;
    int *ioff, *count_array, **bonds;            
    double *coord, *distance;     

    count_array = init_int_array(4);

    if (intco_given) {
      num=0;
      if (ip_exist("STRE",0)) ip_count("STRE",&num,0);
      count_array[0]=num;

      num=0;
      if (ip_exist("BEND",0)) ip_count("BEND",&num,0);
      count_array[1]=num;

      num=0;
      if (ip_exist("TORS",0)) ip_count("TORS",&num,0);
      count_array[2]=num;

      num=0;
      if (ip_exist("OUT",0)) ip_count("OUT",&num,0);
      count_array[3]=num;
    }
    else {
      natom = cart_geom.get_natom();
      coord = cart_geom.get_coord();

      ioff = (int *) malloc (32641 * sizeof(int));
      ioff[0]=0;
      for (i=1;i<32641;++i)
        ioff[i] = ioff[i-1] + i;

      /* compute distances */
      distance = init_array( ((natom+1)*natom)/2 );
      count = -1;
      for(i=0;i<natom;++i) {
        for(j=0;j<=i;++j) {
          distance[++count] = sqrt(SQR(coord[3*i+0]-coord[3*j+0]) +
              SQR(coord[3*i+1]-coord[3*j+1]) +
              SQR(coord[3*i+2]-coord[3*j+2]));
        }
      }

      /* determine bonds */
      bonds = init_int_matrix(natom,natom);
      for (i=0;i<natom;++i) {
        for(j=0;j<i;++j) {
          Zmax = MAX((int)cart_geom.get_atomic_num(i),(int)cart_geom.get_atomic_num(j));
          Zmin = MIN((int)cart_geom.get_atomic_num(i),(int)cart_geom.get_atomic_num(j));
          a = ioff[Zmax-1] + (Zmin-1);
          if (bondl[a] != 0.0) {
            if (distance[ioff[i]+j]<(optinfo.scale_connectivity * bondl[a])) {
              bonds[i][j] = 1;
              bonds[j][i] = 1;
            }
          }
          else {
            fprintf(outfile,"WARNING! Optking does not know what bond lengths");
            fprintf(outfile,"to expect for all the atoms.\n");
            fprintf(outfile,"You may have to specify connectivity in input.");
          }
        }
      }

      /* check input for user specified bonds or nobonds */
      if (ip_exist("BONDS",0)) {
        ip_count("BONDS",&num,0);
        for(i=0;i<num;++i) {
          ip_data("BONDS","%d",&a,2,i,0);
          ip_data("BONDS","%d",&b,2,i,0);
          a -= 1;  b -= 1;
          bonds[a][b] = 1;
          bonds[b][a] = 1;
        }
      }

      if (ip_exist("NOBONDS",0)) {
        ip_count("NOBONDS",&num,0);
        for(i=0;i<num;++i) {
          ip_data("NOBONDS","%d",&a,2,i,0);
          ip_data("NOBONDS","%d",&b,2,i,0);
          a -= 1;  b -= 1;
          bonds[a][b] = 0;
          bonds[b][a] = 0;
        }
      }

      /* count number of bonds */
      num=0;
      for(i=0; i<natom; ++i)
        for(j=i+1; j<natom; ++j)
          if(bonds[i][j] == 1) ++num;
      count_array[0]=num;

      /* count number of bends */
      num=0;
      for(i=0;i<natom;++i) {
        for(j=0;j<natom;++j) {
          if(i!=j) {
            for(k=i+1; k<natom; ++k) {
              if(j!=k) {
                if (bonds[i][j] && bonds[j][k]) ++num;
              }
            }
          }
        }
      }
      count_array[1]=num;

      /* count number of torsions */
      num=0;
      for(i=0;i<natom;++i) {
        for(j=0;j<natom;++j) {
          if(i!=j) {
            for(k=0; k<natom; ++k) {
              if((i != k) && (j != k)) {
                for(l=i+1; l<natom; ++l) {
                  if((l != j) && (l != k) && bonds[i][j] && bonds[j][k] && 
                     bonds[k][l]) {
                    ++num;
                  }
                }
              }
            }
          }
        }
      }
      count_array[2]=num;
      count_array[3]=0; // no out of planes

      free(ioff);
      free(distance);
      free(coord);
      free_int_matrix(bonds,natom);
    }
    return count_array;
  }


  void free_info(int nsimples) {
    int i,j,nallatom,natom;
    nallatom = optinfo.nallatom;
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
