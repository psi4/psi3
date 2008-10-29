/*! \file
    \ingroup OPTKING
    \brief main optking function
    OPT.CC Rollin King, begun 1999 
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
--freq_grad_irrep  MODE_FREQ_GRAD_IRREP  use grads in chkpt to compute IRREP freqs
--grad_save        MODE_GRAD_SAVE       save the gradient in chkpt to PSIF list
--energy_save      MODE_ENERGY_SAVE     save the energy in chkpt to PSIF list
--reset_prefix     MODE_RESET_PREFIX    only reset the prefix of PSIF
--disp_num_plus    MODE_DISP_NUM_PLUS   only increment the current displacement number PSIF
--delete_binaries  MODE_DELETE_BINARIES just delete binary files
--test_B           MODE_TEST_BMAT       build B matrix and test numerically
--points           optinfo.points       save the energy in chkpt to PSIF list
--disp_num         optinfo.disp_num     displacement index (by default read from PSIF)
--irrep            optinfo.irrep       the irrep (1,2,...) being displaced or computed

--disp_freq_grad_cart   MODE_DISP_FREQ_GRAD_CART   displaces along SA cart. coordinates
--freq_grad_cart        MODE_FREQ_GRAD_CART        use grads in chkpt to compute freqs
--disp_freq_energy_cart MODE_DISP_FREQ_ENERGY_CART displaces along SA cart. coordinates
--freq_energy_cart      MODE_FREQ_ENERGY_CART      use energies in chkpt to compute freqs
--energy_dat       optinfo.energy_dat   read energy points from the file "energy.dat"
--grad_dat         optinfo.grad_dat   read gradients from file11.dat 
*/

#include <cmath>
#include <cstdio>
#include <libchkpt/chkpt.h>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <libciomr/libciomr.h>
#include <libipv1/ip_lib.h>
#include <physconst.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>	  
#include <psifiles.h>
#include <ccfiles.h>

#include "opt.h"
#include "cartesians.h"
#include "internals.h"
#include "salc.h"
#include "bond_lengths.h"

namespace psi { namespace optking {

void intro(int argc, char **argv);
void free_info(int nsimples);
int *count_internals(cartesians &cart_geom, int intco_given);
extern void zmat_to_intco(void);
extern void get_optinfo();
extern void get_syminfo(internals &simples);
extern void delocalize(internals &simples, cartesians &carts);
extern void disp_user(cartesians &carts, internals &simples, 
                      salc_set &all_salcs);
extern int make_disp_irrep(cartesians &carts, internals &simples, 
                    salc_set &all_salcs);
extern int make_disp_nosymm(cartesians &carts, internals &simples, 
                    salc_set &all_salcs);
extern void load_ref(cartesians &carts);
extern void freq_grad_irrep(cartesians &carts, internals &simples, 
    salc_set &all_salcs);
extern void freq_grad_nosymm(cartesians &carts, internals &simples, 
    salc_set &all_salcs);
extern void grad_energy(cartesians &carts, internals &simples, 
                        salc_set &all_salcs);
extern void grad_save(cartesians &carts);
extern void energy_save(cartesians &carts);
extern int opt_step(cartesians &carts, internals &simples, salc_set &symm_salcs);
extern int opt_step_cart(cartesians &carts, internals &simples, salc_set &symm_salcs);
extern int *read_constraints(internals &simples);

extern int disp_freq_grad_cart(cartesians &carts);
extern void freq_grad_cart(cartesians &carts);
extern int disp_freq_energy_cart(cartesians &carts);
extern void freq_energy_cart(cartesians &carts);
int test_B(cartesians &carts, internals &simples, salc_set &symm);

void chkpt_restart(char *new_prefix);

}} // namespace psi::optking

extern "C" { char *gprgid() { char *prgid = "OPTKING"; return(prgid); } }

using namespace psi::optking;

int main(int argc, char **argv) {

    int i,j,a,b,dim,count,dim_carts,user_intcos, *constraints,xyz;
    int parsed=1, num_disps, disp_length;
    double *f, *coord; 
    char aline[MAX_LINELENGTH], *disp_label, **buffer, *err;

    // Set defaults & read command-line arguments
    optinfo.mode = MODE_OPT_STEP;
    optinfo.disp_num = 0;
    optinfo.points = 3;
    optinfo.points_freq_grad_ints = 3; // always 3 for now
    optinfo.energy_dat = 0;
    optinfo.grad_dat = 0;
    for (i=1; i<argc; ++i) {
      if (!strcmp(argv[i],"--disp_nosymm")) {
        optinfo.mode = MODE_DISP_NOSYMM;
        parsed++;
      }
      else if (!strcmp(argv[i],"--disp_freq_grad_cart")) {
        optinfo.mode = MODE_DISP_FREQ_GRAD_CART;
        parsed++;
      }
      else if (!strcmp(argv[i],"--disp_freq_energy_cart")) {
        optinfo.mode = MODE_DISP_FREQ_ENERGY_CART;
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
      else if (!strcmp(argv[i],"--freq_grad_cart")) {
        optinfo.mode = MODE_FREQ_GRAD_CART;
        parsed++;
      }
      else if (!strcmp(argv[i],"--freq_energy_cart")) {
        optinfo.mode = MODE_FREQ_ENERGY_CART;
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
      else if (!strcmp(argv[i],"--reset_prefix")) {
        optinfo.mode = MODE_RESET_PREFIX;
        parsed++;
      }
      else if (!strcmp(argv[i],"--disp_num_plus")) {
        optinfo.mode = MODE_DISP_NUM_PLUS;
        parsed++;
      }
      else if (!strcmp(argv[i],"--delete_binaries")) {
        optinfo.mode = MODE_DELETE_BINARIES;
        parsed++;
      }
      else if (!strcmp(argv[i],"--test_B")) {
        optinfo.mode = MODE_TEST_BMAT;
        parsed++;
      }
      else if (!strcmp(argv[i],"--points")) {
        sscanf(argv[++i], "%d", &optinfo.points);
        parsed+=2;
      }
      else if (!strcmp(argv[i],"--disp_num")) {
        sscanf(argv[++i], "%d", &optinfo.disp_num);
        parsed+=2;
      }
      else if (!strcmp(argv[i],"--irrep")) {
        sscanf(argv[++i], "%d", &optinfo.irrep);
        optinfo.irrep -= 1;
        parsed+=2;
      }
      else if (!strcmp(argv[i],"--energy_dat")) {
        optinfo.energy_dat = 1;
        parsed++;
      }
      else if (!strcmp(argv[i],"--grad_dat")) {
        optinfo.grad_dat = 1;
        parsed++;
      }
      else {
        //printf("command line argument not understood.\n");
        //exit(1);
		break;
      }
    }

    psi_start(&infile,&outfile,&psi_file_prefix,argc-parsed,argv+parsed,0);
    /* init_in_out() sets the value of "infile", so we need to save it */
    fp_input = infile;
    
    intro(argc, argv);

    ip_cwk_add(":OPTKING");

    // determine if simples and salcs are present in intco.dat 
    optinfo.simples_present = 0;
    optinfo.salcs_present = 0;
    ffile_noexit(&fp_intco, "intco.dat", 2);
    if (fp_intco != NULL)
      ip_append(fp_intco, outfile) ;
    if ( ip_exist(":INTCO",0) ) {
      ip_cwk_add(":INTCO");
      optinfo.simples_present = 1;
    }
    if ( ip_exist("SYMM",0) || ip_exist("ASYMM",0) )
      optinfo.salcs_present = 1;
    if (fp_intco != NULL)
      fclose(fp_intco);

    psio_init(); psio_ipv1_config();
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
    else if ( (optinfo.mode != MODE_DISP_LOAD) && (optinfo.mode != MODE_LOAD_REF)
           && (optinfo.mode != MODE_RESET_PREFIX) && (optinfo.mode != MODE_DISP_NUM_PLUS)
           && (optinfo.mode != MODE_DELETE_BINARIES) ) {
      fprintf(outfile,
      "\nCartesian geometry in a.u. with masses\n");
      carts.print(2,outfile,0,disp_label,0);
    }
    delete [] disp_label;
    fflush(outfile);

    internals simples(carts, optinfo.simples_present);

    /* read in constraints */
    ffile_noexit(&fp_fintco, "fintco.dat", 2);
    if (fp_fintco != NULL) ip_append(fp_fintco, outfile) ;
    if (optinfo.mode == MODE_OPT_STEP) optinfo.constraints = read_constraints(simples);
    if (fp_fintco != NULL) fclose(fp_fintco);

    coord = carts.get_coord();
    simples.compute_internals(carts.get_natom(), coord);
    simples.compute_s(carts.get_natom(), coord);
 //simples.print_s();
    free(coord);
    if ( (optinfo.mode != MODE_DISP_LOAD) && (optinfo.mode != MODE_LOAD_REF)
      && (optinfo.mode != MODE_RESET_PREFIX) && (optinfo.mode != MODE_DISP_NUM_PLUS)
      && (optinfo.mode != MODE_DELETE_BINARIES) ) {
      fprintf(outfile,"\nSimple Internal Coordinates and Values\n");
      simples.print(outfile,1);
    }
    fflush(outfile);

    /* obtain symmetry info, including simple transformation matrix */
    get_syminfo(simples);
    fflush(outfile);

    /*** If SYMM is not user given, produce SYMM containing delocalized \
     *** internal coordinates or else use redundant simples          ***/
    if (!(optinfo.salcs_present)) {
      if (optinfo.delocalize) {
        fprintf(outfile,"\nForming delocalized internal coordinates.\n");
        delocalize(simples, carts);
      }
      else {
        fprintf(outfile,"\nPutting simple, possibly redundant, internal coordinates in intco.dat.\n");
        ffile(&fp_intco, "intco.dat", 2);
        count = 0;
        for( ; ; ) {
          err = fgets(aline, MAX_LINELENGTH, fp_intco);
          if (err == NULL) break;
          ++count;
        }
        rewind(fp_intco);

        /* read all but the last line of the file into memory... */
        buffer = (char **) malloc((count-1) * sizeof(char *));
        for(i=0; i < count-1; i++) {
          buffer[i] = (char *) malloc(MAX_LINELENGTH * sizeof(char));
          err = fgets(buffer[i], MAX_LINELENGTH, fp_intco);
        }
        rewind(fp_intco);

        /* ...and overwite */
        for(i=0; i < count-1; i++) {
          fprintf(fp_intco, "%s", buffer[i]);
          free(buffer[i]);
        }
        free(buffer);

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
    fflush(outfile);

    // do optimization step by gradients
    if (optinfo.mode == MODE_OPT_STEP) {
      fprintf(outfile," \n ** Taking normal optimization step. **\n");
      salc_set symm_salcs("SYMM");
      if (!optinfo.redundant)
        symm_salcs.print();
      if (optinfo.test_B) {
        test_B(carts,simples,symm_salcs);
        coord = carts.get_coord(); // restore internals to undisplaced state
        simples.compute_internals(carts.get_natom(), coord);
        simples.compute_s(carts.get_natom(), coord);
      }
      if (optinfo.cartesian)
        a = opt_step_cart(carts, simples, symm_salcs);
      else
        a = opt_step(carts, simples, symm_salcs);
      free_info(simples.get_num());
      exit_io();
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

    if (optinfo.mode == MODE_DISP_FREQ_GRAD_CART) {
      fprintf(outfile,
      "\n * Performing displacements along symmetry adapted cartesian coordinates *\n\n");
      i = disp_freq_grad_cart(carts);
      free_info(simples.get_num());
      exit_io();
      return i;
    }
    if (optinfo.mode == MODE_DISP_FREQ_ENERGY_CART) {
      fprintf(outfile,
      "\n * Performing displacements along symmetry adapted cartesian coordinates *\n\n");
      i = disp_freq_energy_cart(carts);
      free_info(simples.get_num());
      exit_io();
      return i;
    }


    if ((optinfo.mode == MODE_DISP_NOSYMM) || (optinfo.mode == MODE_DISP_IRREP)) {
      // generate unique displacements
      salc_set all_salcs;
      all_salcs.print();
      if (optinfo.mode == MODE_DISP_IRREP) {
        num_disps = make_disp_irrep(carts, simples, all_salcs);
      }
      else {
        num_disps = make_disp_nosymm(carts, simples, all_salcs);
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
      int *irrep_per_disp;
      psio_address next;

      open_PSIF();
      psio_read_entry(PSIF_OPTKING, "OPT: Current disp_num",
          (char *) &(optinfo.disp_num), sizeof(int));
      psio_read_entry(PSIF_OPTKING, "OPT: Num. of disp.",
          (char *) &(total_num_disps), sizeof(int));

      micro_geoms = block_matrix(total_num_disps, dim_carts);
      psio_read_entry(PSIF_OPTKING, "OPT: Displaced geometries",
          (char *) &(micro_geoms[0][0]), total_num_disps*dim_carts* sizeof(double));

      irrep_per_disp = (int *) malloc(total_num_disps * sizeof(int));
      psio_read_entry(PSIF_OPTKING, "OPT: Irrep per disp",
          (char *) &(irrep_per_disp[0]), total_num_disps*sizeof(int));
      close_PSIF();

      chkpt_restart( syminfo.irrep_lbls[ irrep_per_disp[optinfo.disp_num]] );

      geom2D = block_matrix(carts.get_natom(),3);
      for (i=0; i<carts.get_natom(); ++i)
        for (j=0; j<3; ++j)
          geom2D[i][j] = micro_geoms[optinfo.disp_num][3*i+j];

      chkpt_init(PSIO_OPEN_OLD);

      chkpt_wt_disp_irrep( irrep_per_disp[optinfo.disp_num] );

      /* set flag to tell cscf to use checkpoint occupations (ignoring DOCC, etc.);
       these are determined and written by input */
      chkpt_wt_override_occ(1);

      chkpt_wt_geom(geom2D);
      chkpt_close();
      fprintf(outfile,"\n ** Geometry for displacement %d sent to chkpt. **\n", 
              optinfo.disp_num+1);
      free_block(micro_geoms);

      print_mat(geom2D,carts.get_natom(),3,outfile);

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
        freq_grad_irrep(carts, simples, all_salcs);
      else
        freq_grad_nosymm(carts, simples, all_salcs);
      free_info(simples.get_num());
      exit_io();
      return(0);
    }

    if (optinfo.mode==MODE_FREQ_GRAD_CART) {
      fprintf(outfile,"\n ** Calculating frequencies from gradients. **\n");
      freq_grad_cart(carts);
      free_info(simples.get_num());
      exit_io();
      return(0);
    }
    if (optinfo.mode==MODE_FREQ_ENERGY_CART) {
      fprintf(outfile,"\n ** Calculating frequencies from energies. **\n");
      freq_energy_cart(carts);
      free_info(simples.get_num());
      exit_io();
      return(0);
    }

    if (optinfo.mode == MODE_TEST_BMAT) {
      fprintf(outfile," \n ** Testing B matrix **\n");
      salc_set symm_salcs("SYMM");
      test_B(carts,simples,symm_salcs);
      free_info(simples.get_num());
      exit_io();
      return (0);
    }

    if (optinfo.mode == MODE_RESET_PREFIX) {
      /* delete CC temporary files */
      if((strcmp(optinfo.wfn, "MP2")==0)       || (strcmp(optinfo.wfn, "CCSD")==0)  ||
         (strcmp(optinfo.wfn, "CCSD_T")==0)    || (strcmp(optinfo.wfn, "EOM_CCSD")==0)  ||
         (strcmp(optinfo.wfn, "LEOM_CCSD")==0) || (strcmp(optinfo.wfn, "BCCD")==0)  ||
         (strcmp(optinfo.wfn,"BCCD_T")==0)     || (strcmp(optinfo.wfn, "SCF")==0)  ||
         (strcmp(optinfo.wfn,"CIS")==0)        || (strcmp(optinfo.wfn,"RPA")==0)  ||
         (strcmp(optinfo.wfn,"CC2")==0)        || (strcmp(optinfo.wfn,"CC3")==0)  ||
         (strcmp(optinfo.wfn,"EOM_CC3")==0) ) {
          fprintf(outfile, "Deleting binary files\n");
          for(i=CC_MIN; i <= CC_MAX; i++) {
            psio_open(i,1); psio_close(i,0);
          }
          psio_open(35,1); psio_close(35,0);
          psio_open(72,1); psio_close(72,0);
      }
      else if ( (strcmp(optinfo.wfn, "DETCI")==0) ) {
        fprintf(outfile, "Deleting CI binary files\n");
        psio_open(35,1); psio_close(35,0);
        psio_open(72,1); psio_close(72,0);
        psio_open(50,1); psio_close(50,0);
        psio_open(51,1); psio_close(51,0);
        psio_open(52,1); psio_close(52,0);
        psio_open(53,1); psio_close(53,0);
      }
      fprintf(outfile,"Resetting checkpoint prefix.\n");
      chkpt_init(PSIO_OPEN_OLD);
      chkpt_reset_prefix();
      chkpt_commit_prefix();
      chkpt_close();
    }

    if (optinfo.mode == MODE_DELETE_BINARIES) {
      if((strcmp(optinfo.wfn, "MP2")==0)       || (strcmp(optinfo.wfn, "CCSD")==0)  ||
         (strcmp(optinfo.wfn, "CCSD_T")==0)    || (strcmp(optinfo.wfn, "EOM_CCSD")==0)  ||
         (strcmp(optinfo.wfn, "LEOM_CCSD")==0) || (strcmp(optinfo.wfn, "BCCD")==0)  ||
         (strcmp(optinfo.wfn,"BCCD_T")==0)     || (strcmp(optinfo.wfn, "SCF")==0)  ||
         (strcmp(optinfo.wfn,"CIS")==0)        || (strcmp(optinfo.wfn,"RPA")==0)  ||
         (strcmp(optinfo.wfn,"CC2")==0)        || (strcmp(optinfo.wfn,"CC3")==0)  ||
         (strcmp(optinfo.wfn,"EOM_CC3")==0) ) {
          fprintf(outfile, "Deleting binary files\n");
          for(i=CC_MIN; i <= CC_MAX; i++) {
            psio_open(i,1); psio_close(i,0);
          }
          psio_open(35,1); psio_close(35,0);
          psio_open(72,1); psio_close(72,0);
      }
      else if ( (strcmp(optinfo.wfn, "DETCI")==0) ) {
        fprintf(outfile, "Deleting CI binary files\n");
        psio_open(35,1); psio_close(35,0);
        psio_open(72,1); psio_close(72,0);
        psio_open(50,1); psio_close(50,0);
        psio_open(51,1); psio_close(51,0);
        psio_open(52,1); psio_close(52,0);
        psio_open(53,1); psio_close(53,0);
      }
    }

    /* increment current displacement number */
    if (optinfo.mode == MODE_DISP_NUM_PLUS) {
      open_PSIF();
      psio_read_entry(PSIF_OPTKING, "OPT: Current disp_num",
                (char *) &(optinfo.disp_num), sizeof(int));
      optinfo.disp_num += 1;
      psio_write_entry(PSIF_OPTKING, "OPT: Current disp_num",
                (char *) &(optinfo.disp_num), sizeof(int));
      close_PSIF();
    }
    return(0);
}

namespace psi { namespace optking {

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
  free_int_matrix(syminfo.ct);
  free_int_matrix(syminfo.ict);
  free_int_matrix(syminfo.fict);
  for (i=0; i<syminfo.nirreps; ++i) {
    free(syminfo.irrep_lbls[i]);
    free(syminfo.clean_irrep_lbls[i]);
  }
  free(syminfo.irrep_lbls);
  free(syminfo.clean_irrep_lbls);
  
  free_int_matrix(syminfo.ict_ops);
  free_int_matrix(syminfo.ict_ops_sign);

  // free optinfo
  free(optinfo.to_dummy);
  free(optinfo.to_nodummy);
  if (optinfo.nfragment > 1) {
    free(optinfo.natom_per_fragment);
    free(optinfo.nallatom_per_fragment);
    free(optinfo.nref_per_fragment);
    free(optinfo.fragment_coeff);
  }
  return;
}


void load_ref(cartesians &carts) {
  int i,j,cnt;
  double *geom, **geom2D,*grad, energy;

  /*
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
  fprintf(outfile,"Geometry written.\n");
  print_mat(geom2D,carts.get_natom(),3,outfile);
  */

  chkpt_reset_prefix();

  /*
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_geom(geom2D);
  chkpt_wt_grad(grad);
  chkpt_wt_etot(energy);
  chkpt_close();
  */
}

void chkpt_restart(char *new_prefix) {
  int nallatom, natom, *atom_dummy;
  double **fgeom, *zvals;
  char **felement;

  fprintf(outfile,"\nSetting chkpt prefix to irrep %s.\n",new_prefix);

  chkpt_init(PSIO_OPEN_OLD);
  nallatom = chkpt_rd_nallatom();
  natom = chkpt_rd_natom();
  fgeom = chkpt_rd_fgeom();
  atom_dummy = chkpt_rd_atom_dummy();
  zvals = chkpt_rd_zvals();
  felement = chkpt_rd_felement();

  chkpt_set_prefix(new_prefix);
  chkpt_commit_prefix();
  chkpt_close();
                                                                                                                
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_nallatom(nallatom);
  chkpt_wt_natom(natom);
  chkpt_wt_fgeom(fgeom);
  chkpt_wt_atom_dummy(atom_dummy);
  chkpt_wt_zvals(zvals);
  chkpt_wt_felement(felement);
  chkpt_close();
  return;
}

}} /* namespace psi::optking */
