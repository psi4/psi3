#include<stdio.h>
#include<math.h>
#include<libipv1/ip_lib.h>
#include<libint/libint.h>

#include"defines.h"
#define EXTERN
#include"global.h"
#ifdef INCLUDE_Fock
 #include"scf_parsing.h"
#endif

void parsing()
{
  int errcod;
  int cutoff_exp;
  long int max_bytes;
  
  UserOptions.print_lvl = 1;
  errcod = ip_data("PRINT","%d",&(UserOptions.print_lvl),0);

  /*--- This piece of code from CPHF by Ed Seidl ---*/
  if (ip_exist("MEMORY", 0)) {
    fndcor(&max_bytes, infile, outfile);
    UserOptions.max_memory = max_bytes / sizeof(double);
  }
  else
    UserOptions.max_memory = MAX_NUM_DOUBLES;
  UserOptions.memory = UserOptions.max_memory;
  
  cutoff_exp = CUTOFF;
  errcod = ip_data("CUTOFF","%d",&cutoff_exp,0);
  UserOptions.cutoff = 1.0/pow(10.0,(double)cutoff_exp);

  UserOptions.make_oei = 1;
  UserOptions.make_fock = 0;
  UserOptions.make_eri = 1;
  UserOptions.symm_ints = 1;
  errcod = ip_boolean("MAKE_ERI",&(UserOptions.make_eri),0);

  errcod = ip_data("S_FILE","%d",&(IOUnits.itapS),0);
  errcod = ip_data("T_FILE","%d",&(IOUnits.itapT),0);
  errcod = ip_data("V_FILE","%d",&(IOUnits.itapV),0);
  if (UserOptions.make_eri)
    errcod = ip_data("ERI_FILE","%d",&(IOUnits.itap33),0);

  UserOptions.scf_only = 0;
  errcod = ip_string("WFN",&UserOptions.wfn,0);
  if (UserOptions.wfn == NULL)
    punt("Keyword WFN is missing");
  if (!strcmp("SCF",UserOptions.wfn))
    UserOptions.scf_only = 1;

  UserOptions.num_threads = 1;
  errcod = ip_data("NUM_THREADS","%d",&(UserOptions.num_threads),0);
  if (UserOptions.num_threads < 1)
    UserOptions.num_threads = 1;

  UserOptions.restart = 0;
  errcod = ip_boolean("RESTART",&UserOptions.restart,0);
  if (UserOptions.restart) {
    errcod = ip_data("RESTART_TASK","%d",&UserOptions.restart_task,0);
    if (UserOptions.restart_task < 0)
      punt("RESTART_TASK < 0");
  }

  return;

}

void parsing_cmdline(int argc, char *argv[])
{
   int i, errcod;
   char *refstring;

   /* try to detect if we've called CINTS by itself but don't 
      really need it, as in a direct SCF calculation
    */
   if (argc < 2) {
     ip_cwk_add(":SCF");
     ip_cwk_add(":CSCF");
     errcod = ip_boolean("DIRECT_SCF",&i,0);
     if (errcod == IPE_OK && i == 1) {
        UserOptions.print_lvl = 0;
        stop_io();
        exit(0);
     }
   }

   for (i=1; i<argc; i++) {
       
       /*--- build Fock option ---*/
       if (strcmp(argv[i], "--fock") == 0) {
#ifdef INCLUDE_Fock
	   scf_parsing();
           /*UserOptions.make_oei = 0;
	   UserOptions.make_eri = 0;
	   UserOptions.make_fock = 1;
	   UserOptions.print_lvl = 0;
	   UserOptions.symm_ints = 0;
	   UserOptions.make_dft = 0;
	   errcod = ip_string("REFERENCE",&refstring,0);
	   if (errcod != IPE_OK)
	     punt("REFERENCE keyword is missing");
	   else if (!strcmp(refstring,"RHF") || !strcmp(refstring,""))
	     UserOptions.reftype = rhf;
	   else if (!strcmp(refstring,"ROHF"))
	     UserOptions.reftype = rohf;
	   else if (!strcmp(refstring,"UHF"))
	     UserOptions.reftype = uhf;
	   else if (!strcmp(refstring,"RKS")){
	       UserOptions.reftype = rhf;
	       UserOptions.make_dft = 1;
	   }
	   else if (!strcmp(refstring,"UKS")){
	       UserOptions.reftype = uhf;
	       UserOptions.make_dft = 1;
	   }
	   else
	   punt("The specified REFERENCE not implemented");*/
#else
	   punt("--fock option is not supported by your CINTS executable.\nRecompile the code including files in Fock subdirectory.");
#endif
	   return;
       }

       /*--- build oeints option ---*/
       if (strcmp(argv[i], "--oeints") == 0) {
#ifdef INCLUDE_Default_Ints
           UserOptions.make_oei = 1;
	   UserOptions.make_eri = 0;
	   UserOptions.make_fock = 0;
	   UserOptions.print_lvl = 0;
	   UserOptions.symm_ints = 1;
	   UserOptions.num_threads = 1;
#else
	   punt("--oeints option is not supported by your CINTS executable.\nRecompile the code including files in Default_Ints subdirectory.");
#endif
	   return;
       }

       /*--- build ERIs option ---*/
       if (strcmp(argv[i], "--teints") == 0) {
#ifdef INCLUDE_Default_Ints
           UserOptions.make_oei = 0;
	   UserOptions.make_eri = 1;
	   UserOptions.make_fock = 0;
	   UserOptions.print_lvl = 0;
	   UserOptions.symm_ints = 1;
	   UserOptions.num_threads = 1;
#else
	   punt("--teints option is not supported by your CINTS executable.\nRecompile the code including files in Default_Ints subdirectory.");
#endif
	   return;
       }

       /*--- compute 1st derivatives option ---*/
       if (strcmp(argv[i], "--deriv1") == 0) {
#ifdef INCLUDE_Default_Deriv1
           UserOptions.make_oei = 0;
	   UserOptions.make_eri = 0;
	   UserOptions.make_fock = 0;
	   UserOptions.make_deriv1 = 1;
	   UserOptions.symm_ints = 0;
	   errcod = ip_string("DERTYPE",&UserOptions.dertype,0);
	   if (UserOptions.dertype == NULL)
	     punt("Keyword DERTYPE is missing");
	   if (strcmp(UserOptions.dertype,"FIRST"))
	     punt("Only DERTYPE=FIRST can be presently handled");
	   if (!strcmp("SCF",UserOptions.wfn)) {
	     errcod = ip_string("REFERENCE",&refstring,0);
	     if (errcod != IPE_OK)
	       punt("REFERENCE keyword is missing");
	     else if (!strcmp(refstring,"RHF") || !strcmp(refstring,""))
	       UserOptions.reftype = rhf;
	     else if (!strcmp(refstring,"ROHF"))
	       UserOptions.reftype = rohf;
	     else if (!strcmp(refstring,"UHF"))
	       UserOptions.reftype = uhf;
	     else if (!strcmp(refstring,"TWOCON"))
	       UserOptions.reftype = twocon;
	     else
	       punt("SCF gradients with specified REFERENCE not implemented");
	   }
	   else
	     UserOptions.num_threads = 1;
#else
	   punt("--deriv1 option is not supported by your CINTS executable.\nRecompile the code including files in Default_Deriv1 subdirectory.");
#endif
	   return;
       }

       /*--- compute 2nd derivatives ---*/
       if(!strcmp(argv[i], "--deriv2")) {
	 UserOptions.make_oei = 0;
	 UserOptions.make_eri = 0;
	 UserOptions.make_fock = 0;
	 UserOptions.symm_ints = 0;
	 UserOptions.make_deriv2 = 1;
	 return;
       }

       /*--- compute one-electron property integrals option ---*/
       if (strcmp(argv[i], "--oeprop") == 0) {
#ifdef INCLUDE_OEProp_Ints
           UserOptions.make_oei = 0;
	   UserOptions.make_eri = 0;
	   UserOptions.make_fock = 0;
	   UserOptions.make_oeprop = 1;
	   UserOptions.symm_ints = 0;
#else
	   punt("--oeprop option is not supported by your CINTS executable.\nRecompile the code including files in OEProp_Ints subdirectory.");
#endif
	   return;
       }
       
       /*--- compute MP2 energy option ---*/
       if (strcmp(argv[i], "--mp2") == 0) {
#ifdef INCLUDE_MP2
           UserOptions.make_oei = 0;
	   UserOptions.make_eri = 0;
	   UserOptions.make_fock = 0;
	   UserOptions.make_deriv1 = 0;
	   UserOptions.make_mp2 = 1;
	   UserOptions.symm_ints = 0;
	   if (!strcmp("MP2",UserOptions.wfn)) {
	     errcod = ip_string("REFERENCE",&refstring,0);
	     if (errcod != IPE_OK)
	       punt("REFERENCE keyword is missing");
	     else if (!strcmp(refstring,"RHF") || !strcmp(refstring,""))
	       UserOptions.reftype = rhf;
	     else if (!strcmp(refstring,"UHF")) {
	       UserOptions.reftype = uhf;
	       punt("UMP2 energy evaluation is not yet implemented");
	     }
	     else
	       punt("MP2 energy evaluation with specified REFERENCE not implemented");
	   }
	   else
	     punt("MP2 energy is requested, but WFN != MP2 in your input file");
#else
	   punt("--mp2 option is not supported by your CINTS executable.\nRecompile the code including files in MP2 subdirectory.");
#endif
	   return;
       }

       /*--- build te integrals for R12 methods option ---*/
       if (strcmp(argv[i], "--r12ints") == 0) {
#ifdef INCLUDE_R12_Ints
           UserOptions.make_oei = 0;
	   UserOptions.make_eri = 0;
	   UserOptions.make_fock = 0;
	   UserOptions.make_deriv1 = 0;
	   UserOptions.make_mp2 = 0;
	   UserOptions.make_r12ints = 1;
	   UserOptions.symm_ints = 1;
	   UserOptions.num_threads = 1;
#else
	   punt("--r12ints option is not supported by your CINTS executable.\nRecompile the code including files in R12_Ints subdirectory.");
#endif
	   return;
       }

       /*--- compute MP2-R12 energy option ---*/
       if (strcmp(argv[i], "--mp2r12") == 0) {
#ifdef INCLUDE_MP2R12
           UserOptions.make_oei = 0;
	   UserOptions.make_eri = 0;
	   UserOptions.make_fock = 0;
	   UserOptions.make_deriv1 = 0;
	   UserOptions.make_mp2 = 0;
	   UserOptions.make_mp2r12 = 1;
	   UserOptions.symm_ints = 0;
	   if (!strcmp("MP2",UserOptions.wfn)) {
	     errcod = ip_string("REFERENCE",&refstring,0);
	     if (errcod != IPE_OK)
	       punt("REFERENCE keyword is missing");
	     else if (!strcmp(refstring,"RHF") || !strcmp(refstring,""))
	       UserOptions.reftype = rhf;
	     else
	       punt("Direct MP2-R12/A integrals transformation with specified REFERENCE not implemented");
	   }
	   else
	     punt("Direct MP2-R12/A integrals transformation is requested, but WFN != MP2 in your input file");
#else
	   punt("--mp2r12 option is not supported by your CINTS executable.\nRecompile the code including files in MP2R12 subdirectory.");
#endif
	   return;
       }
   }

   return;
}


