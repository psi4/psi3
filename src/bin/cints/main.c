#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<ip_libv1.h>
#include<libciomr.h>
#include<file30.h>
#include<psio.h>
#include<psifiles.h>
#include<libint.h>

#include"defines.h"
#include"global.h"

#include"small_fns.h"
#include"parsing.h"
#include"prints.h"
#include"molecule.h"
#include"basisset.h"
#include"symm.h"
#include"dcr.h"
#include"gto.h"
#ifdef INCLUDE_Default_Ints
 #include"enuc.h"
 #include"oe_ints.h"
 #include"te_ints.h"
#endif
#ifdef INCLUDE_HF_Fock
 #include"fock.h"
#endif
#ifdef INCLUDE_Default_Deriv1
 #include"deriv1.h"
#endif
#ifdef INCLUDE_MP2
 #include"mp2.h"
#endif
#ifdef INCLUDE_R12_Ints
 #include"r12_te_ints.h"
#endif
#ifdef INCLUDE_MP2R12
 #include"mp2r12.h"
#endif

/*-------------------------------
  External functions declaration
 -------------------------------*/
char *gprgid();
void init_globals();

int main(int argc, char *argv[])
{
   /*--- Local variables ---*/
   int i,j,k,l,m,count;

   init_globals();
   start_io();
   parsing();
   /*--- Parse the command line ---*/
   parsing_cmdline(argc,argv);
   print_intro();
   setup();
   
   /*--- Prepare data ---*/
   init_molecule();
   init_symmetry();
   init_basisset();
   init_dcr();
   init_gto();

   /*--- Print out some stuff ---*/
   print_scalars();
   print_basisset();

   /*--- Compute the integrals ---*/
#ifdef INCLUDE_Default_Ints
   if (UserOptions.make_oei) {
     /* Molecule.Enuc = */ compute_enuc();
     file30_wt_enuc(Molecule.Enuc);
     oe_ints();
   }
   if (UserOptions.make_eri)
     te_ints();
#endif
#ifdef INCLUDE_HF_Fock
   if (UserOptions.make_fock)
     fock();
#endif
   /*--- Compute the derivative integrals ---*/
#ifdef INCLUDE_Default_Deriv1
   if (UserOptions.make_deriv1)
     deriv1();
#endif
#ifdef INCLUDE_MP2
   if (UserOptions.make_mp2)
     mp2();
#endif
#ifdef INCLUDE_R12_Ints
   if (UserOptions.make_r12ints)
     r12_te_ints();
#endif
#ifdef INCLUDE_MP2R12
   if (UserOptions.make_mp2r12)
     mp2r12();
#endif
   
   /*--- Cleanup ---*/
   cleanup_gto();
   cleanup_symmetry();
   cleanup_basisset();
   cleanup_molecule();

   stop_io();
   exit(0);
}
