/*
/*************************************************************************/
/*                                                                        */
/*   CSCF:                                                                */
/*      Written by Edward Seidl (a NON-hog)                               */
/*      September 1990                                                    */
/*      Parts liberally ripped off from HFS group SCF                     */
/*        code written in FORTRAN (Boo!!!)                                */
/*                                                                        */
/*      modified April 17, 1991 to use new input format developed by      */
/*        Curtis Janssen                                                  */
/**************************************************************************/
/*                                                                        */
/*   Description of input                                                 */
/*                                                                        */
/*   LABEL = string                                                       */
/*         This is a character string to be included in the output.       */
/*         There is no default.                                           */
/*                                                                        */
/*   WFN = string                                                         */
/*         This is the type of wavefunction which is ultimately desired.  */
/*         The default is SCF.                                            */
/*                                                                        */
/*   DIRECT_SCF = boolean                                                 */
/*         Flag to request the direct formation of the Fock matrix        */
/*                                                                        */
/*   OPENTYPE = string                                                    */
/*         This specifies the state desired.  It can be one of NONE       */
/*         (for a closed shell singlet), SINGLET (for an open shell       */
/*         singlet), HIGHSPIN (for any high spin open shell system),      */
/*         TWOCON (for a two configuration singlet), or SPECIAL.          */
/*         If SPECIAL is given, then alpha and beta coupling              */
/*         coefficients must be given with the ALPHA and BETA keywords.   */
/*         The default is NONE.                                           */
/*                                                                        */
/*   DOCC = integer_vector                                                */
/*         This gives the number of doubly occupied orbitals in each      */
/*         irreducible representation.  There is no default.              */
/*                                                                        */
/*   SOCC = integer_vector                                                */
/*         This gives the number of singly occupied orbitals in each      */
/*         irreducible representation.  If OPENTYPE = NONE this defaults  */
/*         to the zero vector.  Otherwise, there is no default.           */
/*                                                                        */
/*   DERTYPE = string                                                     */
/*         This specifies the order of derivative that is to even-        */
/*         tually  be  done.   It  is  used  by the scf program to        */
/*         determine if certain files are to be written and it  is        */
/*         also  used  to determine the default convergence of the        */
/*         wavefunction.  The default is FIRST.                           */
/*                                                                        */
/*   MAXITER = integer                                                    */
/*         This gives the maximum number of iterations.  The default      */
/*         is 40.                                                         */
/*                                                                        */
/*   CONVERGENCE = integer                                                */
/*         The convergence criterion is 10**(-integer).  The default is   */
/*         7 if both DERTYPE = NONE and WFN = SCF are given and 10        */
/*         otherwise.                                                     */
/*                                                                        */
/*   LEVELSHIFT = real                                                    */
/*        This specifies the level shift. The default is 1.0.             */
/*                                                                        */
/*                                                                        */
/*   There are also a large number of less  commonly  used  input         */
/*   parameters.   If  you  do  not understand what the following         */
/*   options mean, then make sure that they do not appear in your         */
/*   input.   The defaults will work in the overwhelming majority         */
/*   of cases.  These are specified with the following keywords:          */
/*                                                                        */
/*                                                                        */
/*   REORDER = boolean                                                    */
/*        The molecular orbitals will be  reordered  if  this  is         */
/*        true,  in  which  case,  the  MOORDER parameter must be         */
/*        present.  The default is false.                                 */
/*                                                                        */
/*   MOORDER = integer_vector                                             */
/*        This specifies a molecular orbital  reordering  vector.         */
/*        It  will  only  be  used if REORDER = YES.  This vector         */
/*        contains first the ordering for  the  orbitals  in  the         */
/*        first  irreducible  representation  and then the second         */
/*        and so on.   The  first  orbital  of  each  irreducible         */
/*        representation is numbered 1.  There is no default.             */
/*                                                                        */
/*   ALPHA = real_vector                                                  */
/*        If OPENTYPE = SPECIAL, then this  parameter  gives  the         */
/*        alpha coupling coefficients.  The number of elements in         */
/*        this vector is MM(MM+1)/2, where MM is  the  number  of         */
/*        irreducible  representations containing singly occupied         */
/*        molecular orbitals.  There is no default.                       */
/*                                                                        */
/*   BETA = real_vector                                                   */
/*        If OPENTYPE = SPECIAL, then this  parameter  gives  the         */
/*        beta  coupling coefficients.  The number of elements in         */
/*        this vector is MM(MM+1)/2, where MM is  the  number  of         */
/*        irreducible  representations containing singly occupied         */
/*        molecular orbitals.  There is no default.                       */
/*                                                                        */
/*   RESTART = boolean                                                    */
/*        The calculation will restart from the old  wavefunction         */
/*        if RESTART is true.  If the old wavefunction does not           */
/*        exist, then the cscf program will generate its own  ini-        */
/*        tial  guess  automatically.   Possible  values for this         */
/*        parameter are TRUE, YES, 1,  FALSE,  NO,  and  0.   The         */
/*        default is true.                                                */
/*                                                                        */
/*   IPRINT = integer                                                     */
/*        This is a print option.  The default is 0.                      */
/*                                                                        */
/*   ROTATE = boolean                                                     */
/*        The molecular orbitals will not be rotated if  this  is         */
/*        false.   The rotation only affects the virtual orbitals         */
/*        for open shell systems.  This parameter  must  be  true         */
/*        for  correlated  gradients  and  it  must  be false for         */
/*        second and higher derivatives.  The default is false if         */
/*        WFN = SCF and true otherwise.                                   */
/*                                                                        */
/*   DIIS = boolean                                                       */
/*        This determines whether diis will be used.  The default is      */
/*        false for OPENTYPE = TWOCON and true otherwise.                 */
/*                                                                        */
/*   NDIIS = integer                                                      */
/*        This gives the number of error matrices to use in the diis      */
/*        procedure.  The default is 6 for closed shell, 4 for open       */
/*        shell, and 3 for tcscf.                                         */
/*                                                                        */
/*   DIISSTART = integer                                                  */
/*        This gives the first iteration for which DIIS  will  be         */
/*        used.  The default is 0.                                        */
/*                                                                        */
/*   DIISDAMP = real                                                      */
/*        This gives the damping factor for the diis procedure.  The      */
/*        default is 0.0 for closed shell, 0.02 for open shell, and       */
/*        0.01 for tcscf.                                                 */
/*                                                                        */
/*   INCR = real                                                          */
/*        This is used in tcscf to determine how often the ci             */
/*        coefficients are recalculated.  A small number (~0.25)          */
/*        will cause them to be recalculated nearly every scf             */
/*        iteration.  The default is 0.5.                                 */
/*                                                                        */
/*                                                                        */
/*   FOCK_TYPE = integer                                                  */
/*        Only used for tcscf and open shell calculations.                */
/*        If FOCK_TYPE = 0 use a simple form for fock_eff                 */
/*         "     "     = 1 use form of fock_eff suitable for              */
/*                         high-spin cases                                */
/*         "     "     > 1 experimental fock matrices                     */
/*                                                                        */
/**************************************************************************/



/*-------------------------------------------------------------------------
  READ THIS FIRST: This code is a hack of the original CSCF which employed
  symmetric orthogonalization procedure. This code uses the canonical
  orthogonalization procedure. To decrease the amount of
  rewriting I had to do I left the initialization part (init_scf and
  init_scf2) intact. Thus, although I have to allocate more space than
  I need, I should be fine otherwise.

                                             - Edward Valeev, August'99
 -------------------------------------------------------------------------*/


static char *rcsid = "$Id$";

#include "includes.h"
#include "common.h"
#include <ip_libv1.h>
#include <psio.h>

int main(argc,argv)
   int argc;
   char *argv[];
{
   int i,nn;
   char *prog_name="CSCF3.0: An SCF program written in C";
   char *output="APPEND  ";
   struct symm *s;
   ip_value_t *ipvalue=NULL;

   ffile(&infile,"input.dat",2);
   ffile(&outfile,"output.dat",1);
   /*ffile(&JK,"jandk.dat",0);
   ffile(&gmat,"gmat.dat",0);
   ffile(&diis_out,"diis_out.dat",0);*/

   ip_set_uppercase(1);
   ip_initialize(infile,outfile);

#if 0
   ip_cwk_clear();
#endif
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":SCF");

   ip_string("OUTPUT",&output,0);
   if(!strcmp(output,"TERMINAL")) {
     outfile = stdout;
     }
   else if(!strcmp(output,"WRITE")) {
     fclose(outfile);
     ffile(&outfile,"output.dat",0);
     }

   if(argc > 1) {
     ip_print_tree(stdout,NULL);
     if(ipvalue != NULL) ip_print_value(stdout,ipvalue);
     }

   tstart(outfile);
   
   fprintf(outfile,"\n%13c------------------------------------------\n",' ');
   fprintf(outfile,"\n%16c%s\n",' ',prog_name);
   fprintf(outfile,"\n%14cWritten by too many people to mention here\n",' ');
   fprintf(outfile,"\n%13c------------------------------------------\n\n\n",' ');
   
   
   itap30 = 30;
   itap33 = PSIF_AO_TEI;
   itap34 = PSIF_SO_INTS;
   itapS  = PSIF_AO_S;
   itapT  = PSIF_AO_T;
   itapV  = PSIF_AO_V;
   itapDSCF = PSIF_DSCF;
   itap92 = 92;
   itap93 = 93;
   file30_init();

   /* EFV 10/24/98 Get the integral format: IWL = true */
   use_iwl = 1;
   ip_boolean("USE_IWL",&use_iwl,0);

/* JPK 6/1/00 integral accuracy: dynamic(default)=1, static=0 */
   dyn_acc = 1;
   eri_cutoff = 1.0E-14;
   ip_boolean("DYN_ACC",&dyn_acc,0);
   tight_ints=0;
   delta = 1.0;
                                                      
/* open integrals file(s) */

   if (use_iwl)
     psio_init();
   else
     rfile(itap34);
   
   /* STB (6/30/99) - Function added because in order to initialize things
      one must know whether you are doing UHF or restricted */   
   
   occ_init();
   
/* initialize some constants and arrays */
   if(uhf)
       init_uhf();
   else
       init_scf();

/* read input.dat, get occupations, flags, etc. */

   scf_input(ipvalue);

/* set up other useful arrays */

   init_scf2();

/* get one electron integrals */

   rdone_iwl();
   
   if(print & 1) {
      for (i=0; i < num_ir; i++) {
         s = &scf_info[i];
         if (nn=s->num_so) {
            fprintf(outfile,"\nsmat for irrep %s\n",s->irrep_label);
            print_array(s->smat,nn,outfile);
            fprintf(outfile,"\ntmat for irrep %s\n",s->irrep_label);
            print_array(s->tmat,nn,outfile);
            fprintf(outfile,"\nhmat for irrep %s\n",s->irrep_label);
            print_array(s->hmat,nn,outfile);
            }
         }
      }

/* form S-1/2 matrix sahalf */

   shalf();

   if (print & 1) {
      for (i=0; i < num_ir ; i++) {
         s = &scf_info[i];
         if (nn=s->num_so) {
            fprintf(outfile,"\nsahalf for irrep %s\n",s->irrep_label);
            print_mat(s->sahalf,nn,s->num_mo,outfile);
            }
         }
      }

/* if no initial vector, form one from core hamiltonian */

   if (inflg < 0) form_vec();
   fflush(outfile);

/* guess or designate orbital occupations*/
   guess();

/* form first density matrix */
   if(uhf && inflg)
       schmit_uhf(0);
   else
       schmit(0);
   
   if (!twocon){
       if(!uhf)
	   dmat();
       else{
	   if(!inflg)
	       cmatsplit();
	   dmatuhf();
	   
       }
   }
   
/* if there are no two-electron integrals, exit */

/* Decide how to form the Fock matrix */

   if (!direct_scf) {

/* read in tei's and form supermatrix */

     num_ints = 0;
     num_bufs = 0;
     rfile(itap92);
     rfile(itap93);
     rdtwo();
     
     if (!use_iwl)
       if(delete_ints) rclose(itap34,4);
       else rclose(itap34,3);
   }
   else {
/* form the Fock matrix directly */

   /*check for rohf singlet...doesn't work direct*/
 
      if(!uhf && singlet) {
  
         fprintf(outfile,"\n  rohf open shell singlet doesn't work direct\n");
         fprintf(outfile,"  remove 'direct_scf = true' from input\n");
         fprintf(stdout,"rohf open shell singlet doesn't work direct\n");
         fprintf(stdout,"remove 'direct_scf = true' from input\n");
 
         psio_done();
         file30_close();
         exit(1);
      }
                                               

     formg_direct();
     if(dyn_acc)  fprintf(outfile,"\n  Using inexpensive integrals");   
   }

/* iterate */

   iter = 0;
   converged = 0;
      
   if(twocon) scf_iter_2();
   else if(uhf) uhf_iter();
   else scf_iter();

   psio_done();

   cleanup();
   }
