/**************************************************************************/
/*                                                                        */
/*   CPHF:                                                                */
/*      Written by Edward Seidl (a NON-hog)                               */
/*      October 1990                                                      */
/*      Parts liberally ripped off from HFS group CPXXAOS and DIPDER      */
/*        codes written in FORTRAN (Boo!!!)                               */
/*      Uses PK  supermatrix on file37                                    */
/*                                                                        */
/*      Modified                                                          */
/*      ETS                                                               */
/*      March 21, 1991 to do pure angular momentum functions              */
/*                                                                        */
/*      Modified                                                          */
/*      ETS                                                               */
/*      April 17, 1991 to use new input format developed by               */
/*        Curtis Janssen                                                  */
/*                                                                        */
/*      References: 1. Osamura,Yamaguchi,Saxe,Fox,Vincent,HFS             */
/*                     J.Mol.Struct. 103(1983) p.183 (THEOCHEM)           */
/*                  2. Yamaguchi,Osamura,HFS                              */
/*                     JACS 105(1983) p.7506                              */
/*                                                                        */
/**************************************************************************/
/*                                                                        */
/*   Description of input                                                 */
/*                                                                        */
/*   IR_INT = boolean                                                     */
/*      If true, dipole derivatives are calculated.  The default is true. */
/*                                                                        */
/*   CONVERGENCE = integer                                                */
/*      Convergence of the cphf equations. The default is 12              */
/*                                                                        */
/*   IPRINT = integer                                                     */
/*      A printing flag.                                                  */
/*       IPRINT =    0 minimum output possible                            */
/*          "   =    1 print 2nd derivs, dipole derivs                    */
/*          "   =    2 print eigenvectors and dependent pairs             */
/*          "   =    4 print results of derivs()                          */
/*          "   =    8 print results of famat()                           */
/*          "   =   16 print results of bamat()                           */
/*          "   =   32 print results of bfmat() and dipole_derivs()       */
/*          "   =   64 print results of cphf()                            */
/*          "   =  128 print results of uxmat()                           */
/*          "   =  256 print results of wamat()                           */
/*          "   =  512 print results of hamat()                           */
/*          "   = 1024 print results of bafmat()                          */
/*          "   = 2048 print results of hemat()                           */
/*                                                                        */
/*       IPRINT values can be added together, eg. if you want output      */
/*       from derivs and uxmat, set IPRINT = 132                          */
/*                                                                        */
/*   MAXWORDS = integer                                                   */
/*      The maximum number of real words of core memory which can be      */
/*      allocated.                                                        */
/*                                                                        */
/**************************************************************************/

/* $Log$
 * Revision 1.1  2000/02/04 22:50:47  evaleev
 * Initial revision
 *
/* Revision 1.6  1997/09/12 13:54:43  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 1.5  1997/08/25  21:53:37  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 1.4  1996/06/24  23:11:37  sherrill
 * Make cphf read memory keyword instead of maxwords keyword.
 *
 * Revision 1.3  1995/01/19  19:48:38  seidl
 * replace some nulls with spaces
 *
 * Revision 1.2  1991/07/30  04:51:43  seidl
 * fix for f and g functions
 * move to psi area
 *
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";


#include "includes.h"
#include "common.h"
#include <ip_libv1.h>


void main()
{
   int i,j,ij;
   int errcod;
   int ierr,max_words,max_bytes;
   int *i30t;
   PSI_FPTR junk;
   double *temp;

   itap37 = 37;
   itap42 = 42;
   itap43 = 43;
   itap44 = 44;

   ffile(&infile,"input.dat",2);
   ffile(&outfile,"output.dat",1);

   ip_set_uppercase(1);
   ip_initialize(infile,outfile);

   rfile(itap37);
   rfile(itap42);
   rfile(itap44);

   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":CPHF");

   tstart(outfile);

   fprintf(outfile,"\n%8cCPHF: A General cphf code in ao basis\n",' ');
   fprintf(outfile,"%14cOooh.  And now it does 5d's too.\n\n",' ');

   fprintf(outfile,"%8cRCSID: %s\n\n",' ',rcsid);

   /* previously, we read maxwords from input.  Now use fndcor to
      read the memory keyword, if it's available.                 */

   if (ip_exist("MEMORY", 0)) {
      fndcor(&max_bytes, infile, outfile);
      max_words = max_bytes / sizeof(double);
      }
   else {
      max_words = MAX_WORDS;
      }

   print=0;
   errcod = ip_data("IPRINT","%d",&print,0);

   conv=12;
   errcod = ip_data("CONVERGENCE","%d",&conv,0);

   dipole=1;
   errcod = ip_boolean("IR_INT",&dipole,0);

   ioff[0] = 0;
   for (i = 1; i < 1024 ; i++) {
      ioff[i] = ioff[i-1] + i;
      }

   twocon=0;
   i30t = (int *) init_array(100);
   rfile(30);
   wreadw(30,(char *) i30t,sizeof(int)*200,(PSI_FPTR) sizeof(int)*100,&junk);
   i=i30t[42];
   if(i < 0) twocon=1;
   free(i30t);
   rclose(30,3);

 /* get some stuff from file40 */
   init_master();

   fprintf(outfile,"%8cCONVERGENCE   = %7d\n",' ',conv);
   fprintf(outfile,"%8cIR_INT        = %7d\n",' ',dipole);
   if(print) fprintf(outfile,"%8cIPRINT        = %7d\n",' ',print);
   fprintf(outfile,"%8cMAXWORDS      = %7d\n",' ',max_words);
   fprintf(outfile,"%8cnbfso         = %7d\n",' ',nbfso);
   fprintf(outfile,"%8cnbfao         = %7d\n",' ',nbfao);
   fprintf(outfile,"%8cnbatri        = %7d\n",' ',nbatri);
   fprintf(outfile,"%8cnbstri        = %7d\n",' ',nbstri);
   fprintf(outfile,"%8cnatom         = %7d\n",' ',natom);
   fprintf(outfile,"%8ciopen         = %7d\n",' ',iopen);
   fprintf(outfile,"%8cnind          = %7d\n",' ',nind);
   fprintf(outfile,"%8cndep          = %7d\n",' ',ndep);

   fprintf(outfile,"\n%8cnuclear repulsion energy = %20.10f\n",' ',enuc);
   fprintf(outfile,"%8cscf energy               = %20.10f\n\n",' ',escf);

   if(twocon) fprintf(outfile,"%8ccalculation for tcscf wavefunction\n",' ');
   else if(iopen) fprintf(outfile,"%8ccalculation for open shell system\n",' ');
   else fprintf(outfile,"%8ccalculation for closed shell system\n",' ');

/* read in ao and so eigenvectors and eigenvalues and occ numbers */

   temp = (double *) init_array(nbfao*nbfao);
   e_vals = (double *) init_array(nbfso);
   occ_num = (double *) init_array(nbfso);
   e_vecs_ao = (double **) init_matrix(nbfao,nbfao);
   e_vecs_so = (double **) init_matrix(nbfso,nbfso);

   mread(e_vals,16);
   mread(occ_num,17);

   mread(temp,18);
   for(i=ij=0; i < nbfso ; i++) 
      for(j=0; j < nbfso ; j++,ij++)
         e_vecs_so[j][i] = temp[ij];

   mread(temp,19);
   for(i=ij=0; i < nbfao ; i++) 
      for(j=0; j < nbfao ; j++,ij++)
         e_vecs_ao[j][i] = temp[ij];

   free(temp);

   if(print & 2) {
      fprintf(outfile,"\n ao eigenvector, eigenvalues, and occupations\n");
      eigout(e_vecs_ao,e_vals,occ_num,nbfao,nbfao,outfile);

      fprintf(outfile,"\n so eigenvector, eigenvalues, and occupations\n");
      eigout(e_vecs_so,e_vals,occ_num,nbfso,nbfso,outfile);

      if(iopen) {
         fprintf(outfile,"\n alpa matrix\n");
         print_mat(alpa,ntype1,ntype1,outfile);

         fprintf(outfile,"\n beta matrix\n");
         print_mat(beta,ntype1,ntype1,outfile);
         }
      }

   make_vec();

/* read in derivatives from file42 and transform to mo basis */

   derivs();

/* if twocon, form ha matrices (eqns 33-35 ref 2) */

   if(twocon) hamat();

/* form derivative lagrangian matrix  epsilon_a (eqn 7) */
/* (eqn 32 ref 2)      */

   if(iopen) famat_o();
   else famat_c();

/* if twocon, form h matrix (eqns 5-7 ref 2) */

   if(twocon) hemat();

/* calculate ba matrix for independent pairs (eqn 13) */

   if(iopen) bamat_o(max_words);
   else bamat(max_words);
   fflush(outfile);

/* bf matrix for dipole stuff */

   if(dipole) bfmat();

/* if twocon, form ci part of a and ba matrix (eqns 25 & 27 ref 2) */

   if(twocon) bafmat();

/* do the cphf iterations */

   cphf_iter(max_words);

/* finish up u matrices */

   uxmat();

/* form wa matrix (eqn 6) */

   if(iopen) wamat_o(max_words);
   else wamat_c(max_words);

/* calculate scf second derivatives (eqn 5) */

   if(!iopen) scf2nd();
   else scf2nd_o();

/* dipole stuff */

   if(dipole) {
      dipmo();
      polar();
      }

   rclose(itap37,3);
   rclose(itap40,3);
   rclose(itap42,3);
   rclose(itap44,3);
   rclose(work,4);
   tstop(outfile);

   ip_done();
   exit(0);
   }
