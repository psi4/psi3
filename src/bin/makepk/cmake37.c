
/* $Log$
 * Revision 1.3  2002/12/17 20:29:28  mabrams
 * Arbitrary input and output filenames
 *
/* Revision 1.2  2002/03/25 02:17:36  janssen
/* Get rid of tmpl.  Use new naming scheme for libipv1 includes.
/*
/* Revision 1.1.1.1  2000/02/04 22:51:32  evaleev
/* Started PSI 3 repository
/*
/* Revision 1.3  1997/09/12 13:54:49  crawdad
/* Changing marco name from ULL to PSI_FPTR.
/*
 * Revision 1.2  1997/08/25  21:52:24  crawdad
 * Making changes for extension of PSI file size limit.
 *
 * Revision 1.1  1991/06/15  22:06:29  seidl
 * Initial revision
 * */

/**************************************************************************/
/*                                                                        */
/*   MAKEPK:                                                              */
/*      Written by Edward Seidl (a NON-hog)                               */
/*      September 1990                                                    */
/*      Parts liberally ripped off from HFS group MAKE37 and MASTER       */
/*        codes written in FORTRAN (Boo!!!)                               */
/*      Creates PKI supermatrix on file37                                 */
/*      Creates master file file40 as well                                */
/*      Uses a more efficient algorithm (I hope)                          */
/*                                                                        */
/*      modified march 21, 1991 by ets to do 5ds                          */
/*      modified april 17, 1991 by ets to use new input format developed  */
/*       by Curtis Janssen                                                */
/*                                                                        */
/**************************************************************************/
/*                                                                        */
/*   Description of input                                                 */
/*                                                                        */
/*   WFN = string                                                         */
/*      This is the type of wavefunction which is ultimately desired.     */
/*      The default is SCF.                                               */
/*                                                                        */
/*   DELETE34 = boolean                                                   */
/*      If true, file34 is removed.  The default is true if WFN = SCF     */
/*      and false otherwise.                                              */
/*                                                                        */
/*   CUTOFF = integer                                                     */
/*      Threshold for eliminating integrals is 10**(-CUTOFF)              */
/*                                                                        */
/*   IPRINT = integer                                                     */
/*      A printing flag.                                                  */
/*                                                                        */
/**************************************************************************/

static char *rcsid = "$Id$";

#include "includes.h"
#include "common.h"
#include <libipv1/ip_lib.h>
#include <libqt/qt.h>


void main(int argc, char *argv[])
   {
      int i,isadr;
      int errcod;
      int ns,lapu,iclos,nkind;
      PSI_FPTR next;
      int nn,ierr;
      int delete34;
      int *i30;
      char *wfn="SCF";
      char *bool;
      double *dum_buf,*smat,*tmat,*vmat;

      itap30 = 30;
      itap34 = 34;
      itap37 = 37;
      itap40 = 40;
      
      /* Arbitrary input and output filenames */
      init_in_out(argc-1,argv+1);
      
      ip_set_uppercase(1);
      ip_initialize(infile,outfile);

      rfile(itap30);
      rfile(itap34);
      rfile(itap37);

      ip_cwk_clear();
      ip_cwk_add(":DEFAULT");
      ip_cwk_add(":MAKEPK");

      tstart(outfile);

      fprintf(outfile,"\n%8cMAKEPK: A Program to write a supermatrix file\n\n",
                                                               ' ');

      toler = 12;
      errcod = ip_data("CUTOFF","%d",&toler,0);

      print=0;
      errcod = ip_data("IPRINT","%d",&print,0);

      ci_calc=0;
      errcod = ip_string("WFN",&wfn,0);
      if(strcmp(wfn,"SCF")) ci_calc=1;

      delete34 = (ci_calc) ? 0 : 1;
      errcod = ip_string("DELETE34",&bool,0);
      if(errcod == IPE_OK) {
         if(!strcmp(bool,"NO") || !strcmp(bool,"FALSE") || !strcmp(bool,"0"))
           delete34=0;
         else delete34=1;
         }

      fprintf(outfile,"%8ccutoff for keeping integrals = %e\n",
                                              ' ',pow(10.0,(double) -toler));
      fprintf(outfile,"%8cprint  = %d\n",' ',print);
      if(ci_calc) fprintf(outfile,"%8cci calculation so no file37 produced\n",
                                                                   ' ');
      if(delete34) fprintf(outfile,"%8cfile34 will be deleted\n",' ');

      i30 = (int *) init_array(100);
      dum_buf = (double *) init_array(MAX_BASIS);

      ioff[0] = 0;
      for (i = 1; i < 1024 ; i++) {
         ioff[i] = ioff[i-1] + i;
         }

      wreadw(itap30,(char *) i30,sizeof(int)*200,(PSI_FPTR) sizeof(int)*100,&next);

      nbfso = i30[17];
      natom = i30[18];
      nbfao = i30[21];
      mxcoef = i30[41];
      iopen = i30[42];

      free(i30);

      if(iopen < 0) {
         iopen = -iopen;
         twocon = 1;
         }

/* read header information from integral tape */

      wreadw(itap34,(char *) &nkind,(PSI_FPTR) sizeof(int)*1,0,&next);
      wreadw(itap34,(char *) &iclos,sizeof(int)*1,next,&next);
      wreadw(itap34,(char *) blabel,sizeof(char)*80,next,&next);
      wreadw(itap34,(char *) &repnuc,sizeof(double)*1,next,&next);
      wreadw(itap34,(char *) &num_ir,sizeof(int)*1,next,&next);
      wreadw(itap34,(char *) degen,sizeof(int)*num_ir,next,&next);
      wreadw(itap34,(char *) dum_buf,sizeof(int)*num_ir,next,&next);
      wreadw(itap34,(char *) num_so,sizeof(int)*num_ir,next,&next);
      wreadw(itap34,(char *) &ns,sizeof(int)*1,next,&next);
      wreadw(itap34,(char *) dum_buf,sizeof(int)*2*ns,next,&next);
      wreadw(itap34,(char *) &lapu,sizeof(int)*1,next,&next);
      wreadw(itap34,(char *) dum_buf,sizeof(int)*lapu,next,&next);
      wreadw(itap34,(char *) dum_buf,sizeof(int)*lapu,next,&next);

      free(dum_buf);
/* set integral file pointer to sector boundary */

      isadr = i2sec(next) + 1;
      rsetsa(itap34,isadr);
      pos34 = sec2i(isadr-1);

      if (nkind != 1 && nkind != 2) {
         fprintf(outfile,"integral file screwed up, fix ints somebody!!!\n");
         exit(1);
         }

      n_so_typs=0;

      for (i=0; i < num_ir ; i++)
         if (nn=num_so[i]) {
            block_num[i] = n_so_typs;
            n_so_typs++;
            }

      nbatri = nbfao*(nbfao+1)/2;
      nbstri = nbfso*(nbfso+1)/2;
      nbasis = nbfso;
      ntri = nbstri;

      ideg[0]=0;
      if (n_so_typs != 1) {
         int ii = 0;
         for (i=1; i < num_ir; i++) {
            if (num_so[i] <= 0) {
               ideg[i]=ideg[i-1];
               }
            else {
               do {
                  nn=num_so[ii];
                  ii++;
                  } while(!nn);
                ideg[i]=ideg[i-1]+nn;
                }
             }
         }

      maxbuf=8192;

/* initialize (ise?) master file file40 */

      init_master();

/* get one electron integrals from file34 */
      
      smat = (double *) init_array(nbstri);
      tmat = (double *) init_array(nbstri);
      vmat = (double *) init_array(nbstri);

      rdone(smat);
      mwrit(smat,13);
      mwrit(smat,20);

      rdone(tmat);
      rdone(vmat);

      for (i=0; i < nbstri ; i++) smat[i] = tmat[i]+vmat[i];

      mwrit(smat,14);
      mwrit(smat,21);

      free(smat);
      free(tmat);
      free(vmat);

/* now get information about the scf calculation from file30 */

      scf_stuff();

/* and sort the eigenvector matrices */

      vec_sort();

/* form lagrangian matrix and alpb and betb for open shell */

      if(iopen || ci_calc) {
         abmat();
         density();
         }
      else clscf();

/* for now read in tei's */

      rdtwo();

      if(iopen || ci_calc) zetmat();

      rclose(itap30,3);

      if(delete34) rclose(itap34,4);
      else rclose(itap34,3);

      if(ci_calc) rclose(itap37,4);
      else rclose(itap37,3);

      rclose(itap40,3);
      tstop(outfile);

      ip_done();
      exit(0);
      }
