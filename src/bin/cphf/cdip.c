/**************************************************************************/
/*                                                                        */
/*   CDIP:                                                                */
/*      Written by Edward Seidl (a NON-hog)                               */
/*      May 1991                                                          */
/*      Parts liberally ripped off from HFS group DIPDER                  */
/*        code  written in FORTRAN (Boo!!!)                               */
/**************************************************************************/

/* $Log$
 * Revision 1.1  2000/02/04 22:50:46  evaleev
 * Initial revision
 *
/* Revision 1.2  1995/01/19 19:48:36  seidl
/* replace some nulls with spaces
/*
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
   int ierr,max_words;
   int *i30t;
   double *temp;

   itap43 = 43;

   ffile(&infile,"input.dat",2);
   ffile(&outfile,"output.dat",1);

   ip_set_uppercase(1);
   ip_initialize(infile,outfile);

   rfile(itap43);

   ip_cwk_clear();
   ip_cwk_add(":DEFAULT");
   ip_cwk_add(":CDIP");

   tstart(outfile);

   fprintf(outfile,"\n%8cCDIP: Calculates dipole derivatives\n",' ');

   print=1;
   errcod = ip_data("IPRINT","%d",&print,0);

   ioff[0] = 0;
   for (i = 1; i < 1024 ; i++) {
      ioff[i] = ioff[i-1] + i;
      }

   twocon=0;
   i30t = (int *) init_array(100);
   rfile(30);
   wreadw(30,i30t,sizeof(int)*200,sizeof(int)*100,&i);
   i=i30t[42];
   if(i < 0) twocon=1;
   free(i30t);
   rclose(30,3);

 /* get some stuff from file40 */
   init_master();

   if(print) fprintf(outfile,"%8cIPRINT        = %7d\n",' ',print);
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

   link();

   tstop(outfile);

   ip_done();
   exit(0);
   }
