
/* $Log$
 * Revision 1.1  2000/02/04 22:50:47  evaleev
 * Initial revision
 *
/* Revision 1.2  1997/08/25 21:53:40  crawdad
/* Making changes for extension of PSI file size limit.
/*
 * Revision 1.1  1991/06/15  22:45:28  seidl
 * Initial revision
 * */

static char *rcsid = "$Id$";

#define EXTERN
#include "includes.h"
#include "common.h"

derivs()
   {
      int i,j,ij,how_many;
      double *temp,**temp2;
      double **deriv1,**deriv2,**onet;

      temp = (double *) init_array(nbatri*natom3);
      deriv1 = (double **) init_matrix(3,natom);
      deriv2 = (double **) init_matrix(natom3,natom3);
      onet = (double **) init_matrix(nbfao,nbfao);
      temp2 = (double **) init_matrix(nbfao,nbfao);

      work=91;
      rfile(work);

      how_many = (twocon) ? 6 : ntypes;

/* read in scf first derivs */

      sread(itap42,(char *) temp,sizeof(double)*natom3);
      for(i=ij=0; i < natom ; i++)
         for(j=0; j < 3 ; j++,ij++)
            deriv1[j][i] = temp[ij];

      if(print & 1) {
         fprintf(outfile,"\n scf first derivatives\n");
         print_mat(deriv1,3,natom,outfile);
         }

      free_matrix(deriv1,3);

/* read in scf second derivs */

      sread(itap42,(char *) temp,sizeof(double)*natom3*natom3);
         for(i=ij=0; i < natom3 ; i++)
            for(j=0; j < natom3 ; j++,ij++)
               deriv2[j][i] = temp[ij];

      if(print & 1) {
         fprintf(outfile,"\n scf second derivatives\n");
         print_mat(deriv2,natom3,natom3,outfile);
         }

      free_matrix(deriv2,natom3);

/* read in overlap derivative integrals and transform to mo basis */
      
      for (i=0; i < natom3 ; i++) {
         sread(itap42,(char *) temp,sizeof(double)*nbatri);
         if(print & 4) {
            fprintf(outfile,"\n sa (ao basis) i = %5d\n",i);
            print_array(temp,nbfao,outfile);
            }
         ao_to_mo(temp,onet,e_vecs_ao,temp2,nbfao,nbfso);
         if(print & 4) {
            fprintf(outfile,"\n sa (mo basis) i = %5d\n",i);
            print_array(temp,nbfso,outfile);
            }
         rwrit(itap44,(char *) temp,sizeof(double)*nbstri,sa_loc[i]);
         }

/* read in one electron derivative integrals and transform to mo basis */

      for (i=0; i < natom3 ; i++) {
         sread(itap42,(char *) temp,sizeof(double)*nbatri);
         if(print & 4) {
            fprintf(outfile,"\n ha (ao basis) i = %5d\n",i);
            print_array(temp,nbfao,outfile);
            }
         ao_to_mo(temp,onet,e_vecs_ao,temp2,nbfao,nbfso);
         if(print & 4) {
            fprintf(outfile,"\n ha (mo basis) i = %5d\n",i);
            print_array(temp,nbfso,outfile);
            }
         rwrit(itap44,(char *) temp,sizeof(double)*nbstri,ha_loc[i]);
         }

/* read in two electron derivative integrals and transform to mo basis */

      for (j=0; j < how_many ; j++) {
         for (i=0; i < natom3 ; i++) {
            sread(itap42,(char *) temp,sizeof(double)*nbatri);
            if(print & 4) {
               fprintf(outfile,"\n ta (ao basis) j = %5d i = %5d\n",j,i);
               print_array(temp,nbfao,outfile);
               }
            ao_to_mo(temp,onet,e_vecs_ao,temp2,nbfao,nbfso);
            if(print & 4) {
               fprintf(outfile,"\n ta (mo basis) j = %5d i = %5d\n",j,i);
               print_array(temp,nbfso,outfile);
               }
            swrit(work,(char *) temp,sizeof(double)*nbstri);
            }
         }

      free(temp);
      free_matrix(onet,nbfao);
      free_matrix(temp2,nbfao);
      }
