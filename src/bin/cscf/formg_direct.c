static char *rcsid = "$Id$";

#define EXTERN
#include <psio.h>
#include "includes.h"
#include "common.h"

void formg_direct()
{
   static double *gfull, *gfull_o;
   int stat;
   int i,j,k,jj,kk,l,off,joff,nn,max;      
   int ntri;

   ntri=ioff[nbasis];

   if(gfull == NULL) {
      gfull = (double *) init_array(ntri);
      if (iopen || uhf)
	gfull_o = (double *) init_array(ntri);
   }

      if(acc_switch) {
         if (uhf) {
            for(k=0; k < num_ir ; k++) {
               if(nn=scf_info[k].num_so) {
                  for(i=0; i < nn ; i++) {
                     for(j=0; j <= i ; j++) {
                        spin_info[0].scf_spin[k].gmat[ioff[i]+j] = 0.0;
                        spin_info[1].scf_spin[k].gmat[ioff[i]+j] = 0.0;
                     }
                  }
               }
            }
         }
         else if (iopen) {
            for(k=joff=0; k < num_ir ; k++) {
               if(nn=scf_info[k].num_so) {
                  for(i=0; i < nn ; i++) {
                     for(j=0; j <= i ; j++) {
                        scf_info[k].gmat[ioff[i]+j] = 0.0;
                        scf_info[k].gmato[ioff[i]+j] = 0.0;
                     }
                  }
               }
            }
         }
         else {
            for(i=0;i<num_ir;i++) {
               max = scf_info[i].num_so;
               for(j=0;j<max;j++) {
                  for(k=0;k<=j;k++) {
                     scf_info[i].gmat[ioff[j]+k] = 0.0;
                  }
               }
            }
         }
      }
                                                             
   stat = system("cints --fock");
   switch (stat) {
   case 0:
       /* everything is OK */
       break;

   default:
       /* Something went wrong */
       fprintf(outfile,"  formg_direct: System call to CINTS failed. Check to see if it's in your PATH\n");
       fprintf(stderr,"System call to CINTS failed. Check to see if it's in your PATH.\n");
       exit(1);
   }
   
   psio_open(itapDSCF, PSIO_OPEN_OLD);
   if (!uhf) {
     psio_read_entry(itapDSCF, "Total G-matrix", (char *) gfull, sizeof(double)*ntri);
     if (iopen) {
       psio_read_entry(itapDSCF, "Open-shell G-matrix", (char *) gfull_o, sizeof(double)*ntri);
     }
   }
   else {
     psio_read_entry(itapDSCF, "Alpha G-matrix", (char *) gfull, sizeof(double)*ntri);
     psio_read_entry(itapDSCF, "Beta G-matrix", (char *) gfull_o, sizeof(double)*ntri);
   }
   psio_close(itapDSCF, 0);

   if (uhf) {
     for(k=joff=0; k < num_ir ; k++) {
       if(nn=scf_info[k].num_so) {
	 for(i=0; i < nn ; i++) {
	   for(j=0; j <= i ; j++) {
	     spin_info[0].scf_spin[k].gmat[ioff[i]+j] += gfull[ioff[i+joff]+j+joff];
	     spin_info[1].scf_spin[k].gmat[ioff[i]+j] += gfull_o[ioff[i+joff]+j+joff];
	   }
	 }
       }
       joff += nn;
     }
   }
   else if (iopen) {
     for(k=joff=0; k < num_ir ; k++) {
       if(nn=scf_info[k].num_so) {
	 for(i=0; i < nn ; i++) {
	   for(j=0; j <= i ; j++) {
	     scf_info[k].gmat[ioff[i]+j] += gfull[ioff[i+joff]+j+joff];
	     scf_info[k].gmato[ioff[i]+j] += gfull_o[ioff[i+joff]+j+joff];
	   }
	 }
       }
       joff += nn;
     }
   }
   else {
     for(k=joff=0; k < num_ir ; k++) {
       if(nn=scf_info[k].num_so) {
	 for(i=0; i < nn ; i++) {
	   for(j=0; j <= i ; j++) {
	     scf_info[k].gmat[ioff[i]+j] += gfull[ioff[i+joff]+j+joff];
	   }
	 }
       }
       joff += nn;
     }
   }

   return;
}

