#include <math.h>

#if defined(AIX)
# define C_MMULT mmult
#else
# define C_MMULT mmult_
#endif

/***************************************************************/
/*                                                             */
/* a reasonably fast matrix multiply (at least on the DEC3100) */
/* written by ETS                                              */
/*                                                             */
/* AF,BF,and CF are fortran arrays                             */
/*                                                             */
/* ta,tb and tc indicate whether the corresponding arrays are  */
/*              to be converted to their transpose             */
/*                                                             */
/* nr,nl,nc are the number of rows,links,and columns in the    */
/*          final matrices to be multiplied together           */
/*          if ta=0 AF should have the dimensions nr x nl      */
/*          if ta=1 AF should have the dimensions nl x nr      */
/*          if tb=0 BF should have the dimensions nl x nc      */
/*          if tb=1 BF should have the dimensions nc x nl      */
/*          if tc=0 CF should have the dimensions nr x nc      */
/*          if tc=1 CF should have the dimensions nc x nr      */
/*                                                             */
/* add is 1 if this matrix is to be added to the one passed    */
/*        in as CF, 0 otherwise                                */
/***************************************************************/

static int keep_nr=0;
static int keep_nl=0;
static int keep_nc=0;
static double **a,**b;

C_MMULT(AF,da,ta,BF,db,tb,CF,dc,tc,nr,nl,nc,add)
   double *AF,*BF,*CF;
   int *ta,*tb,*tc,*nr,*nl,*nc,*add;
   int *da,*db,*dc;

{
   int odd_nr,odd_nc,odd_nl;
   int i,j,k,ij;
   double t00,t01,t10,t11;
   double *att,*bt;
   double *at1,*bt1;

   if(!a) {
      a = (double **) init_matrix(*nr,*nl);
      b = (double **) init_matrix(*nc,*nl);
      keep_nr = *nr;
      keep_nl = *nl;
      keep_nc = *nc;
      }

   if(*nl > keep_nl) {
      free_matrix(a,keep_nr);
      free_matrix(b,keep_nc);
      keep_nl = *nl;
      keep_nr = (*nr > keep_nr) ? *nr : keep_nr;
      keep_nc = (*nc > keep_nc) ? *nc : keep_nc;
      a = (double **) init_matrix(keep_nr,keep_nl);
      b = (double **) init_matrix(keep_nc,keep_nl);
      }
   if(*nr > keep_nr) {
      free_matrix(a,keep_nr);
      keep_nr = *nr;
      a = (double **) init_matrix(keep_nr,keep_nl);
      }
   if(*nc > keep_nc) {
      free_matrix(b,keep_nc);
      keep_nc = *nc;
      b = (double **) init_matrix(keep_nc,keep_nl);
      }

   odd_nr = (*nr)%2;
   odd_nc = (*nc)%2;
   odd_nl = (*nl)%2;

   if(*ta)
      for(i=0; i < *nr ; i++) {
         ij= *da*i;
         for(j=0; j < *nl ; j++,ij++)
            a[i][j] = AF[ij];
         }
   else
      for(j=0; j < *nl ; j++) {
         ij = *da*j;
         for(i=0; i < *nr ; i++,ij++)
            a[i][j] = AF[ij];
         }

   if(*tb)
      for(j=0; j < *nl ; j++) {
         ij = *db*j;
         for(i=0; i < *nc ; i++,ij++)
            b[i][j] = BF[ij];
         }
   else
      for(i=0; i < *nc ; i++) {
         ij = *db*i;
         for(j=0; j < *nl ; j++,ij++)
            b[i][j] = BF[ij];
         }
      
   if(*tc) {
      for(i=ij=0; i < (*nr)-1 ; i+=2) {
         ij = i*(*dc);
         for(j=0; j < (*nc)-1 ; j+=2,ij+=2) {
            att=a[i]; bt=b[j];
            at1=a[i+1]; bt1=b[j+1];
            if(*add) {
               t00 = CF[ij];
               t01 = CF[ij+(*dc)];
               t10 = CF[ij+1];
               t11 = CF[ij+1+(*dc)];
               }
            else t00=t01=t10=t11=0.0;
            for(k= *nl; k ; k--,att++,bt++,at1++,bt1++) {
               t00 += *att * *bt;
               t10 += *att * *bt1;
               t01 += *at1 * *bt;
               t11 += *at1 * *bt1;
               }
            CF[ij]=t00;
            CF[ij+(*dc)]=t01;
            CF[ij+1]=t10;
            CF[ij+1+(*dc)]=t11;
            }
         if(odd_nc) {
            att=a[i]; bt=b[j];
            at1=a[i+1];
            if(*add) {
               t00 = CF[ij];
               t01 = CF[ij+(*dc)];
               }
            else t00=t01=0.0;
            for(k= *nl; k ; k--,att++,bt++,at1++) {
               t00 += *att * *bt;
               t01 += *at1 * *bt;
               }
            CF[ij]=t00;
            CF[ij+(*dc)]=t01;
            ij++;
            }
         }
      if(odd_nr) {
         ij=i*(*dc);
         for(j=0; j < (*nc)-1 ; j+=2,ij+=2) {
            att=a[i]; bt=b[j];
            bt1=b[j+1];
            if(*add) {
               t00 = CF[ij];
               t10 = CF[ij+1];
               }
            else t00=t10=0.0;
            for(k= *nl; k ; k--,att++,bt++,bt1++) {
               t00 += *att * *bt;
               t10 += *att * *bt1;
               }
            CF[ij]=t00;
            CF[ij+1]=t10;
            }
         if(odd_nc) {
            att=a[i]; bt=b[j];
            t00 = (*add) ? CF[ij] : 0.0;
            for(k= *nl; k ; k--,att++,bt++)
               t00 += *att * *bt;
            CF[ij]=t00;
            }
         }
      }
   else {
      for(j=ij=0; j < (*nc)-1 ; j+=2) {
         ij = j*(*dc);
         for(i=0; i < (*nr)-1 ; i+=2,ij+=2) {
            att=a[i]; bt=b[j];
            at1=a[i+1]; bt1=b[j+1];
            if(*add) {
               t00 = CF[ij];
               t01 = CF[ij+(*dc)];
               t10 = CF[ij+1];
               t11 = CF[ij+1+(*dc)];
               }
            else t00=t01=t10=t11=0.0;
            for(k= *nl; k ; k--,att++,bt++,at1++,bt1++) {
               t00 += *att * *bt;
               t01 += *att * *bt1;
               t10 += *at1 * *bt;
               t11 += *at1 * *bt1;
               }
            CF[ij]=t00;
            CF[ij+(*dc)]=t01;
            CF[ij+1]=t10;
            CF[ij+1+(*dc)]=t11;
            }
         if(odd_nr) {
            att=a[i]; bt=b[j];
            bt1=b[j+1];
            if(*add) {
               t00 = CF[ij];
               t01 = CF[ij+(*dc)];
               }
            else t00=t01=0.0;
            for(k= *nl; k ; k--,att++,bt++,bt1++) {
               t00 += *att * *bt;
               t01 += *att * *bt1;
               }
            CF[ij]=t00;
            CF[ij+(*dc)]=t01;
            ij++;
            }
         }
      if(odd_nc) {
         ij = j*(*dc);
         for(i=0; i < (*nr)-1 ; i+=2,ij+=2) {
            att=a[i]; bt=b[j];
            at1=a[i+1];
            if(*add) {
               t00 = CF[ij];
               t10 = CF[ij+1];
               }
            else t00=t10=0.0;
            for(k= *nl; k ; k--,att++,bt++,at1++) {
               t00 += *att * *bt;
               t10 += *at1 * *bt;
               }
            CF[ij]=t00;
            CF[ij+1]=t10;
            }
         if(odd_nr) {
            att=a[i]; bt=b[j];
            t00 = (*add) ? CF[ij] : 0.0;
            for(k= *nl; k ; k--,att++,bt++)
               t00 += *att * *bt;
            CF[ij]=t00;
            }
         }  
      }
   }
