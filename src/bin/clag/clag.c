/****************************************************************************/
/* clag: the main controlling program for calculating the lagrangian and CI */ 
/*      energy. The lagrangian is written to file 75 and the Ci energy is   */
/*      printed in the output as a simple check                             */
/****************************************************************************/

#include <stdio.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libiwl/iwl.h>
#include <libipv1/ip_lib.h>
#include <libfile30/file30.h>
#include <psifiles.h>
#include <libpsio/psio.h>
#include <math.h>
#include "clag.h"
#define INDEX(i,j) ( (i>j) ? (ioff[(i)] + (j)): (ioff[(j)] + (i)) )

/*
** define input parsing files and ioff array 
*/

FILE *outfile;    /* pointer to the file output.dat */
FILE *infile;     /* pointer to the file input.dat  */
int *ioff;        /* the ioff array                 */
int print_lvl=1;  /* diagnostic info flag           */

/***************************************************************************/
/* The main procedure                                                      */
/***************************************************************************/
main(int argc, char **argv) 

{

  double **opdm;                       /* the one particle density matrix */
  double *tpdm;                        /* the two particle density matrix */
  double **lag;                        /* the lagrangian we are finding   */
  int i,j,ij;                          /* a simple running variable       */
  int errcod;                          /* error flag for input parsing    */ 
  int ntri, ntri2;                     /* number of one and two e ints    */  
  int nmo;                             /* number of molecular orbitals    */
  int nfzv;                            /* number of frozen virtual orbs   */
  int npop;                            /* number of populated orbitals;
                                          or nmo - nfzv                   */
  int *orbspi;                         /* orbitals per irrep array        */
  int *docc;                           /* doubly occupied orbs per irrep  */
  int *socc;                           /* singly occupied orbs per irrep  */
  int *frdocc;                         /* frozen doubly occupied array    */
  int *fruocc;                         /* frozen unoccupied orb array     */
  int **ras_opi;                       /* orbs per [ras_space][irrep]     */
  int *pitz_to_corr;                   /* map orbs Pitzer->correlated ord */
  int *corr_to_pitz;                   /* map orbs correlated order->Pitz */
  int nirreps;                         /* number of irreps                */
  double efzc;                         /* frozen core energy              */
  int oei_file = PSIF_OEI;             /* where 1e mo ints are stored     */
  int oei_erase = 0;                   /* 0=do not erase 1e ints 1=do     */
  int tei_file = PSIF_MO_TEI;          /* where 2e mo ints are stored     */
  int opdm_file = PSIF_MO_OPDM;        /* file number for one-pdm         */ 
  int tpdm_file = PSIF_MO_TPDM;        /* file number for two-pdm         */ 
  int lag_file = PSIF_MO_LAG;          /* file number for largrangian     */

  /* the following was for Yukio Yamaguchi's CAS code, I think 
   * normally we would be using DETCAS now
   */

  int cas_onel_file = 81;              /* CAS interface one-elec ints file*/ 
  int cas_twoel_file = 82;             /* CAS interface two-elec ints file*/
  int cas_opdm_file = 83;              /* CAS interface onepdm file       */
  int cas_tpdm_file = 84;              /* CAS interface twopdm file       */
  int cas_lag_file = 85;               /* CAS interface lagrangian file   */
  int write_cas_files = 0;             /* write out a files for CASSCF?   */

  double *onel_ints, *twoel_ints;      /* 1e and 2e ints                  */
  double enuc = 0.0;                   /* nuclear repulsion energy        */ 
  double eci_30;                       /* ci energy from file 30          */
  double lagtr;                        /* trace of lagrangian             */

  for (i=1; i<argc; i++) 
     if (strcmp("-quiet", argv[i]) == 0) print_lvl = 0; 
  
  /*
  ** initialize the io parser
  */
  init_io(); 
  ip_set_uppercase(1);
  ip_initialize(infile, outfile);
  ip_cwk_clear();
  ip_cwk_add(":DEFAULT");  
  ip_cwk_add(":CLAG");  
  errcod = ip_data("PRINT","%d", &print_lvl,0);  
  errcod = ip_boolean("WRITE_CAS_FILES",&write_cas_files,0);

  /*
  ** print out header information
  */
  fprintf(outfile,"CLAG: PROGRAM TO FORM LAGRANGIAN AND CALCULATE CI ENERGY\n");
  fprintf(outfile,"WRITTEN BY BRIAN HOFFMAN, MATT LEININGER");
  fprintf(outfile," AND DAVID SHERRILL\n\n"); 

  /*
  ** calculate some needed numbers 
  */
  file30_init();
  nmo = file30_rd_nmo();
  enuc = file30_rd_enuc();
  eci_30 = file30_rd_ecorr(); 
  nirreps = file30_rd_nirreps();
  orbspi = file30_rd_orbspi();
  file30_close();

  docc = init_int_array(nirreps);
  socc = init_int_array(nirreps);
  frdocc = init_int_array(nirreps);
  fruocc = init_int_array(nirreps);
  ras_opi = init_int_matrix(4,nirreps);
  pitz_to_corr = init_int_array(nmo);

  /* get orbital information */
  ras_set(nirreps, nmo, 0, orbspi, docc, socc, frdocc, fruocc, ras_opi,
          pitz_to_corr, 1);

  /* get the array which maps correlated orbitals back to pitzer order */
  corr_to_pitz = init_int_array(nmo);
  for (i=0; i<nmo; i++) {
    j = pitz_to_corr[i];
    corr_to_pitz[j] = i;
  }

  for (i=0,nfzv=0; i<nirreps; i++) nfzv += fruocc[i];
  npop = nmo - nfzv; 
  ntri = (nmo*(nmo+1))/2;   
  ntri2 = (ntri*(ntri+1))/2; 

  /*
  ** set up the ioff array
  */
  ioff = init_int_array(nmo*nmo+1);  
  for (i=1; i<nmo*nmo+1; i++)
     ioff[i] = ioff[i-1] + i;

  /*
  ** read in the integral and the density matricies
  */
  opdm = rdopdm(npop, print_lvl, opdm_file); 
  tpdm = rdtpdm(npop, print_lvl, tpdm_file);   

  onel_ints = init_array(ntri); 
  twoel_ints = init_array(ntri2); 

  if (print_lvl>4) {
    fprintf(outfile, "\nOne-electron integrals\n");
  }

  if (!iwl_rdone(oei_file, PSIF_MO_OEI, onel_ints, ntri, oei_erase,
            (print_lvl>4), outfile)) {
    fprintf(outfile, "Failed to read one-electron integrals\n");
    exit(1);
  }

  if (print_lvl>4) {
    fprintf(outfile, "\nTwo-electron integrals\n");
  }

  iwl_rdtwo(tei_file, twoel_ints, ioff, nmo, 0, 0, (print_lvl>4), outfile); 


  /* 
  ** test the trace of the pdms
  */
  trace_opdm(opdm, npop);
  trace_tpdm(tpdm, npop);


  /*
  ** form lagrangian matrix and write to file 75
  */
  lag = block_matrix(nmo,nmo);
  lagtr = lagcalc(opdm,tpdm,onel_ints,twoel_ints,lag,nmo,npop,
                  print_lvl,lag_file); 
  ci_energy(opdm, tpdm, onel_ints, twoel_ints, npop, enuc, eci_30, lagtr); 
  

  /*
  ** write out the two-pdm in a form that the CAS program will like
  */
  if (write_cas_files) {
    onel_to_cas(onel_ints, corr_to_pitz, nmo, print_lvl, cas_onel_file);
    twoel_to_cas(twoel_ints, corr_to_pitz, nmo, print_lvl, cas_twoel_file);
    onepdm_to_cas(opdm, corr_to_pitz, nmo, npop, print_lvl, cas_opdm_file);
    twopdm_to_cas(tpdm, corr_to_pitz, nmo, npop, print_lvl, cas_tpdm_file);
    lag_to_cas(lag, corr_to_pitz, nmo, print_lvl, cas_lag_file);
  }

  /*
  ** free memory
  */
  free(docc); free(socc);
  free(frdocc); free(fruocc);
  free_int_matrix(ras_opi, 4);
  free(onel_ints);
  free(twoel_ints);
  free_block(opdm);
  free(tpdm);
  free_block(lag);


  /*
  ** close files and end the program
  */
  close_io(); 
  return(0);
}


/****************************************************************************/
/* init_io(): Function opens input and output files                         */
/****************************************************************************/
void init_io(void)
{
   ffile(&infile,"input.dat",2) ;
   ffile(&outfile,"output.dat",1);
   if (print_lvl > 0) tstart(outfile) ;
   psio_init();
}

/****************************************************************************/
/* close_io(): Function closes down I/O and exits                           */
/****************************************************************************/

void close_io(void)
{
   psio_done();
   fclose(infile);
   if (print_lvl > 0) tstop(outfile);
   fclose(outfile);
}

/****************************************************************************/
/* rdopdm: reads the one particle density matrix from opdm_file             */
/*         and returns opdm as an in core matrix                            */
/****************************************************************************/
double **rdopdm(int nbf, int print_lvl, int opdm_file)
{
 
 int i, root, errcod;
 double **opdm; 
 PSI_FPTR index = 0; 

 rfile(opdm_file); 
 if (flen(opdm_file)==0) {
   fprintf(outfile,"\tCan't find the one-particle density matrix file %d\n",
     opdm_file);
   rclose(opdm_file, 4);
   exit(1);
  }

 opdm = block_matrix(nbf, nbf); 
 root = 1;
 errcod = ip_data("ROOT","%d",&root,0);
/* fprintf(outfile,"GOH Root = %d\n",root); */
/* index = 0; */
 for (i =0; i<root; i++) {
   wreadw(opdm_file, (char *) opdm[0], 
             (sizeof(double)*nbf*nbf), index, &index);
 }

 if (print_lvl > 2) {
   fprintf(outfile,"One-Particle Density Matrix\n");
   print_mat(opdm, nbf, nbf, outfile); 
   fprintf(outfile,"\n\n"); 
   }

 rclose (opdm_file,3);
 return (opdm);
}

/****************************************************************************/
/* rdtpdm: reads the two particle density matrix from tpdm_file             */
/*         and returns the tpdm as an array                                 */
/****************************************************************************/
double *rdtpdm(int nbf, int print_lvl, int tpdm_file)
{

 double *tpdm;
 int numslots, sqnbf;   
 int *ioff_lt, i;                    /* offsets for left (or right) indices */
 PSI_FPTR index = 0 ; 
 struct iwlbuf TBuff;

 iwl_buf_init(&TBuff, tpdm_file, 0.0, 1, 1);

 sqnbf = nbf*nbf ; 
 numslots = (sqnbf*(sqnbf+1))/2 ;  
 tpdm = init_array(numslots);  

 /* Construct the ioff_lt array (same here as ioff_rt) : different than
  * regular ioff because there is no perm symmetry between left indices
  * or right indices.
  */
 ioff_lt = init_int_array(nbf);
 for (i=0; i<nbf; i++) {
   ioff_lt[i] = i * nbf;
 }
 
 iwl_buf_rd_all(&TBuff, tpdm, ioff_lt, ioff_lt, 1, ioff, 
                (print_lvl>5), outfile);
  
 if (print_lvl > 3) {
   fprintf(outfile,"Two-Particle Density Matrix\n");
   print_array(tpdm, sqnbf, outfile);
   fprintf(outfile,"\n\n"); 
   }

 iwl_buf_close(&TBuff, 1);
 free(ioff_lt);
 return (tpdm);
}  


/***************************************************************************/
/* trace_opdm: test the trace of the one-particle density matrix           */
/***************************************************************************/
void trace_opdm(double **opdm, int nbf)
{
  int i;
  double sum;

  for (sum=0.0,i=0; i<nbf; i++) {
    sum += opdm[i][i];
  }

  fprintf(outfile, "\n\tTrace of one-pdm = %16.12lf\n", sum);

}


/***************************************************************************/
/* trace_tpdm: test the trace of the two-particle density matrix           */
/***************************************************************************/
void trace_tpdm(double *tpdm, int nbf)
{
  int i,j,ii,jj,iijj;
  double sum;

  for (sum=0.0,i=0; i<nbf; i++) {
    ii = i*nbf + i;
    for (j=0; j<nbf; j++) {
      jj = j*nbf + j;
      iijj = INDEX(ii,jj);
      sum += tpdm[iijj];
    }
  }

  fprintf(outfile, "\tTrace of two-pdm = %16.12lf\n", sum);

}


/***************************************************************************/
/* gpgrid: program id                                                      */
/***************************************************************************/
char *gprgid()
{
   char *prgid = "CLAG";

   return(prgid);
}

