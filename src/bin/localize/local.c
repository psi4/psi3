#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libipv1/ip_lib.h>
#include <libciomr/libciomr.h>
#include <libint/libint.h>
#include <libchkpt/chkpt.h>
#include <libiwl/iwl.h>
#include <psifiles.h>
#include <libqt/qt.h>

FILE *infile, *outfile;
char *psi_file_prefix;

void init_io(int argc, char *argv[]);
void exit_io(void);
void title(void);
void pipek_mezey(void);
void boys(int NAO, int NMO, int first_mo, int last_mo, double **CC);


#define METHOD_PIPEK_MEZEY 0
#define METHOD_BOYS 1


int main(int argc, char *argv[])
{
  int method_code = METHOD_PIPEK_MEZEY;
  int localize_virtuals = 0;
  char *method_string;
  int errcod, i, nirreps, nmo, nao, nocc, nfzc;
  int *clsdpi, *openpi, *frdocc;
  double **C;

  init_io(argc, argv);
  title();

  errcod = ip_string("LOCALIZATION_METHOD",&method_string,0);
  if (errcod == IPE_OK) {
    if (strcmp(method_string,"PIPEK_MEZEY")==0) {
      method_code = METHOD_PIPEK_MEZEY;
    }
    else if (strcmp(method_string,"BOYS")==0) {
      method_code = METHOD_BOYS;
    }
    else {
      fprintf(outfile, "Error: unrecognized LOCALIZATION_METHOD %s\n",
        method_string);
      exit(PSI_RETURN_FAILURE);
    }
  }

  /* this option only works with Boys at the moment */
  errcod = ip_boolean("LOCALIZE_VIRTUALS",&localize_virtuals,0);

  /* Pipek-Mezey localization */
  if (method_code == METHOD_PIPEK_MEZEY)
    pipek_mezey();

  /* Boys localization */
  else if (method_code == METHOD_BOYS) {
    chkpt_init(PSIO_OPEN_OLD);
    nao = chkpt_rd_nao();
    nmo = chkpt_rd_nmo();
    nirreps = chkpt_rd_nirreps();
    if(nirreps != 1) {
      fprintf(outfile, 
        "\n\tError: localization is only valid in C1 symmetry!\n");
      exit(PSI_RETURN_FAILURE);
    }
    clsdpi = chkpt_rd_clsdpi();
    openpi = chkpt_rd_openpi();
    for (i=0,nocc=0; i<nirreps; i++) nocc += clsdpi[i] + openpi[i];
    free(clsdpi);
    free(openpi);
    C = chkpt_rd_scf();
    
    /* localize the occupied orbitals (- frozen core) */
    frdocc = get_frzcpi();
    nfzc = frdocc[0];
    fprintf(outfile, "%d core orbitals\n", nfzc);
    free(frdocc);
    boys(nao, nmo, nfzc, nocc, C);

    /* if requested, attempt to localize the virtuals */
    if (localize_virtuals)
      boys(nao, nmo, nocc, nmo, C);

    chkpt_wt_scf(C);
    chkpt_close();
  }

  /* unrecognized localization option */
  else {
    fprintf(outfile, "Error: unrecognized METHOD_CODE %d\n", method_code);  
    exit(PSI_RETURN_FAILURE);
  } 

  free(method_string);
  exit_io();
  exit(PSI_RETURN_SUCCESS);
}


void pipek_mezey(void)
{
  int iter, s, t, A, k, l, m, p, q, inew, iold, max, max_col;
  int phase_ok, phase_chk;
  int i, j, ij, am, atom, shell_length, offset, stat;
  int nirreps, nao, nmo, nso, natom, nshell, noei, nocc, errcod, nfzc;
  int *stype, *snuc, *aostart, *aostop, *ao2atom, *l_length;
  int *clsdpi, *openpi, *orbspi, *dummy, *order, *frdocc;
  double *ss, **S, **scf, **u, **Ctmp, **C, *evals, **MO_S, **scf_old, **X;

  int *orb_order, *orb_boolean, puream;
  double P, PiiA, Pst, Pss, Ptt, Ast, Bst, AB;
  double Uss, Utt, Ust, Uts, Cks, Ckt, **U, **V, **VV, **F;
  double cos4a, alpha, alphamax, alphalast, conv, norm;
  int print;

  alphalast = 1.0;

  chkpt_init(PSIO_OPEN_OLD);
  nao = chkpt_rd_nao();
  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();
  natom = chkpt_rd_natom();
  nshell = chkpt_rd_nshell();
  stype = chkpt_rd_stype();
  snuc = chkpt_rd_snuc();
  u = chkpt_rd_usotao();
  nirreps = chkpt_rd_nirreps();
  clsdpi = chkpt_rd_clsdpi();
  openpi = chkpt_rd_openpi();
  orbspi = chkpt_rd_orbspi();
  C = chkpt_rd_scf();
  evals = chkpt_rd_evals();
  chkpt_close();

  /* A couple of error traps */
  if(nirreps != 1) {
    fprintf(outfile, "\n\tError: localization is only valid in C1 symmetry!\n");
    exit(PSI_RETURN_FAILURE);
  }
  if(openpi[0]) {
    fprintf(outfile, "\n\tError: localization available for closed-shells only!\n");
    exit(PSI_RETURN_FAILURE);
  }

  /* Frozen orbital info */
  frdocc = get_frzcpi();
  nfzc = frdocc[0];
  free(frdocc);

  /* Compute the length of each AM block */
  ip_boolean("PUREAM", &(puream), 0);
  l_length = init_int_array(LIBINT_MAX_AM);
  l_length[0] = 1;
  for(l=1; l < (LIBINT_MAX_AM); l++) {
    if(puream) l_length[l] = 2 * l + 1;
    else l_length[l] = l_length[l-1] + l + 1;
  }

  /* Set up the atom->AO and AO->atom lookup arrays */
  aostart = init_int_array(natom);
  aostop = init_int_array(natom);
  for(i=0,atom=-1,offset=0; i < nshell; i++) {
    am = stype[i] - 1;
    shell_length = l_length[am];

    if(atom != snuc[i]-1) {
      if(atom != -1) aostop[atom] = offset-1;
      atom = snuc[i]-1;
      aostart[atom] = offset;
    }

    offset += shell_length;
  }
  aostop[atom] = offset-1;

  ao2atom = init_int_array(nso);
  for(i=0; i < natom; i++)
    for(j=aostart[i]; j <= aostop[i]; j++) ao2atom[j] = i;

  /* Get the overlap integrals -- these should be identical to AO S */
  noei = nso*(nso+1)/2;
  ss = init_array(noei);
  stat = iwl_rdone(PSIF_OEI,PSIF_SO_S,ss,noei,0,0,outfile);
  S = block_matrix(nso,nso);
  for(i=0,ij=0; i < nso; i++)
    for(j=0; j <= i; j++,ij++) {
      S[i][j] = S[j][i] = ss[ij];
    }
  free(ss);

  /* Compute nocc --- closed-shells only */
  for(i=0,nocc=0; i < nirreps; i++) nocc += clsdpi[i];

  fprintf(outfile, "\tNumber of doubly occupied orbitals: %d\n\n", nocc);

  fprintf(outfile, "\tIter     Pop. Localization   Max. Rotation Angle       Conv\n");
  fprintf(outfile, "\t------------------------------------------------------------\n");

  U = block_matrix(nocc, nocc);
  V = block_matrix(nocc, nocc);
  for(i=0; i < nocc; i++) V[i][i] = 1.0;
  VV = block_matrix(nocc, nocc);

  for(iter=0; iter < 100; iter++) {

    P = 0.0;
    for(i=nfzc; i < nocc; i++) {
      for(A=0; A < natom; A++) {
	PiiA = 0.0;

	for(l=aostart[A]; l <= aostop[A]; l++)
	  for(k=0; k < nso; k++) 
	    PiiA += C[k][i] * C[l][i] * S[k][l];

	P += PiiA * PiiA;
      }
    }

    /* Compute 2x2 rotations for Pipek-Mezey localization */
    alphamax = 0.0;
    for(s=nfzc; s < nocc; s++) {
      for(t=nfzc; t < s; t++) {

	Ast = Bst = 0.0;

	for(A=0; A < natom; A++) {

	  Pst = Pss = Ptt = 0.0;
	      
	  for(l=aostart[A]; l <= aostop[A]; l++) {
	    for(k=0; k < nso; k++) {
	      Pst += 0.5 * (C[k][s] * C[l][t] +
			    C[l][s] * C[k][t]) * S[k][l];

	      Pss += C[k][s] * C[l][s] * S[k][l];

	      Ptt += C[k][t] * C[l][t] * S[k][l];
	    }
	  }

	  Ast += Pst * Pst - 0.25 * (Pss - Ptt) * (Pss - Ptt);
	  Bst += Pst * (Pss - Ptt);

	} /* A-loop */

	/* Compute the rotation angle */
	AB = Ast * Ast + Bst * Bst;
	if(fabs(AB) > 0.0) {
	  cos4a = -Ast/sqrt(AB);
	  alpha = 0.25 * acos(cos4a) * (Bst > 0 ? 1 : -1);
	}
	else alpha = 0.0;

	/* Keep up with the maximum 2x2 rotation angle */
	alphamax = (fabs(alpha) > alphamax ? alpha : alphamax);

	Uss = cos(alpha);
	Utt = cos(alpha);
	Ust = sin(alpha);
	Uts = -Ust;

	/* Now do the rotation */
	for(k=0; k < nso; k++) {
	  Cks = C[k][s];
	  Ckt = C[k][t];
	  C[k][s] = Uss * Cks + Ust * Ckt;
	  C[k][t] = Uts * Cks + Utt * Ckt;
	}

	zero_mat(U, nocc, nocc);
	for(i=0; i < nocc; i++) U[i][i] = 1.0;

	U[s][s] = Uss;
	U[t][t] = Utt;
	U[s][t] = Ust;
	U[t][s] = Uts;

	zero_mat(VV, nocc, nocc);
	for(i=0; i < nocc; i++) {
	  for(j=0; j < nocc; j++) {
	    for(k=0; k < nocc; k++) {
	      VV[i][j] += V[i][k] * U[j][k];
	    }
	  }
	}

	for(i=0; i < nocc; i++)
	  for(j=0; j < nocc; j++)
	    V[i][j] = VV[i][j];
	      
      } /* t-loop */
    } /* s-loop */

    conv = fabs(alphamax) - fabs(alphalast);
    fprintf(outfile, "\t%4d  %20.10f  %20.10f  %4.3e\n", iter, P, alphamax, conv);
    if((iter > 2) && ((fabs(conv) < 1e-12) || alphamax == 0.0)) break;
    alphalast = alphamax;

    fflush(outfile);
      
  } /* iter-loop */

  /*  print_mat(V, nocc, nocc, outfile);  */

  /* Transform occupied orbital eigenvalues */
  F = block_matrix(nocc, nocc);
  for(i=0; i < nocc; i++)
    for(j=0; j < nocc; j++)
      for(k=0; k < nocc; k++) 
	F[i][j] += V[k][i] * evals[k] * V[k][j];

  /*
    fprintf(outfile, "\nTransformed Orbital Energies:\n");
    print_mat(F, nocc, nocc, outfile);
  */

  /* Compute a reordering array based on the diagonal elements of F */
  orb_order = init_int_array(nocc);
  orb_boolean = init_int_array(nocc);
  for(i=0; i < nocc; i++) { orb_order[i] = 0;  orb_boolean[i] = 0; }

  for(i=0,max=0; i < nocc; i++) /* First, find the overall maximum */
    if(fabs(F[i][i]) > fabs(F[max][max])) max = i;

  orb_order[0] = max;  orb_boolean[max] = 1;

  for(i=1; i < nocc; i++) {
    max = 0;
    while(orb_boolean[max]) max++; /* Find an unused max */
    for(j=0; j < nocc; j++) 
      if((fabs(F[j][j]) >= fabs(F[max][max])) && !orb_boolean[j]) max = j;
    orb_order[i] = max; orb_boolean[max] = 1;
  }

  /*
    for(i=0; i < nocc; i++) fprintf(outfile, "%d %d\n", i, orb_order[i]);
  */

  /*
    fprintf(outfile, "\n\tPipek-Mezey Localized MO's (before sort):\n");
    print_mat(C, nso, nmo, outfile);
  */

  /* Now reorder the localized MO's according to F */
  Ctmp = block_matrix(nso,nocc);
  for(i=0; i < nocc; i++)
    for(j=0; j < nso; j++) Ctmp[j][i] = C[j][i];

  for(i=0; i < nocc; i++) {
    iold = orb_order[i];
    for(j=0; j < nso; j++) C[j][i] = Ctmp[j][iold];
    evals[i] = F[iold][iold];
  }
  free_block(Ctmp);

  print = 0;
  errcod = ip_boolean("PRINT_MOS", &(print), 0);
  if(print) {
    fprintf(outfile, "\n\tPipek-Mezey Localized MO's (after sort):\n");
    print_mat(C, nso, nmo, outfile);
  }

  /* Check MO normalization */
  /*
    for(i=0; i < nmo; i++) {
    norm = 0.0;
    for(j=0; j < nso; j++) 
    for(k=0; k < nso; k++) {
    norm += C[j][i] * C[k][i] * S[j][k];
    }

    fprintf(outfile, "norm[%d] = %20.10f\n", i, norm);
    }
  */

  /* correct orbital phases for amplitude restarts */
  chkpt_init(PSIO_OPEN_OLD);
  scf_old = chkpt_rd_local_scf();
  chkpt_close();
  if (scf_old != NULL) {
    MO_S = block_matrix(nmo, nmo);
    X = block_matrix(nso, nmo);
    C_DGEMM('n','n',nso, nmo, nso, 1, &(S[0][0]), nso, &(C[0][0]), nmo,
	    0, &(X[0][0]), nmo);
    C_DGEMM('t','n',nmo, nmo, nso, 1, &(scf_old[0][0]), nmo, &(X[0][0]), nmo,
	    0, &(MO_S[0][0]), nmo);
    free_block(X);

    /*
    fprintf(outfile, "Approximate Overlap Matrix\n");
    print_mat(MO_S, nmo, nmo, outfile);
    */

    for(p=0; p < nmo; p++) {
      max = 0.0;
      for(q=0; q < nmo; q++) {
	if(fabs(MO_S[p][q]) > max) {
	  max = fabs(MO_S[p][q]); max_col = q;
	}
      }
      if(max_col != p) phase_ok = 0;
    }

    chkpt_init(PSIO_OPEN_OLD);
    chkpt_wt_phase_check(phase_ok);
    chkpt_close();
    if(phase_ok) {
      for(p=0; p < nmo; p++) {
	if(MO_S[p][p] < 0.0) {
	  for(q=0; q < nso; q++)
	    C[q][p] *= -1.0;
	}
      }
    }

    free_block(MO_S);
    free_block(scf_old);
  }
  free_block(S);

  /* Write the new MO's to chkpt */
  chkpt_init(PSIO_OPEN_OLD);
  chkpt_wt_scf(C);
  chkpt_wt_local_scf(C);
  chkpt_close();

  free_block(C);
  free(evals);
  free(clsdpi);
  free(openpi);
  free(orbspi);

  fprintf(outfile, "\n\tLocalization of occupied orbitals complete.\n");

}


void boys(int NAO,int NMO,int first,int last,double **CC)
{
  int mu,i,j,ij,k,l,count;
  int iDmax,jDmax;
  double *dmxAOarray;
  double *dmyAOarray;
  double *dmzAOarray;
  double **dmxAO;
  double **dmyAO;
  double **dmzAO;
  double **dmxMO;
  double **dmyMO;
  double **dmzMO;
  double *dmMO_diag;
  double **Ax;
  double **Ay;
  double **Az;
  double **Bx;
  double **By;
  double **Bz;
  double **CC_sorted;
  double D, DOld=0.0;
  double Dmax;
  double nom,denom;
  double cosa,sina,cos4a,sin4a;
  double x2plus,x2minus;
  double x1,x2,y1,y2,sin4acandidate1,sin4acandidate2;
  double AngleEpsilon,DChangeEpsilon;
  int boyscountMax = 200;
  double tval, tval1, tval2;
  double minD;
  int minj, *sort_done;

  AngleEpsilon=1.0E-6;
  DChangeEpsilon=1.0E-6;
                                                                                
  fprintf(outfile, "\nPerforming Boys localization for %d AO's and %d MO's\n",
    NAO, NMO);
  fprintf(outfile, "Localizing orbitals %d through %d\n", first, last-1);

  fprintf(outfile, "\tIter     D           Conv\n");
  fprintf(outfile, "\t----------------------------------\n");

  dmxAOarray=init_array(NAO*(NAO+1)/2);
  dmyAOarray=init_array(NAO*(NAO+1)/2);
  dmzAOarray=init_array(NAO*(NAO+1)/2);

  iwl_rdone(PSIF_OEI,PSIF_AO_MX,dmxAOarray,NAO*(NAO+1)/2,0,0,outfile);
  iwl_rdone(PSIF_OEI,PSIF_AO_MY,dmyAOarray,NAO*(NAO+1)/2,0,0,outfile);
  iwl_rdone(PSIF_OEI,PSIF_AO_MZ,dmzAOarray,NAO*(NAO+1)/2,0,0,outfile);
                                                                                
  dmxAO=block_matrix(NAO,NAO);
  dmyAO=block_matrix(NAO,NAO);
  dmzAO=block_matrix(NAO,NAO);
                                                                                
  for(i=0,ij=0;i<NAO;i++){
    for(j=0;j<=i;j++,ij++){
      dmxAO[i][j]=dmxAOarray[ij];
      dmyAO[i][j]=dmyAOarray[ij];
      dmzAO[i][j]=dmzAOarray[ij];
      dmxAO[j][i]=dmxAO[i][j];
      dmyAO[j][i]=dmyAO[i][j];
      dmzAO[j][i]=dmzAO[i][j];
    }
  }
                                                                                
  free(dmxAOarray);
  free(dmyAOarray);
  free(dmzAOarray);
                                                                                
  dmxMO=block_matrix(NMO,NMO);
  dmyMO=block_matrix(NMO,NMO);
  dmzMO=block_matrix(NMO,NMO);
                                                                                
  Ax=block_matrix(NMO,NMO);
  Ay=block_matrix(NMO,NMO);
  Az=block_matrix(NMO,NMO);
  Bx=block_matrix(NMO,NMO);
  By=block_matrix(NMO,NMO);
  Bz=block_matrix(NMO,NMO);

  if (first >= last) return;
  
  for(count=0;count<boyscountMax;count++){
    for(i=first;i<last;i++){
      for(j=first;j<last;j++){
        dmxMO[i][j]=0.0;
        dmyMO[i][j]=0.0;
        dmzMO[i][j]=0.0;
        for(k=0;k<NAO;k++){
          for(l=0;l<NAO;l++){
            dmxMO[i][j]+=CC[k][i]*CC[l][j]*dmxAO[k][l];
            dmyMO[i][j]+=CC[k][i]*CC[l][j]*dmxAO[k][l];
            dmzMO[i][j]+=CC[k][i]*CC[l][j]*dmxAO[k][l];
          }
        }
      }
    }
                                                                                
    D=0.0;
    for(i=first;i<last;i++){
      D+=dmxMO[i][i]*dmxMO[i][i]+dmyMO[i][i]*dmyMO[i][i]+
        dmzMO[i][i]*dmzMO[i][i];
    }
                                                                                
    fprintf(outfile,"\t%3d  %12.6lf    %2.2E\n",count,D,fabs(D-DOld));
    if(count!=0){
      if(fabs(D-DOld)<DChangeEpsilon){
        fprintf(outfile,"\nConvergence to %2.2E is achieved.\n",
          count,DChangeEpsilon);
        break;
      }
    }
                                                                                
    DOld=D;
                                                                                
    for(i=first;i<last;i++){
      for(j=first;j<last;j++){
        tval = dmxMO[i][i]-dmxMO[j][j];
        Ax[i][j]=dmxMO[i][j]*dmxMO[i][j]-0.25*tval*tval;
        Bx[i][j]=dmxMO[i][j]*tval;
        tval = dmyMO[i][i]-dmyMO[j][j];
        Ay[i][j]=dmyMO[i][j]*dmyMO[i][j]-0.25*tval*tval;
        By[i][j]=dmyMO[i][j]*tval;
        tval = dmzMO[i][i]-dmzMO[j][j];
        Az[i][j]=dmzMO[i][j]*dmzMO[i][j]-0.25*tval*tval;
        Bz[i][j]=dmzMO[i][j]*tval;
      }
    }
                                                                                
    Dmax=0.0;
    iDmax=1;
    jDmax=1;
    for(i=first;i<last;i++){
      for(j=first;j<i;j++){
        D=Ax[i][j]+Ay[i][j]+Az[i][j]+
          sqrt(Ax[i][j]*Ax[i][j]+Bx[i][j]*Bx[i][j])+
          sqrt(Ay[i][j]*Ay[i][j]+By[i][j]*By[i][j])+
          sqrt(Az[i][j]*Az[i][j]+Bz[i][j]*Bz[i][j]);
        if(D>Dmax){
          Dmax=D;
          iDmax=i;
          jDmax=j;
        }
      };
    };
                                                                                
    /*printf("Dmax occupied  = %lf\n",Dmax);*/
                                                                                
    nom=Ax[iDmax][jDmax]+Ay[iDmax][jDmax]+Az[iDmax][jDmax];
    denom=sqrt(Ax[iDmax][jDmax]*Ax[iDmax][jDmax]+
               Bx[iDmax][jDmax]*Bx[iDmax][jDmax])
      +sqrt(Ay[iDmax][jDmax]*Ay[iDmax][jDmax]+
           By[iDmax][jDmax]*By[iDmax][jDmax])+
      +sqrt(Az[iDmax][jDmax]*Az[iDmax][jDmax]+
            Bz[iDmax][jDmax]*Bz[iDmax][jDmax]);
    cos4a=-nom/denom;
                                                                                
    nom=Bx[iDmax][jDmax]+By[iDmax][jDmax]+Bz[iDmax][jDmax];
    sin4a=nom/denom;
                                                                                
    x2plus=0.5*(1.0+sqrt((1.0-0.5*(1.0-cos4a))));
    x2minus=0.5*(1.0-sqrt((1.0-0.5*(1.0-cos4a))));
                                                                                
    x1=sqrt(x2plus);
    x2=sqrt(x2minus);
                                                                                
    y1=sqrt(1-x2plus);
    y2=sqrt(1-x2minus);
                                                                                
    sin4acandidate1=4.0*x1*y1*(x1*x1-y1*y1);
    sin4acandidate2=4.0*x2*y2*(x2*x2-y2*y2);
                                                                                
    /*printf("sin4a occupied = %lf\n",sin4a);
    printf("sin4acandidate1 occupied = %lf\n",sin4acandidate1);
    printf("sin4acandidate2 occupied = %lf\n",sin4acandidate2);*/
                                                                                
    if (fabs(sin4acandidate1-sin4a)<AngleEpsilon){
      cosa=x1;
      sina=y1;
      /*printf("cos(a) occupied = %lf\n",cosa);
        printf("sin(a) occupied = %lf\n",sina);*/
    }
    else{
      if (fabs(sin4acandidate2-sin4a)<AngleEpsilon){
        cosa=x2;
        sina=y2;
        /*printf("cos(a) occupied = %lf\n",cosa);
          printf("sin(a) occupied = %lf\n",sina);*/
      }
      else{
        printf("Cannot find the rotation angle in occupieds!\n");
      }
    }
                                                                                
    for(mu=0;mu<NAO;mu++){
      tval1 =  cosa*CC[mu][iDmax]+sina*CC[mu][jDmax];
      tval2 = -sina*CC[mu][iDmax]+cosa*CC[mu][jDmax];
      CC[mu][iDmax] = tval1;
      CC[mu][jDmax] = tval2;
    }
                                                                                
  } /* end boys iterations */

  /* now sort the MO's according to the size of dm_MO[] */
  CC_sorted = block_matrix(NAO, last-first+1);
  dmMO_diag = init_array(last-first);
  for (i=first; i<last; i++) {
    dmMO_diag[i-first] = dmxMO[i][i]*dmxMO[i][i] 
                        + dmyMO[i][i]*dmyMO[i][i]
                        + dmzMO[i][i]*dmzMO[i][i];

  }

  sort_done = init_int_array(last-first);

  for (i=0; i<last-first; i++) {

    /* find the next biggest element of dmMO_diag */
    minD = 100000.0;
    minj = 0;
    for (j=0; j<last-first; j++) {
      if ((dmMO_diag[j] <= minD) && !sort_done[j]) {
        minD = dmMO_diag[j];
        minj = j;
      }
    }
    sort_done[minj] = 1;
    for (k=0; k<NAO; k++) {
      CC_sorted[k][i] = CC[k][minj+first]; 
    }
  }

  for (i=0; i<NAO; i++) {
    for (j=first,k=0; j<last; j++,k++) {
      CC[i][j] = CC_sorted[i][k];
    }
  }

  free_block(dmxAO);  free_block(dmyAO);  free_block(dmzAO);
  free_block(dmxMO);  free_block(dmyMO);  free_block(dmzMO);
  free_block(Ax);     free_block(Ay);     free_block(Az);
  free_block(Bx);     free_block(By);     free_block(Bz);
  free(dmMO_diag);    free(sort_done);
  free_block(CC_sorted);
}




void init_io(int argc, char * argv[])
{
  extern char *gprgid();
  char *progid;

  progid = (char *) malloc(strlen(gprgid())+2);
  sprintf(progid, ":%s",gprgid());

  psi_start(argc-1,argv+1,0);
  ip_cwk_add(":INPUT"); /* for checking puream keyword */
  ip_cwk_add(progid);
  free(progid);
  tstart(outfile);

  psio_init();
}

void title(void)
{
  fprintf(outfile, "\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t*         LOCALIZE       *\n");
  fprintf(outfile, "\t\t\t*                        *\n");
  fprintf(outfile, "\t\t\t**************************\n");
  fprintf(outfile, "\n");
}

void exit_io(void)
{
  psio_done();
  tstop(outfile);
  psi_stop();
}

char *gprgid()
{
   char *prgid = "LOCALIZE";

   return(prgid);
}
