#include "scf.h"
#include <liboptions/liboptions.h>
#include <libmoinfo/libmoinfo.h>
#include "memory_manager.h"

#include <iostream>
#include <cmath>

extern FILE* outfile;

namespace psi{ namespace mcscf{

SCF::SCF()
{
  startup();
}

SCF::~SCF()
{
  cleanup();
}

void SCF::startup()
{
  ioff    = moinfo_scf->get_ioff();
  nirreps = moinfo_scf->get_nirreps();
  nso     = moinfo_scf->get_nso();
  sopi    = moinfo_scf->get_sopi();
  docc    = moinfo_scf->get_docc();
  actv    = moinfo_scf->get_actv();

  // Compute the block_offset for the so basis
  allocate1(double,block_offset,nirreps);
  block_offset[0] = 0;
  for(int h=1; h < nirreps; ++h){
    block_offset[h] = block_offset[h-1] + sopi[h-1];
  }

  // Allocate the pairs
  allocate1(int,pairpi,nirreps);
  allocate1(int,pair_offset,nirreps);
  allocate2(int,pair,nso,nso);
  allocate2(int,pair_sym,nso,nso);
  pairs = 0;
  PK = 0;
  K = 0;

  if(options_get_str("REFERENCE")=="RHF"){
    reference = rhf;
  }
  else if(options_get_str("REFERENCE")=="ROHF"){
    reference = rohf;
  }
  else if(options_get_str("REFERENCE")=="UHF"){
    reference = uhf;
  }
  else if(options_get_str("REFERENCE")=="TCSCF"){
    reference = tcscf;
  }

  root = options_get_int("ROOT") - 1;

  // OUT OF CORE ALGORITHM
  out_of_core  = false;

  // DIIS
  use_diis = options_get_bool("USE_DIIS");
  ndiis    = options_get_int("NDIIS");
  current_diis = 0;

  turn_on_actv = options_get_int("TURN_ON_ACTV");

  epsilon     .allocate("epsilon",nirreps,sopi);

  e           .allocate("e",nirreps,sopi,sopi);
  C           .allocate("C",nirreps,sopi,sopi);
  C_t         .allocate("C_t",nirreps,sopi,sopi);
  C_T         .allocate("C_T",nirreps,sopi,sopi);
  Dc          .allocate("Dc",nirreps,sopi,sopi);
  Feff_t      .allocate("Feff_t",nirreps,sopi,sopi);
  Feff_oAO    .allocate("Feff_oAO",nirreps,sopi,sopi);
  Fc          .allocate("Fc",nirreps,sopi,sopi);
  Fc_t        .allocate("Fc_t",nirreps,sopi,sopi);
  G           .allocate("G",nirreps,sopi,sopi);
  H           .allocate("H",nirreps,sopi,sopi);
  O           .allocate("O",nirreps,sopi,sopi);
  S           .allocate("S",nirreps,sopi,sopi);
  S_sqrt_inv  .allocate("S^-1/2",nirreps,sopi,sopi);
  S_sqrt      .allocate("S^1/2",nirreps,sopi,sopi);
  T           .allocate("T",nirreps,sopi,sopi);

  for(int i = 0; i < ndiis; ++i){
    diis_F[i] .allocate("diis_F[" + to_string(i) + "]",nirreps,sopi,sopi);
    diis_e[i] .allocate("diis_e[" + to_string(i) + "]",nirreps,sopi,sopi);
  }

  if(reference == rohf){
    Do        .allocate("Do",nirreps,sopi,sopi);
    Fo        .allocate("Fo",nirreps,sopi,sopi);
    Fo_t      .allocate("Fo_t",nirreps,sopi,sopi);
  }
  if(reference == tcscf){
    int count = 0;
    for(int h = 0; h < nirreps; ++h){
      for(int n = 0; n < actv[h]; ++n){
        tcscf_sym[count] = h;
        tcscf_mos[count] = docc[h] + n;
        count++;
      }
    }

    nci = count;
    fprintf(outfile,"\n  TCSCF MOs = [");
    for(int I = 0; I < nci; ++I)
      fprintf(outfile,"%d (%s)%s",tcscf_mos[I] + block_offset[tcscf_sym[I]],
                                 moinfo_scf->get_irr_labs(tcscf_sym[I]),
                                 I != nci - 1 ? "," : "");
    fprintf(outfile,"]");

    Favg      .allocate("Favg",nirreps,sopi,sopi);
    Favg_t    .allocate("Favg_t",nirreps,sopi,sopi);

    allocate1(double,ci,nci);
    allocate1(double,ci_grad,nci);
    allocate2(double,H_tcscf,nci,nci);
    for(int I = 0; I < nci; ++I){
      Dtc[I]    .allocate("Dtc[" + to_string(I) + "]",nirreps,sopi,sopi);
      Dsum[I]   .allocate("Dsum[" + to_string(I) + "]",nirreps,sopi,sopi);
      Ftc[I]    .allocate("Ftc[" + to_string(I) + "]",nirreps,sopi,sopi);
      Ftc_t[I]  .allocate("Ftc_t[" + to_string(I) + "]",nirreps,sopi,sopi);
      ci[I] = 0.0;// (I == 0 ? 0.7071067811865475244 : -0.7071067811865475244) ;
    }
  }
}

void SCF::cleanup()
{
  release1(block_offset);
  release1(pairpi);
  release1(pair_offset);
  release1(pairs);
  release2(pair);
  release2(pair_sym);

  if(reference == tcscf){
    release1(ci);
    release1(ci_grad);
    release2(H_tcscf);
  }

  release1(PK);
  release1(K);
}

void SCF::transform(SBlockMatrix& Initial, SBlockMatrix& Final, SBlockMatrix& Transformation){
  T.multiply(false,false,Initial,Transformation);
  Final.multiply(true,false,Transformation,T);
}

}} /* End Namespaces */

/*





  maxdiis = options_get_int("MAXDIIS");
  nuclear_energy = moinfo_scf->get_nuclear_energy();



  
  rohf  = false;
  
    rohf  = true;

  tcscf = false;
  if(options_get_str("CORR_REFERENCE")=="TCSCF"){
    tcscf = true;
    
    int count = 0;
    for(int h = 0; h < nirreps; ++h){
      for(int n = 0; n < actv[h]; ++n){
        tcscf_sym[count] = h;
        tcscf_mos[count] = docc[h] + n;
        count++;
      }
    }
    if(count > 2){
      fprintf(outfile,"\n  Too many electrons in corr_actv (%d)",count*2);
      fflush(outfile);
      exit(1);
    }else{
      fprintf(outfile,"\n  TCSCF MOs = %d (%s), %d (%s)",
                           tcscf_mos[0],moinfo_scf->get_irr_labs(tcscf_sym[0]),
                           tcscf_mos[1],moinfo_scf->get_irr_labs(tcscf_sym[1]));
    }
  }

  if(rohf)
    fprintf(outfile,"\n  Running a ROHF computation");


 

  // Compute the so symmetry
  allocate1(int,so_sym,nso);
  int k=0;
  for(int h=0; h < nirreps; ++h){
    for(int i=0; i < sopi[h]; ++i){
      so_sym[k] = h;
      ++k;
    }
  }

  allocate1(double,epsilon,nso);

  H = NULL;
  S = NULL;
  S_sqrt_inv = NULL;
  C    = new SBlockMatrix("C",nirreps,sopi,sopi);
  C_t  = new SBlockMatrix("C transformed",nirreps,sopi,sopi);
  D    = new SBlockMatrix("D",nirreps,sopi,sopi);
  Dc   = new SBlockMatrix("Dc",nirreps,sopi,sopi);
  Do   = new SBlockMatrix("Do",nirreps,sopi,sopi);
  F    = new SBlockMatrix("Fock",nirreps,sopi,sopi);
  F_t  = new SBlockMatrix("Fock transformed",nirreps,sopi,sopi);
  Fc   = new SBlockMatrix("Fock closed-shell",nirreps,sopi,sopi);
  Fo   = new SBlockMatrix("Fock open-shell",nirreps,sopi,sopi);
  Fc_t = new SBlockMatrix("Fock closed-shell transformed",nirreps,sopi,sopi);
  Fo_t = new SBlockMatrix("Fock open-shell   transformed",nirreps,sopi,sopi);
  G    = new SBlockMatrix("G",nirreps,sopi,sopi);
  Temp = new SBlockMatrix("Temp",nirreps,sopi,sopi);
  SDF  = new SBlockMatrix("SDF",nirreps,sopi,sopi);
  FDS  = new SBlockMatrix("FDS",nirreps,sopi,sopi);
  CSC  = new SBlockMatrix("CSC",nirreps,sopi,sopi);
  C_T  = new SBlockMatrix("C Transposed",nirreps,sopi,sopi);
  S_sqrt  = new SBlockMatrix("S^1/2",nirreps,sopi,sopi);
  Feff_t  = new SBlockMatrix("effective Fock transformed",nirreps,sopi,sopi);
  for(int n = 0; n < maxdiis; ++n){
    diis_F[n]  = new SBlockMatrix("diis_F",nirreps,sopi,sopi);
    diis_e[n]  = new SBlockMatrix("diis_e",nirreps,sopi,sopi);
  }
  if(tcscf){
    Dtc[0] = new SBlockMatrix("D TCSCF 0",nirreps,sopi,sopi);
    Dtc[1] = new SBlockMatrix("D TCSCF 1",nirreps,sopi,sopi);
    Ex     = new SBlockMatrix("Exchange",nirreps,sopi,sopi);
    Ftc[0] = new SBlockMatrix("Fock TCSCF 0",nirreps,sopi,sopi);
    Ftc[1] = new SBlockMatrix("Fock TCSCF 1",nirreps,sopi,sopi);
    Ftc_t[0] = new SBlockMatrix("Fock TCSCF 0 transformed",nirreps,sopi,sopi);
    Ftc_t[1] = new SBlockMatrix("Fock TCSCF 1 transformed",nirreps,sopi,sopi);

  }
}



void SCF::read_C(){
  for(int h =0; h < nirreps; ++h){
    double** c = C->get_block(h);
    double** psi_C = moinfo_scf->get_scf_mos(h);
    // Set all blocks to Fc
    for(int i = 0; i < sopi[h]; ++i){
      for(int j = 0; j < sopi[h]; ++j){             
        c[i][j] = psi_C[i][j];
      }
    }
  }
  
  C->print();
}



void SCF::cleanup()
{
#if STORE_TEI
  for(int h=0; h < nirreps; ++h){
    if(pairpi[h] > 0){
      release1(tei_so[h]);
    }
  }
  release1(tei_so);
#endif



  if(CSC!=NULL)
    delete CSC;
  if(C!=NULL)
    delete C;
  if(C_t!=NULL)
    delete C_t;
  if(C_T!=NULL)
    delete C_T;
  if(D!=NULL)
    delete D;
  if(Dc!=NULL)
    delete Dc;
  if(Do!=NULL)
    delete Do;
  if(tcscf){
    delete Dtc[0];
    delete Dtc[1];
    delete Ex;
    delete Ftc[0];
    delete Ftc[1];
    delete Ftc_t[0];
    delete Ftc_t[1];

  }
  if(F!=NULL)
    delete F;
  if(F_t!=NULL)
    delete F_t;
  if(Fc!=NULL)
    delete Fc;
  if(Fo!=NULL)
    delete Fo;
  if(Fc_t!=NULL)
    delete Fc_t;
  if(Fo_t!=NULL)
    delete Fo_t;
  if(Feff_t!=NULL)
    delete Feff_t;
  if(G!=NULL)
    delete G;
  if(H!=NULL)
    delete H;
  if(S!=NULL)
    delete S;
  if(S_sqrt_inv!=NULL)
    delete S_sqrt_inv;
  if(S_sqrt!=NULL)
    delete S_sqrt;
  if(Temp!=NULL)
    delete Temp;
  if(SDF!=NULL)
    delete SDF;
  if(FDS!=NULL)
    delete FDS;
  for(int n = 0; n < maxdiis; ++n){
    delete diis_F[n];
    delete diis_e[n];
  }

}
*/


