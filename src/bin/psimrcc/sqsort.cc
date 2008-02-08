/***************************************************************************
 *  SQ
 *  Copyright (C) 2007 by Francesco Evangelista (frank@ccc.uga.edu)
 *  A second quantization code
 ***************************************************************************/

#include "sqsort.h"
#include "moinfo.h"
#include "utilities.h"
#include "debugging.h"
#include "memory_manager.h"

#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.h>
#include <psifiles.h>

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

SQSort::SQSort()
{
  init();
}

SQSort::~SQSort()
{
  cleanup();
}

void SQSort::init()
{
  frozen_core_energy = 0.0;
  ioff = moinfo->get_ioff();
  // Allocate memory for the one electron integrals
  allocate2(double,f,moinfo->get_nmo(),moinfo->get_nmo());
  allocate2(double,h,moinfo->get_nmo(),moinfo->get_nmo());
  // Allocate memory for the two electron integrals
  int ntei = four(moinfo->get_nmo()-1,moinfo->get_nmo()-1,moinfo->get_nmo()-1,moinfo->get_nmo()-1)+1;
  allocate1(double,tei,ntei);
  // Read the integrals
  read_one_electron_integrals();
  read_two_electron_integrals();
  fprintf(outfile,"\n  SQSort: frozen core energy = %20.12f", frozen_core_energy);
}

void SQSort::cleanup()
{
  // Allocate memory for the one electron integrals
  release2(f);
  release2(h);
  // Allocate memory for the two electron integrals
  int ntei = four(moinfo->get_norbs()-1,moinfo->get_norbs()-1,moinfo->get_norbs()-1,moinfo->get_norbs()-1)+1;
  release1(tei);
}


void SQSort::read_one_electron_integrals()
{
  // One electron integrals
  int ntri = ioff[moinfo->get_norbs()];
  double* htri;
  double** h_pitzer;
  allocate1(double,htri,ntri);
  allocate2(double,h_pitzer,moinfo->get_norbs(),moinfo->get_norbs());

  // Read all the one electorn integrals from psi assuming Pitzer ordering
  iwl_rdone(PSIF_OEI,PSIF_MO_OEI,htri,ntri,0,0,outfile);

  // Save the one electron integrals sorted by energy
  for(int p = 0;p<moinfo->get_norbs();p++)
    for(int q=p;q<moinfo->get_norbs();q++){
      h_pitzer[p][q]=htri[p+ioff[q]];
      h_pitzer[q][p]=h_pitzer[p][q];
      int p_orbs = moinfo->get_all_to_nonfrozen(p);
      int q_orbs = moinfo->get_all_to_nonfrozen(q);
      if((p_orbs>=0) && (q_orbs>=0)){
        h[p_orbs][q_orbs] = h_pitzer[p][q];
        h[q_orbs][p_orbs] = h_pitzer[p][q];
      }
    }

  // Compute the frozen core energy contribution of the one electron integrals
  frozen_core_energy_oei_contribution(h_pitzer);

  release1(htri);
  release2(h_pitzer);
  fprintf(outfile,"\n  SQSort: read %d one electron integrals", moinfo->get_norbs()*moinfo->get_norbs());
}

void SQSort::read_two_electron_integrals()
{
  // Read all the (frozen + non-frozen) TEI in Pitzer order
  // and store them in a in-core block-matrix

  // Allocate memory for the two electron integrals
  double* tei_pitzer;
  int ntei_pitzer = four(moinfo->get_norbs()-1,moinfo->get_norbs()-1,moinfo->get_norbs()-1,moinfo->get_norbs()-1)+1;
  allocate1(double,tei_pitzer,ntei_pitzer);

  int elements = 0;
  int ilsti;
  struct iwlbuf ERIIN;
  iwl_buf_init(&ERIIN,PSIF_MO_TEI,0.0,1,1);
    do{
      ilsti     = ERIIN.lastbuf;
      int nbuf  = ERIIN.inbuf;
      int fi    = 0;
      for(int index = 0; index < nbuf; index++){
        int p  = abs(ERIIN.labels[fi]);
        int q  = ERIIN.labels[fi+1];
        int r  = ERIIN.labels[fi+2];
        int s  = ERIIN.labels[fi+3];
        tei_pitzer[four(p,q,r,s)] = ERIIN.values[index];

        // Eliminate the frozen core and frozen virtual orbitals
        int p_orbs = moinfo->get_all_to_nonfrozen(p);
        int q_orbs = moinfo->get_all_to_nonfrozen(q);
        int r_orbs = moinfo->get_all_to_nonfrozen(r);
        int s_orbs = moinfo->get_all_to_nonfrozen(s);
        if((p_orbs>=0) && (q_orbs>=0) && (r_orbs>=0) && (s_orbs>=0))
          tei[four(p_orbs,q_orbs,r_orbs,s_orbs)] = ERIIN.values[index];

        fi += 4;
        elements++;
      }
      if(!ilsti)
        iwl_buf_fetch(&ERIIN);
    } while(!ilsti);
  fflush(outfile);
  iwl_buf_close(&ERIIN,1);

  frozen_core_energy_tei_contribution(tei_pitzer);

  release1(tei_pitzer);
}



void SQSort::frozen_core_energy_oei_contribution(double**& h_pitzer)
{
  // One electron contribution to the frozen core energy from each irrep
  for(int h=0;h<moinfo->get_nirreps();h++){
    for(int i=0;i<moinfo->get_focc(h);i++){
      int ii = i + moinfo->get_first_orbs_pitzer(h);
      frozen_core_energy += 2.0 * h_pitzer[ii][ii];
    }
  }
}

void SQSort::frozen_core_energy_tei_contribution(double*& tei_pitzer)
{
  // Two electron contribution to the frozen core energy
  for(int hi=0;hi<moinfo->get_nirreps();hi++){
    for(int i=0;i<moinfo->get_focc(hi);i++){
      for(int hj=0;hj<moinfo->get_nirreps();hj++){
        for(int j=0;j<moinfo->get_focc(hj);j++){
          int ii = i + moinfo->get_first_orbs_pitzer(hi);
          int jj = j + moinfo->get_first_orbs_pitzer(hj);
          frozen_core_energy += 2.0 * tei_pitzer[four(ii,ii,jj,jj)];
          frozen_core_energy -=       tei_pitzer[four(ii,jj,ii,jj)];
        }
      }
    }
  }
  //Modify the active part of H to include the core effects;
  for(int p = 0;p<moinfo->get_nmo();p++){
    for(int q = 0;q<moinfo->get_nmo();q++){
      for(int hi=0;hi<moinfo->get_nirreps();hi++){
        for(int i=0;i<moinfo->get_focc(hi);i++){
          int ii = i + moinfo->get_first_orbs_pitzer(hi);
          int pp = moinfo->get_nonfrozen_to_all(p);
          int qq = moinfo->get_nonfrozen_to_all(q);
          if((pp>=0) && (qq>=0))
            h[p][q] += 2.0 * tei_pitzer[four(ii,ii,pp,qq)]-tei_pitzer[four(ii,pp,ii,qq)];
        }      
      }
    }
  }
  DEBUGGING(4,
    print_dmatrix106(h,moinfo->get_nmo(),moinfo->get_nmo(),outfile,"H matrix (frozen core)");
  );
}


#if 0




void Integrals::init(MOInfo& moinfo)
{
    long int l,lnbfso;

    energyNuclear = moinfo.energyNuclear;

    nbfso   = moinfo.nbfso;
    nfocc   = moinfo.nfocc;
    nfvir   = moinfo.nfvir;
    nmo     = nbfso-nfocc-nfvir;
    ndocc   = moinfo.ndocc;
    nelp    = moinfo.nel/2;
    moSort  = moinfo.moSort;
    if(nfocc+nfvir>0)
        freeze = true;
 
    lnbfso   = moinfo.nbfso;   
    ioff = new long int[(1+lnbfso*(lnbfso+1)/2+lnbfso)];
    for(l=0;l<(1+lnbfso*(lnbfso+1)/2+lnbfso);l++)
        ioff[l]=l*(l+1)/2;
  
    // Two electron integrals
    ntot        = four(lnbfso-1,lnbfso-1,lnbfso-1,lnbfso-1)+1;
    ntot_store  = four(nmo+nfvir-1,nmo+nfvir-1,nmo+nfvir-1,nmo+nfvir-1)+1;
    ntot_active = four(nmo+nfvir-1,nmo+nfvir-1,nmo+nfvir-1,nmo+nfvir-1)+1;
    h           = initMatrix(nbfso,nbfso);
    h1          = initMatrix(nmo,nmo);
    g           = initMatrix(nmo,nmo);
    h1vec       = new double[INDEX(nmo,nmo)+1];
    gvec        = new double[INDEX(nmo,nmo)+1];
    tei         = new double[ntot];
    for(int i=0;i<ntot;i++) tei[i]=0.0;
    if(freeze){
    tei_active  = new double[ntot_store];
    }
    fockD       = new double[nbfso];

}

void Integrals::free()
{
    delete[] ioff;
    delete[] tei;
    delete[] h1vec;
    delete[] gvec;
    delete[] fockD;
    freeMatrix(h1,nmo,nmo);
    freeMatrix(g,nmo,nmo);
    if(freeze){
        freeMatrix(h,nmo,nmo);
    }else{
        freeMatrix(h,nbfso,nbfso);
    }
}

void Integrals::readPsi3(MOInfo& moinfo)
{
    int i,j,k,l;
    int ii,jj,kk,ll;
    int ilsti,nbuf,fi,ntei;
    int index,elements;
    long int index4;
    struct iwlbuf ERIIN;
    double value;

    // One electron integrals
    ntri        = ioff[nbfso];
    htri        = new double[ntri];
    iwl_rdone(PSIF_OEI,PSIF_MO_OEI,htri,ntri,0,0,outfile);    

    // Save the one electron integrals sorted by energy
    for(i=0;i<nbfso;i++)
        for(j=i;j<nbfso;j++){
        ii=moSort[i];
        jj=moSort[j];
        h[ii][jj]=htri[ioff[j]+i];
        h[jj][ii]=h[ii][jj];
    }
    // Free htri array
    delete[] htri;

     if(debug>5){
        cout << "\n\n\tH matrix\n\n";
        printMatrix(h,nbfso,nbfso);
     }

    // Two electron integrals
    elements=0;
    iwl_buf_init(&ERIIN,PSIF_MO_TEI,0.0,1,1);
    do{
        ilsti = ERIIN.lastbuf;
        nbuf = ERIIN.inbuf;
        fi=0;
        for(index=0;index<nbuf;index++){
            i=ERIIN.labels[fi];
            i=abs(i);
            j=ERIIN.labels[fi+1];
            k=ERIIN.labels[fi+2];
            l=ERIIN.labels[fi+3];
            ii=moSort[i];
            jj=moSort[j];
            kk=moSort[k];
            ll=moSort[l];
            value=ERIIN.values[index];
            
            index4=four(ii,jj,kk,ll); //Store as (ij|kl)
//              cout << "\n\t(" << i << " " << j<< "|" << k<< " " << l << ") = "<< value;
            tei[index4]=value;
            if(freeze)
                if(((ii>=nfocc)&&(jj>=nfocc))&&((kk>=nfocc)&&(ll>=nfocc))){
                index4=four(ii-nfocc,jj-nfocc,kk-nfocc,ll-nfocc);  //Store active part of (ij|kl)
                tei_active[index4]=value;
                }
                fi+=4;
                elements++;
        }
        if(!ilsti)
            iwl_buf_fetch(&ERIIN);
    } while(!ilsti);
    iwl_buf_close(&ERIIN,1);
    ntei=elements;
    
    cout << "\n\n\tElectron Integrals\n";
    cout << "\n\tRead\t" <<  ntri << "\tone electron integrals";
    cout << "\n\tRead\t" <<  ntei << "\ttwo electron integrals" << endl;
    cout << "\n\tStored\t" <<  ntot << "\ttwo electron integrals" << endl;
    
    
    // Compute the Fock matrix diagonal elements
    if(debug>4) cout << "\n\tFock matrix diagonal elements [For Closed Shell Reference]" << endl;
    for(i=0;i<nbfso;i++){
        fockD[i]=h[i][i];
        for(j=0;j<(nelp);j++)
            fockD[i]+=2.0*tei[four(i,i,j,j)]-tei[four(i,j,i,j)];
        if(debug>4) cout << "\n\tf[" << i << "] = " << fockD[i]; 
    }

    /*// Compute the Fock matrix diagonal elements
    if(debug>4) cout << "\n\tFock Core matrix diagonal elements [For Closed Shell Reference]" << endl;
    for(i=0;i<nbfso;i++){
        fockCore[i]=0.0;
        for(j=0;j<nfocc;j++)
            fockCore[i]+=2.0*tei[four(i,i,j,j)]-tei[four(i,j,i,j)];
        if(debug>4) cout << "\n\tf[" << i << "] = " << fockCore[i]; 
}*/

    for(i=0;i<(nfocc+ndocc);i++){
        energyZero+=2.0*fockD[i];
    }
    // SCF Energy computed from MO integrals
    for(i=0;i<(nelp);i++){
        energySCF+=2.0*h[i][i];
        for(j=0;j<(nelp);j++){
            energySCF+=2.0*tei[four(i,i,j,j)] - tei[four(i,j,i,j)];
        }
    }
    energySCF+=energyNuclear;

    cout << "\n\tSCF Energy computed from MO integrals = " << energySCF <<endl;
   
    if(freeze)
    freezeCoreVir();

    //Compute h' matrix
    for(int p=0;p<nmo;p++){
        for(int q=0;q<nmo;q++){
            h1[p][q]=h[p][q];
            for(i=0;i<nmo;i++){
                h1[p][q]-=0.5*tei[four(p,i,i,q)];
            }
            h1vec[INDEX(p,q)]=h1[p][q];
        }
    }
    
    //Compute the g matrix [See Olsen J. Chem. Phys. 2185, 89 (4) 1988]
    for(int k=0;k<nmo;k++){
        for(int l=0;l<nmo;l++){
            g[k][l]=h[k][l];
            for(j=0;j<k;j++){
                g[k][l]-=tei[four(k,j,j,l)];
            }
            if(k==l)
                g[k][l]-=0.5*tei[four(k,k,k,l)];
            if(k>l)
                g[k][l]-=tei[four(k,k,k,l)];
            gvec[INDEX(k,l)]=g[k][l];
        }
    }

    if(debug>0){
        //Print the h' matrix
        cout << "\n\n\tH' Matrix";
        for(int p=0;p<nmo;p++){
            cout << "\n";
            for(int q=0;q<nmo;q++){
                cout << "\t" << h1[p][q];
            }
        }
    }

    if((moinfo.multp==1)&&(moinfo.sigmaMs0))
    {
        double value;
        for(int i=0;i<INDEX(nmo-1,nmo-1)+1;i++){
            value=0.5*tei[INDEX(i,i)];
            tei[INDEX(i,i)]=value;
        }
    }
}

void Integrals::freezeCoreVir()
{ 
    int i,j,p,q;
    double  **tempH, *temp;

    tempH           = initMatrix(nbfso,nbfso);
  
    coreEnergy=0.0;
    for(i=0;i<nfocc;i++){
        coreEnergy+=2*h[i][i];
        for(j=0;j<nfocc;j++)
            coreEnergy+=2.0*tei[four(i,i,j,j)]-tei[four(i,j,i,j)];
    }
    cout << "\n\tFrozen-core energy =\t" << coreEnergy << endl;
  
    //Modify the active part of H to include the core effects;
    for(p=nfocc;p<nbfso;p++){
        for(q=nfocc;q<nbfso;q++){
            tempH[p][q]=h[p][q];
            for(i=0;i<nfocc;i++){
                tempH[p][q]+=2.0*tei[four(i,i,p,q)]-tei[four(i,p,i,q)];
            }
        }
    }
    
    freeMatrix(h,nbfso,nbfso);
    h = initMatrix(nmo,nmo);
    //Save the active part of H;
    for(p=0;p<nmo;p++){
        for(q=0;q<nmo;q++){
            h[p][q]=tempH[p+nfocc][q+nfocc];
        }
    }
    
    if(debug>5){
        cout << "\n\n\tH matrix (frozen core)\n\n";
        printMatrix(h,nmo,nmo);
    }

    freeMatrix(tempH,nbfso,nbfso);
     

    //Trow away the useless (core part) tei;
    delete[] tei;
    tei  = new double[ntot_active];
    for(i=0;i<ntot_active;i++){
        tei[i]=tei_active[i];  
    }
    delete[] tei_active;
    ntot=ntot_active;

    //Trow away the useless (core part) fockD;
    coreFock = 0.0;
    for(p=0;p<nfocc;p++)
        coreFock += 2.0*fockD[p];
    if(debug>4) cout << "\n\tFrozen-core Fock Energy = " << coreFock << endl;

    temp = new double[nmo];
    for(p=0;p<nmo;p++)
        temp[p] = fockD[p+nfocc];
    delete[] fockD;
    fockD = new double[nmo];
    if(debug>4) cout << "\n\tFock matrix diagonal elements (frozen-core)" << endl;
    for(p=0;p<nmo;p++){
        fockD[p] = temp[p];    
        if(debug>4) cout << "\n\tf[" << p << "] = " << fockD[p]; 
    }
    delete[] temp; 
}
#endif





}} /* End Namespaces */