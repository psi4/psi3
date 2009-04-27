/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include <iostream>
#include <algorithm>

#include <libmoinfo/libmoinfo.h>
#include <libutil/libutil.h>

#include "index.h"

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

int                       CCIndex::nirreps=-1;

CCIndex::CCIndex(std::string str):
nelements(0),ntuples(0), greater_than_or_equal(false), greater_than(false), memory(0.0), label(str),tuples(0),
one_index_to_tuple(0),one_index_to_irrep(0),
two_index_to_tuple(0),two_index_to_irrep(0),
three_index_to_tuple(0),three_index_to_irrep(0)
{
  if(nirreps<0) nirreps = moinfo->get_nirreps();
  init();
}

CCIndex::~CCIndex()
{
  cleanup();
}

void CCIndex::init()
{
  // New orbital spaces must be added here
  for(int i =0;i<label.size();i++)
    if( label[i]=='o' || label[i]=='a' || label[i]=='v' || label[i]=='s' || label[i]=='n' || label[i]=='f')
      nelements++;

  // Get the orbital spaces data pointers
  for(int i=0;i<label.size();i++)
  {
    if(label[i]=='o'){
      mospi.push_back(moinfo->get_occ());
      indices_to_pitzer.push_back(moinfo->get_occ_to_mo());
    }else if(label[i]=='v'){
      mospi.push_back(moinfo->get_vir());
      indices_to_pitzer.push_back(moinfo->get_vir_to_mo());
    }else if(label[i]=='a'){
      mospi.push_back(moinfo->get_actv());
      indices_to_pitzer.push_back(moinfo->get_actv_to_mo());
    }else if(label[i]=='f'){
      mospi.push_back(moinfo->get_fvir());
      indices_to_pitzer.push_back(moinfo->get_fvir_to_mo());
    }else if(label[i]=='s'){
      mospi.push_back(moinfo->get_sopi());
    }else if(label[i]=='n'){
      mospi.push_back(moinfo->get_mopi());
    }
  }
  for(int i=0;i<nelements;i++){
    first_mos.push_back(vector<int>(nirreps,0));
    dimension.push_back(0);
  }
  for(int i=0;i<nelements;i++){
    for(int h=0;h<nirreps;h++){
      first_mos[i][h]=dimension[i];
      dimension[i]+=mospi[i][h];
    }
  }

  switch(nelements){
    case 0:
      make_zero_index();
      break;
    case 1:
      make_one_index();
      break;
    case 2:
      make_two_index();
      break;
    case 3:
      make_three_index();
      break;
    default:{
      fprintf(outfile,"\n\n\tThe CCIndex class cannot handle %s because there are more than three indices!!!\n\n",label.c_str());
      fflush(outfile);
      exit(1);
    }
  }
  _memory_manager_->add_allocated_memory(memory);
}

void CCIndex::cleanup()
{
  if(tuples!=0)
    release2(tuples);
//    free_smatrix(tuples,ntuples,dimension.size());
  if(one_index_to_tuple!=0)
    delete[] one_index_to_tuple;
  if(one_index_to_irrep!=0)
    delete[] one_index_to_irrep;
  if(two_index_to_tuple!=0)
    release2(two_index_to_tuple);
  if(two_index_to_irrep!=0)
    release2(two_index_to_irrep);
  if(three_index_to_tuple!=0)
    release3(three_index_to_tuple);
  if(three_index_to_irrep!=0)
    release3(three_index_to_irrep);
}

void CCIndex::make_zero_index()
{
  std::vector<std::vector<short> >  pairs;                  // The pairs ordered as a vector
  ntuples = 0;
  for(int h=0;h<nirreps;h++){
    first.push_back(ntuples);
    if(h==0){
      std::vector<short> pair;
      pairs.push_back(pair);
      ntuples++;
    }
    last.push_back(ntuples);
    pairpi.push_back(last[h]-first[h]);
  }
  // Allocate the memory for the tuples and store them
  allocate2(short,tuples,1,1);
//  memory+=(double)init_smatrix(tuples,1,1)/1048576.0;
  tuples[0][0] = 0;
}

void CCIndex::make_one_index()
{
  // The pairs ordered as a vector
  std::vector<std::vector<short> >  pairs;

  // Allocate the 1->tuple mapping array and set them to -1
  one_index_to_tuple = new size_t[dimension[0]];
  one_index_to_irrep = new int[dimension[0]];
  for(int i=0;i<dimension[0];i++){
    one_index_to_tuple[i] =  0;
    one_index_to_irrep[i] = -1;
  }
  memory+=(double)dimension[0]*2.0*sizeof(int)/1048576.0;

  ntuples = 0;
  for(int h=0;h<nirreps;h++){
    first.push_back(ntuples);
    for(int p=0;p<mospi[0][h];p++){
      one_index_to_tuple[ntuples]=p;
      one_index_to_irrep[ntuples]=h;
      std::vector<short> pair;
      pair.push_back(ntuples);
      pairs.push_back(pair);
      ntuples++;
    }
    last.push_back(ntuples);
    pairpi.push_back(last[h]-first[h]);
  }

  // Allocate the memory for the tuples and store them
//  memory+=(double)init_smatrix(tuples,ntuples,1)/1048576.0;
  allocate2(short,tuples,ntuples,1);
  for(int n=0;n<pairs.size();n++)
    tuples[n][0] = pairs[n][0];
}

void CCIndex::make_two_index()
{
  std::vector<std::vector<short> >  pairs;                  // The pairs ordered as a vector

  // Allocate the 2->tuple mapping array and set them to -1
  allocate2(size_t,two_index_to_tuple,dimension[0],dimension[1]);
  allocate2(int,two_index_to_irrep,dimension[0],dimension[1]);

  for(int i=0;i<dimension[0];i++){
    for(int j=0;j<dimension[1];j++){
      two_index_to_tuple[i][j] = 0;
      two_index_to_irrep[i][j] = -1;
    }
  }

  // [X>=Y]
  if(label.find(">=")!=string::npos){
    greater_than_or_equal = true;
    ntuples               = 0;
    for(int h=0;h<nirreps;h++){
      first.push_back(ntuples);
      for(int p_sym=0;p_sym<nirreps;p_sym++){
        int q_sym = h^p_sym;
        int p = first_mos[0][p_sym];
        for(int p_rel=0;p_rel<mospi[0][p_sym];p_rel++){
          int q = first_mos[1][q_sym];
          for(int q_rel=0;q_rel<mospi[1][q_sym];q_rel++){
            if(p>=q){
              two_index_to_tuple[p][q]=ntuples-first[h];
              two_index_to_irrep[p][q]=h;
              std::vector<short> pair;
              pair.push_back(p);
              pair.push_back(q);
              pairs.push_back(pair);
              ntuples++;
            }
            q++;
          }
          p++;
        }
      }
      last.push_back(ntuples);
      pairpi.push_back(last[h]-first[h]);
    }
  }else if(label.find(">")!=string::npos){
    greater_than          = true;
    ntuples               = 0;
    for(int h=0;h<nirreps;h++){
      first.push_back(ntuples);
      for(int p_sym=0;p_sym<nirreps;p_sym++){
        int q_sym = h^p_sym;
        int p = first_mos[0][p_sym];
        for(int p_rel=0;p_rel<mospi[0][p_sym];p_rel++){
          int q = first_mos[1][q_sym];
          for(int q_rel=0;q_rel<mospi[1][q_sym];q_rel++){
            if(p>q){
              two_index_to_tuple[p][q]=ntuples-first[h];
              two_index_to_irrep[p][q]=h;
              std::vector<short> pair;
              pair.push_back(p);
              pair.push_back(q);
              pairs.push_back(pair);
              ntuples++;
            }
            q++;
          }
          p++;
        }
      }
      last.push_back(ntuples);
      pairpi.push_back(last[h]-first[h]);
    }
  }else{
    ntuples               = 0;
    for(int h=0;h<nirreps;h++){
      first.push_back(ntuples);
      for(int p_sym=0;p_sym<nirreps;p_sym++){
        int q_sym = h^p_sym;
        int p = first_mos[0][p_sym];
        for(int p_rel=0;p_rel<mospi[0][p_sym];p_rel++){
          int q = first_mos[1][q_sym];
          for(int q_rel=0;q_rel<mospi[1][q_sym];q_rel++){
            two_index_to_tuple[p][q]=ntuples-first[h];
            two_index_to_irrep[p][q]=h;
            std::vector<short> pair;
            pair.push_back(p);
            pair.push_back(q);
            pairs.push_back(pair);
            ntuples++;
            q++;
          }
          p++;
        }
      }
      last.push_back(ntuples);
      pairpi.push_back(last[h]-first[h]);
    }
  }

  // Allocate the memory for the tuples and store them
//  memory+=(double)init_smatrix(tuples,ntuples,2)/1048576.0;
  allocate2(short,tuples,ntuples,2);
  for(int n=0;n<pairs.size();n++){
    tuples[n][0] = pairs[n][0];
    tuples[n][1] = pairs[n][1];
  }
}


void CCIndex::make_three_index()
{
  if(label.find(">")!=string::npos){
    fprintf(outfile,"\n\n\tThe CCIndex class cannot handle restricted loops for triplets!!!\n\n");
    fflush(outfile);
    exit(1);
  }

  std::vector<std::vector<short> >  pairs;                  // The pairs ordered as a vector

  // Allocate the 3->tuple mapping array and set them to -1
  allocate3(size_t,three_index_to_tuple,dimension[0],dimension[1],dimension[2]);
  allocate3(int,three_index_to_irrep,dimension[0],dimension[1],dimension[2]);
  for(int i=0;i<dimension[0];i++){
    for(int j=0;j<dimension[1];j++){
      for(int k=0;k<dimension[2];k++){
        three_index_to_tuple[i][j][k] =  0;
        three_index_to_irrep[i][j][k] = -1;
      }
    }
  }


  ntuples  = 0;
  for(int h=0;h<nirreps;h++){
    first.push_back(ntuples);
    for(int p_sym=0;p_sym<nirreps;p_sym++){
      for(int q_sym=0;q_sym<nirreps;q_sym++){
        int r_sym = h^p_sym^q_sym;
        int p = first_mos[0][p_sym];
        for(int p_rel=0;p_rel<mospi[0][p_sym];p_rel++){
          int q = first_mos[1][q_sym];
          for(int q_rel=0;q_rel<mospi[1][q_sym];q_rel++){
            int r = first_mos[2][r_sym];
            for(int r_rel=0;r_rel<mospi[2][r_sym];r_rel++){
              three_index_to_tuple[p][q][r]=ntuples-first[h];
              three_index_to_irrep[p][q][r]=h;
              std::vector<short> pair;
              pair.push_back(p);
              pair.push_back(q);
              pair.push_back(r);
              pairs.push_back(pair);
              ntuples++;
              r++;
            }
          q++;
          }
        p++;
        }
      }
    }
    last.push_back(ntuples);
    pairpi.push_back(last[h]-first[h]);
  }

  // Allocate the memory for the tuples and store them
//  memory+=(double)init_smatrix(tuples,ntuples,3)/1048576.0;
  allocate2(short,tuples,ntuples,3);
  for(int n=0;n<pairs.size();n++){
    tuples[n][0] = pairs[n][0];
    tuples[n][1] = pairs[n][1];
    tuples[n][2] = pairs[n][2];
  }
}

void CCIndex::print()
{
  fprintf(outfile,"\n\n---------------------------------");
  fprintf(outfile,"\n\tPair Type %s has %d elements",label.c_str(),ntuples);
  fprintf(outfile,"\n---------------------------------");
  int index=0;
  for(int h=0;h<nirreps;h++){
    if(pairpi[h]>0)
      fprintf(outfile,"\n\t%s",moinfo->get_irr_labs(h));
    for(int tuple=0;tuple<pairpi[h];tuple++){
      fprintf(outfile,"\n\t\t( ");
      for(int k=0;k<nelements;k++)
        fprintf(outfile,"%d ",tuples[index][k]);
      fprintf(outfile,")");
      index++;
    }
  }
  fprintf(outfile,"\n---------------------------------");
}

}} /* End Namespaces */
