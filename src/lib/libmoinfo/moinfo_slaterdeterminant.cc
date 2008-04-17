#include <iostream>
#include "moinfo.h"

extern FILE *infile, *outfile;

using namespace std;

MOInfo::SlaterDeterminant::SlaterDeterminant()
{
}


MOInfo::SlaterDeterminant::~SlaterDeterminant()
{
}


/**
 * @fn MOInfo::SlaterDeterminant::is_closed_shell()
 */
bool MOInfo::SlaterDeterminant::is_closed_shell()
{
  for(int i=0;i<moinfo->get_nmo();i++)
    if(bits[i]!=bits[i+moinfo->get_nmo()])
      return(false);
  return(true);
}


/**
 * @fn MOInfo::SlaterDeterminant::is_spin_flip(SlaterDeterminant& det)
 */
bool MOInfo::SlaterDeterminant::is_spin_flipped(SlaterDeterminant& det)
{
  for(int i=0;i<moinfo->get_nmo();i++){
    if(bits[i]!=det.test(i+moinfo->get_nmo()))
      return(false);
    if(bits[i+moinfo->get_nmo()]!=det.test(i))
      return(false);
  }
  return(true);
}


/**
 * @fn MOInfo::SlaterDeterminant::print(int n)
 */
void MOInfo::SlaterDeterminant::print()
{
  fprintf(outfile,"|");
  for(int i=0;i<moinfo->get_nmo();i++){
    fprintf(outfile,"%s",get_occupation_symbol(i));
  }
  fprintf(outfile,">");
}

/**
 * @fn MOInfo::SlaterDeterminant::print(int n)
 */
void MOInfo::SlaterDeterminant::print_occ()
{
  fprintf(outfile,"|");
  int counter = 0;
  for(int h=0;h<moinfo->get_nirreps();h++){
    fprintf(outfile,"[");
    for(int i=0;i<moinfo->get_docc(h);i++){
      fprintf(outfile,"%c",get_occupation_symbol(counter));
      counter++;
    }
    for(int i=0;i<moinfo->get_actv(h);i++){
      fprintf(outfile,"%c",get_occupation_symbol(counter));
      counter++;
    }
    counter+=moinfo->get_avir(h);
    fprintf(outfile,"]");
  }
  fprintf(outfile,">");
}

/**
 * @fn MOInfo::SlaterDeterminant::get_internal_excitations(...)
 */
void MOInfo::SlaterDeterminant::get_internal_excitations(SlaterDeterminant& det,double& sign,
                                                 vector<pair<int,int> >& alpha_operators,
                                                 vector<pair<int,int> >& beta_operators)
{
  int ann, cre;
  int nmo = moinfo->get_nmo();
  bitdet bits_exc = det.get_bits();
  bitdet bits_tmp = bits;
  sign = 1.0;
  // Find one set of excitations at a time
  ann = -1; cre = -1;
  while(cre<nmo){
    while(++ann<nmo)
      if(bits[ann] && !bits_exc[ann]) break;
    while(++cre<nmo)
      if(!bits[cre] && bits_exc[cre]) break;
    if(cre<nmo){
      alpha_operators.push_back(make_pair(moinfo->get_all_to_occ(ann),moinfo->get_all_to_vir(cre)));
      sign *= annihilate(bits_tmp,ann);
      sign *=     create(bits_tmp,cre);
    }
  }
  ann = -1; cre = -1;
  while(cre<nmo){
    while(++ann<nmo)
      if(bits[ann+nmo] && !bits_exc[ann+nmo]) break;
    while(++cre<nmo)
      if(!bits[cre+nmo] && bits_exc[cre+nmo]) break;
    if(cre<nmo){
      beta_operators.push_back(make_pair(moinfo->get_all_to_occ(ann),moinfo->get_all_to_vir(cre)));
      sign *= annihilate(bits_tmp,ann+nmo);
      sign *=     create(bits_tmp,cre+nmo);
    }
  }
}

double MOInfo::SlaterDeterminant::annihilate(bitdet& bits_det,int so)
{
  if(bits_det.test(so)){
    bits_det.flip(so);
    double sign = 1.0;
    for(int i=0;i<so;i++)
      if(bits_det.test(i)) sign *= -1.0;
    return(sign);
  }
  else
    return(0.0);
}

double MOInfo::SlaterDeterminant::create(bitdet& bits_det,int so)
{
  if(!bits_det.test(so)){
    bits_det.flip(so);
    double sign = 1.0;
    for(int i=0;i<so;i++)
      if(bits_det.test(i)) sign *= -1.0;
    return(sign);
  }
  else
    return(0.0);
}

vector<int> MOInfo::SlaterDeterminant::get_aocc()
{
  vector<int> aocc;
  for(int i=0;i<moinfo->get_nmo();i++)
    if(bits[i])
      aocc.push_back(moinfo->get_all_to_occ(i));
  return(aocc);
}

vector<int> MOInfo::SlaterDeterminant::get_bocc()
{
  vector<int> bocc;
  for(int i=0;i<moinfo->get_nmo();i++)
    if(bits[i+moinfo->get_nmo()])
      bocc.push_back(moinfo->get_all_to_occ(i));
  return(bocc);
}

vector<int> MOInfo::SlaterDeterminant::get_avir()
{
  vector<int> avir;
  for(int i=0;i<moinfo->get_nmo();i++)
    if(!bits[i])
      avir.push_back(moinfo->get_all_to_vir(i));
  return(avir);
}

vector<int> MOInfo::SlaterDeterminant::get_bvir()
{
  vector<int> bvir;
  for(int i=0;i<moinfo->get_nmo();i++)
    if(!bits[i+moinfo->get_nmo()])
      bvir.push_back(moinfo->get_all_to_vir(i));
  return(bvir);
}

char MOInfo::SlaterDeterminant::get_occupation_symbol(int i)
{
  int occupation=bits[i]+bits[i+moinfo->get_nmo()];
  if(occupation==0)                     return('0');
  if(occupation==2)                     return('2');
  if((occupation==1) && bits.test(i))   return('+');
  if((occupation==1) && bits.test(i+moinfo->get_nmo())) return('-');
  return(' ');
}
