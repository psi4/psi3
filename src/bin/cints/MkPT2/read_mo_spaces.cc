#include <cstring>
#include <cstdio>
#include <cstdlib>
#include<libipv1/ip_lib.h>
#include<libciomr/libciomr.h>
#include<libchkpt/chkpt.h>

#include <algorithm>
#include <vector>
#include <string>

#include <mkpt2_ints.h>

extern FILE* outfile;

namespace psi{ namespace CINTS{ namespace mkpt2{


void correlate(char *ptgrp, int irrep, int& nirreps_old, int& nirreps_new,int*& arr)
{ /* This is a hack from input!!! (ACS) */
  int  i;

  if (strcmp(ptgrp,"C1 ") == 0)
    nirreps_old = 1;
  else if (strcmp(ptgrp,"Cs ") == 0)
    nirreps_old = 2;
  else if (strcmp(ptgrp,"Ci ") == 0)
    nirreps_old = 2;
  else if (strcmp(ptgrp,"C2 ") == 0)
    nirreps_old = 2;
  else if (strcmp(ptgrp,"C2v") == 0)
    nirreps_old = 4;
  else if (strcmp(ptgrp,"D2 ") == 0)
    nirreps_old = 4;
  else if (strcmp(ptgrp,"C2h") == 0)
    nirreps_old = 4;
  else if (strcmp(ptgrp,"D2h") == 0)
    nirreps_old = 8;
  else {
    fprintf(outfile,"point group %s unknown.\n",ptgrp);
    exit(1);
  }

  arr = new int[nirreps_old];

  if (irrep == 0) { /* return identity */
    nirreps_new = nirreps_old;
    for (i=0; i<nirreps_old; ++i)
      arr[i] = i;
    return;
  }

  nirreps_new = nirreps_old / 2;
  if ((strcmp(ptgrp,"C1 ") == 0) || (strcmp(ptgrp,"Cs ") == 0) ||
      (strcmp(ptgrp,"Ci ") == 0) || (strcmp(ptgrp,"C2 ") == 0) ) {
        arr[0] = 0; arr[1] = 0;
  }
  else if ( (strcmp(ptgrp,"C2v") == 0) || (strcmp(ptgrp,"D2 ") == 0) ||
            (strcmp(ptgrp,"C2h") == 0) ) {
    if (irrep == 1) {
      arr[0] = 0;  arr[1] = 0; arr[2] = 1;  arr[3] = 1;
    }
    else if (irrep == 2) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 0;  arr[3] = 1;
    }
    else if (irrep == 3) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 1;  arr[3] = 0;
    }
  }
  else if (strcmp(ptgrp,"D2h") == 0) {
    /* 1,2,3 give C2h displaced geometries */
    if (irrep == 1) {
      arr[0] = 0;  arr[1] = 0; arr[2] = 1;  arr[3] = 1;
      arr[4] = 2;  arr[5] = 2; arr[6] = 3;  arr[7] = 3;
    }
    else if (irrep == 2) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 0;  arr[3] = 1;
      arr[4] = 2;  arr[5] = 3; arr[6] = 2;  arr[7] = 3;
    }
    else if (irrep == 3) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 1;  arr[3] = 0;
      arr[4] = 2;  arr[5] = 3; arr[6] = 3;  arr[7] = 2;
    }
    /* 4 gives D2 displaced geometries */
    else if (irrep == 4) { /* D2 */
      arr[0] = 0;  arr[1] = 1; arr[2] = 2;  arr[3] = 3;
      arr[4] = 0;  arr[5] = 1; arr[6] = 2;  arr[7] = 3;
    }
    /* displacements along irreps 5,6,7 make C2v structures */
    /* care is taken to make sure definition of b1 and b2 will
       match those that input will generate - the following seems to work:
       b1u disp: has C2(z), b2 irrep symmetric wrt sigma(yz)
       b2u disp: has C2(y), b2 irrep symmetric wrt sigma(xy)
       b3u disp: has C2(x), b2 irrep symmetric wrt sigma(xz) */
    else if (irrep == 5) { /* b1u */
      arr[0] = 0;  arr[1] = 1; arr[2] = 2;  arr[3] = 3;
      arr[4] = 1;  arr[5] = 0; arr[6] = 3;  arr[7] = 2;
    }
    else if (irrep == 6) { /* b2u */
      arr[0] = 0;  arr[1] = 3; arr[2] = 1;  arr[3] = 2;
      arr[4] = 1;  arr[5] = 2; arr[6] = 0;  arr[7] = 3;
    }
    else if (irrep == 7) { /* b3u */
      arr[0] = 0;  arr[1] = 2; arr[2] = 3;  arr[3] = 1;
      arr[4] = 1;  arr[5] = 3; arr[6] = 2;  arr[7] = 0;
    }
  }
  else {
    fprintf(outfile,"Point group unknown for correlation table.\n");
  }

  return;
}

void read_mo_space(int nirreps_ref, int& n, std::vector<int>& mo, std::string labels)
{
  bool read = false;

  std::vector<std::string> label_vec = split(labels);
  for(int k = 0; k < label_vec.size(); ++k){
    int size;
    int stat = ip_count(const_cast<char *>(label_vec[k].c_str()),&size,0); // Cast to avoid warnings
    if(stat == IPE_OK){
      // Defaults is to set all to zero
      mo.assign(nirreps_ref,0);
      n = 0;
      if(read){
        fprintf(outfile,"\n\n  libmoinfo has found a redundancy in the input keywords %s , please fix it!",labels.c_str());
        fflush(outfile);
        exit(1);
      }else{
        read = true;
      }
      if(size==nirreps_ref){
        for(int i=0;i<size;i++){
          stat=ip_data(const_cast<char *>(label_vec[k].c_str()),const_cast<char *>("%d"),(&(mo)[i]),1,i); // Cast to avoid warnings
          n += mo[i];
        }
      }else{
        fprintf(outfile,"\n\n  The size of the %s array (%d) does not match the number of irreps (%d), please fix the input file",label_vec[k].c_str(),size,nirreps_ref);
        fflush(outfile);
        exit(1);
      }
    }
  }
}

std::vector<int> read_chkpt_intvec(int n, int* array)
{
  // Read an array from chkpt and save as a STL vector
  std::vector<int> stl_vector(&array[0],&array[n]);
  free(array);
  return(stl_vector);
}

std::vector<std::string> split(const std::string& str)
{
  // Split a string
  typedef std::string::const_iterator iter;
  std::vector<std::string> splitted_string;
  iter i = str.begin();
  while(i != str.end()){
    // Ignore leading blanks
    i = std::find_if(i,str.end(), not_space);
    // Find the end of next word
    iter j = std::find_if(i,str.end(),space);
    // Copy the characters in [i,j)
    if(i!=str.end())
      splitted_string.push_back(std::string(i,j));
    i = j;
  }
  return(splitted_string);
}

bool space(char c)
{
  return isspace(c);
}

bool not_space(char c)
{
  return !isspace(c);
}
//void read_mo_space(int nirreps_ref,int& n, int* mo, const char* label)
//{
//  int size;
//  int stat=ip_count(label,&size,0);
//  n=0;
//  if(stat == IPE_OK){
//    if(size==nirreps_ref){
//      for(int i=0;i<size;i++){
//        stat=ip_data(label,"%d",(&(mo)[i]),1,i);
//        n+=mo[i];
//      }
//    }else{
//      fprintf(outfile,"\n\nThe size of the %s array (%d) does not match the number of irreps (%d), please fix the input file",label,size,nirreps_ref);
//      fflush(outfile);
//      exit(1);
//    }
//  }
//}

}}} /* End namespaces */
