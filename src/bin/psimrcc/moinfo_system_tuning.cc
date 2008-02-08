/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <iostream>
#include <cmath>
#include "error.h"
#include "moinfo.h"
#include "utilities.h"
#include "squlli.h"

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

void sq_operations();
void psimrcc_blas();
void psimrcc_disk();

void MOInfo::tuning()
{
  fprintf(outfile,"\n\n  Tuning PSIMRCC to work best with the systems");
  fprintf(outfile,"\n  sizeof(short)  = %d bytes      sizeof(size_t) = %d bytes",sizeof(short),sizeof(size_t));
  fprintf(outfile,"\n  sizeof(int)    = %d bytes      sizeof(ull)    = %d bytes",sizeof(int),sizeof(unsigned long long));
  fprintf(outfile,"\n  sizeof(long)   = %d bytes      sizeof(double) = %d bytes",sizeof(long),sizeof(double));
  psimrcc_blas();
  psimrcc_disk();
  sq_operations();
}

void psimrcc_blas()
{
}

void psimrcc_disk()
{
}

void sq_operations()
{
  typedef std::bitset<32> bit32;
  typedef std::bitset<64> bit64;
  
  int ntests = 10000;
  Timer bitset32;
  fprintf(outfile,"\n\n  Testing bitset<32> performance:");
  // Clear and assign one bit
  bit32 test32;
  for(int i=0;i<ntests;++i){
    test32.reset();
  }
  double time_clear = bitset32.get();
  fprintf(outfile," %6.3f ms",time_clear*1000.0);

  // Assign one bit
  for(int i=0;i<ntests;++i){
    test32.set(7);
    test32.set(2);
    test32.set(18);
    test32.set(4);
    test32.set(9);
    test32.reset(7);
    test32.reset(2);
    test32.reset(18);
    test32.reset(4);
    test32.reset(9);
  }
  double time_assign = bitset32.get();
  fprintf(outfile," %6.3f ms",(time_assign-time_clear)*1000.0);

  Timer bitset64;
  fprintf(outfile,"\n  Testing bitset<64> performance:");
  // Clear and assign one bit
  bit64 test64;
  for(int i=0;i<ntests;++i){
    test64.reset();
  }
  time_clear = bitset64.get();
  fprintf(outfile," %6.3f ms",time_clear*1000.0);
  // Assign one bit
  for(int i=0;i<ntests;++i){
    test64.set(7);
    test64.set(2);
    test64.set(18);
    test64.set(4);
    test64.set(9);
    test64.reset(7);
    test64.reset(2);
    test64.reset(18);
    test64.reset(4);
    test64.reset(9);
  }
  time_assign = bitset64.get();
  fprintf(outfile," %6.3f ms",(time_assign-time_clear)*1000.0);


  SQULL testull;

  fprintf(outfile,"\n  Testing SQULL      performance:");
  // Clear and assign one bit
  Timer ull;
  for(int i=0;i<ntests;++i){
    testull.reset();
  }
  time_clear = ull.get();
  fprintf(outfile," %6.3f ms",time_clear*1000.0);
  // Assign one bit
  for(int i=0;i<ntests;++i){
    testull.set(7);
    testull.set(2);
    testull.set(18);
    testull.set(4);
    testull.set(9);
    testull.reset(7);
    testull.reset(2);
    testull.reset(18);
    testull.reset(4);
    testull.reset(9);
  }
  time_assign = ull.get();
  fprintf(outfile," %6.3f ms",(time_assign-time_clear)*1000.0);
}

}} /* End Namespaces */