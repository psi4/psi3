/***************************************************************************
 *  SQ
 *  Copyright (C) 2007 by Francesco Evangelista (frank@ccc.uga.edu)
 *  A second quantization code
 ***************************************************************************/

#include "squlli.h"
#include <libciomr/libciomr.h>

extern FILE *infile, *outfile;

namespace psi{ namespace psimrcc{

using namespace std;

int  SQULL::maxbits = 0;
bool SQULL::initialized = false;
unsigned long long SQULL::flags[128];
unsigned long long SQULL::flods[128];

void binary(unsigned long long number);

SQULL::SQULL():bits(0)
{
  if(!initialized)
    init();
}

SQULL::~SQULL()
{
}

void SQULL::init()
{
//   fprintf(outfile, "\nInitializing the SQULL class");
  maxbits = sizeof(ull)*8;
//   fprintf(outfile, "\nmaxbits = %d",maxbits);
  ull one  = 1;
  ull zero = 0;
  flags[0]=one;
  flods[0]=zero;
  for(int i=0;i<maxbits;i++)
    flags[i] = one << i;
  for(int i=1;i<maxbits;i++){
    flods[i]=0;
    for(int j=0;j<i;j++)
      flods[i] = flods[i] | flags[j];
  }
  
  for(int i=0;i<maxbits;i++){
//     fprintf(outfile, "\nflags[%2d] = ",i);
//     binary(flags[i]);
  }
  
  for(int i=0;i<maxbits;i++){
//     fprintf(outfile, "\nflods[%2d] = ",i);
//     binary(flods[i]);
  }
  initialized = true;
}

/*
void binary(unsigned long long number)
{
  int remainder;

  if(number <= 1){
    fprintf(outfile, "%d",number);
    return;
  }
  remainder = number%2;
  binary(number >> 1);
  fprintf(outfile, "%d",remainder);
}*/

void SQULL::binary(ull number)
{
  for (int i=0; i<maxbits;i++) {
    int bit = ((number >> i) & 1);
    fprintf(outfile, "%d",bit);
  }
}



void SQULL::print()
{
  for (int i=0; i<maxbits;i++) {
    int bit = ((bits >> i) & 1);
    fprintf(outfile, "%d",bit);
  }
//   //fprintf(outfile,"The binary for %20qu is \n",n);
//   for(int i=0;i<maxbits;i++)
//     fprintf(outfile,"%c",test(i) ? '1' : '0');
}




// char SQULL::kthBit(bits i, int k){
//   if (isKthBitOn(i, k))
//     return '1';
//   else
//     return '0';
// }

/*
int SQULL::isKthBitOn(bits n, int k){  // k = 1-32
  bits one=1;
  bits i = n & (one << (k-1));
  return (i != 0);
}





void SQULL::print_bits_nice(bits n)
{
  bits i;
   //fprintf(outfile,"The binary for %20qu is \n",n);
  fprintf(outfile,"\n");
  for(i=1;i<=nmo;i++)
    fprintf(outfile,"%c",kthBit(n,i));
  fprintf(outfile,"\n");
  fprintf(outfile,"  number of 1's = %d\n",bitcount(n));
}



inline int SQULL::bitcount(bits n)
{
     return bits_in_16bits [n         & 0xffffu]
         +  bits_in_16bits [(n >> 16) & 0xffffu]
         +  bits_in_16bits [(n >> 32) & 0xffffu]
         +  bits_in_16bits [(n >> 48) & 0xffffu];
}*/



}} /* End Namespaces */
