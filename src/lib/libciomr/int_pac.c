#include "includes.h"

void int_pac(i,ib,j,jb,ipak,jpak)
   int *i,*ib,*j,*jb;
   unsigned int *ipak,*jpak;

   {
      *ipak = *ib;
      *ipak <<= 4;
      *ipak += *i;


      *jpak = *jb;
      *jpak <<= 4;
      *jpak += *j;
      }

void int_unpac(i,ib,j,jb,ipak,jpak)
   int *i,*ib,*j,*jb;
   unsigned int *ipak, *jpak;

   {
      *ib = *ipak >> 4;
      *i = *ipak & 15;

      *jb = *jpak >> 4;
      *j = *jpak & 15;
      }
