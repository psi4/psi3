c----------------------------------------------------------------------------
c  The formatted read integer triangular array of dimension 2 routine.
      integer function frdit2(keywrd,result,n1,ntri)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer n1,ntri
      character*(*) keywrd
      integer result(ntri)

      integer i,j,ij,tmp

      integer frdpos,frdlpn,frdint,frdrpn

      frdit2 = frdpos(keywrd)
      if (frdit2 .ne. EOK) return

      frdit2 = frdlpn()
      if (frdit2 .ne. EOK) return

      ij = 0
      do 101 i = 1, n1
        frdit2 = frdlpn()
        if (frdit2 .ne. EOK) return
        do 201 j = 1, i
          ij = ij + 1
          frdit2 = frdint(tmp)
          if (frdit2.ne.EOK) return
          result(ij) = tmp
 201      continue
        frdit2 = frdrpn()
        if (frdit2 .ne. EOK) return
 101    continue

      frdit2 = frdrpn()
      end
