c-----------------------------------------------------------------------------
c The master integer datum reader.
      integer function frdfm(keywrd,nres,res,nidx,idx)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) keywrd
      integer nres
      real*8 res(nres)
      integer nidx
      integer idx(nidx)

      integer i
      integer errcod
      real*8 tmp
      character*80 emess
      common /lparer/ errcod, emess

      integer frdpos, frdflt, arpos

c The user may want to ignore errors, so make the default zero.
      do 101 i = 1, nres
        res(i) = 0.0
  101   continue

c Position the input
      frdfm = frdpos(keywrd)
      if (frdfm .ne. EOK) return

c Position to the starting datum.
      frdfm = arpos(nidx,idx)
      if (frdfm .ne. EOK) return

c Begin reading in the result vector.
      do 103 i = 1, nres
c       Get the next token (the result)
        frdfm = frdflt(tmp)
        if (frdfm .ne. EOK) return

c Copy the token to the result
        res(i) = tmp
  103   continue

      frdfm = EOK
      emess = ' '
      return
      end
