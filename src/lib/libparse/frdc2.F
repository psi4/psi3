c----------------------------------------------------------------------------
c This reads a character datum from the file unit and places it in result.
c An error code is returned. The datum is the i2th element of a vector found
c in the i1th position of the top level vector
      integer function frdc2(keywrd,result,i1,i2)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer i1,i2
      character*(*) keywrd, result

      integer idx(2)
      character*(lentok) res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdcm

      idx(1) = i1
      idx(2) = i2

      frdc2 = frdcm(keywrd,1,res,2,idx)
      if (frdc2 .ne. EOK) return
      result = res(1)

      if (result .ne. res(1)) then
        errcod = Etrunc
        frdc2 = errcod
        emess = 'Character data for '''
        call scat(emess,keywrd)
        call scat(emess,''' truncated (frdc2)')
        endif

      return
      end
