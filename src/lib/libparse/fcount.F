c This routine counts the number of elements in the top level vector.
      integer function fcnt(keywrd,result)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) keywrd
      integer result

      integer idx(1)

      integer fcntm

      fcnt = fcntm(keywrd,result,0,idx)
      return
      end
c-----------------------------------------------------------------------------
c This routine counts the number of elements in a next to top level vector.
      integer function fcnt1(keywrd,result,i1)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) keywrd
      integer result, i1

      integer idx(1)

      integer fcntm

      idx(1) = i1
      fcnt1 = fcntm(keywrd,result,1,idx)
      return
      end
c-----------------------------------------------------------------------------
c This routine counts the number of elements in a next to next to
c top level vector.
      integer function fcnt2(keywrd,result,i1,i2)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) keywrd
      integer result, i1, i2

      integer idx(2)

      integer fcntm

      idx(1) = i1
      idx(2) = i2
      fcnt2 = fcntm(keywrd,result,2,idx)
      return
      end
c-----------------------------------------------------------------------------
c See fnct2
      integer function fcnt3(keywrd,result,i1,i2,i3)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) keywrd
      integer result, i1, i2, i3

      integer idx(3)

      integer fcntm

      idx(1) = i1
      idx(2) = i2
      idx(3) = i3
      fcnt3 = fcntm(keywrd,result,3,idx)
      return
      end
c-----------------------------------------------------------------------------
c See fnct2
      integer function fcnt4(keywrd,result,i1,i2,i3,i4)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) keywrd
      integer result, i1, i2, i3, i4

      integer idx(4)

      integer fcntm

      idx(1) = i1
      idx(2) = i2
      idx(3) = i3
      idx(4) = i4
      fcnt4 = fcntm(keywrd,result,4,idx)
      return
      end
c-----------------------------------------------------------------------------
c This is the master routine that counts the number of items in a vector.
      integer function fcntm(keywrd,result,nidx,idx)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer result,nidx
      character*(*) keywrd
      integer idx(nidx)

      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdpos, arpos, arcnt

c Position the input
      errcod = frdpos(keywrd)
      if (errcod .ne. EOK) then
        fcntm = errcod
        return
        endif

c Position to the starting datum.
      errcod = arpos(nidx,idx)
      if (errcod .ne. EOK) then
        fcntm = errcod
        return
        endif

c Count the number of elements
      fcntm = arcnt(result)
      return
      end
c-----------------------------------------------------------------------------
c This counts the number of items in a vector.  If the next token is not a '('
c then the datum is scalar and the result is -1.
      integer function arcnt(result)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer result

      integer errcod
      character*(lentok) token
      character*80 emess
      common /lparer/ errcod, emess

      integer gtoken, matchp

      errcod = gtoken(token)
      if (errcod .ne. EOK) then
        arcnt = errcod
        return
        endif

c Initialize various things
      errcod = EOK
      arcnt = errcod
      emess = ' '
      result = 0

      if (token .ne. '(') then
c       This is a scalar
        result = -1
        return
        endif

 801  continue
        errcod = gtoken(token)
        if (errcod .ne. EOK) then
          arcnt = errcod
          return
          endif

c See if we have reached the end of the array
        if (token .eq. ')') return

c If a new array starts here, then skip it.
        if (token .eq. '(') then
          errcod = matchp()
          if (errcod .ne. EOK) then
            arcnt = errcod
            return
            endif
          endif

        result = result + 1
        goto 801

      end
