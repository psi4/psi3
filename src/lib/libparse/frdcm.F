c----------------------------------------------------------------------------
c This is the master character datum reader. The first index in idx always
c corresponds to the number of the outermost vector, and so on.  The result
c vector in res always corresponds to the innermost vector.
      integer function frdcm(keywrd,nres,res,nidx,idx)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) keywrd
      integer nres
      character*(*) res(nres)
      integer nidx
      integer idx(nidx)

      integer i
      integer srrcod
      character*80 smess
      integer errcod
      character*(lentok) token
      character*80 emess
      common /lparer/ errcod, emess

      integer frdpos, frdchr, arpos

c The user may want to ignore errors, so make the default blank.
      do 101 i = 1, nres
        res(i) = ' '
  101   continue

c Position the input
      frdcm = frdpos(keywrd)
      if (frdcm .ne. EOK) return

c Position to the starting datum.
      frdcm = arpos(nidx,idx)
      if (frdcm .ne. EOK) return

c Initialize the return code to EOK
      srrcod = EOK
      smess = ' '

c Begin reading in the result vector.
      do 103 i = 1, nres
c       Get the next token (the result)
        frdcm = frdchr(token)
        if (frdcm .ne. EOK) return

c Copy the token to the result
        res(i) = token

c If the token was truncated set the appropriate error flag, but continue
c with the processing.
        if (res(i).ne.token) then
          srrcod = Etrunc
          smess = keywrd
          call scatb(smess,': result has been truncated (frdcm)')
          endif
  103   continue

      frdcm = srrcod
      emess = smess
      return
      end
