      subroutine cochec(c)
#if defined(UNDEFINED)
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "share.in"
#if defined(sun)
      integer c(c(MAXCOR))
#else
      integer c(1)
#endif
      integer locat

c the starting location
      locat = c(COSTRT)

 801  continue
      if (     (c(locat) .ne. stfree)
     &    .and.(c(locat) .ne. stres)) then
        write(6,*) 'cochec: the storage pool is corrupted,'
        write(6,*) '       check for overwrites'
        call codump(c)
        call qabort
        endif

      locat = c(locat+2)
c See if we are back where we started.
      if (c(COSTRT) .eq. locat) then
        return
        endif
      goto 801
      end
