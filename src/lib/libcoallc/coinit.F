      subroutine coinit(c,mc)
#if defined(sun)
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
      integer mc
      integer c(mc)
#include "share.in"

c Implementing a circular list of free space for core allocation.

c The initial status is free
      c(NEQUIV+1) = stfree
c The initial size is mc - NEQUIV - 4
      c(NEQUIV+2) = mc - NEQUIV - 4
c The initial forward link is here.
      c(NEQUIV+3) = NEQUIV + 1
c The initial backward link is here.
      c(NEQUIV+4) = NEQUIV + 1

c This is where we start looking for memory:
      c(COSTRT) = NEQUIV + 1

c Initialize the maxcor in the c array:
      c(MAXCOR) = mc

      return
      end

