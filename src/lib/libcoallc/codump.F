      subroutine codump(c)
#if defined(sun)
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
      integer block, start, printc

      block = NEQUIV + 1
      start = block

      printc = 100

 801  continue

      printc = printc + 1
      if (printc .gt. 20) then
        printc = 0
        write(6,1) 'location','status','size','link','backlink'
        endif

      write(6,2) block, c(block),c(block+1),c(block+2),c(block+3)

   1  format(5(2x,a8,2x))
   2  format(5(2x,i9,1x))

      if ((c(block) .ne. stfree) .and. (c(block) .ne. stres)) then
        write(6,*) 'codump: allocate chain cannot be followed'
        return
        endif

      block = c(block+2)
      if (block .eq. start) return
      goto 801
      end
