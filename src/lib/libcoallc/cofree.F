      subroutine cofree(c,data)
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
      integer data
      integer block, link, blink

c If data is zero coallc did not allocate anything (allocation size was
c zero or out of memory)
      if (data .eq. 0) then
        return
        endif

c Convert from the data address to the block address
      block = data - 4

      if (c(block) .ne. stres) then
        write(6,*) 'cofree: data not allocated (or overwritten)'
        call mabort
        endif

      c(block) = stfree

      link = c(block+2)
      blink = c(block+3)

c If the back link is free then join the blocks of memory together.
c (don't follow the list from bottom to top)
      if (     (c(blink) .eq. stfree)
     &    .and.(blink .lt. block)) then
        c(blink+1) = c(blink+1) + c(block+1) + 4
        c(blink+2) = c(block+2)
        block = blink
        endif

c If the forward link is free then join the block together
c (don't follow the list from top to bottom)
      if (     (c(link) .eq. stfree)
     &    .and.(link .gt. block)) then
        c(block+1) = c(block+1) + c(link+1) + 4
        c(block+2) = c(link+2)
        endif

c Set up the back link in the block which this block links to
      c(c(block+2)+3) = block

c At the next allocate put data at this freed area.
      c(COSTRT) = block
      return
      end
