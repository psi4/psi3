      integer function coallc(c,size)
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
      integer size, nsize
      integer ostrt

c Make nsize even so we always stay on a 64 bit boundary on 32 bit machines
      nsize = size
      if ((nsize/2)*2 .ne. nsize) nsize = nsize + 1

      if (nsize .eq. 0) then
        coallc = 0
        return
        endif

c Save the starting location
      ostrt = c(COSTRT)

c     find a bit of free space which is large enough
 801  continue
      if (c(c(COSTRT)) .eq. stfree) then
        if (c(c(COSTRT)+1) .ge. nsize) goto 802
      else if (c(c(COSTRT)) .ne. stres) then
        write(6,*) 'coallc: memory is corrupted,'
        write(6,*) '       check for overwrites'
        call codump(c)
        call mabort
        endif

      c(COSTRT) = c(c(COSTRT) + 2)
      if (c(COSTRT) .eq. ostrt) then
        write(6,*) 'coallc: out of memory'
        call jabort
        call qabort
        endif
      goto 801

c     space found, finish allocation
 802  continue

c Don't leave little chunks of memory around
      if (c(c(COSTRT)+1) - nsize .lt. 100) then
        c(c(COSTRT)) = stres

        coallc = c(COSTRT) + 4
        c(COSTRT) = c(c(COSTRT)+2)
      else
c       Set up the new free chunk
        c(c(COSTRT)+4+nsize) = stfree
        c(c(COSTRT)+4+nsize+1) = c(c(COSTRT)+1) - nsize - 4
        c(c(COSTRT)+4+nsize+2) = c(c(COSTRT)+2)
        c(c(COSTRT)+4+nsize+3) = c(COSTRT)

c       Set up the backlink of what the free chunk links to
        c(c(c(COSTRT)+4+nsize+2)+3) = c(COSTRT)+4+nsize

c       Set up the reserved block
        c(c(COSTRT)) = stres
        c(c(COSTRT)+1) = nsize
        c(c(COSTRT)+2) = c(COSTRT)+4+nsize

        coallc = c(COSTRT) + 4
        c(COSTRT) = c(COSTRT)+4+nsize
        endif

      return
      end
