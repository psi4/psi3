c This reallocates storage and copies the data (all that exists or all that
c will fit, whichever is appropiate)
      subroutine corllc(c,pntr,size,type)
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
      integer pntr, size
      character*(*) type

      integer csize,osize,npntr,coallc,isizof

c Get the old size of the array
      osize = c(pntr-3)/isizof(type)
      npntr = coallc(c,size*isizof(type))

c The number of items to be copied
      csize = min0(osize,size)

c Copy the data
      if ((type .eq. 'REAL').or.(type.eq.'real')) then
        call copyr(c(pntr),c(npntr),csize)
      else if ((type .eq. 'INT').or.(type.eq.'int')) then
        call copyi(c(pntr),c(npntr),csize)
        endif

      call cofree(c,pntr)
      pntr = npntr
      return
      end
