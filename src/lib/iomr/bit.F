
      integer function iand(i,j)
      implicit integer*4 (a-z)
      iand=and(i,j)
      return
      end

      integer function ishft(i,j)
      implicit integer*4 (a-z)
      if (j .le. 0) then
          ishft=rshift(i,-j)
      else
          ishft=lshift(i, j)
          endif
      return
      end

      integer function ior(i,j)
      implicit integer*4 (a-z)
      ior = or(i,j)
      return
      end

      integer function ieor(i,j)
      implicit integer*4 (a-z)
      ieor = xor(i,j)
      return
      end

#ifdef APOLLO
      integer function rshift(i,j)
      implicit integer*4 (a-z)
      rshift = rshft(i,j)
      return
      end

      integer function lshift(i,j)
      implicit integer*4 (a-z)
      lshift = lshft(i,j)
      return
      end
#endif
