c Copy real arrays. a1 is copied into a2
      subroutine copyr(a1,a2,n)
#if defined(sun)
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
      integer n
      real*8 a1(n), a2(n)

      integer i

      do 101 i = 1, n
        a2(i) = a1(i)
 101    continue

      return
      end
c---------------------------------------------------------------------
c Copy integer arrays. a1 is copied into a2
      subroutine copyi(a1,a2,n)
#if defined(sun)
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
      integer n
      integer a1(n), a2(n)

      integer i

      do 101 i = 1, n
        a2(i) = a1(i)
 101    continue

      return
      end

