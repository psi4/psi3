
      subroutine mxmtrn(a,la,b,lb,i1,i2,j1,j2)
      implicit real*8 (a-h,o-z)
      real*8 a(la,1),b(lb,1)

      ii = 1
      do 101 i = i1, i2
        jj = 1
        do 201 j = j1, j2
          b(jj,ii) = a(i,j)
          jj = jj + 1
  201     continue
        ii = ii + 1
  101   continue

      return
      end
