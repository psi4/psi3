
      subroutine mxmtnc(a,ia,b,ib,c,ic,i1,i2,j1,j2,k1,k2)
      implicit real*8 (a-h, o-z)
      dimension a(1), b(1), c(1)
      integer aoff, aoffs, boff, boffs, boffst, coff, coffs
c c(i,j) = a(k,i)*b(k,j)

      aoffs = (i1-1)*ia + k1
      boffs = (j1-1)*ib + k1
      coffs = (j1-1)*ic + i1

      kinc = k2 - k1

      do 102 i = i1, i2
        boffst = boffs
        coff = coffs
        do 202 j = j1, j2
          t = c(coff)
          aoff = aoffs
          do 302 boff = boffst, boffst + kinc
            t      = t      + a(aoff)*b(boff)
            aoff = aoff + 1
  302       continue
          c(coff) = t
          boffst = boffst + ib
          coff = coff + ic
  202     continue
        aoffs = aoffs + ia
        coffs = coffs + 1
  102   continue
      end
