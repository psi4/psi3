
      subroutine mxmttc(a,ia,b,ib,c,ic,i1,i2,j1,j2,k1,k2)
      implicit real*8 (a-h, o-z)
      dimension a(1), b(1), c(1)
      integer aoff, aoffs, boff, boffs, boffst, coff, coffs
      integer boffl, boffls
c c(i,j) = a(k,i)*b(j,k)

      aoffs = (i1-1)*ia + k1
      boffs = (k1-1)*ib + j1
      coffs = (j1-1)*ic + i1
      boffls = (k2-1)*ib + j1

      kinc = k2 - k1

      do 102 i = i1, i2
        boffst = boffs
        boffl  = boffls
        coff = coffs
        do 202 j = j1, j2
          t = c(coff)
          aoff = aoffs
          do 302 boff = boffst, boffl, ib
            t      = t      + a(aoff)*b(boff)
            aoff = aoff + 1
  302       continue
          c(coff) = t
          boffst = boffst + 1
          boffl  = boffl + 1
          coff = coff + ic
  202     continue
        aoffs = aoffs + ia
        coffs = coffs + 1
  102   continue
      end
