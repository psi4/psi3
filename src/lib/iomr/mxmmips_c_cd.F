
      block data mxmini
      integer ivec,ivec1,kvec,kvec1
      real*8 ascr
      common /mxmprm/ ivec,ivec1,kvec,kvec1
c     The dim of ascr is ivec*kvec
      common /mxmscr/ ascr(6144)
      data ivec/128/
      data ivec1/127/
      data kvec/48/
      data kvec1/47/
      end

      subroutine mxmnn(a,ia, b,ib, c,ic, l ,m ,n )
      implicit real*8 (a-h, o-z)
      dimension a(ia,1), b(ib,1), c(ic,1)
      common/mxmprm/ ivec,ivec1,kvec,kvec1
      common/mxmscr/ ascr(6144)
      do 100 ii = 1, l, ivec
          do 200 kk = 1, m, kvec
            i2 = min0(ii+ivec1,l)
            k2 = min0(kk+kvec1,m)
            call mxmtrn(a,ia,ascr,kvec,ii,i2,kk,k2)
      call mxmtnd(ascr,kvec,b,ib,c,ic,ii,i2,1,n,kk,k2)
  200          continue
  100      continue
      return
      end

      subroutine mxmtn(a,ia, b,ib, c,ic, l ,m ,n )
      implicit real*8 (a-h, o-z)
      dimension a(ia,1), b(ib,1), c(ic,1)
      common/mxmprm/ ivec,ivec1,kvec,kvec1
      do 100 ii = 1, l, ivec
          do 200 kk = 1, m, kvec
      call mxmtnc(a,ia,b,ib,c,ic,ii,min0(ii+ivec1,l),1,n,kk,
     &            min0(kk+kvec1,m))
  200          continue
  100      continue
      return
      end

      subroutine mxmnt(a,ia, b,ib, c,ic, l ,m ,n )
      implicit real*8 (a-h, o-z)
      dimension a(ia,1), b(ib,1), c(ic,1)
      common/mxmprm/ ivec,ivec1,kvec,kvec1
      common/mxmscr/ ascr(6144)
      do 100 ii = 1, l, ivec
          do 200 kk = 1, m, kvec
            i2 = min0(ii+ivec1,l)
            k2 = min0(kk+kvec1,m)
            call mxmtrn(a,ia,ascr,kvec,ii,i2,kk,k2)
      call mxmttd(ascr,kvec,b,ib,c,ic,ii,i2,1,n,kk,k2)
  200          continue
  100      continue
      return
      end

      subroutine mxmtt(a,ia, b,ib, c,ic, l ,m ,n )
      implicit real*8 (a-h, o-z)
      dimension a(ia,1), b(ib,1), c(ic,1)
      common/mxmprm/ ivec,ivec1,kvec,kvec1
      do 100 ii = 1, l, ivec
          do 200 kk = 1, m, kvec
      call mxmttc(a,ia,b,ib,c,ic,ii,min0(ii+ivec1,l),1,n,kk,
     &            min0(kk+kvec1,m))
  200          continue
  100      continue
      return
      end
