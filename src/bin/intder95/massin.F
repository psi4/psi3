C//////////////////////////////////////////////////////////////////////
      SUBROUTINE MASSIN (XMASS,NA,IFLAG,n,maskey,flag)
      IMPLICIT REAL*8 (A-H,O-Z)
#include <libparse/error.h>
      CHARACTER A*12,B*72,junk*80
      character maslbl*12,maskey*8
      LOGICAL P,Q,flag
      integer nmass, getmas, errcod, frdpos, frdlpn, fcnt, frdchr
      integer matchp
      integer n, i, j
      DIMENSION XMASS(NA)
      PARAMETER(ZERO=0.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
   1  FORMAT(A72)
   2  FORMAT(F12.6)
      IFLAG=0
      if (flag) then
      MN=(NA-1)/6+1
      MI=-5
      DO 1000 MM=1,MN
      MI=MI+6
      MF=MI+5
      IF(MF.GT.NA) MF=NA
      READ(IIN1,1) B
      NI=-11
      DO 500 M=MI,MF
      NI=NI+12
      NF=NI+11
      A(1:12)=B(NI:NF)
      DO 100 I=1,12
         P=LGE(A(I:I),'A')
         Q=LLE(A(I:I),'Z')
         IF(P.AND.Q) then
           iflag = getmas(na,m,a,i,xmass)
           goto 500
         endif
 100  CONTINUE
      READ(A,2) XMASS(M)
      IF(XMASS(M).LE.ZERO) IFLAG=1
 500  CONTINUE
 1000 CONTINUE
      RETURN

c do PSI 2.0 input
      else 
        nmass = 0
        errcod = frdpos(maskey)
        if (maskey.eq.'ISOTOPES') then
          errcod = frdlpn()
          errcod = frdlpn()
        elseif (maskey.eq.'MASSES  ') then
          n=1
          errcod = frdlpn()
        else 
          write(iout,*) 'error: keyword ', maskey, ' unrecognized'
        endif

        do 90 i=1,n-1
          errcod = matchp()
          errcod = frdlpn()
 90     continue
        do 2000 m=1,na
          xmass(m)=-1.0D0
          errcod = frdchr(maslbl)
          if(errcod.ne.EOK) goto 9001          

c         convert to upper case and left justify the mass label
          call leftju(maslbl) 
          call touppr(maslbl)
          p=lge(maslbl(1:1),'A')
          q=lle(maslbl(1:1),'Z')
          if(p.and.q) then
            iflag = getmas(na,m,maslbl,1,xmass)
          else
            call chrflt(maslbl,xmass(m))
            if (xmass(m).le.ZERO) iflag=1
          endif
 2000   continue
      endif
        
      return

 9001 write(IOUT,*) 'error parsing MASS vector'
      call pabort
      call qabort

 9002 write(IOUT,*) 'MASS vector has wrong number of elements'
      call pabort
      call qabort

      END
