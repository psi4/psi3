c_*BEGIN_FILE psi.F
      FUNCTION PSI(X,Y,Z,orbtyp,orbnum)
      IMPLICIT REAL*8 (A-H,O-Z)
      common/coord/chgs(112),GEOM(3,112),ict(112),ns(112),ks(112)
      common/funcs/nat,nbfao,msfru
      common/gauss/ex(200),CSPD(5,200)
      common/nshel/knuc(200),ktype(200),kstart(200),kloc(200),
     1 kmin(200),kmax(200),kng(200),nsh,npb
      COMMON/MO/C(200),DC(50,200)
      character*4 orbtyp
      integer orbnum
c
      PI=3.1415926535D0
      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0
      HALF=0.5D0
      PD=ZERO
      inum1=orbnum
      inum2=orbnum
      if (orbtyp.eq.'DENS') inum1=1
      do 99 jp=inum1,inum2
      do 98 i=1,nbfao
 98   c(i)=dc(jp,i)
      NSJ=0
      NSJJ=0
      P=ZERO
      DO 100 I=1,NAT
      XI=X-GEOM(1,I)
      YI=Y-GEOM(2,I)
      ZI=Z-GEOM(3,I)
      RR=XI*XI+YI*YI+ZI*ZI
      NSI=NS(I)
      DO 110 J=1,NSI
      NSJ=NSJ+1
      KTY=KTYPE(NSJ)
      KST=KSTART(NSJ)
      KNGI=KST+KNG(NSJ)-1
c     write (6,*) 'nsj,kty,kst,kngi',nsj,kty,kst,kngi
      GO TO (310,320,330,340,350) KTY
c
c     ...... s orbital
310   NSJJ=NSJJ+1
      PP=ZERO
      DO 120 K=KST,KNGI
120   PP=PP+CSPD(1,K)*DEXP(-EX(K)*RR)
      P=P+C(NSJJ)*PP
c     write (6,*) 'kty,c(nsjj),pp,p',kty,c(nsjj),pp,p
      GO TO 110
c
c     ...... p orbital
320   NSJJ=NSJJ+1
      PP=ZERO
      DO 130 K=KST,KNGI
      CSEX=CSPD(2,K)*DEXP(-EX(K)*RR)
130   PP=PP+CSEX*XI
      P=P+C(NSJJ)*PP
c     write (6,*) 'kty,c(nsjj),pp,p',kty,c(nsjj),pp,p
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 140 K=KST,KNGI
      CSEX=CSPD(2,K)*DEXP(-EX(K)*RR)
140   PP=PP+CSEX*YI
      P=P+C(NSJJ)*PP
c     write (6,*) 'kty,c(nsjj),pp,p',kty,c(nsjj),pp,p
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 150 K=KST,KNGI
      CSEX=CSPD(2,K)*DEXP(-EX(K)*RR)
150   PP=PP+CSEX*ZI
      P=P+C(NSJJ)*PP
c     write (6,*) 'kty,c(nsjj),pp,p',kty,c(nsjj),pp,p
      GO TO 110
c
c     ...... d orbital
330   NSJJ=NSJJ+1
      PP=ZERO
      DO 160 K=KST,KNGI
      CSEX=CSPD(3,K)*DEXP(-EX(K)*RR)
160   PP=PP+CSEX*XI*XI
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 170 K=KST,KNGI
      CSEX=CSPD(3,K)*DEXP(-EX(K)*RR)
170   PP=PP+CSEX*YI*YI
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 180 K=KST,KNGI
      CSEX=CSPD(3,K)*DEXP(-EX(K)*RR)
180   PP=PP+CSEX*ZI*ZI
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 190 K=KST,KNGI
      CSEX=CSPD(3,K)*DEXP(-EX(K)*RR)
190   PP=PP+CSEX*XI*YI*1.732051D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 200 K=KST,KNGI
      CSEX=CSPD(3,K)*DEXP(-EX(K)*RR)
200   PP=PP+CSEX*XI*ZI*1.732051D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 210 K=KST,KNGI
      CSEX=CSPD(3,K)*DEXP(-EX(K)*RR)
210   PP=PP+CSEX*YI*ZI*1.732051D0
      P=P+C(NSJJ)*PP
      GOTO 110
c
c     ...... f orbital
340   NSJJ=NSJJ+1
      PP=ZERO
      DO 221 K=KST,KNGI
      CSEX=CSPD(4,K)*DEXP(-EX(K)*RR)
221   PP=PP+CSEX*XI*XI*XI
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 222 K=KST,KNGI
      CSEX=CSPD(4,K)*DEXP(-EX(K)*RR)
222   PP=PP+CSEX*YI*YI*YI
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 223 K=KST,KNGI
      CSEX=CSPD(4,K)*DEXP(-EX(K)*RR)
223   PP=PP+CSEX*ZI*ZI*ZI
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 224 K=KST,KNGI
      CSEX=CSPD(4,K)*DEXP(-EX(K)*RR)
224   PP=PP+CSEX*XI*XI*YI*2.236068D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 225 K=KST,KNGI
      CSEX=CSPD(4,K)*DEXP(-EX(K)*RR)
225   PP=PP+CSEX*XI*XI*ZI*2.236068D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 226 K=KST,KNGI
      CSEX=CSPD(4,K)*DEXP(-EX(K)*RR)
226   PP=PP+CSEX*YI*YI*XI*2.236068D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 227 K=KST,KNGI
      CSEX=CSPD(4,K)*DEXP(-EX(K)*RR)
227   PP=PP+CSEX*YI*YI*ZI*2.236068D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 228 K=KST,KNGI
      CSEX=CSPD(4,K)*DEXP(-EX(K)*RR)
228   PP=PP+CSEX*ZI*ZI*XI*2.236068D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 229 K=KST,KNGI
      CSEX=CSPD(4,K)*DEXP(-EX(K)*RR)
229   PP=PP+CSEX*ZI*ZI*YI*2.236068D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 230 K=KST,KNGI
      CSEX=CSPD(4,K)*DEXP(-EX(K)*RR)
230   PP=PP+CSEX*XI*YI*ZI*3.8729833D0
      P=P+C(NSJJ)*PP
      GOTO 110
c
c     ...... g orbital
350   NSJJ=NSJJ+1
      PP=ZERO
      DO 231 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
231   PP=PP+CSEX*XI*XI*XI*XI
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 232 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
232   PP=PP+CSEX*YI*YI*YI*YI
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 233 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
233   PP=PP+CSEX*ZI*ZI*ZI*ZI
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 234 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
234   PP=PP+CSEX*XI*XI*XI*YI*2.6457513D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 235 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
235   PP=PP+CSEX*XI*XI*XI*ZI*2.6457513D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 236 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
236   PP=PP+CSEX*YI*YI*YI*XI*2.6457513D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 237 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
237   PP=PP+CSEX*YI*YI*YI*ZI*2.6457513D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 238 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
238   PP=PP+CSEX*ZI*ZI*ZI*XI*2.6457513D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 239 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
239   PP=PP+CSEX*ZI*ZI*ZI*YI*2.6457513D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 240 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
240   PP=PP+CSEX*XI*XI*YI*YI*3.4156503D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 241 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
241   PP=PP+CSEX*XI*XI*ZI*ZI*3.4156503D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 242 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
242   PP=PP+CSEX*YI*YI*ZI*ZI*3.4156503D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 243 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
243   PP=PP+CSEX*XI*XI*YI*ZI*5.9160798D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 244 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
244   PP=PP+CSEX*YI*YI*ZI*XI*5.9160798D0
      P=P+C(NSJJ)*PP
c
      NSJJ=NSJJ+1
      PP=ZERO
      DO 245 K=KST,KNGI
      CSEX=CSPD(5,K)*DEXP(-EX(K)*RR)
245   PP=PP+CSEX*ZI*ZI*YI*XI*5.9160798D0
      P=P+C(NSJJ)*PP
C  
110   CONTINUE
100   CONTINUE
      if (orbtyp.eq.'DENS') pd=pd+p*p
99    continue
      if (orbtyp.eq.'DENS') then
      psi=pd
      else
      PSI=P
      endif
      RETURN
C
      END
