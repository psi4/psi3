c_*BEGIN_FILE psi.F
      FUNCTION PSI(X,Y,Z,orbtyp,orbnum)
      IMPLICIT REAL*8 (A-H,O-Z)
      common/funcs/nat,num
      COMMON/MO/AA(500),DC(200,500)
      dimension bas0(200)
      character*4 orbtyp
      integer orbnum
c
      call phi(x,y,z,bas0)
      inum1=orbnum
      inum2=orbnum
      if (orbtyp.eq.'DENS') inum1=1
c
      pd=0
      do 99 jp=inum1,inum2
c
      p=0
      do 97 jq=1,num
      p=p+dc(jp,jq)*bas0(jq)
 97   continue
c
      if (orbtyp.eq.'DENS') pd=pd+p*p
99    continue
c
       psi=p
       if (orbtyp.eq.'DENS') psi=pd
c
      RETURN
      END

c
c
      SUBROUTINE PHI(X,Y,Z,BAS0)
      IMPLICIT REAL*8 (A-H,O-Z)
c
      COMMON/coord/ZAN(112),GEOM(3,112)
      common/funcs/nat,num
      common/gauss/ex(500),CSPD(500,5)
      common/nshel/katom(500),ktype(500),kstart(500),kloc(500),
     1 kmin(500),kmax(500),kng(500),nshell,npb
      COMMON/MO/C(500),DC(200,500)
      dimension BAS0(200)
c
      zero=0.0d0
      IBAS=0
      DO 100 I=1,NSHELL
      KTY=KTYPE(I)
      KST=KSTART(I)
      KNGI=KST+KNG(I)-1
      KAT=KATOM(I)
      if(i.ne.1.and.kat.eq.katom(i-1)) goto 300
      XI=X-GEOM(1,KAT)
      YI=Y-GEOM(2,KAT)
      ZI=Z-GEOM(3,KAT)
      RR=XI*XI+YI*YI+ZI*ZI
300   GO TO (310,320,330,340,350) KTY
c
c     ...... s orbital
310   IBAS=IBAS+1
      PP=ZERO
      DO 120 K=KST,KNGI
120   PP=PP+CSPD(K,1)*DEXP(-EX(K)*RR)
      BAS0(IBAS)=PP
      GO TO 100
c
c     ...... p orbital
320   IBAS=IBAS+1
      PP=ZERO
      DO 130 K=KST,KNGI
130   PP=PP+CSPD(K,2)*DEXP(-EX(K)*RR)
      BAS0(IBAS)=PP*XI
      CONTINUE
      IBAS=IBAS+1
      BAS0(IBAS)=PP*YI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*ZI
      GO TO 100
c
c     ...... d orbital
330   IBAS=IBAS+1
      PP=ZERO
      DO 160 K=KST,KNGI
160   PP=PP+CSPD(K,3)*DEXP(-EX(K)*RR)
      BAS0(IBAS)=PP*XI*XI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*YI*YI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*ZI*ZI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*XI*YI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*XI*ZI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*YI*ZI
      GOTO 100
c
c     ...... f orbital
340   IBAS=IBAS+1
      PP=ZERO
      DO 221 K=KST,KNGI
221   PP=PP+CSPD(K,4)*DEXP(-EX(K)*RR)
      BAS0(IBAS)=PP*XI*XI*XI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*YI*YI*YI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*ZI*ZI*ZI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*XI*XI*YI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*XI*XI*ZI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*YI*YI*XI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*YI*YI*ZI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*ZI*ZI*XI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*ZI*ZI*YI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*XI*YI*ZI
      GOTO 100
c
c     ...... g orbital
350   IBAS=IBAS+1
      PP=ZERO
      DO 231 K=KST,KNGI
231   PP=PP+CSPD(K,5)*DEXP(-EX(K)*RR)
      BAS0(IBAS)=PP*XI*XI*XI*XI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*YI*YI*YI*YI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*ZI*ZI*ZI*ZI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*XI*XI*XI*YI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*XI*XI*XI*ZI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*YI*YI*YI*XI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*YI*YI*YI*ZI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*ZI*ZI*ZI*XI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*ZI*ZI*ZI*YI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*XI*XI*YI*YI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*XI*XI*ZI*ZI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*YI*YI*ZI*ZI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*XI*XI*YI*ZI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*YI*YI*ZI*XI
      IBAS=IBAS+1
      BAS0(IBAS)=PP*ZI*ZI*YI*XI
c
100   CONTINUE
      RETURN
      END
