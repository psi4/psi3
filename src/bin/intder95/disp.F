C     ////////////////////////////////////////////////////////////
C  NOTE THE OVERLAP IN MEMORY OF THE XR AND SRX ARRAYS.
      SUBROUTINE DISP(NAD,NA,NC,NS,NDISP,IOPT,
     $        XA,XMASS,S,SS,B,BS,A,TYPE,U,IA,IU,
     $        XD,XE,SD,XR,SRX,XT,IFLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
#include <error.h>
      integer fcnt, frdpos, frdlpn, frdint, frdflt
      integer errcod, mdisp, idispo, gmfile
      logical duminp,dumsqm
      CHARACTER TYPE*5,LABEL*4,opr*1,cpr*1
      DIMENSION XECK(6)
      DIMENSION IOPT(30),TYPE(NS),IA(NS,6),IU(NS,0:1),S(NS),SS(NS)
      DIMENSION U(NS,1),XA(NAD,3),XMASS(1)
      DIMENSION B(NS,NC),A(NC,NC),BS(NC,NC)
      DIMENSION SRX(NC,NC),XR(NC,1),XT(NC,NC)
      DIMENSION SD(NS),XD(NAD,3),XE(NA,3)
      DIMENSION IDISPV(6),XDISPV(5)
      PARAMETER(RAD=57.29577951308232D0)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,HALF=0.50D0,SIX=6.0D0)
      PARAMETER(TOLINV=1.0D-10,TOLDISP=1.0D-14,MXITER=20)
      common /input/ duminp, dumsqm
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
  1   FORMAT(//'NUCLEAR MASSES'/)
  2   FORMAT(4(I4,F12.6,4X))
  3   FORMAT(A4,I4)
  4   FORMAT(/,1X,'EXECUTION STOPPED DUE TO INPUT ERROR IN DISP.')
  5   FORMAT(///,1X,'VALUES OF SIMPLE INTERNAL COORDINATES ',
     $      '(ANG. OR DEG.) FOR REFERENCE GEOMETRY'/)
  6   FORMAT(4(I4,F16.10))
  7   FORMAT(///,1X,'VALUES OF SYMMETRY INTERNAL COORDINATES ',
     $      '(ANG. OR RAD.) FOR REFERENCE GEOMETRY'/)
  8   FORMAT(I5,F20.10)
  9   FORMAT(//,1X,'INTERNAL (SYMMETRY) COORDINATE FINAL VALUES'/)
 10   FORMAT(I5,F20.10)
 11   FORMAT(/,1X,' ITER = ',I5,'   MAX INTERNAL DEVIATION = ',E10.3)
 12   FORMAT(1X,'B MATRIX DETERMINANT = 0 .')
 13   FORMAT(/' B MATRIX INVERSION UNSUCCESSFUL FOR GIVEN TOLERANCE: ',
     $      E10.3)
 14   FORMAT(/' NEW CARTESIAN GEOMETRY FOUND IN ',I5,' ITERATIONS.')
 15   FORMAT(/' CARTESIAN GEOMETRY NOT FOUND IN ',I5,' ITERATIONS.')
 16   FORMAT(/' NEW CARTESIAN GEOMETRY (BOHR)'/)
 17   FORMAT(3F20.10)
 20   FORMAT(//' INTERNAL DISPLACEMENT SUBROUTINE'/)
 21   FORMAT(/' INTERNAL DISPLACEMENTS'/)
 22   FORMAT(I5,F20.10)
 23   FORMAT(15X,'  XDISP=',E10.3)
 25   FORMAT(/' DISPLACEMENT',I5)
 26   FORMAT('# GEOMUP #################')
 27   FORMAT(3F20.10)
 28   FORMAT(/'ECKART CONDITION VECTORS'/)
 29   FORMAT(3F20.10)
 30   format( 'geometry = (')
 31   format( ' (',3f20.10,')')
 32   format( ')')
 33   format( '% displacement = ',5(a1,i3,1x,f7.4,a1))
C
Cwa 4/28/1999
      DO 50  I=1,5
      IDISPV(I)=0
 50   XDISPV(I)=ZERO
      opr='('
      cpr=')'
C
      gmfile = 8
      IFLAG=0
      NSYM=IOPT(3)
      NSX=NS
      IF(NSYM.GT.0) NSX=NSYM
      WRITE(IOUT,20)
      IF(NDISP.LT.0) THEN
           WRITE(IOUT,1)
           WRITE(IOUT,2) (I,XMASS(I),I=1,NA)
      END IF
      if (duminp) then
        READ(IIN1,3) LABEL,MDISP
        IF(LABEL.NE.'DISP') THEN
           WRITE(IOUT,4)
           IFLAG=1
           RETURN
        END IF
      else
        errcod = fcnt('DISPLACEMENTS',mdisp)
c       if(errcod.ne.eok) goto 901

        errcod = frdpos('DISPLACEMENTS')
        if(errcod.ne.EOK) goto 901
        errcod = frdlpn()

      endif
C
      if (duminp) then
         idispo = iin2
         call ffile(idispo,'file07',0)
      else 
         idispo = gmfile
         call ffile(idispo,'geom',1)
      endif

      DO 1000 ID=1,MDISP
      ITER=0
      CALL MACHB(NAD,NC,NS,XA,XMASS,TYPE,IA,B,S)
      IF(NSYM.LE.0) THEN
         DO 100  I=1,NS
 100     SD(I)=S(I)
      ELSE
         DO 110  I=1,NSYM
         SD(I)=ZERO
         DO 110  L=1,IU(I,0)
         K=IU(I,L)
 110     SD(I)=SD(I)+U(I,L)*S(K)
      END IF
      DO 115 J=1,3
      DO 115 I=1,NAD
 115  XD(I,J)=XA(I,J)
      IF(ID.NE.1) GO TO 150
CPRINT
      DO 120  I=1,NS
      XR(I,1)=S(I)
      IF(TYPE(I).EQ.' BEND'.OR.TYPE(I).EQ.' LIN1'.OR.TYPE(I).EQ.
     $   ' TORS'.OR.TYPE(I).EQ.'  OUT') THEN
      XR(I,1)=XR(I,1)*RAD
      END IF
 120  CONTINUE
      WRITE(IOUT,5)
      WRITE(IOUT,6) (I,XR(I,1),I=1,NS)
      IF(NSYM.GT.0) THEN
        WRITE(IOUT,7)
        WRITE(IOUT,6) (I,SD(I),I=1,NSYM)
      END IF
C
 150  CONTINUE
      WRITE(IOUT,25) ID
      WRITE(IOUT,21)
Cwa
      kk=0
      if (duminp) then
 130    READ(IIN1,8) IC,X
        IF(IC.EQ.0) GO TO 135
        WRITE(IOUT,22) IC,X
        SD(IC)=SD(IC)+X
Cwa 4/28/99
        if(kk.lt.5) then
           kk=kk+1
           idispv(kk)=ic
           xdispv(kk)=x
           idispv(6)=kk
        end if 
c
        GO TO 130
      else
        errcod = frdlpn()
        if(errcod.ne.EOK) goto 901
c     now at first displacement vector
 201    errcod = frdint(ic)
        if(errcod.ne.Eevect) then
          errcod = frdflt(x)
          if(errcod.ne.EOK) goto 901
          WRITE(IOUT,22) IC,X
          sd(ic) = sd(ic)+x
Cwa 4/28/99
        if(kk.lt.5) then
           kk=kk+1
           idispv(kk)=ic
           xdispv(kk)=x
           idispv(6)=kk
        end if
c
          goto 201
        endif
      endif


 135  WRITE(IOUT,9)
      WRITE(IOUT,10) (I,SD(I),I=1,NSX)
C
 155  ITER=ITER+1
      CALL MACHB(NAD,NC,NS,XD,XMASS,TYPE,IA,B,S)
      IF(NSYM.LE.0) THEN
         DO 160  I=1,NS
         SS(I)=S(I)
         DO 160  J=1,NC
 160     BS(I,J)=B(I,J)
      ELSE
         DO 162  J=1,NC
         DO 162  I=1,NSYM
           BS(I,J)=ZERO
           DO 164  L=1,IU(I,0)
           K=IU(I,L)
 164       BS(I,J)=BS(I,J)+U(I,L)*B(K,J)
 162     CONTINUE
         DO 166  I=1,NSYM
           SS(I)=ZERO
           DO 166  L=1,IU(I,0)
           K=IU(I,L)
 166       SS(I)=SS(I)+U(I,L)*S(K)
      END IF
C
      ICONV=1
      SDISP=ZERO
      DO 170  I=1,NSX
      X=DABS(SD(I)-SS(I))
 170  IF(X.GT.SDISP) SDISP=X
      IF(SDISP.LE.TOLDISP) ICONV=0
      WRITE(IOUT,11) ITER,SDISP
      IF(ICONV.EQ.0) GO TO 500
C
      IF(NDISP.LT.0) THEN
      DO 180 I=1,NA
 180  XT(I,1)=XMASS(I)
      ELSE
      DO 182 I=1,NA
 182  XT(I,1)=ONE
      END IF
      CALL BINVRT(NSX,NC,XT(1,1),BS,XR,A,SRX,
     $     XT(1,2),IFLAG,TOLINV,1)
      IF(IFLAG.EQ.1) THEN
         WRITE(IOUT,12)
         RETURN
      ELSE IF (IFLAG.EQ.2) THEN
         WRITE(IOUT,13) TOLINV
         RETURN
      END IF
C
      II=0
      DO 200 I=1,NA
      DO 200 J=1,3
      II=II+1
           XE(I,J)=ZERO
           DO 205 K=1,NSX
 205       XE(I,J)=XE(I,J)+A(II,K)*(SD(K)-SS(K))
 200  CONTINUE
      IF(ABS(NDISP).EQ.1) GO TO 400
C
       CALL MACHX(NAD,NC,NS,NSX,IOPT,XD,XMASS,TYPE,IA,A,S,U,IU,XR,SRX)
       DO 300 K=1,NSX
       CALL XIN(NC,NC,SRX,-K,ISCR1)
          XR(K,1)=ZERO
          DO 305 I=1,NC
          II=(I-1)/3+1
          JI=I-3*(II-1)
          DO 305 J=1,NC
          IJ=(J-1)/3+1
          JJ=J-3*(IJ-1)
 305      XR(K,1)=XR(K,1)+SRX(I,J)*XE(IJ,JJ)*XE(II,JI)
 300      XR(K,1)=SD(K)-SS(K)-HALF*XR(K,1)
C
 350   II=0
       XDISP=ZERO
       DO 360 I=1,NA
       DO 360 J=1,3
       II=II+1
           X=ZERO
           DO 365 K=1,NSX
 365       X=X+A(II,K)*XR(K,1)
           W=DABS(XE(I,J)-X)
           IF(W.GT.XDISP) XDISP=W
           XE(I,J)=X
 360   CONTINUE
       WRITE(IOUT,23) XDISP
C
 400  DO 405 J=1,3
      DO 405 I=1,NA
 405   XD(I,J)=XD(I,J)+XE(I,J)
      IF(ITER.LT.MXITER) GO TO 155
C   ROTATE COORDINATES, PRINT NEW CARTESIANS, WRITE OUT GEOM TO INPUT.
 500  IF(ICONV.EQ.0) THEN
         WRITE(IOUT,14) ITER
      ELSE
         WRITE(IOUT,15) MXITER
         RETURN
      END IF
C
      IF(NDISP.GT.0) GO TO 510
      CALL ROTC(NA,NAD,XMASS,XA,XD,IFLAG)
C
 510  DO 525 I=1,6
 525  XECK(I)=ZERO
      DO 530 I=1,NA
      DO 535 J=1,3
 535  XECK(J)=XECK(J)+XMASS(I)*(XD(I,J)-XA(I,J))
      XECK(4)=XECK(4)+XMASS(I)*(XD(I,2)*XA(I,3)-XD(I,3)*XA(I,2))
      XECK(5)=XECK(5)+XMASS(I)*(XD(I,3)*XA(I,1)-XD(I,1)*XA(I,3))
 530  XECK(6)=XECK(6)+XMASS(I)*(XD(I,1)*XA(I,2)-XD(I,2)*XA(I,1))
      WRITE(IOUT,28)
      WRITE(IOUT,29) (XECK(I),I=1,6)
C
      
      DO 605  J=1,3
      DO 605  I=1,NA
 605  XD(I,J)=XD(I,J)/BOHR
      WRITE(IOUT,16)
      DO 610 I=1,NA
 610  WRITE(IOUT,17) (XD(I,J),J=1,3)
      if(duminp) then
        WRITE(idispo,26)
        DO 620 I=1,NA
 620    WRITE(idispo,27) (XD(I,J),J=1,3)
      else
Cwa 4/28/99
        kw=idispv(6)
        write(idispo,33) (opr,idispv(k),xdispv(k),cpr,k=1,kw)
        write(idispo,30)
        do 621 i = 1,na
 621      write(idispo,31) (XD(I,J),J=1,3)
        write(idispo,32)
      endif
C

 1000 CONTINUE

c      if (.not.duminp) close(idispo)
       close(idispo)
      
      RETURN

 901  write(6,*) 'error parsing DISPLACEMENTS vector'
      call pabort
      call qabort
 902  write(6,*) 'can only handle 4 coordinates per displacement'
      call pabort
      call qabort

      END

