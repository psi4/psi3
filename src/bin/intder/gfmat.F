C//////////////////////////////////////////////////////////////////////
      SUBROUTINE GFMAT(NA,NAD,NC,NSX,NFREQ,IRINT,IA,
     $                       XMASS,XA,W,BS,F,XR,XS,XT,XU,IFLAG,lbls)
      IMPLICIT REAL*8 (A-H,O-Z)
      character*80 lbls(nsx)
      integer max
      integer frdfa2
      real*8 btmp(600)
      LOGICAL ITEST
      DIMENSION XA(NAD,3),BS(NC,NC),F(NC,NC),XS(NC,NC),XT(NC,NC)
      DIMENSION XMASS(NA),W(NSX),XR(NC,NC),XU(1),IA(5,NSX)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,WAVE=1302.781D0,CINT=42.25472D0)
      PARAMETER(PEDCUT=4.0D0,TWO=2.0D0)
croag
      common/zpve/nvib, izvlim
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /
      DATA CL / 2.99792458D+00 /
croag
C     WAVE=SQRT(10*N0)/(2*PI*C)      C IN (M/S), N0=AVOGADRO'S NUMBER
C     CINT=
   1  FORMAT(//,1X,'NUCLEAR CARTESIAN COORDINATES (ANG.) AND MASSES'/)
   2  FORMAT(2X,'ATOM',4X,'MASS',13X,'X',13X,'Y',13X,'Z'/)
   3  FORMAT(1X,I4,F12.6,3X,3F14.10)
   4  FORMAT(2I5)
   5  FORMAT(/,1X,'FILE16 HEADER INCONSISTENT WITH NA')
   6  FORMAT(3F20.10)
   7  FORMAT(//,1X,'QUADRATIC FORCE CONSTANTS IN MDYN/A, MDYN/RAD, ',
     $   'OR MDYN*A/RAD**2'/)
   8  FORMAT(//,1X,'G MATRIX')
   9  FORMAT(//,1X,'EIGENVALUES OF THE G MATRIX'/)
  10  FORMAT(7F12.6)
  11  FORMAT(7F10.6)
  12  FORMAT(///,1X,'LOWEST EIGENVALUES OF THE G MATRIX'/)
  13  FORMAT(///,1X,'VIBRATIONAL FREQUENCIES (CM-1) AND EIGENVECTORS'//)
  14  FORMAT(//,1X,'POTENTIAL ENERGY DISTRIBUTIONS (%) ',
     $             'AMONG DIAGONAL ELEMENTS'//)
  15  FORMAT(//,1X,'VIBRATIONAL ASSIGNMENTS')
  16  FORMAT(//,3X,'MODE',35X,'CM-1',17X,'DOMINANT COMPONENTS OF PED'/)
  17  FORMAT(1X,I5,1x,a30,F10.1,2X,5(2X,I3,1X,'(',F5.1,')'))
  18  FORMAT(3F20.10)
  19  FORMAT(//,24X,'IR INTENSITY'/
     $          3X,'MODE',35X,'CM-1',4X,'(KM/MOLE AND RELATIVE)',
     $          15X,'DOMINANT COMPONENTS OF PED'/)
  20  FORMAT(1X,I5,1x,a30,
     +  F10.1,1X,F12.5,F10.5,3X,5(2X,I3,1X,'(',F5.1,')'))
  21  FORMAT(//,1X,'DIPOLE MOMENT DERIVATIVES IN (SYMMETRY) ',
     $  'INTERNAL COORDINATES'//4X,'I',7X,'M(X)',6X,'M(Y)',6X,'M(Z)'/)
  22  FORMAT(1X,I4,3X,3F10.6)
  23  FORMAT(//,1X,'DIPOLE MOMENT DERIVATIVES IN NORMAL COORDINATES'
     $      //4X,'I',7X,'M(X)',6X,'M(Y)',6X,'M(Z)'/)
  24  FORMAT(1X,I4,3X,3F10.6)
  25  FORMAT(2I5)
  26  FORMAT(/,1X,'FILE18 HEADER INCONSISTENT WITH NA')
  27  FORMAT(3F20.10)
  29  FORMAT(//,2X,' ZERO-POINT VIBRATIONAL ENERGY'//
     &          2X,F20.5,' IN CM-1',5X,F20.5,' IN KCAL/MOL'/)
      WRITE(IOUT,1)
      WRITE(IOUT,2)
      DO 100  I=1,NA
 100  WRITE(IOUT,3) I,XMASS(I),(XA(I,J),J=1,3)
      IF(NFREQ.EQ.4) THEN
cets050491          READ(IIN1,11) ((F(M,N),N=M,NSX),M=1,NSX)
cets050491          DO 101  M=1,NSX
cets050491          DO 101  N=M,NSX
cets050491 101      F(N,M)=F(M,N)
          errcod = frdfa2('Q_F_CONSTS',f,nsx,nsx)
          if(errcod.ne.EOK) then
            write(6,*) 'trouble reading in quadratic force constants'
            write(6,*) 'need Q_F_CONSTS'
            call pabort
            call qabort
          endif
          GO TO 102
      END IF
cbfy  REWIND I16
      READ(I16,4) M,N
      IF(M.NE.NA) THEN
          WRITE(IOUT,5)
          IFLAG=1
          RETURN
      END IF
      READ(I16,6) ((F(M,N),N=1,NSX),M=1,NSX)
 102  WRITE(IOUT,7)
      CALL TABLE1(NC,NSX,F)
      NV=NSX*(NSX+1)/2
      DO 105  N=1,NSX
      DO 105  M=1,NSX
 105  XR(M,N)=ZERO
      DO 110  I=1,NC
      II=(I-1)/3+1
      XX=ONE/XMASS(II)
         DO 115  N=1,NSX
         DO 115  M=1,NSX
 115     XR(M,N)=XR(M,N)+BS(M,I)*BS(N,I)*XX
 110  CONTINUE
      IF(LPRT(3,NPRT).EQ.1) THEN
         WRITE(IOUT,8)
         CALL TABLE1(NC,NSX,XR)
      END IF
C     NOTE THAT XU AND XS ARE IN THE SAME MEMORY LOCATION.
      II=0
      DO 120  M=1,NSX
      DO 120  N=1,M
      II=II+1
 120  XU(II)=XR(M,N)
      CALL RSP(NC,NSX,NV,XU,W,1,XT,BS(1,1),BS(1,2))
      IF(LPRT(3,NPRT).EQ.1) THEN
         WRITE(IOUT,9)
         WRITE(IOUT,10) (W(I),I=1,NSX)
         GO TO 126
      END IF
      WRITE(IOUT,12)
      IGE=MIN0(5,NSX)
      WRITE(IOUT,10) (W(I),I=1,IGE)
 126  DO 128  I=1,NSX
 128  W(I)=DSQRT(W(I))
      DO 130  N=1,NSX
      DO 130  M=1,NSX
         XR(M,N)=ZERO
         DO 140  K=1,NSX
 140     XR(M,N)=XR(M,N)+XT(M,K)*XT(N,K)*W(K)
 130  CONTINUE
      DO 150  N=1,NSX
      DO 150  M=1,NSX
           XX=ZERO
           DO 152  J=1,NSX
 152       XX=XX+(F(M,J)+F(J,M))*XR(J,N)/TWO
 150  XT(M,N)=XX
      II=0
      DO 154  M=1,NSX
      DO 154  N=1,M
      II=II+1
           XX=ZERO
           DO 156  J=1,NSX
 156       XX=XX+XR(M,J)*XT(J,N)
 154  XU(II)=XX
      CALL RSP(NC,NSX,NV,XU,W,1,XT,BS(1,1),BS(1,2))
      DO 160  I=1,NSX
           XX=ONE
           IF(W(I).LT.ZERO) XX=-XX
 160       W(I)=DSQRT(DABS(W(I)))*XX*WAVE
      DO 165  N=1,NSX
      DO 165  M=1,NSX
           XS(M,N)=ZERO
           DO 170  I=1,NSX
 170       XS(M,N)=XS(M,N)+XR(M,I)*XT(I,N)
 165  CONTINUE
      DO 175  N=1,NSX
           XX=ZERO
           DO 180  M=1,NSX
           XT(M,N)=XS(M,N)*F(M,M)*XS(M,N)
 180       XX=XX+DABS(XT(M,N))
           IF(XX.EQ.ZERO) GO TO 175
           DO 185  M=1,NSX
 185       XT(M,N)=1.0D2*XT(M,N)/XX
 175  CONTINUE
      WRITE(IOUT,13)
      CALL TABLE2(NC,NSX,NSX,W,XS)
      WRITE(IOUT,14)
      CALL TABLE3(NC,NSX,W,XT)
croag 
croag calculate the zpve and print it out.
croag
C
C  LOWER LIMIT TO THE ZPVE CORRECTION
      ZLIMIT=FLOAT(IZVLIM)
C
C  A CONVERSION FACTOR FROM CM-1 TO KCAL.MOLE-1
      CMKM=(CL*PH*AVN*1.0D-04)/4.184D+00
C
C   EVALUATE VIBRATIONAL ZERO-POINT ENERGY
C   THE SMALLEST FREQUENCY IS LIMITED TO ZLIMIT
      VZERO=0.0D+00
      DO 103 I=1,NVIB
         IF(W(I).LT.ZLIMIT) GO TO 103
         VZERO=VZERO+W(I)
  103 CONTINUE
      VZERO=VZERO*0.5D+00
      VKM=VZERO*CMKM
      WRITE(iout,29) VZERO,VKM
croag
croag end zpve
croag
      DO 200  J=1,NSX
      IA(5,J)=4
      DO 200  I=1,4
      IA(I,J)=0
 200  XR(I,J)=ZERO
      DO 205  J=1,NSX
          DO 210  I=1,NSX
          IF(DABS(XT(I,J)).GT.DABS(XR(1,J))) THEN
             XR(1,J)=XT(I,J)
             IA(1,J)=I
          END IF
 210      CONTINUE
          DO 212  I=1,NSX
          ITEST=I.NE.IA(1,J)
          IF(DABS(XT(I,J)).GT.DABS(XR(2,J)).AND.ITEST) THEN
             XR(2,J)=XT(I,J)
             IA(2,J)=I
          END IF
 212      CONTINUE
          DO 214  I=1,NSX
          ITEST=I.NE.IA(1,J).AND.I.NE.IA(2,J)
          IF(DABS(XT(I,J)).GT.DABS(XR(3,J)).AND.ITEST) THEN
             XR(3,J)=XT(I,J)
             IA(3,J)=I
          END IF
 214      CONTINUE
          DO 216  I=1,NSX
          ITEST=I.NE.IA(1,J).AND.I.NE.IA(2,J).AND.I.NE.IA(3,J)
          IF(DABS(XT(I,J)).GT.DABS(XR(4,J)).AND.ITEST) THEN
             XR(4,J)=XT(I,J)
             IA(4,J)=I
          END IF
 216      CONTINUE
          DO 218  K=1,4
             IF(DABS(XR(K,J)).LT.PEDCUT) THEN
             IA(5,J)=K-1
             GO TO 205
             END IF
 218      CONTINUE
 205  CONTINUE
      DO 230  J=1,NSX
          DO 235  K=1,5
          XX=XS(IA(K,J),J)
          IF(XX.LT.ZERO) XR(K,J)=-XR(K,J)
 235      CONTINUE
          XX=XR(1,J)
          DO 240  K=1,5
          IF(XX.LT.ZERO) XR(K,J)=-XR(K,J)
 240      CONTINUE
 230      CONTINUE
      IF(IRINT.EQ.0) GO TO 290
      IF(IRINT.EQ.1) THEN
         READ(I18,25) I,J
         IF(I.NE.NA) THEN
         WRITE(IOUT,26)
         IFLAG=2
         RETURN
         END IF
         READ(I18,27) ((BS(I,J),I=1,NSX),J=1,3)
      ELSE
cets050491      DO 250  I=1,NSX
cets050491 250  READ(IIN1,18) (BS(I,J),J=1,3)
        errcod = frdfa2('INT_DIPDER',btmp,3,nsx)
        if(errcod.ne.EOK) then
          write(6,*) 'trouble reading in dipole derivatives'
          call pabort
          call qabort
        endif
        ij=0
        do 250 i=1,nsx
          do 250 j=1,3
            ij=ij+1
            bs(i,j) = btmp(ij)
 250    continue
      END IF
      WRITE(IOUT,21)
      DO 255  I=1,NSX
 255  WRITE(IOUT,22) I,(BS(I,J),J=1,3)
      XZ=ZERO
      IZ=0
      DO 265  I=1,NSX
          XY=ZERO
          DO 270  J=1,3
          XX=ZERO
          DO 275  K=1,NSX
 275      XX=XX+XS(K,I)*BS(K,J)
          BS(I,J+3)=XX
 270      XY=XY+XX*XX
          IF(XY.GT.XZ) THEN
          IZ=I
          XZ=XY
          END IF
 265  XT(I,1)=XY*CINT
      XZ=XZ*CINT
      IF(XZ.EQ.ZERO) XZ=ONE
      DO 280  I=1,NSX
 280  XT(I,2)=XT(I,1)/XZ
      IF(LPRT(3,NPRT).EQ.2) THEN
          WRITE(IOUT,23)
          DO 282  I=1,NSX
 282      WRITE(IOUT,24) I,(BS(I,J),J=4,6)
      END IF
 290  WRITE(IOUT,15)
      IF(IRINT.EQ.0) THEN
          WRITE(IOUT,16)
          DO 300  J=1,NSX
          max=ia(1,j)
 300      WRITE(IOUT,17) J,lbls(max),W(J),(IA(K,J),XR(K,J),K=1,IA(5,J))
      ELSE
          WRITE(IOUT,19)
          DO 310  J=1,NSX
          max=ia(1,j)
          WRITE(IOUT,20) J,lbls(max),W(J),XT(J,1),XT(J,2),
     $                  (IA(K,J),XR(K,J),K=1,IA(5,J))
 310      CONTINUE
      END IF
      RETURN
      END
