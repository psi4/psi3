C//////////////////////////////////////////////////////////////////////
      SUBROUTINE NORMCO(NA,NAD,NC,NSX,NFREQ,IRINT,
     $                       XMASS,XA,BS,F,XR,XS,XT,XU,IFLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),BS(NC,NC),F(NC,NC),XS(NC,NC),XT(NC,NC)
      DIMENSION XMASS(NA),W(6),XR(NC,NC),XU(1),IZ(6)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,WAVE=5140.4886D0,CINT=42.25472D0)
      PARAMETER(FRQCUT=5.0D0,TWO=2.0D0)
croag
      common/zpve/nvib, izvlim
      DATA PH,AVN / 6.626176D+00 , 6.022045D+00 /
      DATA CL / 2.99792458D+00 /
croag
   1  FORMAT(//,1X,'NUCLEAR CARTESIAN COORDINATES (ANG.) AND MASSES'/)
   2  FORMAT(2X,'ATOM',4X,'MASS',13X,'X',13X,'Y',13X,'Z'/)
   3  FORMAT(1X,I4,F12.6,3X,3F14.10)
   4  FORMAT(2I5)
   5  FORMAT(/,1X,'FILE15 HEADER INCONSISTENT WITH NA')
   6  FORMAT(3F20.10)
   7  FORMAT(//,1X,'QUADRATIC FORCE CONSTANTS IN HARTREE/BOHR**2'/)
   8  FORMAT(///,1X,'ZERO FREQUENCIES (CM-1) AND EIGENVECTORS'//)
   9  FORMAT(///,1X,'ZERO FREQUENCIES (CM-1)'/)
  10  FORMAT(6F12.4)
  11  FORMAT(///,1X,'VIBRATIONAL FREQUENCIES (CM-1) AND EIGENVECTORS'//)
  12  FORMAT(2I5)
  13  FORMAT(/,1X,'FILE17 HEADER INCONSISTENT WITH NA')
  14  FORMAT(3F20.10)
  15  FORMAT(//,1X,'DIPOLE MOMENT DERIVATIVES IN NORMAL COORDINATES'
     $      //4X,'I',7X,'M(X)',6X,'M(Y)',6X,'M(Z)'/)
  16  FORMAT(1X,I4,3X,3F10.6)
  17  FORMAT(//,24X,'IR INTENSITY'/
     $          3X,'MODE',4X,'CM-1',4X,'(KM/MOLE AND RELATIVE)')
  18  FORMAT(1X,I5,F10.1,1X,F12.5,F10.5)
  19  FORMAT(//,3X,'MODE',4X,'CM-1'/)
  20  FORMAT(1X,I5,F10.1)
  29  FORMAT(//,2X,' ZERO-POINT VIBRATIONAL ENERGY'//
     &          2X,F20.5,' IN CM-1',5X,F20.5,' IN KCAL/MOL'/)
      WRITE(IOUT,1)
      WRITE(IOUT,2)
      DO 100  I=1,NA
 100  WRITE(IOUT,3) I,XMASS(I),(XA(I,J),J=1,3)
cbfy  REWIND I15
      READ(I15,4) I,J
      IF(I.NE.NA) THEN
          WRITE(IOUT,5)
          IFLAG=1
          RETURN
      END IF
      READ(I15,6) ((F(I,J),J=1,NC),I=1,NC)
      WRITE(IOUT,7)
      CALL TABLE1(NC,NC,F)
      NV=NC*(NC+1)/2
      DO 105  J=1,NC
      JK=(J-1)/3+1
      DO 105  I=1,NC
      IK=(I-1)/3+1
 105  XR(I,J)=(F(J,I)+F(I,J))/(TWO*DSQRT(XMASS(IK)*XMASS(JK)))
      II=0
      DO 110  I=1,NC
      DO 110  J=1,I
      II=II+1
 110  XU(II)=XR(I,J)
C     NOTE THAT XU AND XS ARE IN THE SAME MEMORY LOCATION.
      CALL RSP(NC,NC,NV,XU,BS(1,1),1,XT,BS(1,2),BS(1,3))
      DO 120  I=1,NC
      IK=(I-1)/3+1
      XX=ONE/DSQRT(XMASS(IK))
      DO 120  J=1,NC
 120  XR(I,J)=XT(I,J)*XX
      DO 130  I=1,NC
           XX=ONE
           IF(BS(I,1).LT.ZERO) XX=-XX
 130       XU(I)=DSQRT(DABS(BS(I,1)))*XX*WAVE
      II=0
      IJ=0
      DO 135  J=1,NC
           XX=DABS(XU(J))
           IF(XX.LT.FRQCUT) THEN
           IJ=IJ+1
           IZ(IJ)=J
           W(IJ)=XU(J)
           DO 138  I=1,NC
 138       BS(I,IJ)=XR(I,J)
           GO TO 135
           END IF
      II=II+1
      XU(II)=XU(J)
      DO 140  I=1,NC
 140  XR(I,II)=XR(I,J)
 135  CONTINUE
      NZ=IJ
      NF=II
      IF(LPRT(3,NPRT).EQ.3) THEN
           WRITE(IOUT,8)
           CALL TABLE2(NC,NC,NZ,W,BS)
      ELSE
           WRITE(IOUT,9)
           WRITE(IOUT,10) (W(I),I=1,NZ)
      END IF
      WRITE(IOUT,11)
      CALL TABLE2(NC,NC,NF,XU,XR)
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
croag 3.30.94
c	ok, xu is in ascending order rather than descending, and
c	there are nonzero numbers in the 3n-nvib spaces.  the
c	obvious solution is to sum to NF, as determined above
croag 3.30.94
c      DO 102 I=1,NVIB
      DO 102 I=1,NF
         IF(xu(I).LT.ZLIMIT) GO TO 102
         VZERO=VZERO+xu(I)
  102 CONTINUE
      VZERO=VZERO*0.5D+00
      VKM=VZERO*CMKM
      WRITE(iout,29) VZERO,VKM
croag
croag end zpve
croag
      IF(IRINT.EQ.0) GO TO 184
      READ(I17,12) I,J
           IF(I.NE.NA) THEN
           WRITE(IOUT,13)
           IFLAG=2
           RETURN
           END IF
      READ(I17,14) ((BS(I,J),I=1,NC),J=1,3)
      XZ=ZERO
      IY=0
      DO 165  I=1,NF
          XY=ZERO
          DO 170  J=1,3
          XX=ZERO
          DO 175  K=1,NC
 175      XX=XX+XR(K,I)*BS(K,J)
          BS(I,J+3)=XX
 170      XY=XY+XX*XX
          IF(XY.GT.XZ) THEN
          IY=I
          XZ=XY
          END IF
 165  XT(I,1)=XY*CINT
      XZ=XZ*CINT
      IF(XZ.EQ.ZERO) XZ=ONE
      DO 180  I=1,NF
 180  XT(I,2)=XT(I,1)/XZ
      IF(LPRT(3,NPRT).EQ.2) THEN
          WRITE(IOUT,15)
          DO 182  I=1,NF
 182      WRITE(IOUT,16) I,(BS(I,J),J=4,6)
      END IF
 184  IF(IRINT.EQ.0) THEN
      WRITE(IOUT,19)
      DO 186  I=1,NF
 186  WRITE(IOUT,20) I,XU(I)
      ELSE
      WRITE(IOUT,17)
      DO 190  I=1,NF
 190  WRITE(IOUT,18) I,XU(I),XT(I,1),XT(I,2)
      END IF
      RETURN
      END
