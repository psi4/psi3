C//////////////////////////////////////////////////////////////////////
      SUBROUTINE SQMFC(NA,NAD,NC,NSX,NSF,NISO,NMODE,XA,BS,XM,V,WT,
     $    GH,D,G,G0,F0,DF,HF,FA,W,S1,S2,S3,F,XR,XS,XT,XU,
     $    C,TA,JA,IFLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
#include <error.h>
      REAL*8 JA
      LOGICAL ITEST
      logical duminp,dumsqm
      integer dum
      integer fwkclr, fwkadd, frdboo
      integer frdpos, fcnt, frdint, frdflt, frdlpn
      integer fcnt1, fcnt2, fcnt3, frdi, frdi1
      integer frdi3, frdf3, frdi2, frdf2, frdfa1, errcod
      integer dumint
      character*8 sfhess
      DIMENSION HF(NSF,NSF),G0(NSF),F0(NSF),DF(NSF)
      DIMENSION NSC(100,0:100),IA(100,100)
      DIMENSION XA(NAD,3),BS(NC,NC),XM(NA,NISO),V(NSX,NISO),G(NSF)
      DIMENSION GH(NSX,NSX,NISO),FA(NSF),F(NSX,NSX),D(NSX,NSX)
      DIMENSION S1(NSX),S2(NSX),S3(NSX),W(NSX),XR(NC,2),XU(1),XT(NC,NC)
      DIMENSION XS(NC,NC),WT(NSX,NISO),TA(NSX,NSX,NSF)
      DIMENSION C(NSF,NSF),JA(NSX,NSF,NISO)
C     PARAMETER(CINT=42.25472D0)
C     WAVE=SQRT(10*N0)/(2*PI*C)      C IN (M/S), N0=AVOGADRO'S NUMBER
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0)
      PARAMETER(VCUT=1.0D-2,EDCUT=4.0D0,MAXIT=100,XL00=1.0D0)
      PARAMETER(XLMAX=10.0D0,ONEH=1.0D2,FDISP=1.0D-3)
      PARAMETER(SCONV=1.0D-10,GCONV=1.0D-5,XLCUT=1.0D-6,F2CUT=1.0D-6)
      PARAMETER(SCONV2=1.0D-7,GCONV2=1.0D-5)
      common /input/ duminp,dumsqm
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
   1  FORMAT(5I5)
   2  FORMAT(F10.5,I5)
   3  FORMAT(16I5)
   4  FORMAT(7F10.7)
   5  FORMAT(I5)
   6  FORMAT(I5,2F10.7)
   7  FORMAT(//,1X,'NUCLEAR CARTESIAN COORDINATES (ANG.)'/)
   8  FORMAT(2X,'ATOM',9X,'X',13X,'Y',13X,'Z'/)
   9  FORMAT(1X,I4,3X,3F14.10)
  10  FORMAT(//,1X,'MASSES FOR ISOTOPOMERS'/)
  11  FORMAT(1X,'ISO=',5(I10,4X)/)
  12  FORMAT(5X,5F14.5)
  13  FORMAT(//,1X,'EXPERIMENTAL FREQUENCIES FOR ISOTOPOMERS',
     $       ' (RELATIVE WEIGHTS GIVEN IN PARENTHESES)')
  14  FORMAT(1X,'ISO=',8I16/)
  15  FORMAT(5X,8(F9.2,' (',F7.5,')'))
  16  FORMAT(2I5)
  17  FORMAT(/,1X,'FILE16 HEADER INCONSISTENT WITH NA')
  18  FORMAT(3F20.10)
  19  FORMAT(//,1X,'UNSCALED QUADRATIC FORCE CONSTANTS IN MDYN/A, ',
     $   'MDYN/RAD, OR MDYN*A/RAD**2'/)
  20  FORMAT(//' ISOTOPOMER ',I4,/,1X,
     $    'LOWEST EIGENVALUES OF THE G MATRIX'/)
  21  FORMAT(7F12.6)
  22  FORMAT(//,1X,'SCALED QUADRATIC FORCE CONSTANTS IN MDYN/A, ',
     $   'MDYN/RAD, OR MDYN*A/RAD**2'/)
  23  FORMAT(//,1X,'OPTIMIZED SCALE FACTORS'/)
  24  FORMAT(4(I5,1X,F10.7))
  25  FORMAT(///,1X,'ISOTOPOMER ',I4,/,1X,
     $           'VIBRATIONAL FREQUENCIES (CM-1) AND EIGENVECTORS'//)
  26  FORMAT(//,1X,'ISOTOPOMER ',I4,/,1X,'POTENTIAL ENERGY ',
     $         'DISTRIBUTIONS (%) AMONG DIAGONAL ELEMENTS'//)
  27  FORMAT(//,1X,'ISOTOPOMER ',I4,/,1X,'VIBRATIONAL ASSIGNMENTS')
  28  FORMAT(//,3X,'MODE',4X,'EXPT',5X,'THEORY',3X,'RESIDUAL',3X,
     $       '% ERROR',15X,'DOMINANT COMPONENTS OF PED'/)
  29  FORMAT(1X,I5,2F10.1,2F10.2,2X,5(2X,I3,1X,'(',F5.1,')'))
  30  FORMAT(/' ITER','  LEAST-SQUARES SUM ','  DIRECTIONAL DERIVATIVE',
     $     '   DISPLACEMENT NORM')
  31  FORMAT(I5,3F20.10)
  32  FORMAT(1X,'GRADIENTS (RELATIVE)')
  33  FORMAT(5F12.5)
  34  FORMAT(1X,'LAMBDA UPDATE AVERTED   ',E12.5,5X,E12.5)
  35  FORMAT(1X,'LAMBDA=',E12.5)
  36  FORMAT(//,1X,'OPTIMIZATION OF SCALE FACTORS'/)
  38  FORMAT(/,1X,' DIAGONAL ELEMENTS OF INVERSE HESSIAN'/)
  39  FORMAT(5E12.5)
  40  FORMAT(1X,I5,10X,F10.1,22X,5(2X,I3,1X,'(',F5.1,')'))
  41  FORMAT(//,'  NSF',' NISO',' NOPT','   NH','  NWT'/,5I5)
  42  FORMAT(//,1X,'INITIAL DIAGONAL ELEMENTS OF INVERSE HESSIAN'/)
  43  FORMAT(5E10.3)
  44  FORMAT(2I5,10X,F20.10)
  45  FORMAT(//,1X,'ISOTOPOMER ',I4,/,1X,'TOTAL ENERGY ',
     $         'DISTRIBUTIONS (%)'//)
  46  FORMAT(//,3X,'MODE',4X,'EXPT',5X,'THEORY',3X,'RESIDUAL',3X,
     $       '% ERROR',15X,'DOMINANT COMPONENTS OF TED'/)
  47  FORMAT(//,1X,'SCALE FACTOR',7X,'CONNECTING COORDINATES'/)
  48  FORMAT(3X,I5,8X,16I5)
  49  FORMAT(//,1X,'OPTIMIZED SCALE FACTORS (STD. ERRORS)'/)
  50  FORMAT(2(I5,1X,F10.7,' (',F8.5,')'))
  51  FORMAT(/' STD. DEV. OF OBSERVATION OF UNIT WEIGHT',F12.7,
     $       /' EFFECTIVE DEGREES OF FREEDOM',F12.7,
     $       /' INTEGER DEGREES OF FREEDOM',I5)
  52  FORMAT(/'CORRELATION COEFFICIENTS'/)
  53  FORMAT(//,1X,'EIGENVALUES AND EIGENVECTORS OF SCALE FACTOR ',
     $      'HESSIAN'/)
  56  FORMAT(//,3X,'MODE',4X,'EXPT',5X,'THEORY',3X,'RESIDUAL',3X,
     $  '% ERROR',2X,'STD. DEV.',14X,'DOMINANT COMPONENTS OF PED'/)
  57  FORMAT(//,3X,'MODE',4X,'EXPT',5X,'THEORY',3X,'RESIDUAL',3X,
     $  '% ERROR',2X,'STD. DEV.',14X,'DOMINANT COMPONENTS OF TED'/)
  58  FORMAT(1X,I5,2F10.1,2F10.2,F10.3,2X,5(2X,I3,1X,'(',F5.1,')'))
  59  FORMAT(1X,I5,10X,F10.1,20X,F10.3,2X,5(2X,I3,1X,'(',F5.1,')'))
  60  FORMAT(/,1X,'HESSIAN FOR ITERATION',I5)
  61  FORMAT(3F20.10)
  62  format(//,'Zero-Point Vibrational Energy (cm-1) = ',f14.4/)
C*****
C   Note that no memory is allocated for JA and TA unless ABS(NH) > 1.
C   Thus, TA and JA cannot be used in first-order algorithms.
C*****
CSA
      WAVE=WAVE0
CSA
      IFLAG=0
      IF(NSF.GT.100) THEN
         IFLAG=1
         RETURN
      END IF
      DO 100  N=1,NSF
      DO 100  M=1,NSF
 100  HF(M,N)=ZERO
      if(dumsqm) then
        CALL LOCATE(IIN1,'# SQMFC ##',IERR)
        READ(IIN1,1) M,N,NOPT,NH,NWT
        DO 102 I=1,NSF
        READ(IIN1,2) FA(I),NN
        NSC(I,0)=NN
 102    READ(IIN1,3) (NSC(I,J),J=1,NN)
      else
c this is grim.  and ugly.
        errcod = fwkclr()
        errcod = fwkadd(iin1,'SQMFC')
        errcod = fwkadd(iin1,'INTDER:SQMFC')
        errcod = fwkadd(iin1,'INTDER')
        errcod = fwkadd(iin1,'DEFAULT')
        nh = 0
        errcod = frdc('SF_HESS',sfhess)
        if(sfhess.eq.'GUESS   ') nh = 0
        if(sfhess.eq.'INPUT   ') nh = 1
        if(sfhess.eq.'ANAL    ') nh = 2
        if(sfhess.eq.'NUM_FD  ') nh = -1
        if(sfhess.eq.'LLS     ') nh = -2

        errcod = frdboo('OPT_SF',dum)
        nopt = 0
        if(dum.eq.1) nopt = 1
        errcod = frdint('WEIGHT',nwt)
        errcod = frdpos('SCALE_FACTORS')
        errcod = frdlpn()
        do 103 i=1,nsf
          errcod = frdlpn()
          j=0
          errcod = frdflt(fa(i))
 105      errcod = frdint(nn)
          if(errcod.ne.Eevect) then
            j=j+1
            nsc(i,j) = nn
            goto 105
          endif
          nsc(i,0) = j
 103    continue
      endif
      DO 104 N=1,NISO
      DO 104 I=1,NSX
      V(I,N)=ZERO
 104  WT(I,N)=ZERO
      NFT=0
      DO 106 N=1,NISO
      CALL MASSIN(XM(1,N),NA,IFLAG,n,'ISOTOPES',dumsqm)
      IF(NOPT.EQ.0) GO TO 106
      if(dumsqm) then
        READ(IIN1,5) NFUND
      else
        nfund = 0
        errcod = fcnt1('FUNDAMENTALS',nfund,N)
        if (errcod.ne.EOK) goto 901
        if (nfund.eq.0) goto 901
      endif

      NFT=NFT+NFUND
      DO 108 J=1,NFUND
      if (dumsqm) then
        READ(IIN1,6) II,V(II,N),WT(II,N)
      else
        errcod = fcnt2('FUNDAMENTALS',dumint,n,j)
        if (errcod.ne.EOK) goto 901
        errcod = frdi3('FUNDAMENTALS',II,n,j,1)
        if (errcod.ne.EOK) goto 901
        errcod = frdf3('FUNDAMENTALS',V(II,N),n,j,2)
        if (errcod.ne.EOK) goto 901
        if (dumint.gt.2) then
          errcod = frdf3('FUNDAMENTALS',WT(II,N),n,j,3)
        endif
      endif

      IF(WT(II,N).EQ.ZERO) THEN
        IF(NWT.EQ.0) THEN
           WT(II,N)=ONE/V(II,N)
        ELSE IF(NWT.EQ.1) THEN
           WT(II,N)=ONE
        ELSE IF(NWT.EQ.2) THEN
           WT(II,N)=ONE/(V(II,N)*V(II,N))
        END IF
      END IF
 108  CONTINUE
 106  CONTINUE
C  Normalize the weights.
      IF(NOPT.EQ.1) THEN
        WTSUM=ZERO
        DO 110  N=1,NISO
        DO 110  I=1,NSX
 110    WTSUM=WTSUM+WT(I,N)
        DO 112  N=1,NISO
        DO 112  I=1,NSX
 112    WT(I,N)=WT(I,N)/WTSUM
      END IF
C
      IF(NH.EQ.1.AND.NOPT.EQ.1) THEN
        if (dumsqm) then
          READ(IIN1,43) (HF(I,I),I=1,NSF)
        else
          errcod = frdfa1('SF_HESS_DIAG',HF,NSF)
          if (errcod.ne.EOK) goto 902
        endif
      END IF
C
      WRITE(IOUT,41) NSF,NISO,NOPT,NH,NWT
      WRITE(IOUT,7)
      WRITE(IOUT,8)
      DO 114  I=1,NA
 114  WRITE(IOUT,9) I,(XA(I,J),J=1,3)
      WRITE(IOUT,47)
      DO 116  I=1,NSF
      NJ=NSC(I,0)
 116  WRITE(IOUT,48) I,(NSC(I,J),J=1,NJ)
      WRITE(IOUT,10)
      WRITE(IOUT,11) (I,I=1,NISO)
      DO 120  I=1,NA
 120  WRITE(IOUT,12) (XM(I,J),J=1,NISO)
      IF(NOPT.EQ.0) GO TO 124
      WRITE(IOUT,13)
      WRITE(IOUT,14) (I,I=1,NISO)
      DO 122  I=1,NSX
 122  WRITE(IOUT,15) (V(I,J),WT(I,J),J=1,NISO)
      IF(NH.EQ.1.AND.NOPT.EQ.1) THEN
      WRITE(IOUT,42)
      WRITE(IOUT,43) (HF(I,I),I=1,NSF)
      END IF
C
 124  REWIND I16
      READ(I16,16) M,N
      IF(M.NE.NA) THEN
          WRITE(IOUT,17)
          IFLAG=2
          RETURN
      END IF
      READ(I16,18) ((F(M,N),N=1,NSX),M=1,NSX)
      DO 126  M=1,NSX
      DO 126  N=1,M-1
      XX=(F(M,N)+F(N,M))/TWO
      F(M,N)=XX
 126  F(N,M)=XX
      WRITE(IOUT,19)
      CALL TABLE1(NSX,NSX,F)
C
      DO 128  J=1,NSF
      DO 130  I=1,NSX
 130  IA(I,J)=0
      DO 132  I=1,NSC(J,0)
 132  IA(NSC(J,I),J)=1
 128  CONTINUE
C  FORM G(1/2) MATRICES FOR ISOTOPOMERS.
      NV=NSX*(NSX+1)/2
      DO 150  NN=1,NISO
        II=0
        DO 152  M=1,NSX
        DO 152  N=1,M
        II=II+1
        XU(II)=ZERO
          DO 154  I=1,NC
          IJ=(I-1)/3+1
 154      XU(II)=XU(II)+BS(M,I)*BS(N,I)/XM(IJ,NN)
 152    CONTINUE
C    (NOTE THAT XU AND XS ARE IN THE SAME MEMORY LOCATION.)
      CALL RSP(NC,NSX,NV,XU,W,1,XT,S1,S2)
      WRITE(IOUT,20) NN
      IGE=MIN0(5,NSX)
      WRITE(IOUT,21) (W(I),I=1,IGE)
      DO 156  I=1,NSX
 156  W(I)=DSQRT(W(I))
      DO 158  N=1,NSX
      DO 158  M=1,NSX
         XX=ZERO
         DO 160  K=1,NSX
 160     XX=XX+XT(M,K)*XT(N,K)*W(K)
 158  GH(M,N,NN)=XX
 150  CONTINUE
      IF(NOPT.EQ.0) GO TO 400
C  BEGIN ITERATIONS FOR SCALE FACTORS
      NTER=1
      IFD=1
      IFS=1
      XL=XL00
      IOP=1
      SUM0=ZERO
      WRITE(IOUT,36)
 200  CONTINUE
      IF(NH.EQ.-1.AND.IFS.NE.0) THEN
         IF(IFS.EQ.1) THEN
            FA(IFD)=FA(IFD)+FDISP
         ELSE IF (IFS.EQ.-1) THEN
            FA(IFD)=FA(IFD)-TWO*FDISP
         END IF
      END IF
      DO 202  I=1,NSF
      DO 202  K=1,NSC(I,0)
      J=NSC(I,K)
 202  S1(J)=DSQRT(FA(I))
      SUM=ZERO
      SUM1=ZERO
      DO 204  N=1,NSF
      IF(ABS(NH).GT.1) THEN
      DO 206  M=1,NSF
      C(M,N)=ZERO
 206  HF(M,N)=ZERO
      END IF
 204  G(N)=ZERO
C
C Form gradients and Hessian (if NH = 2).
C
      DO 210  NN=1,NISO
      DO 212  N=1,NSX
      DO 212  M=1,NSX
           XX=ZERO
           DO 214  J=1,NSX
 214       XX=XX+F(M,J)*GH(J,N,NN)*S1(J)
 212  XT(M,N)=XX
      II=0
      DO 216  M=1,NSX
      DO 216  N=1,M
      II=II+1
           XX=ZERO
           DO 220  J=1,NSX
 220       XX=XX+GH(M,J,NN)*S1(J)*XT(J,N)
 216  XU(II)=XX
      CALL RSP(NC,NSX,NV,XU,S2,1,XT,W,S3)
      DO 230  I=1,NSX
           W(I)=DSIGN(DSQRT(DABS(S2(I))),S2(I))*WAVE
           S3(I)=WT(I,NN)*(ONE-V(I,NN)/W(I))
           SUM=SUM+WT(I,NN)*(W(I)-V(I,NN))**2
 230       SUM1=SUM1+WT(I,NN)*W(I)*(W(I)-V(I,NN)/TWO)
      DO 236  M=1,NSX
      DO 236  N=1,NSX
           XS(M,N)=ZERO
           DO 236 L=1,NSX
 236       XS(M,N)=XS(M,N)+GH(M,L,NN)*XT(L,N)
      IF(ABS(NH).GT.1) GO TO 260
C  NH = 0, 1, or -1 section.
      DO 240  M=1,NSX
      DO 240  N=1,M
        XR(M,N)=ZERO
        DO 242  I=1,NSX
 242    XR(M,N)=XR(M,N)+S3(I)*XS(M,I)*XS(N,I)
        XR(M,N)=XR(M,N)*F(M,N)
 240    XR(N,M)=XR(M,N)
      DO 244  N=1,NSX
      DO 244  M=1,NSX
 244    XS(M,N)=S1(N)/(TWO*S1(M))
      DO 246  MA=1,NSF
        DO 248  M=1,NSX
        DO 248  N=1,NSX
        IF(IA(M,MA)+IA(N,MA).EQ.2) THEN
           G(MA)=G(MA)+XR(M,N)
        ELSE IF(IA(M,MA).EQ.1.AND.IA(N,MA).EQ.0) THEN
           G(MA)=G(MA)+XR(M,N)*XS(M,N)
        ELSE IF(IA(M,MA).EQ.0.AND.IA(N,MA).EQ.1) THEN
           G(MA)=G(MA)+XR(M,N)*XS(N,M)
        END IF
 248    CONTINUE
 246  CONTINUE
      GO TO 210
C       NH=2 and NH=-2 section
 260  CONTINUE
      DO 262  N=1,NSX
      DO 262  M=1,NSX
 262    XT(M,N)=S1(N)/(TWO*S1(M))
      DO 264  MA=1,NSF
          DO 268  M=1,NSX
          DO 268  N=1,NSX
            D(N,M)=ZERO
            IF(IA(M,MA)+IA(N,MA).EQ.2) THEN
               D(N,M)=F(M,N)
            ELSE IF(IA(M,MA).EQ.1.AND.IA(N,MA).EQ.0) THEN
               D(N,M)=F(M,N)*XT(M,N)
            ELSE IF(IA(M,MA).EQ.0.AND.IA(N,MA).EQ.1) THEN
               D(N,M)=F(M,N)*XT(N,M)
            END IF
 268        CONTINUE
        DO 270  J=1,NSX
        DO 270  M=1,NSX
        XX=ZERO
        DO 272  N=1,NSX
 272      XX=XX+D(N,M)*XS(N,J)
 270    XR(M,J)=XX
        DO 276  I=1,NSX
        DO 276  J=1,I
        XX=ZERO
        DO 274  M=1,NSX
 274      XX=XX+XS(M,I)*XR(M,J)
 276    TA(I,J,MA)=XX
 264  CONTINUE
C
      XX=WAVE*WAVE/TWO
      DO 278  MA=1,NSF
      DO 278  I=1,NSX
      JA(I,MA,NN)=TA(I,I,MA)*XX/W(I)
 278  G(MA)=G(MA)+S3(I)*TA(I,I,MA)
      IF(NH.EQ.-2) GO TO 289
C
      DO 280  M=1,NSX
      DO 280  N=1,M
        XR(M,N)=ZERO
        DO 282  I=1,NSX
 282    XR(M,N)=XR(M,N)+S3(I)*XS(M,I)*XS(N,I)
        XR(M,N)=XR(M,N)*F(M,N)
 280    XR(N,M)=XR(M,N)
      DO 284  M=1,NSX
      DO 284  N=1,NSX
        XS(M,N)=S1(N)/(FOUR*S1(M)**3)
 284    XT(M,N)=ONE/(FOUR*S1(M)*S1(N))
C
      DO 286  MA=1,NSF
      DO 286  MB=1,MA
        DO 288  M=1,NSX
        DO 288  N=1,M
        IR=IA(M,MA)+IA(N,MA)
        IS=IA(M,MB)+IA(N,MB)
        IT=IR*IS
        XX=ZERO
        IF(MA.EQ.MB) THEN
           IF(IA(M,MA).EQ.1.AND.IA(N,MA).EQ.0) THEN
             XX=-XR(M,N)*XS(M,N)
           ELSE IF(IA(M,MA).EQ.0.AND.IA(N,MA).EQ.1) THEN
             XX=-XR(M,N)*XS(N,M)
           END IF
        ELSE
           IF(IT.NE.0) THEN
             XX=XR(M,N)*XT(M,N)
           END IF
        END IF
        IF(M.NE.N) XX=TWO*XX
 288    HF(MA,MB)=HF(MA,MB)+XX
 286    CONTINUE
C
 289  DO 290  I=1,NSX
      XU(I+NSX)=WT(I,NN)/(TWO*S2(I))
 290  XU(I)=XU(I+NSX)*V(I,NN)/W(I)
      IF(NH.EQ.-2) GO TO 293
C  Note that forming XU overwrites XS
      DO 292  I=2,NSX
      DO 292  J=1,I-1
        IF(DABS(W(I)-W(J)).LT.VCUT) THEN
          XR(I,J)=ZERO
        ELSE
          XR(I,J)=TWO*(S3(I)-S3(J))/(S2(I)-S2(J))
        END IF
 292  CONTINUE
 293  DO 294  MA=1,NSF
      DO 294  MB=1,MA
        XX=ZERO
        XY=ZERO
        DO 296  I=1,NSX
        XZ=TA(I,I,MA)*TA(I,I,MB)
        XX=XX+XZ*XU(I)
 296    XY=XY+XZ*XU(I+NSX)
        IF(NH.GE.2) THEN
          HF(MA,MB)=HF(MA,MB)+XX
        ELSE IF (NH.EQ.-2) THEN
          HF(MA,MB)=HF(MA,MB)+XY
        END IF
        C(MA,MB)=C(MA,MB)+XY/TWO
 294  CONTINUE
      IF(NH.EQ.-2) GO TO 210
      DO 297  MA=1,NSF
      DO 297  MB=1,MA
        XX=ZERO
        DO 298  J=1,NSX
        DO 298  I=J+1,NSX
 298    XX=XX+TA(I,J,MA)*TA(I,J,MB)*XR(I,J)
        HF(MA,MB)=HF(MA,MB)+XX
 297  CONTINUE
 210  CONTINUE
C
C  Multiply by conversion factors to obtain final derivatives.
C
      XX=WAVE*WAVE
      DO 304  MA=1,NSF
 304  G(MA)=G(MA)*XX
C
      IF(NH.EQ.-1.AND.IFS.NE.0) THEN
         IF(IFS.EQ.1) THEN
            DO 308  MA=1,NSF
 308        HF(MA,IFD)=G(MA)
            IFS=-1
         ELSE IF (IFS.EQ.-1) THEN
            DO 310  MA=1,NSF
 310        HF(MA,IFD)=(HF(MA,IFD)-G(MA))/(TWO*FDISP)
            IFS=1
            FA(IFD)=FA(IFD)+FDISP
            IFD=IFD+1
         END IF
         IF(IFD.GT.NSF) THEN
            IFS=0
         END IF
         GO TO 200
      END IF
C
      IF(ABS(NH).GT.1) THEN
      DO 306  MA=1,NSF
      DO 306  MB=1,MA
      C(MA,MB)=C(MA,MB)*XX
      HF(MA,MB)=HF(MA,MB)*XX
      C(MB,MA)=C(MA,MB)
 306  HF(MB,MA)=HF(MA,MB)
      END IF
C
C  Scale factor update section
C
      IF(ABS(NH).GT.1.OR.NH.EQ.-1) GO TO 350
      IF(IOP.EQ.1) THEN
        IF(NTER.EQ.1.AND.NH.EQ.0) THEN
          DO 314  I=1,NSF
  314     HF(I,I)=ONE/SUM1
        ELSE IF(NTER.GT.1) THEN
          DO 320  I=1,NSF
  320     W(I)=G(I)-G0(I)
          A0=ZERO
          A1=ZERO
          DO 324  I=1,NSF
          S2(I)=ZERO
          DO 322  J=1,NSF
  322     S2(I)=S2(I)+HF(I,J)*W(J)
          A0=A0+S2(I)*W(I)
  324     A1=A1+DF(I)*W(I)
          DO 326  J=1,NSF
          DO 326  I=1,NSF
  326     HF(I,J)=HF(I,J)+XL*DF(I)*DF(J)/A1-S2(I)*S2(J)/A0
        END IF
        A0=SUM
        A1=ZERO
        SNORM=ZERO
        FNORM=ZERO
        DO 330  I=1,NSF
        S2(I)=G(I)*FA(I)/SUM
        SNORM=SNORM+S2(I)*S2(I)
        G0(I)=G(I)
        F0(I)=FA(I)
        DF(I)=ZERO
        DO 332  J=1,NSF
  332   DF(I)=DF(I)-HF(I,J)*G(J)
        A1=A1+DF(I)*G(I)
        XX=XL*DF(I)
        FNORM=FNORM+XX*XX
  330   FA(I)=FA(I)+XX
        IOP=2
        SNORM=DSQRT(SNORM)
        FNORM=DSQRT(FNORM)
        WRITE(IOUT,30)
        WRITE(IOUT,31) NTER,SUM,A1,FNORM
        WRITE(IOUT,32)
        WRITE(IOUT,33) (S2(I),I=1,NSF)
        IF(SNORM.GT.GCONV.OR.DABS(ONE-SUM0/SUM).GT.SCONV) THEN
           SUM0=SUM
           GO TO 200
        END IF
      ELSE IF (IOP.EQ.2) THEN
        XX=ZERO
        DO 334  I=1,NSF
  334   XX=XX+DF(I)*G(I)
        A3=XL*(XX-A1)-TWO*(SUM-A0-XL*A1)
        A3=A3/XL**3
        A2=(XX-A1-THREE*A3*XL*XL)/(TWO*XL)
        DISC=A2*A2-THREE*A1*A3
        XLN=(-A2+DSQRT(DABS(DISC)))/(THREE*A3)
        IF(DISC.LT.ZERO.OR.XLN.LT.XLCUT) THEN
          WRITE(IOUT,34) DISC,XLN
          GO TO 336
        END IF
        XL=XLN
  336   WRITE(IOUT,35) XL
        DO 338  I=1,NSF
  338   FA(I)=F0(I)+XL*DF(I)
        IOP=1
        NTER=NTER+1
        IF(NTER.LT.MAXIT) GO TO 200
      END IF
      GO TO 380
C   NH = 2 and NH = -2 section
  350 CONTINUE
      WRITE(ICHECK,60) NTER
      DO 351  I=1,NSF
  351 WRITE(ICHECK,61) (HF(I,J),J=1,NSF)
      DO 352  J=1,NSF
      DO 352  I=1,NSF
      XR(I,J)=HF(I,J)
  352 XR(I,J+NSF)=ZERO
      DO 354  I=1,NSF
  354 XR(I,I+NSF)=ONE
C  Note that the above filling of XR overwrites XT.
C  Save HF before inversion for later use:
      II=0
      DO 353  I=1,NSF
      DO 353  J=1,I
      II=II+1
  353 XU(II)=HF(I,J)
      CALL FLIN(XR,NC,NSF,NSF,DET)
      DO 356  J=1,NSF
      DO 356  I=1,NSF
  356 HF(I,J)=XR(I,J+NSF)
      SNORM=ZERO
      FNORM=ZERO
      A1=ZERO
      DO 360  I=1,NSF
        S2(I)=G(I)*FA(I)/SUM
        SNORM=SNORM+S2(I)*S2(I)
        DF(I)=ZERO
        DO 358  J=1,NSF
  358   DF(I)=DF(I)-HF(I,J)*G(J)
        A1=A1+DF(I)*G(I)
        FA(I)=FA(I)+DF(I)
  360   FNORM=FNORM+DF(I)*DF(I)
        SNORM=DSQRT(SNORM)
        FNORM=DSQRT(FNORM)
        WRITE(IOUT,30)
        WRITE(IOUT,31) NTER,SUM,A1,FNORM
        WRITE(IOUT,32)
        WRITE(IOUT,33) (S2(I),I=1,NSF)
        IF(SNORM.GT.GCONV2.OR.DABS(ONE-SUM0/SUM).GT.SCONV2) THEN
           IF(NTER.LT.MAXIT) THEN
             SUM0=SUM
             NTER=NTER+1
             IF(NH.EQ.-1) THEN
               IFS=1
               IFD=1
             END IF
             GO TO 200
           END IF
        END IF
  380 IF(NTER.GE.MAXIT) THEN
        IFLAG=3
        RETURN
      END IF
      WRITE(IOUT,38)
      WRITE(IOUT,39) (HF(I,I),I=1,NSF)
 400  CONTINUE
C  Statistics section
      IF(NOPT.EQ.1.AND.ABS(NH).GE.2) THEN
        DO 502  MA=1,NSF
        DO 502  NA=1,NSF
           D(MA,NA)=ZERO
           DO 502  LA=1,NSF
 502       D(MA,NA)=D(MA,NA)+HF(MA,LA)*C(LA,NA)
        XX=ZERO
        XY=ZERO
        DO 504  MA=1,NSF
        XX=XX+D(MA,NA)
        DO 504  NA=1,NSF
 504    XY=XY+D(MA,NA)*D(MA,MA)
        DO 506  MA=1,NSF
        DO 506  NA=1,NSF
           C(MA,NA)=ZERO
           DO 506  LA=1,NSF
 506       C(MA,NA)=C(MA,NA)+D(MA,LA)*HF(LA,NA)
        XZ=FOUR*(XY-XX)
        STX=DSQRT(SUM/(NFT-XZ))
        DO 508  MA=1,NSF
 508    S2(MA)=DSQRT(C(MA,MA))*STX*TWO
        DO 510  MA=2,NSF
        DO 510  NA=1,MA-1
        XR(MA,NA)=C(MA,NA)/DSQRT(C(MA,MA)*C(NA,NA))
 510    XR(NA,MA)=XR(MA,NA)
        DO 512  MA=1,NSF
 512    XR(MA,MA)=ONE
        WRITE(IOUT,49)
        WRITE(IOUT,50) (I,FA(I),S2(I),I=1,NSF)
        WRITE(IOUT,51) STX,NFT-XZ,NFT-NSF
        WRITE(IOUT,52)
        CALL TABLE1(NC,NSF,XR)
        NVF=NSF*(NSF+1)/2
        CALL RSP(NC,NSF,NVF,XU,W,1,XT,S2,S3)
        WRITE(IOUT,53)
        CALL TABLE4(NC,NSF,NSF,W,XT)
      ELSE
        WRITE(IOUT,23)
        WRITE(IOUT,24) (I,FA(I),I=1,NSF)
      END IF
C  Scale force constants; Perform vibrational analysis for isotopomers
      DO 404  I=1,NSF
      DO 404  K=1,NSC(I,0)
      J=NSC(I,K)
 404  S1(J)=DSQRT(FA(I))
      DO 406  J=1,NSX
      DO 406  I=1,NSX
 406  F(I,J)=F(I,J)*S1(I)*S1(J)
      WRITE(IOUT,22)
      CALL TABLE1(NSX,NSX,F)
C
      DO 408 M=1,NSX
      DO 408 N=1,M
       IF(DABS(F(M,N)).GT.F2CUT) WRITE(ICHECK,44) M,N,F(M,N)
 408  CONTINUE
C
      DO 410  NN=1,NISO
      DO 414  N=1,NSX
      DO 414  M=1,NSX
           XX=ZERO
           DO 416  J=1,NSX
 416       XX=XX+F(M,J)*GH(J,N,NN)
 414  XT(M,N)=XX
      II=0
      DO 418  M=1,NSX
      DO 418  N=1,M
      II=II+1
           XX=ZERO
           DO 420  J=1,NSX
 420       XX=XX+GH(M,J,NN)*XT(J,N)
 418  XU(II)=XX
      CALL RSP(NC,NSX,NV,XU,W,1,XT,XR(1,1),XR(1,2))
      DO 430  I=1,NSX
           XX=DSIGN(DSQRT(DABS(W(I))),W(I))
 430       W(I)=XX*WAVE
      DO 436  M=1,NSX
      DO 436  N=1,NSX
           XS(M,N)=ZERO
           DO 440 L=1,NSX
 440       XS(M,N)=XS(M,N)+GH(M,L,NN)*XT(L,N)
 436  CONTINUE
      WRITE(IOUT,25) NN
      CALL TABLE2(NC,NSX,NSX,W,XS)
      IF(NMODE.EQ.1) THEN
        DO 445  N=1,NSX
             XX=ZERO
             DO 446  M=1,NSX
             XR(M,N)=XS(M,N)*F(M,M)*XS(M,N)
 446         XX=XX+DABS(XR(M,N))
             IF(XX.EQ.ZERO) GO TO 445
             DO 447  M=1,NSX
 447         XR(M,N)=ONEH*XR(M,N)/XX
 445    CONTINUE
      ELSE
        DO 450  N=1,NSX
        DO 450  M=1,NSX
 450         XR(M,N)=XS(M,N)
C   Note that the following step overwrites the XT array in memory.
        DO 451  N=1,NSX
        DO 451  M=1,NSX
 451         XR(M,N+NSX)=ZERO
        DO 452  M=1,NSX
 452         XR(M,M+NSX)=ONE
        CALL FLIN(XR,NC,NSX,NSX,DET)
        DO 453  N=1,NSX
        DO 453  M=1,NSX
 453         XR(M,N)=XS(M,N)*XR(N,M+NSX)*ONEH
      END IF
      IF(NMODE.EQ.1) THEN
        WRITE(IOUT,26) NN
      ELSE
        WRITE(IOUT,45) NN
      END IF
      CALL TABLE3(NC,NSX,W,XR)
      DO 460  J=1,NSX
      IA(5,J)=4
      DO 460  I=1,4
      IA(I,J)=0
 460  XT(I,J)=ZERO
      DO 465  J=1,NSX
          DO 470  I=1,NSX
          IF(DABS(XR(I,J)).GT.DABS(XT(1,J))) THEN
             XT(1,J)=XR(I,J)
             IA(1,J)=I
          END IF
 470      CONTINUE
          DO 472  I=1,NSX
          ITEST=I.NE.IA(1,J)
          IF(DABS(XR(I,J)).GT.DABS(XT(2,J)).AND.ITEST) THEN
             XT(2,J)=XR(I,J)
             IA(2,J)=I
          END IF
 472      CONTINUE
          DO 474  I=1,NSX
          ITEST=I.NE.IA(1,J).AND.I.NE.IA(2,J)
          IF(DABS(XR(I,J)).GT.DABS(XT(3,J)).AND.ITEST) THEN
             XT(3,J)=XR(I,J)
             IA(3,J)=I
          END IF
 474      CONTINUE
          DO 476  I=1,NSX
          ITEST=I.NE.IA(1,J).AND.I.NE.IA(2,J).AND.I.NE.IA(3,J)
          IF(DABS(XR(I,J)).GT.DABS(XT(4,J)).AND.ITEST) THEN
             XT(4,J)=XR(I,J)
             IA(4,J)=I
          END IF
 476      CONTINUE
          DO 478  K=1,4
             IF(DABS(XT(K,J)).LT.EDCUT) THEN
             IA(5,J)=K-1
             GO TO 465
             END IF
 478      CONTINUE
 465  CONTINUE
      DO 480  J=1,NSX
          DO 485  K=1,4
          XX=XS(IA(K,J),J)
          IF(XX.LT.ZERO) IA(K,J)=-IA(K,J)
 485      CONTINUE
          IF(IA(1,J).LT.0) THEN
          DO 490  K=1,4
 490      IA(K,J)=-IA(K,J)
          END IF
 480      CONTINUE
      WRITE(IOUT,27) NN
      IF(NOPT.EQ.1.AND.ABS(NH).GE.2) THEN
        DO  550  I=1,NSX
        WT(I,NN)=ZERO
            DO 552  MA=1,NSF
            DO 552  NA=1,NSF
 552        WT(I,NN)=WT(I,NN)+JA(I,MA,NN)*JA(I,NA,NN)*C(MA,NA)
 550    WT(I,NN)=DSQRT(WT(I,NN))*TWO*STX
        IF(NMODE.EQ.1) THEN
          WRITE(IOUT,56)
        ELSE
          WRITE(IOUT,57)
        END IF
        DO 560  J=1,NSX
        IF(V(J,NN).NE.ZERO) THEN
        XX=V(J,NN)-W(J)
        XY=XX/V(J,NN)*ONEH
        WRITE(IOUT,58) J,V(J,NN),W(J),XX,XY,WT(J,NN),
     $                  (IA(K,J),XT(K,J),K=1,IA(5,J))
        ELSE
        WRITE(IOUT,59) J,W(J),WT(J,NN),(IA(K,J),XT(K,J),K=1,IA(5,J))
        END IF
 560    CONTINUE
      ELSE
        IF(NMODE.EQ.1) THEN
          WRITE(IOUT,28)
        ELSE
          WRITE(IOUT,46)
        END IF
        DO 495  J=1,NSX
        IF(V(J,NN).NE.ZERO) THEN
        XX=V(J,NN)-W(J)
        XY=XX/V(J,NN)*ONEH
      WRITE(IOUT,29) J,V(J,NN),W(J),XX,XY,(IA(K,J),XT(K,J),K=1,IA(5,J))
        ELSE
          WRITE(IOUT,40) J,W(J),(IA(K,J),XT(K,J),K=1,IA(5,J))
        END IF
 495  CONTINUE
      END IF
C
      zpve = ZERO
      do 505 j=1,nsx
        if(w(j).gt.ZERO) zpve = zpve+w(j)
 505  continue
      zpve = zpve/TWO
      write(iout,62) zpve

 410  CONTINUE
      RETURN

 901  write(IOUT,*) 'error parsing FUNDAMENTALS vector'
      call pabort
      call qabort

 902  write(IOUT,*) 'error parsing SF_HESS_DIAG vector'
      call pabort
      call qabort

      END
