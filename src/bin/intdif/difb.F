C     //////////////////////////////////////////////////////////////////
      SUBROUTINE DIFB(NPOINT,NDER,NGUESS,S,E,oldinp)
      IMPLICIT REAL*8 (A-H,O-Z)
#include <error.h>
      integer frdf2,frdflt,frdpos,frdlpn
      integer errcod
      logical oldinp
cjtf7.24.95
      logical done
      real*8 grad, avgdsp
cdc7.24.95
      REAL*8 WE,MU,BE,DE,MASS(2),AE,WEXE
      DIMENSION E(NPOINT),S(NPOINT),H(5,6),G(5),X(5),DIFN(5),CUTN(5)
      dimension disp(4)
      COMMON/IO/IIN,IOUT,I12,I16,I21,IDER,ISCR1,ISCR2
      PARAMETER(AMU2KG=1.6605402D-27,C1=1.05052D-3,C2=1.24155D-6)
      PARAMETER(HART=4.359814445D0,MAXIT=30)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0)
      PARAMETER(FIVE=5.0D0,SIX=6.0D0)
    1 FORMAT(1X,'NPOINT MUST EXCEED 2.  EXECUTION TERMINATED.')
    2 FORMAT(/1X,'ENERGY POINTS'/)
    3 FORMAT(2F20.12)
    4 FORMAT(/1X,'INPUT GUESS FORCE CONSTANTS',//1X,'F3=',F20.8,5X,
     $    ' F4=',F20.8)
    5 FORMAT(//1X,'ITERATIONS FOR PARAMETERS',//2X,'ITER',10X,'E0',
     $     16X,'R0',16X,'F2')
    6 FORMAT(I5,F20.12,F16.8,F16.7)
    7 FORMAT(1X,'CONVERGENCE NOT REACHED IN',I5,' ITERATIONS.')
    8 FORMAT(//1X,'EQUILIBRIUM ENERGY       =',F20.12,' HARTREES')
    9 FORMAT(1X,'EQUILIBRIUM DISTANCE     =',F16.8,' A')
   10 FORMAT(1X,'QUADRATIC FORCE CONSTANT =',F16.7,' MDYNE/A')
   11 FORMAT(/1X,'INPUT GUESS FORCE CONSTANT',//1X,'F4=',F20.8)
   12 FORMAT(//1X,'ITERATIONS FOR PARAMETERS',//2X,'ITER',10X,'E0',
     $     16X,'R0',16X,'F2',14X,'F3')
   13 FORMAT(I5,F20.12,F16.8,F16.7,F16.6)
   14 FORMAT(1X,'CUBIC FORCE CONSTANT     =',F16.6,' MDYNE/A**2')
   15 FORMAT(//1X,'ITERATIONS FOR PARAMETERS',//2X,'ITER',10X,'E0',
     $     16X,'R0',16X,'F2',14X,'F3',14X,'F4')
   16 FORMAT(I5,F20.12,F16.8,F16.7,F16.6,F16.5)
   17 FORMAT(1X,'QUARTIC FORCE CONSTANT   =',F16.5,' MDYNE/A**3')
cdc.7.24.95
 18   FORMAT(//1X,'MASSES OF ATOMS =',F12.8,' A.M.U.',F12.8,' A.M.U.')
 19   FORMAT(1X,'HARMONIC VIBRATIONAL FREQUENCY,    We =',F12.5,' CM-1')
 20   FORMAT(1X,'ROTATIONAL CONSTANT,               Be =',F12.5,' CM-1')
 21   FORMAT(1X,'CENTRIFUGAL DISTORTION CONSTANT,   De =',D12.5,' CM-1')
 22   FORMAT(//1X,'ERROR READING MASSES! NO SPEC. CONST. WILL BE GIVEN')
 23   FORMAT(1X,'VIBRATION-ROTATION COUPLING CONST, Ae =',F12.5,' CM-1')
 24   FORMAT(1X,'ANHARMONIC CONSTANT,             WeXe =',F12.5,' CM-1')
cjtf.7.24.95
   25 format(//1x,'GRADIENT BY FINITE DIFFERENCE AT X=',f16.9)
   26 format(1x,'dE/dr = ',f16.9,',  BUT DANIEL IS STILL A GOOB')
   27 format(1x,'avgdsp = ',f16.9)

      CUTN(1)=1.0D-12
      CUTN(2)=1.0D-8
      CUTN(3)=1.0D-7
      CUTN(4)=1.0D-6
      CUTN(5)=1.0D-5
      IF(NDER.GE.1) GO TO 500
      IF(NPOINT.LT.3) THEN
         WRITE(IOUT,1)
         RETURN
      END IF
      WRITE(IOUT,2)
      DO 100  I=1,NPOINT
      if(oldinp) then
        READ(IIN,3) S(I),E(I)
      else
        errcod = frdf2('POINTS',s(i),i,1)
        errcod = frdf2('POINTS',e(i),i,2)
      endif
      WRITE(IOUT,3) S(I),E(I)
  100 E(I)=E(I)*HART
      IF(NPOINT.NE.3) GO TO 120
      IF(NGUESS.NE.0) THEN
          if(oldinp) then
            READ(IIN,3) F3,F4
          else
            errcod = frdf('CUBE_F',f3)
            errcod = frdf('QUART_F',f4)
          endif
          WRITE(IOUT,4) F3,F4
      ELSE
          F3=ZERO
          F4=ZERO
      END IF
cjtf7.24.95
  107 continue
      done = .true.
      do 109 i=1,2
         if (s(i).gt.s(i+1)) then
            t1 = s(i)
            s(i) = s(i+1)
            s(i+1) = t1
            t1 = e(i)
            e(i) = e(i+1)
            e(i+1) = t1
            done = .false.
         end if
  109 continue
      if (.not.done) goto 107

      disp(1) = s(2)-s(1)
      disp(2) = s(3)-s(2)
      if ((disp(2)-disp(1))**2.gt.cutn(1)) goto 118
      avgdsp = (disp(1)+disp(2))/2.0
      grad = -(e(1)-e(3))/(hart*2.0*avgdsp)
      write(iout,25) s(2)
      write(iout,26) grad
  118 continue

      R=(E(1)-E(2))/(E(1)-E(3))
      XX=(S(1)**2)*(ONE-R)-S(2)**2+R*S(3)**2
      XX=XX/(S(1)*(ONE-R)-S(2)+R*S(3))
      X(2)=XX/TWO
      XX=(S(1)-X(2))**2-(S(2)-X(2))**2
      X(3)=TWO*(E(1)-E(2))/XX
      X(1)=E(1)-X(3)*(S(1)-X(2))**2/TWO
      WRITE(IOUT,5)
      ITER=0
  106 ITER=ITER+1
      DO 102  I=1,3
      H(I,1)=ONE
      XX=S(I)-X(2)
      H(I,2)=-X(3)*XX-F3*(XX**2)/TWO-F4*(XX**3)/SIX
      H(I,3)=XX**2/TWO
      G(I)=X(1)+X(3)*(XX**2)/TWO+F3*(XX**3)/SIX+F4*(XX**4)/(SIX*FOUR)
  102 H(I,4)=E(I)-G(I)
      CALL FLIN (H,5,3,1,DET)
      ICONV=0
      DO 104  I=1,3
      X(I)=X(I)+H(I,4)
      DIFN(I)=DABS(H(I,4))
      IF(DIFN(I).GT.CUTN(I)) ICONV=1
  104 CONTINUE
      WRITE(IOUT,6) ITER,X(1)/HART,X(2),X(3)
      IF(ICONV.EQ.1.AND.ITER.LT.MAXIT) GO TO 106
      IF(ICONV.EQ.1) THEN
          WRITE(IOUT,7) MAXIT
      END IF
      X(1)=X(1)/HART
      WRITE(IOUT,8) X(1)
      WRITE(IOUT,9) X(2)
      WRITE(IOUT,10) X(3)
cdc7.24.95
C Check for Masses vector in input.dat.  If found, calculate We and output
      errcod = frdpos('MASSES')
      if(errcod.ne.EOK) go to 600
      errcod = frdlpn()
      errcod = frdflt(MASS(1))
      if(errcod.ne.EOK) go to 600
      errcod = frdflt(MASS(2))
      if(errcod.ne.EOK) go to 600
      MU = (MASS(1) * MASS(2))/(MASS(1) + MASS(2))
      MU = MU * AMU2KG
      write(iout, 18) MASS(1), MASS(2)
      WE = 5.3088375D-11*DSQRT(X(3)/MU)
      write(iout, 19) WE
      BE = 2.7992775D-26 / (MU * X(2) * X(2))
      write(iout, 20) BE
      DE = (4.0 * BE * BE * BE)/(WE * WE)
      write(iout, 21) DE
      GO TO 1000
  120 IF(NPOINT.NE.4) GO TO 140
cjtf7.24.95
  127 continue
      done = .true.
      do 129 i=1,3
         if (s(i).gt.s(i+1)) then
            t1 = s(i)
            s(i) = s(i+1)
            s(i+1) = t1
            t1 = e(i)
            e(i) = e(i+1)
            e(i+1) = t1
            done = .false.
         end if
  129 continue
      if (.not.done) goto 127

      disp(1) = s(2)-s(1)
      disp(2) = s(3)-s(2)
      disp(3) = s(4)-s(3)
      if ((disp(2)-2.0*disp(1))**2.gt.cutn(1)) goto 138
      if ((2.0*disp(3)-disp(2))**2.gt.cutn(1)) goto 138
      avgdsp = (disp(1)+disp(2)/2.0+disp(3))/3.0
      grad = (e(1)-e(4)-8.0*e(2)+8.0*e(3))/(hart*12.0*avgdsp)
      write(iout,25) (s(3)+s(2))/2.0
      write(iout,26) grad
  138 continue

      IF(NGUESS.NE.0) THEN
          if(oldinp) then
            READ(IIN,3) F4
          else
            errcod = frdf('QUART_F',f4)
          endif
          WRITE(IOUT,11) F4
      ELSE
          F4=ZERO
      END IF
      R=(E(1)-E(2))/(E(1)-E(3))
      XX=(S(1)**2)*(ONE-R)-S(2)**2+R*S(3)**2
      XX=XX/(S(1)*(ONE-R)-S(2)+R*S(3))
      X(2)=XX/TWO
      XX=(S(1)-X(2))**2-(S(2)-X(2))**2
      X(3)=TWO*(E(1)-E(2))/XX
      X(1)=E(1)-X(3)*(S(1)-X(2))**2/TWO
      X(4)=ZERO
      WRITE(IOUT,12)
      ITER=0
  122 ITER=ITER+1
      DO 124  I=1,4
      H(I,1)=ONE
      XX=S(I)-X(2)
      H(I,2)=-X(3)*XX-X(4)*(XX**2)/TWO-F4*(XX**3)/SIX
      H(I,3)=XX**2/TWO
      H(I,4)=XX**3/SIX
      G(I)=X(1)+X(3)*(XX**2)/TWO+X(4)*(XX**3)/SIX+F4*(XX**4)/(SIX*FOUR)
  124 H(I,5)=E(I)-G(I)
      CALL FLIN (H,5,4,1,DET)
      ICONV=0
      DO 126  I=1,4
      X(I)=X(I)+H(I,5)
      DIFN(I)=DABS(H(I,5))
      IF(DIFN(I).GT.CUTN(I)) ICONV=1
  126 CONTINUE
      WRITE(IOUT,13) ITER,X(1)/HART,X(2),X(3),X(4)
      IF(ICONV.EQ.1.AND.ITER.LT.MAXIT) GO TO 122
      IF(ICONV.EQ.1) THEN
          WRITE(IOUT,7) MAXIT
      END IF
      X(1)=X(1)/HART
      WRITE(IOUT,8) X(1)
      WRITE(IOUT,9) X(2)
      WRITE(IOUT,10) X(3)
      WRITE(IOUT,14) X(4)
cdc7.24.95
C Check for Masses vector in input.dat.  If found, calculate some stuff
      errcod = frdpos('MASSES')
      if(errcod.ne.EOK) go to 600
      errcod = frdlpn()
      errcod = frdflt(MASS(1))
      if(errcod.ne.EOK) go to 600
      errcod = frdflt(MASS(2))
      if(errcod.ne.EOK) go to 600
      MU = (MASS(1) * MASS(2))/(MASS(1) + MASS(2))
      MU = MU * AMU2KG
      write(iout, 18) MASS(1), MASS(2)
      WE = 5.3088375D-11*DSQRT(X(3)/MU)
      write(iout, 19) WE
      BE = 2.7992775D-26 / (MU * X(2) * X(2))
      write(iout, 20) BE
      DE = (4.0 * BE * BE * BE)/(WE * WE)
      write(iout, 21) DE
      AE = (-6.0*BE*BE/WE)*(1+(C1*WE*X(4))/(DSQRT(BE*X(3)*X(3)*X(3))))
      write(iout,23) AE
      GO TO 1000
  140 IF(NPOINT.NE.5) GO TO 160
cjtf7.24.95
  147 continue
      done = .true.
      do 149 i=1,4
         if (s(i).gt.s(i+1)) then
            t1 = s(i)
            s(i) = s(i+1)
            s(i+1) = t1
            t1 = e(i)
            e(i) = e(i+1)
            e(i+1) = t1
            done = .false.
         end if
  149 continue
      if (.not.done) goto 147

      disp(1) = s(2)-s(1)
      disp(2) = s(3)-s(2)
      disp(3) = s(4)-s(3)
      disp(4) = s(5)-s(4)
      if ((disp(2)-disp(1))**2.gt.cutn(1)) goto 158
      if ((disp(3)-disp(2))**2.gt.cutn(1)) goto 158
      if ((disp(4)-disp(3))**2.gt.cutn(1)) goto 158
      avgdsp = (disp(1)+disp(2)+disp(3)+disp(4))/4.0
      grad = (e(1)-e(5)-8.0*e(2)+8.0*e(4))/(hart*12.0*avgdsp)
      write(iout,25) s(3)
      write(iout,26) grad
  158 continue

      R=(E(1)-E(2))/(E(1)-E(3))
      XX=(S(1)**2)*(ONE-R)-S(2)**2+R*S(3)**2
      XX=XX/(S(1)*(ONE-R)-S(2)+R*S(3))
      X(2)=XX/TWO
      XX=(S(1)-X(2))**2-(S(2)-X(2))**2
      X(3)=TWO*(E(1)-E(2))/XX
      X(1)=E(1)-X(3)*(S(1)-X(2))**2/TWO
      X(4)=ZERO
      X(5)=ZERO
      WRITE(IOUT,15)
      ITER=0
  142 ITER=ITER+1
      DO 144  I=1,5
      H(I,1)=ONE
      XX=S(I)-X(2)
      H(I,2)=-X(3)*XX-X(4)*(XX**2)/TWO-X(5)*(XX**3)/SIX
      H(I,3)=XX**2/TWO
      H(I,4)=XX**3/SIX
      H(I,5)=XX**4/(SIX*FOUR)
      G(I)=X(1)+X(3)*(XX**2)/TWO+X(4)*(XX**3)/SIX
     $       +X(5)*(XX**4)/(SIX*FOUR)
  144 H(I,6)=E(I)-G(I)
      CALL FLIN (H,5,5,1,DET)
      ICONV=0
      DO 146  I=1,5
      X(I)=X(I)+H(I,6)
      DIFN(I)=DABS(H(I,6))
      IF(DIFN(I).GT.CUTN(I)) ICONV=1
  146 CONTINUE
      WRITE(IOUT,16) ITER,X(1)/HART,X(2),X(3),X(4),X(5)
      IF(ICONV.EQ.1.AND.ITER.LT.MAXIT) GO TO 142
      IF(ICONV.EQ.1) THEN
          WRITE(IOUT,7) MAXIT
      END IF
      X(1)=X(1)/HART
      WRITE(IOUT,8) X(1)
      WRITE(IOUT,9) X(2)
      WRITE(IOUT,10) X(3)
      WRITE(IOUT,14) X(4)
      WRITE(IOUT,17) X(5)
cdc7.24.95
C Check for Masses vector in input.dat.  If found, calculate some stuff
      errcod = frdpos('MASSES')
      if(errcod.ne.EOK) go to 600
      errcod = frdlpn()
      errcod = frdflt(MASS(1))
      if(errcod.ne.EOK) go to 600
      errcod = frdflt(MASS(2))
      if(errcod.ne.EOK) go to 600
      MU = (MASS(1) * MASS(2))/(MASS(1) + MASS(2))
      MU = MU * AMU2KG
      write(iout, 18) MASS(1), MASS(2)
      WE = 5.3088375D-11*DSQRT(X(3)/MU)
      write(iout, 19) WE
      BE = 2.7992775D-26 / (MU * X(2) * X(2))
      write(iout, 20) BE
      DE = (4.0 * BE * BE * BE)/(WE * WE)
      write(iout, 21) DE
      AE = (-6.0*BE*BE/WE)*(1+(C1*WE*X(4))/(DSQRT(BE*X(3)*X(3)*X(3))))
      write(iout,23) AE
      WEXE = C2*(WE/X(3))*(WE/X(3)) * ((5.0*X(4)*X(4))/(3.0*X(3))-X(5))
      write(iout,24) WEXE
  160 CONTINUE
  500 GO TO 1000
cdc7.24.95
 600  write(iout,22) 
 1000 CONTINUE
      RETURN
      END
