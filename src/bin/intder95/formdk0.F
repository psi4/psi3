C///////////////////////////////////////////////////////////////////////
        SUBROUTINE FORMDK0 (NA,NC,XA,XMASS,RFLAG,IOPT, DOT,CP1,DOT2,CP2,
     $                     EFLAG,DK1,DK2,DK3)
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  COMPUTE THE DERIVATIVES OF THE EXTERNAL COORDS OF THE MOL IN RESPECT
C  TO THE CARTESIAN COORDINATES OF THE ATOMS AT K=0 (IE ALL EULERIAN ANG
C  ZERO AND THE CENTER OF MASS AT THE ORIGIN.
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  ON ENTRY:
C       NA                      NUMBER OF ATOMS
C       NC                      NA*3
C       XA      (NC)            REFERENCE GEOMETRY (ANG.)
C       XMASS   (NA)            MASS OF ATOMS
C	RFLAG	LOGICAL		SET TRUE IF LINEAR MOLECULE
C       DOT     (NA,NC)         SCRATCH SPACE (USED ONLY IF NDER > 1)
C       CP1     (3,NA,NC)       SCRATCH SPACE (	"		"   )
C	DOT2	(NA,NC,NC)	SCRATCH SPACE (USED ONLY IF NDER > 2)
C	CP2	(3,NA,NC,NC)	SCRATCH SPACE (	"		"   )
C  ON RETURN:
C       EFLAG   LOGICAL         SET TRUE ON FATAL ERROR
C       DK1     (6,NC)          1ST DERVIVATIVES
C       DK2     (3,NC,NC)       2ND DERVIVATIVES ROTATIONAL VARIABLES ON
C       DK3     (3,NC,NC,NC)    3RD DERVIVATIVES ROTATIONAL VARIABLES ON
C  ERROR CONDITIONS:
C	INERTIA TENSOR IS SINGULAR THUS REF GEO OR MASSES NOT VALID
C  NOTES:
C       1) HIGHEST DERIV FORMED DETERMINED BY IOPT(4)
C	2) DERIVAITVE PRINTED GIVEN IOPT(6)
C	3) DERIVATIVES SAVED TO DISK IF IOPT(14)=2
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IMPLICIT REAL*8 (A-H,O-Z)
        LOGICAL EFLAG, RFLAG
        DIMENSION XA(NA,3),XMASS(NA),DOT(NA,NC),CP1(3,NA,NC),IOPT(30)
        DIMENSION DK1(6,NC), DK2(3,NC,NC), DK3(3,NC,NC,NC)
        DIMENSION V(3), V1(3), T1(3,6), IPRMN(3), IPRMX(3)
	DIMENSION DOT2(NA,NC,NC), CP2(3,NA,NC,NC), EU(3), V2(3)
	DIMENSION IDX1(6), IDX2(6), IDX3(6)

        EQUIVALENCE (IL, IPRMN(1)), (IM, IPRMN(2)), (IN, IPRMN(3))
        EQUIVALENCE (IB, IPRMX(1)), (IG, IPRMX(2)), (ID, IPRMX(3))

	PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0)
	PARAMETER (FOUR=4.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14

	DATA (IDX1(J), J=1,6) /1, 1, 2, 2, 3, 3 /
        DATA (IDX2(J), J=1,6) /2, 3, 1, 3, 1, 2 /
        DATA (IDX3(J), J=1,6) /3, 2, 3, 1, 2, 1 /

  1     FORMAT(///'ERROR IN FORMDK0: INERTIA TENSOR SINGULAR:'
     $           /'           REFERENCE GEOMETRY OR MASSES NOT VALID.')
  2	FORMAT (3F20.10)
  3	FORMAT (6F12.6)
  4	FORMAT (/' DK(I) MATRIX FOR TRANSLATIONAL COORDINATE',I5/)
  5	FORMAT (/' DK(I) MATRIX FOR ROTATIONAL COORDINATE',I5/)
  6	FORMAT (/' DK(I,J) MATRIX FOR ROTATIONAL COORDINATE',I5/)
  7	FORMAT (/' DK(I,J,K) MATRIX FOR ROTATIONAL COORDINATE',I5)
  8	FORMAT (/' I=',I5/)

        EFLAG = .FALSE.

C FORM AND INVERT THE INERTIA TENSOR
        DO 62 I = 1, 3
        DO 62 J = 1, 6
  62    T1(I,J) = ZERO

	IF (RFLAG) THEN
	TM = ZERO
	SUM = ZERO
	DO 40 I = 1, NA
	TM = TM + XMASS(I)
  40	SUM = SUM + XMASS(I)*XA(I,3)**2
	SUM = ONE / SUM
	T1(1,4) = SUM
	T1(2,5) = SUM

	ELSE
        DO 63  I = 1, 3
  63    T1(I,I+3) = ONE

	TM = ZERO
        DO 70  K=1,NA
	TM = TM + XMASS(K)
        T1(1,1)=T1(1,1)+(XA(K,2)**2+XA(K,3)**2)*XMASS(K)
        T1(2,2)=T1(2,2)+(XA(K,1)**2+XA(K,3)**2)*XMASS(K)
        T1(3,3)=T1(3,3)+(XA(K,1)**2+XA(K,2)**2)*XMASS(K)
        T1(2,1)=T1(2,1)-XA(K,1)*XA(K,2)*XMASS(K)
        T1(3,1)=T1(3,1)-XA(K,1)*XA(K,3)*XMASS(K)
  70    T1(3,2)=T1(3,2)-XA(K,2)*XA(K,3)*XMASS(K)
        T1(1,2)=T1(2,1)
        T1(1,3)=T1(3,1)
        T1(2,3)=T1(3,2)
        CALL FLIN(T1,3,3,3,DET)
        IF (DET.EQ.ZERO) THEN
                WRITE (IOUT,1)
                EFLAG = .TRUE.
                RETURN
                END IF
	END IF

C  FIRST DERIVATIVES (TRANSLATIONAL VARIABLES)
	DO 49 I = 1, 3
	DO 49 K = 1, NC
  49	DK1(I,K) = ZERO
        DO 50 K = 1, NA
        XXX = XMASS(K) / TM
        DO 50 I = 1, 3
        L = 3*(K-1)+I
  50    DK1(I,L) = XXX

C  FIRST DERIVATIVES (ROTATIONAL VARIABLES)
        DO 10 I = 1, NA
        IX = (I-1)*3 + 1
        IY = IX + 1
        IZ = IY + 1
        DK1(4,IX) = ( T1(1,5)*XA(I,3) - T1(1,6)*XA(I,2)) * XMASS(I)
        DK1(5,IX) = ( T1(2,5)*XA(I,3) - T1(2,6)*XA(I,2)) * XMASS(I)
        DK1(6,IX) = ( T1(3,5)*XA(I,3) - T1(3,6)*XA(I,2)) * XMASS(I)
        DK1(4,IY) = (-T1(1,4)*XA(I,3) + T1(1,6)*XA(I,1)) * XMASS(I)
        DK1(5,IY) = (-T1(2,4)*XA(I,3) + T1(2,6)*XA(I,1)) * XMASS(I)
        DK1(6,IY) = (-T1(3,4)*XA(I,3) + T1(3,6)*XA(I,1)) * XMASS(I)
        DK1(4,IZ) = ( T1(1,4)*XA(I,2) - T1(1,5)*XA(I,1)) * XMASS(I)
        DK1(5,IZ) = ( T1(2,4)*XA(I,2) - T1(2,5)*XA(I,1)) * XMASS(I)
        DK1(6,IZ) = ( T1(3,4)*XA(I,2) - T1(3,5)*XA(I,1)) * XMASS(I)
10      CONTINUE

C	SAVE TO DISK IF NSTOP SET
	IF (IOPT(14).EQ.2) THEN
		REWIND I31
		WRITE (I31,2) ((DK1(I,J), J=1,NC), I=1,6)
		END IF

C	PRINT IF SPECIFIED BY PRINT OPTION
	IF (LPRT(1,IOPT(6)).GE.3) THEN
		DO 20 I = 1, 3
		WRITE (IOUT,4) I
  20		WRITE (IOUT,3) (DK1(I,J), J=1,NC)
		DO 21 I = 4, 6
		WRITE (IOUT,5) I-3
  21		WRITE (IOUT,3) (DK1(I,J), J=1,NC)
		END IF

	IF (IOPT(4).LE.2) RETURN

C  FORM (A . KIB)*MI AND (A X KIB) QUANTITIES
        DO 100 I = 1, NA
	DO 100 J = 1, NC
        DOT(I,J) = (XA(I,1)*DK1(4,J)+XA(I,2)*DK1(5,J)+
     $                     XA(I,3)*DK1(6,J)) * XMASS(I)
        CP1(1,I,J) = XA(I,2)*DK1(6,J) - XA(I,3)*DK1(5,J)
        CP1(2,I,J) = XA(I,3)*DK1(4,J) - XA(I,1)*DK1(6,J)
        CP1(3,I,J) = XA(I,1)*DK1(5,J) - XA(I,2)*DK1(4,J)
  100   CONTINUE

C  SECOND DERIVATIVES
        DO 125 IN = 1, NA
        DO 125 IB = 1, 3
        I = (IN-1)*3+IB
        DO 125 JN = 1, IN
        DO 125 JG = 1, 3
        J = (JN-1)*3+JG

	V(1) = ZERO
	V(2) = ZERO
	V(3) = ZERO
        DO 120 K = 1, NA
        V(1) = V(1) + DOT(K,I)*CP1(1,K,J)+DOT(K,J)*CP1(1,K,I)
        V(2) = V(2) + DOT(K,I)*CP1(2,K,J)+DOT(K,J)*CP1(2,K,I)
        V(3) = V(3) + DOT(K,I)*CP1(3,K,J)+DOT(K,J)*CP1(3,K,I)
  120   CONTINUE

        XXX = -XMASS(IN)*XA(IN,IB)
        YYY = -XMASS(JN)*XA(JN,JG)
        V(1) = V(1)/TWO + XXX*DK1(4,J)+YYY*DK1(4,I)
        V(2) = V(2)/TWO + XXX*DK1(5,J)+YYY*DK1(5,I)
        V(3) = V(3)/TWO + XXX*DK1(6,J)+YYY*DK1(6,I)

        V(IB) = V(IB) + DOT(IN,J)
        V(JG) = V(JG) + DOT(JN,I)

        DK2(1,I,J) = V(1)*T1(1,4)+V(2)*T1(1,5)+V(3)*T1(1,6)
        DK2(2,I,J) = V(1)*T1(2,4)+V(2)*T1(2,5)+V(3)*T1(2,6)
        DK2(3,I,J) = V(1)*T1(3,4)+V(2)*T1(3,5)+V(3)*T1(3,6)
        DK2(1,J,I) = DK2(1,I,J)
        DK2(2,J,I) = DK2(2,I,J)
        DK2(3,J,I) = DK2(3,I,J)

  125   CONTINUE

C       SAVE TO DISK IF NSTOP NONZERO
        IF (IOPT(14).EQ.2) THEN
                REWIND I32
                WRITE (I32,2) (((DK2(I,J,K), K=1,J), J=1,NC), I=1,3)
                END IF

C       PRINT IF SPECIFIED BY THE PRINT OPTION
        IF (LPRT(1,IOPT(6)).GE.3) THEN
                DO 130 I = 1, 3
                WRITE (IOUT,6) I
		DO 130 J = 1, NC
 130            WRITE (IOUT,3) (DK2(I,J,K), K=1,NC)
                END IF

        IF (IOPT(4).LE.3) RETURN

C FORM (A . KQV,RS)*MP AND (A X KQV,RS) TERMS
        DO 400 I = 1, NA
        DO 400 J = 1, NC
        DO 400 K = 1, J

        DOT2(I,J,K) = ( XA(I,1)*DK2(1,J,K) + XA(I,2)*DK2(2,J,K) +
     $                 XA(I,3)*DK2(3,J,K) ) * XMASS(I)
        DOT2(I,K,J) = DOT2(I,J,K)

        CP2(1,I,J,K) = XA(I,2)*DK2(3,J,K) - XA(I,3)*DK2(2,J,K)
        CP2(2,I,J,K) = XA(I,3)*DK2(1,J,K) - XA(I,1)*DK2(3,J,K)
        CP2(3,I,J,K) = XA(I,1)*DK2(2,J,K) - XA(I,2)*DK2(1,J,K)
        CP2(1,I,K,J) = CP2(1,I,J,K)
        CP2(2,I,K,J) = CP2(2,I,J,K)
        CP2(3,I,K,J) = CP2(3,I,J,K)
 400    CONTINUE

C THIRD DERIVATIVES
        DO 500 IL = 1, NA
        DO 500 IB = 1, 3
	ILB = (IL-1)*3+IB
        DO 500 IM = 1, IL
        DO 500 IG = 1, 3
	IMG = (IM-1)*3+IG
        DO 500 IN = 1, IM
        DO 500 ID = 1, 3
        IND = (IN-1)*3+ID

        V1(1) = ZERO
        V1(2) = ZERO
        V1(3) = ZERO

C       LOOP OVER PERMUTATIONS
        DO 600 K = 1, 6
C       NOTE IPRMN(1) AND IL, IPRMX(1) AND IB, ETC HAVE SAME ADDRESS
C	IDX VECTORS ARE PRESET IN A DATA STATEMENT WITH THE PREMUTATIONS
        IP = IPRMN(IDX1(K))
        IU = IPRMX(IDX1(K))
        IPU = (IP-1)*3+IU
        IR = IPRMN(IDX2(K))
        IS = IPRMX(IDX2(K))
        IRS = (IR-1)*3+IS
        IQ = IPRMN(IDX3(K))
        IV = IPRMX(IDX3(K))
        IQV = (IQ-1)*3+IV

	XXX = XMASS(IP)/TWO
	V(1)=(-XA(IP,IU)*DK2(1,IQV,IRS)+DK1(IU+3,IQV)*CP1(1,IP,IRS))*XXX
        V(2)=(-XA(IP,IU)*DK2(2,IQV,IRS)+DK1(IU+3,IQV)*CP1(2,IP,IRS))*XXX
        V(3)=(-XA(IP,IU)*DK2(3,IQV,IRS)+DK1(IU+3,IQV)*CP1(3,IP,IRS))*XXX

        V(IU) = V(IU) + DOT2(IP,IQV,IRS)/TWO

        EU(1) = ZERO
        EU(2) = ZERO
        EU(3) = ZERO
        EU(IU) = ONE
        XXX = ( DK1(4,IQV)*DK1(4,IRS) + DK1(5,IQV)*DK1(5,IRS) +
     $          DK1(6,IQV)*DK1(6,IRS) ) * XMASS(IP) / THREE
        V(1) = V(1) - XXX * (XA(IP,2)*EU(3) - XA(IP,3)*EU(2))
        V(2) = V(2) - XXX * (XA(IP,3)*EU(1) - XA(IP,1)*EU(3))
        V(3) = V(3) - XXX * (XA(IP,1)*EU(2) - XA(IP,2)*EU(1))

        V2(1) = ZERO
        V2(2) = ZERO
        V2(3) = ZERO
        DO 550 J = 1, NA
        V2(1) = V2(1) + DOT2(J,IQV,IRS)*CP1(1,J,IPU) +
     $                  DOT(J,IPU)*CP2(1,J,IQV,IRS)
        V2(2) = V2(2) + DOT2(J,IQV,IRS)*CP1(2,J,IPU) +
     $                  DOT(J,IPU)*CP2(2,J,IQV,IRS)
        V2(3) = V2(3) + DOT2(J,IQV,IRS)*CP1(3,J,IPU) +
     $                  DOT(J,IPU)*CP2(3,J,IQV,IRS)
 550    CONTINUE
        V1(1) = V1(1) + V(1) + V2(1)/FOUR
        V1(2) = V1(2) + V(2) + V2(2)/FOUR
        V1(3) = V1(3) + V(3) + V2(3)/FOUR

 600    CONTINUE

        V2(1) = V1(1)*T1(1,4)+V1(2)*T1(1,5)+V1(3)*T1(1,6)
        V2(2) = V1(1)*T1(2,4)+V1(2)*T1(2,5)+V1(3)*T1(2,6)
        V2(3) = V1(1)*T1(3,4)+V1(2)*T1(3,5)+V1(3)*T1(3,6)
        DK3(1,ILB,IMG,IND) = V2(1)
        DK3(2,ILB,IMG,IND) = V2(2)
        DK3(3,ILB,IMG,IND) = V2(3)
        DK3(1,ILB,IND,IMG) = V2(1)
        DK3(2,ILB,IND,IMG) = V2(2)
        DK3(3,ILB,IND,IMG) = V2(3)
        DK3(1,IMG,ILB,IND) = V2(1)
        DK3(2,IMG,ILB,IND) = V2(2)
        DK3(3,IMG,ILB,IND) = V2(3)
        DK3(1,IMG,IND,ILB) = V2(1)
        DK3(2,IMG,IND,ILB) = V2(2)
        DK3(3,IMG,IND,ILB) = V2(3)
        DK3(1,IND,ILB,IMG) = V2(1)
        DK3(2,IND,ILB,IMG) = V2(2)
        DK3(3,IND,ILB,IMG) = V2(3)
        DK3(1,IND,IMG,ILB) = V2(1)
        DK3(2,IND,IMG,ILB) = V2(2)
        DK3(3,IND,IMG,ILB) = V2(3)

 500    CONTINUE

C       SAVE TO DISK IF NSTOP NONZERO
        IF (IOPT(14).EQ.2) THEN
                REWIND I33
                WRITE (I33,2)
     $			 ((((DK3(I,J,K,L),L=1,K),K=1,J),J=1,NC),I=1,3)
                END IF

C       PRINT IF SPECIFIED BY PRINT OPTION
        IF (LPRT(1,IOPT(6)).GE.3) THEN
                DO 510 I = 1, 3
                WRITE (IOUT,7) I
		DO 510 J = 1, NC
		WRITE (IOUT,8) J
		DO 510 K = 1, NC
  510           WRITE (IOUT,3) (DK3(I,J,K,L), L=1,NC)
                END IF
	RETURN
	END
