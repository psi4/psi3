C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        SUBROUTINE FORMP (IOPT,XMASS,NA,NC,XAR,DK1,DK2,DK3,
     $			  G3,D3,X2,X4,P2,P3,P4)
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  FORM CART PROJECTION MATRICES,
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  ON ENTRY:
C       NA                      NUMBER OF ATOMS
C       NC                      NA*3
C       XAR     (NA,3)	        REFERENCE GEOMETERY
C       DK1     (6,NC)          1ST DERIV'S
C       DK2     (3,NC,NC)       2ND DERIV'S   (USED ONLY IF NDER>1)
C       DK3     (3,NC,NC,NC)    3RD DERIV'S   (USED ONLY IF NDER>2)
C	G3	(NC,NC,NC)	SCRATCH SPACE (USED ONLY IF NDER>1)
C	D3	(NC,NC,NC)	SCRATCH SPACE (USED ONLY IF NDER>1)
C	X2	(NC,NC)		SCRATCH SPACE (USED ONLY IF NDER>1)
C	X4	(NC,NC,NC,NC)	SCRATCH SPACE (USED ONLY IF NDER>2)
C  ON RETURN:
C       P2      (NC,NC)         1ST ORDER PROJECTION MATRIX
C       P3      (NC,NC,NC)      2ND ORDER PROJECTION MATRIX
C       P4      (NC,NC,NC,NC)   3RD ORDER PROJECTION MATRIX
C  \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
        DIMENSION XAR(NA,3), P2(NC,NC), P3(NC,NC,NC), P4(NC,NC,NC,NC)
        DIMENSION DK1(6,NC), DK2(3,NC,NC), DK3(3,NC,NC,NC), IOPT(30)
	DIMENSION XMASS(NA)
        DIMENSION G3(NC,NC,NC), D3(NC,NC,NC), X2(NC,NC), X4(NC,NC,NC,NC)
  1	FORMAT (3F20.10)
  2     FORMAT (/' P(I,J) MATRIX'/)
  3     FORMAT (/' P(I,J,K) MATRIX')
  4     FORMAT (/' P(I,J,K,L) MATRIX')
  5     FORMAT (/' I=',I5/)
  6     FORMAT (6F12.6)
  7     FORMAT (/' I=',I5,' J=',I5/)

        NDER=IOPT(4)
        NSTOP=IOPT(14)

C FORM P2
        DO 10 I = 1, NA
	IX = (I-1)*3+1
        IY = IX + 1
        IZ = IY + 1
        DO 15 J = 1, NC
        P2(IX,J) = -DK1(1,J)-XAR(I,3)*DK1(5,J)+XAR(I,2)*DK1(6,J)
        P2(IY,J) = -DK1(2,J)+XAR(I,3)*DK1(4,J)-XAR(I,1)*DK1(6,J)
  15    P2(IZ,J) = -DK1(3,J)-XAR(I,2)*DK1(4,J)+XAR(I,1)*DK1(5,J)
        P2(IX,IX) = P2(IX,IX) + ONE
        P2(IY,IY) = P2(IY,IY) + ONE
  10    P2(IZ,IZ) = P2(IZ,IZ) + ONE
C SAVE TO DISK IF NSTOP SET
        IF (NSTOP.EQ.3) THEN
                REWIND I35
                WRITE (I35,1) ((P2(I,J), J=1,NC), I=1,NC)
                END IF
C PRINT IF SPECIFIED BY PRINT OPTION
        IF (LPRT(1,IOPT(6)).GE.2) THEN
                WRITE (IOUT,2)
		DO 11 I = 1, NC
  11		WRITE (IOUT,6) (P2(I,J), J=1,NC)
		END IF
        IF (NDER.LE.2) RETURN
C FORM G
        DO 110 J1 = 1, NA
	JX = (J1-1)*3+1
        JY = JX + 1
        JZ = JY + 1
        DO 110 J2 = 1, NC
        DO 110 J3 = 1, J2
        G3(JX,J2,J3) = XAR(J1,3)*DK2(2,J2,J3)-XAR(J1,2)*DK2(3,J2,J3)
        G3(JY,J2,J3) = XAR(J1,1)*DK2(3,J2,J3)-XAR(J1,3)*DK2(1,J2,J3)
        G3(JZ,J2,J3) = XAR(J1,2)*DK2(1,J2,J3)-XAR(J1,1)*DK2(2,J2,J3)
        G3(JX,J3,J2) = G3(JX,J2,J3)
        G3(JY,J3,J2) = G3(JY,J2,J3)
 110    G3(JZ,J3,J2) = G3(JZ,J2,J3)

C FORM D3
	DO 150 JX = 1, NC-1, 3
	JY = JX + 1
	JZ = JY + 1
	DO 150 I1 = 1, NC
	DO 150 I2 = 1, NC
	D3(I1,I2,JX)= P2(JY,I2)*DK1(6,I1)-P2(JZ,I2)*DK1(5,I1)
	D3(I1,I2,JY)=-P2(JX,I2)*DK1(6,I1)+P2(JZ,I2)*DK1(4,I1)
 150	D3(I1,I2,JZ)= P2(JX,I2)*DK1(5,I1)-P2(JY,I2)*DK1(4,I1)
	DO 155 IX = 1, NC-1, 3
	IY = IX + 1
	IZ = IY + 1
	DO 155 J = 1, NC
	D3(IX,J,IY) = D3(IX,J,IY) - DK1(6,J)
	D3(IX,J,IZ) = D3(IX,J,IZ) + DK1(5,J)
	D3(IY,J,IX) = D3(IY,J,IX) + DK1(6,J)
	D3(IY,J,IZ) = D3(IY,J,IZ) - DK1(4,J)
	D3(IZ,J,IX) = D3(IZ,J,IX) - DK1(5,J)
 155	D3(IZ,J,IY) = D3(IZ,J,IY) + DK1(4,J)

C FORM P3
        DO 170 L = 1, NC
        DO 170 I2 = 1, NC
        DO 160 J1 = 1, NC
        DO 160 J2 = 1, NC
 160    X2(J1,J2) = P2(J2,I2)*G3(L,J1,J2)
        DO 170 I1 = 1, NC
        XXX = ZERO
        DO 165 J1 = 1, NC
        XXX = XXX + P2(L,J1)*D3(I1,I2,J1)
        DO 165 J2 = 1, NC
 165    XXX = XXX - P2(J1,I1)*X2(J1,J2)
 170    P3(L,I1,I2) = XXX


C SAVE TO DISK IF NSTOP SET
        IF (IOPT(14).EQ.3) THEN
                REWIND I36
                WRITE (I36,1) (((P3(I,J,K), K=1,NC),J=1,NC),I=1,NC)
                END IF
C PRINT IF SPECIFIED BY PRINT OPTION
        IF (LPRT(1,IOPT(6)).GE.2) THEN
                WRITE (IOUT,3)
		DO 180 I = 1, NC
		WRITE (IOUT,5) I
                DO 180 J = 1, NC
 180          	WRITE (IOUT,6) (P3(I,J,K), K=1,NC)
                END IF
	IF (NDER.LE.3) RETURN

C FORM -G PART OF SUM IN TERM 1 OF P4
        DO 220 J1 = 1, NA
	JX = (J1-1)*3+1
        JY = JX + 1
        JZ = JY + 1
        DO 220 J2 = 1, NC
        DO 220 J3 = 1, J2
	DO 220 J4 = 1, J3
        X4(JX,J2,J3,J4)=-XAR(J1,3)*DK3(2,J2,J3,J4)
     $			+XAR(J1,2)*DK3(3,J2,J3,J4)
        X4(JY,J2,J3,J4)=-XAR(J1,1)*DK3(3,J2,J3,J4)
     $			+XAR(J1,3)*DK3(1,J2,J3,J4)
        X4(JZ,J2,J3,J4)=-XAR(J1,2)*DK3(1,J2,J3,J4)
     $ 			+XAR(J1,1)*DK3(2,J2,J3,J4)
	X4(JX,J2,J4,J3) = X4(JX,J2,J3,J4)
	X4(JX,J3,J2,J4) = X4(JX,J2,J3,J4)
	X4(JX,J3,J4,J2) = X4(JX,J2,J3,J4)
	X4(JX,J4,J2,J3) = X4(JX,J2,J3,J4)
	X4(JX,J4,J3,J2) = X4(JX,J2,J3,J4)
	X4(JY,J2,J4,J3) = X4(JY,J2,J3,J4)
	X4(JY,J3,J2,J4) = X4(JY,J2,J3,J4)
	X4(JY,J3,J4,J2) = X4(JY,J2,J3,J4)
	X4(JY,J4,J2,J3) = X4(JY,J2,J3,J4)
	X4(JY,J4,J3,J2) = X4(JY,J2,J3,J4)
	X4(JZ,J2,J4,J3) = X4(JZ,J2,J3,J4)
	X4(JZ,J3,J2,J4) = X4(JZ,J2,J3,J4)
	X4(JZ,J3,J4,J2) = X4(JZ,J2,J3,J4)
	X4(JZ,J4,J2,J3) = X4(JZ,J2,J3,J4)
 220	X4(JZ,J4,J3,J2) = X4(JZ,J2,J3,J4)

C FORM R PARTS OF SUM IN TERM 1 OF P4
	DO 240 I = 1, NC
	DO 240 JX = 1, NC-1, 3
	JY = JX + 1
	JZ = JY + 1
	DO 240 K = 1, NC
	RX = -P2(I,JY)*DK2(3,K,K)+P2(I,JZ)*DK2(2,K,K)
	RY =  P2(I,JX)*DK2(3,K,K)-P2(I,JZ)*DK2(1,K,K)
	RZ = -P2(I,JX)*DK2(2,K,K)+P2(I,JY)*DK2(1,K,K)
	X4(I,JX,K,K) = X4(I,JX,K,K) + RX
	X4(I,JY,K,K) = X4(I,JY,K,K) + RY
	X4(I,JZ,K,K) = X4(I,JZ,K,K) + RZ
	X4(I,K,JX,K) = X4(I,K,JX,K) + RX
	X4(I,K,JY,K) = X4(I,K,JY,K) + RY
	X4(I,K,JZ,K) = X4(I,K,JZ,K) + RZ
	X4(I,K,K,JX) = X4(I,K,K,JX) + RX
	X4(I,K,K,JY) = X4(I,K,K,JY) + RY
	X4(I,K,K,JZ) = X4(I,K,K,JZ) + RZ
	DO 240 L = 1, K-1
	RX = -P2(I,JY)*DK2(3,K,L)+P2(I,JZ)*DK2(2,K,L)
	RY =  P2(I,JX)*DK2(3,K,L)-P2(I,JZ)*DK2(1,K,L)
	RZ = -P2(I,JX)*DK2(2,K,L)+P2(I,JY)*DK2(1,K,L)
	X4(I,JX,K,L) = X4(I,JX,K,L) + RX
	X4(I,JX,L,K) = X4(I,JX,L,K) + RX
	X4(I,JY,K,L) = X4(I,JY,K,L) + RY
	X4(I,JY,L,K) = X4(I,JY,L,K) + RY
	X4(I,JZ,K,L) = X4(I,JZ,K,L) + RZ
 	X4(I,JZ,L,K) = X4(I,JZ,L,K) + RZ
	X4(I,K,L,JX) = X4(I,K,L,JX) + RX
	X4(I,K,L,JY) = X4(I,K,L,JY) + RY
	X4(I,K,L,JZ) = X4(I,K,L,JZ) + RZ
	X4(I,L,K,JX) = X4(I,L,K,JX) + RX
	X4(I,L,K,JY) = X4(I,L,K,JY) + RY
	X4(I,L,K,JZ) = X4(I,L,K,JZ) + RZ
	X4(I,K,JX,L) = X4(I,K,JX,L) + RX
	X4(I,L,JX,K) = X4(I,L,JX,K) + RX
	X4(I,K,JY,L) = X4(I,K,JY,L) + RY
	X4(I,L,JY,K) = X4(I,L,JY,K) + RY
	X4(I,K,JZ,L) = X4(I,K,JZ,L) + RZ
 240	X4(I,L,JZ,K) = X4(I,L,JZ,K) + RZ

C DO THE SUMMATIONS OF THE FIRST TERMS IN P4

        DO 250  L = 1, NC
        DO 250  I3 = 1, NC
        DO 250  J1 = 1, NC
        DO 250  J2 = 1, NC
        P4(L,J1,J2,I3) = ZERO
        DO 250  J3 = 1, NC
 250    P4(L,J1,J2,I3)=P4(L,J1,J2,I3)+P2(J3,I3)*X4(L,J1,J2,I3)

        DO 252  L = 1, NC
        DO 252  I3 = 1, NC
        DO 252  I2 = 1, NC
        DO 252  J1 = 1, NC
        X4(L,J1,I2,I3) = ZERO
        DO 252  J2 = 1, NC
 252    X4(L,J1,I2,I3)=P4(L,J1,I2,I3)+P2(J2,I2)*P4(L,J1,J2,I3)

        DO 254  L = 1, NC
        DO 254  I3 = 1, NC
        DO 254  I2 = 1, NC
        DO 254  I1 = 1, NC
        P4(L,I1,I2,I3) = ZERO
        DO 254  J1 = 1, NC
 254    P4(L,I1,I2,I3)=P4(L,I1,I2,I3)+P2(J1,I1)*X4(L,J1,I2,I3)

C PUT TERM 3 INTO P4
        DO 230 I1 = 1, NC
        DO 230 I2 = 1, NC
        X44 = DK1(4,I2)*DK1(4,I1)
        X55 = DK1(5,I2)*DK1(5,I1)
        X45 = DK1(4,I2)*DK1(5,I1)
        X46 = DK1(4,I2)*DK1(6,I1)
        X54 = DK1(5,I2)*DK1(4,I1)
        X66 = DK1(6,I2)*DK1(6,I1)
        X56 = DK1(5,I2)*DK1(6,I1)
        X64 = DK1(6,I2)*DK1(4,I1)
        X65 = DK1(6,I2)*DK1(5,I1)
        X6655 = X66 + X55
        X6644 = X66 + X44
        X5544 = X55 + X44
        DO 230 L = 1, NC
        DO 230 I3X = 1, NC-1, 3
        I3Y = I3X + 1
        I3Z = I3Y + 1
        SX = -P2(L,I3X)*X6655 + P2(L,I3Y)*X45 + P2(L,I3Z)*X46
        SY = P2(L,I3X)*X54 - P2(L,I3Y)*X6644 + P2(L,I3Z)*X56
        SZ = P2(L,I3X)*X64 + P2(L,I3Y)*X65 - P2(L,I3Z)*X5544
        DO 235 JX = 1, NC-1, 3
        JY = JX + 1
        JZ = JY + 1
        SX = SX +
     $  P2(L,JX)*(DK1(6,I3X)*D3(I2,I1,JY)-DK1(5,I3X)*D3(I2,I1,JZ)) +
     $  P2(L,JY)*(DK1(4,I3X)*D3(I2,I1,JZ)-DK1(6,I3X)*D3(I2,I1,JX)) +
     $  P2(L,JZ)*(DK1(5,I3X)*D3(I2,I1,JX)-DK1(4,I3X)*D3(I2,I1,JY))
        SY = SY +
     $  P2(L,JX)*(DK1(6,I3Y)*D3(I2,I1,JY)-DK1(5,I3Y)*D3(I2,I1,JZ)) +
     $  P2(L,JY)*(DK1(4,I3Y)*D3(I2,I1,JZ)-DK1(6,I3Y)*D3(I2,I1,JX)) +
     $  P2(L,JZ)*(DK1(5,I3Y)*D3(I2,I1,JX)-DK1(4,I3Y)*D3(I2,I1,JY))
 235    SZ = SZ +
     $  P2(L,JX)*(DK1(6,I3Z)*D3(I2,I1,JY)-DK1(5,I3Z)*D3(I2,I1,JZ)) +
     $  P2(L,JY)*(DK1(4,I3Z)*D3(I2,I1,JZ)-DK1(6,I3Z)*D3(I2,I1,JX)) +
     $  P2(L,JZ)*(DK1(5,I3Z)*D3(I2,I1,JX)-DK1(4,I3Z)*D3(I2,I1,JY))
        P4(L,I1,I2,I3X) = SX
        P4(L,I1,I2,I3Y) = SY
 230    P4(L,I1,I2,I3Z) = SZ
C
C FORM DELTA AND PUT INTO M (HELD IN X4)
	DO 210 I = 1, NC
	DO 211 J = 1, NC
	DO 211 K = 1, NC
 211	X4(J,K,I,I) = ZERO
	DO 210 J = 1, I-1
	VX = DK1(5,I)*DK1(6,J) - DK1(6,I)*DK1(5,J)
	VY = DK1(6,I)*DK1(4,J) - DK1(4,I)*DK1(6,J)
	VZ = DK1(4,I)*DK1(5,J) - DK1(5,I)*DK1(4,J)
	DO 210 IX = 1, NC-1, 3
	IY = IX + 1
	IZ = IY + 1
	DO 210 L = 1, NC
	X4(L,IX,I,J) = P2(L,IY)*VZ - P2(L,IZ)*VY
	X4(L,IY,I,J) = P2(L,IZ)*VX - P2(L,IX)*VZ
	X4(L,IZ,I,J) = P2(L,IX)*VY - P2(L,IY)*VX
	X4(L,IX,J,I) = -X4(L,IX,I,J)
	X4(L,IY,J,I) = -X4(L,IY,I,J)
 210   	X4(L,IZ,J,I) = -X4(L,IZ,I,J)

C SUM 2ND AND 3RD TERMS INTO M
	DO 270 L = 1, NC
	DO 270 I1X = 1, NC-1, 3
	I1Y = I1X + 1
	I1Z = I1Y + 1
	DO 270 I2 = 1, NC
	DO 270 I3 = 1, NC
	XXX = ZERO
	YYY = ZERO
	ZZZ = ZERO
	DO 271 J3 = 1, NC
 	X4(L,I1X,I2,I3) = X4(L,I1X,I2,I3)-
     $		P2(L,I1Y)*DK1(6,J3)*D3(I2,I3,J3)+
     $		P2(L,I1Z)*DK1(5,J3)*D3(I2,I3,J3)
 	X4(L,I1Y,I2,I3) = X4(L,I1Y,I2,I3)+
     $		P2(L,I1X)*DK1(6,J3)*D3(I2,I3,J3)-
     $		P2(L,I1Z)*DK1(4,J3)*D3(I2,I3,J3)
 	X4(L,I1Z,I2,I3) = X4(L,I1Z,I2,I3)-
     $		P2(L,I1X)*DK1(5,J3)*D3(I2,I3,J3)+
     $		P2(L,I1Y)*DK1(4,J3)*D3(I2,I3,J3)
 	XXX = XXX + G3(L,I1X,J3)*P3(J3,I2,I3)
 	YYY = YYY + G3(L,I1Y,J3)*P3(J3,I2,I3)
 271	ZZZ = ZZZ + G3(L,I1Z,J3)*P3(J3,I2,I3)
 	X4(L,I1X,I2,I3) = X4(L,I1X,I2,I3) - XXX
 	X4(L,I1Y,I2,I3) = X4(L,I1Y,I2,I3) - YYY
 270	X4(L,I1Z,I2,I3) = X4(L,I1Z,I2,I3) - ZZZ

C SUM IN TERM 2 INTO P4
	DO 260  L = 1, NC
	DO 260 I1 = 1, NC
	DO 260 I2 = 1, NC
	DO 260 I3 = 1, NC
	XXX = ZERO
	DO 265 J = 1, NC
 265	XXX = XXX + P2(J,I1)*X4(L,J,I2,I3)
 	P4(L,I1,I2,I3) = P4(L,I1,I2,I3) + XXX
  	P4(L,I2,I1,I3) = P4(L,I2,I1,I3) + XXX
 260   	P4(L,I2,I3,I1) = P4(L,I2,I3,I1) + XXX

C SAVE TO DISK IF NSTOP SET
        IF (IOPT(14).EQ.3) THEN
                REWIND I37
                WRITE (I37,1)
     $			((((P4(I,J,K,L),L=1,NC),K=1,NC),J=1,NC),I=1,NC)
                END IF

C PRINT IF SPECIFIED BY PRINT OPTION
        IF (LPRT(1,IOPT(6)).GE.2) THEN
                WRITE (IOUT,4)
                DO 280 I = 1, NC
		DO 280 J = 1, NC
                WRITE (IOUT,7) I, J
                DO 280 K = 1, NC
 280            WRITE (IOUT,6) (P4(I,J,K,L), L=1,NC)
                END IF

	RETURN
        END
