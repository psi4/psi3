C //////////////////////////////////////////////////////////////////////
       SUBROUTINE PRJTST (NC,IOPT,P2,P3,P4,F1,F2,F3,VP2,VP3,X)
C //////////////////////////////////////////////////////////////////////
C VERIFY PROJECTION MATRICIES USING TRANSLATIONAL AND ROTATIONAL
C INVARIANCE RELTIONSHIPS.
C //////////////////////////////////////////////////////////////////////
C ON ENTRY:
C NC       INTEGER          3 * NUMBER OF ATOMS
C P2       (NC,NC)          1ST ORDER PROJECTION MATRIX
C P3       (NC,NC,NC)       2ND    "            "
C P4       (NC,NC,NC,NC)    3RD    "            "
C F1       (NC)             SCRATCH SPACE
C F2       (NC,NC)             "      "
C F3       (NC,NC,NC)          "      "
C VP2      (NC,NC)             "      "
C VP3      (NC,NC,NC)          "      "
C X        (NC,NC)             "      "
C //////////////////////////////////////////////////////////////////////
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION P2(NC,NC), P3(NC,NC,NC), P4(NC,NC,NC,NC)
      DIMENSION F2(NC,NC), F3(NC,NC,NC),X(NC,NC),VP2(NC,NC),F1(NC)
      DIMENSION VP3(NC,NC,NC),IOPT(30)
      LOGICAL FLAG1, FLAG2
      CHARACTER STRING*10
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      PARAMETER (ERTOL=1.0D-5)
      PARAMETER (ZERO=0.0D0,TWO=2.0D0)
  1   FORMAT (2I5)
  2   FORMAT (3F20.10)

      WRITE (IOUT,4)
  4   FORMAT(/'TESTING OF PROJECTION MATRICIES USING ROTATIONAL ',
     $        'INVARIANCE RELATIONS')

      IF (IOPT(4).LT.2 .OR.IOPT(7).NE.3.OR.IOPT(5).EQ.0) THEN
            WRITE (IOUT,5)
  5         FORMAT ('   NO TESTING WAS PERFORMED: INCONSISTENT OPTIONS')
            RETURN
            END IF

      NA = NC/3
C
C TEST 2ND ORDER PROJECTION MATRIX
C THE FOLLOWING ASSUMES THAT THE 1ST ORDER PROJECTION MATRIX IS CORRECT
C
      WRITE(IOUT,9)
  9   FORMAT('   TESTING 2ND ORDER PROJECTION MATRIX, P(I,J,K):')

C READ IN F1 AND F2
      REWIND I11
 1000 READ (I11,1001,END=1002) STRING
 1001 FORMAT (A10)
      GOTO 1000
 1002 DO 1010 K = 1, NA+1
 1010    BACKSPACE I11
      DO 1030 J = 1, NA
 1030   READ (I11,1031) (F1((J-1)*3+K),K=1,3)
 1031 FORMAT (20X,3F20.10)

      REWIND I15
      READ(I15,1) M,N
      READ(I15,2) ((F2(M,N), N=1,NC), M=1,NC)
      DO 20 M=2,NC
      DO 20 N=1,M-1
        XX = (F2(M,N)+F2(N,M))/TWO
        F2(M,N)=XX
  20    F2(N,M)=XX

C COMPUTE Vij*'s

      DO 30 I = 1, NC
      DO 30 J = 1, I
        DO 27 K = 1, NC
        DO 27 L = 1, NC
   27     VP2(I,J) = VP2(I,J) + P2(I,K)*F2(K,L)*P2(L,J)
   30     VP2(J,I) = VP2(I,J)

      FLAG1 = .TRUE.
      DO 50 I = 1, NC
      DO 50 J = 1, NC

        XX = ZERO
        DO 40 K = 1, NC
   40        XX = XX + F1(K)*P3(K,I,J)

        ERR = F2(I,J) - VP2(I,J) - XX
        IF (ABS(ERR).GT.ERRTOL) THEN
             FLAG1 = .FALSE.
             WRITE (IOUT,42) I,J,ERR
   42        FORMAT ('      I= ?  J= ',I3,' K= ',I3,' ERROR= ',F12.8)
             END IF

   50   CONTINUE

        IF (FLAG1) WRITE (IOUT,54)
   54              FORMAT ('      NO ERRORS DETECTED.')



      IF (IOPT(4).LE.2) RETURN
C
C TEST 3RD ORDER PROJECTION MATRIX
C RESULTS WILL BE MEANINGFUL DEPENDING UPON ERROR IN THE
C 2ND ORDER PROJECTION MATRIX.
C
      WRITE (IOUT,109)
 109  FORMAT ('   TESTING 3ND ORDER PROJECTION MATRIX, P(I,J,K,L):')
      IF (FLAG1) WRITE (IOUT,110)
 110  FORMAT ('      RESULTS MAY BE UNRELIABLE DUE TO ERRORS IN '/
     $        'P(I,J,K).')

C READ IN F3

      REWIND I20
      READ(I20, 1) M,N
      READ(I20, 2) (((F3(I,J,K), K=1,J), J=1,I), I=1,NC)
      DO 120 I = 1, NC
      DO 120 J = 1, I-1
        F3(J,I,J) = F3(J,J,I)
        F3(I,I,J) = F3(J,J,I)
        F3(I,J,I) = F3(I,I,J)
        F3(J,I,I) = F3(I,I,J)
        DO 120 K = 1, J-1
          F3(I,K,J)=F3(I,J,K)
          F3(J,I,J)=F3(I,J,K)
          F3(J,K,I)=F3(I,J,K)
          F3(K,I,J)=F3(I,J,K)
          F3(K,J,I)=F3(I,J,K)
  120     CONTINUE

C COMPUTE Vijk*'s

      DO 150 I = 1, NC
      DO 150 M = 1, NC
      X(I,M)=ZERO
      DO 150 L = 1, NC
 150    X(I,M) = X(I,M) + F2(M,L)*P2(L,I)

      DO 160 I = 1, NC
      DO 160 J = 1, NC
      DO 160 K = 1, NC
      DO 160 M = 1, NC
 160    VP3(I,J,K) = VP3(I,J,K) + X(I,M)*P3(M,J,K) + X(J,M)*P3(M,I,K) +
     $              X(K,M)*P3(M,I,J)


      FLAG2 = .TRUE.
      DO 180 I = 1, NC
      DO 180 J = 1, NC
      DO 180 K = 1, NC

        XX = ZERO
        DO 170 L = 1, NC
  170        XX = XX + F1(K)*P4(L,I,J,K)

        ERR = F3(I,J,K) - VP3(I,J,K) - XX
        IF (ABS(ERR).GT.ERRTOL) THEN
             FLAG2= .FALSE.
             WRITE (IOUT,172) I,J,K,ERR
  172        FORMAT ('      I= ?  J= ',I3,' K= ',I3,' L= ',I3,
     $               'ERROR=',F12.8)
             END IF

  180   CONTINUE

C       MESSAGE IF NO ERRORS
        IF (FLAG2) WRITE (IOUT,54)

        RETURN

        END
