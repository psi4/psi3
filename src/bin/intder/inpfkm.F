C////////////////////////////////////////////////////////////////////
      SUBROUTINE INPFKM(NC,NDER,NEQ,NS,F1,F2,F3,F4)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER P,Q,R,S
      DIMENSION F1(NC),F2(NC,NC),F3(NC,NC,NC),F4(NC,NC,NC,NC)
      PARAMETER(ZERO=0.0D0)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
 1    FORMAT(5I5)
 2    FORMAT(3F20.10)
 3    FORMAT(20X,F20.10)
 5    FORMAT(20X,F20.10)
      NA=NC/3
      ider=9
      call ffile(ider,'ider',0)
      rewind ider
      IF(NEQ.EQ.0) GO TO 110
           DO 108  I=1,NS
 108       F1(I)=ZERO
      IF(NDER.LE.1) GO TO 190
 110  DO 112  I=1,NS
      DO 112  J=1,NS
 112  F2(I,J)=ZERO
      IF(NDER.LE.2) GO TO 190
      DO 114  I=1,NS
      DO 114  J=1,I
      DO 114  K=1,J
 114  F3(I,J,K)=ZERO
      IF(NDER.LE.3) GO TO 190
      DO 116  I=1,NS
      DO 116  J=1,I
      DO 116  K=1,J
      DO 116  L=1,K
 116  F4(I,J,K,L)=ZERO
C
 190  IF(NEQ.EQ.0) GO TO 200
 192  READ(ider,1) M
      IF(M.EQ.0) GO TO 200
      BACKSPACE ider
      READ(ider,5) F1(M)
      GO TO 192
 200  READ(ider,1) M,N
      IF(M.EQ.0) GO TO 202
      BACKSPACE ider
      READ(ider,5) F2(M,N)
      GO TO 200
 202  IF(NDER.LE.2) GO TO 212
 204  READ(ider,1) M,N,P
      IF(M.EQ.0) GO TO 206
      BACKSPACE ider
      READ(ider,5) F3(M,N,P)
      GO TO 204
 206  IF(NDER.LE.3) GO TO 212
 208  READ(ider,1) M,N,P,Q
      IF(M.EQ.0) GO TO 212
      BACKSPACE ider
      READ(ider,5) F4(M,N,P,Q)
      GO TO 208
 212  CONTINUE
C
      IF(NEQ.NE.0) THEN
      WRITE(I12,1) NA
      WRITE(I12,3) (F1(I),I=1,NS)
      REWIND I12
      END IF
      DO 214  I=1,NS
      DO 214  J=1,I
 214  F2(J,I)=F2(I,J)
      WRITE (I16,1) NA,NA*6
      WRITE (I16,2) ((F2(I,J),J=1,NS),I=1,NS)
      REWIND I16
      IF(NDER.LE.2) RETURN
      WRITE (I21,1) NA,(NS*(NS+1)*(NS+2))/6
      WRITE (I21,2) (((F3(I,J,K),K=1,J),J=1,I),I=1,NS)
      REWIND I21
      IF(NDER.LE.3) RETURN
      WRITE (I25,1) NA,(NS*(NS+1)*(NS+2)*(NS+3))/24
      WRITE (I25,2) ((((F4(I,J,K,L),L=1,K),K=1,J),J=1,I),I=1,NS)
      REWIND I25
      RETURN
      END
