C///////////////////////////////////////////////////////////////////////
      SUBROUTINE MACHZ(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,U,IU,Z,SR)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER R,P,Q,NSX,RP,IOPT(16)
      CHARACTER TYPE*5
      DIMENSION TYPE(NS),IA(NS,6),S(NS),U(NS,1),IU(NS,0:1)
      DIMENSION XA(NAD,3),A(NC,NC),Z(NC,NC,NC,NC),SR(NC,NC,NC,NC)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
    1 FORMAT(/,1X,'NUMERICAL SR(I,J,K,L) AND Z(M,N,P,Q) MATRICES USED',
     $ /, ' FOR SIMPLE INTERNAL COORDINATE',I5)
      NSYM=IOPT(3)
      IF(NSYM.NE.0) THEN
         ISCR=ISCR10
      ELSE
         ISCR=ISCR9
      END IF
      DO 500  R=1,NS
      CALL ZIN(NC,NC,SR,-R,ISCR8)
      CALL ZIN(NC,NSX,Z,R,ISCR8)
      WRITE(IOUT,1) R
 300  CONTINUE
      CALL FILL4A(NC,NC,SR)
      CALL FILL4A(NC,NSX,Z)
      CALL ZOUT(NC,NC,SR,-R,ISCR)
      CALL ZOUT(NC,NSX,Z,R,ISCR)
 500  CONTINUE
 600  IF(NSYM.EQ.0) RETURN
      DO 650  R=1,NSYM
         DO 610  Q=1,NSYM
         DO 610  P=1,Q
         DO 610  N=1,P
         DO 610  M=1,N
 610     Z(M,N,P,Q)=ZERO
         L=IU(R,0)
         DO 620  I=1,L
         RP=IU(R,I)
         CALL ZIN(NC,NSYM,SR,RP,ISCR10)
         W1=U(R,I)
         DO 630  Q=1,NSYM
         DO 630  P=1,Q
         DO 630  N=1,P
         DO 630  M=1,N
 630     Z(M,N,P,Q)=Z(M,N,P,Q)+W1*SR(M,N,P,Q)
 620     CONTINUE
         CALL FILL4A(NC,NSYM,Z)
         CALL ZOUT(NC,NSYM,Z,R,ISCR9)
 650  CONTINUE
C   OPTION
      DO 750  R=1,NSYM
         DO 710  Q=1,NC
         DO 710  P=1,Q
         DO 710  N=1,P
         DO 710  M=1,N
 710     SR(M,N,P,Q)=ZERO
         L=IU(R,0)
         DO 720  I=1,L
         RP=IU(R,I)
         CALL ZIN(NC,NC,Z,-RP,ISCR10)
         W1=U(R,I)
         DO 730  Q=1,NC
         DO 730  P=1,Q
         DO 730  N=1,P
         DO 730  M=1,N
 730     SR(M,N,P,Q)=SR(M,N,P,Q)+W1*Z(M,N,P,Q)
 720     CONTINUE
         CALL FILL4A(NC,NC,SR)
         CALL ZOUT(NC,NC,SR,-R,ISCR9)
 750  CONTINUE
      RETURN
      END
