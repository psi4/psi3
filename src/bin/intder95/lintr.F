C     ////////////////////////////////////////////////////////////
      SUBROUTINE LINTR(NC,NSX,NSY,IOPT,A,XS,F1,F2,F3,F4,V,
     $            NMM,NRUN,IFLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER P,Q,QQ
      DIMENSION F1(NC),F2(NC,NC),F3(NC,NC,NC),F4(NC,NC,NC,NC)
      DIMENSION IOPT(30),A(NC,NC),XS(NC,NC),V(NMM),DIP(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
      COMMON /DIPOLE/DIP,QQ
 1    FORMAT(20X,3F20.10)
 2    FORMAT(8F10.7)
 3    FORMAT(3F20.10)
 4    FORMAT(F16.12)
 5    FORMAT(6F12.6)
 6    FORMAT(2I5)
 7    FORMAT(/' LINEAR TRANSFORMATION CONTRIBUTIONS TO INTERNAL'/
     $   ' COORDINATE FORCE CONSTANTS:')
 8    FORMAT(/' LINEAR TRANSFORMATION CONTRIBUTIONS TO CARTESIAN'/
     $   ' COORDINATE FORCE CONSTANTS:')
 9    FORMAT(/' QUADRATIC FORCE CONSTANTS')
 10   FORMAT(/' CUBIC FORCE CONSTANTS')
 12   FORMAT(/' QUARTIC FORCE CONSTANTS')
 13   FORMAT(3F20.10)
 14   FORMAT(20X,F20.10)
 15   FORMAT(2I5,3F20.10)
 16   FORMAT(/1X,'TOTAL CHARGE AND DIPOLE MOMENT')
 17   FORMAT(I5,3F20.10)
      IFLAG=0
      NA=IOPT(1)
      NDER=IOPT(4)
      NEQ=IOPT(5)
      NINV=IOPT(7)
CSA
      NNINV=0
      IF (ABS(NINV).NE.3) GO TO 25
         NEQ=0
         NNINV=NINV
         NINV=0
 25   NVEC=IOPT(13)
      INP1=I11
      INP2=I15
      INP3=I20
      INP4=I24
      INP5=I17
      IF(ABS(NNINV).NE.3) THEN
C        CF1=8.238857606D0
C        CF2=15.56918890D0
C        CF3=29.42151140D0
C        CF4=55.59861458D0
         CF1=HART/BOHR
         CF2=CF1/BOHR
         CF3=CF2/BOHR
         CF4=CF3/BOHR
      ELSE
         CF1=ONE
         CF2=ONE
         CF3=ONE
         CF4=ONE
      END IF
CSA
 30   IF(NINV.GE.1) THEN
      INP1=I12
      INP2=I16
      INP3=I21
      INP4=I25
      INP5=I18
      CF1=ONE
      CF2=ONE
      CF3=ONE
      CF4=ONE
      END IF
      IF(NEQ.EQ.0) GO TO 120
      IF(NVEC.EQ.1) THEN
        IF(NRUN.LE.1) THEN
        READ(INP5,15) M,QQ,(DIP(I),I=1,3)
        WRITE(IOUT,16)
        WRITE(IOUT,17) QQ,(DIP(I),I=1,3)
          IF(M.NE.NA) THEN
             IFLAG=5
             RETURN
          END IF
        END IF
        READ(INP5,13) (V(IK),IK=1,NSY)
        GO TO 122
      END IF
      IF(NINV.GE.1) THEN
        READ(INP1,6) M
        IF(M.NE.NA) THEN
             IFLAG=1
             RETURN
        END IF
        READ(INP1,14)(V(IK),IK=1,NSY)
      ELSE
        READ(INP1,1)(V(IK),IK=1,NSY)
      END IF
 122  DO 110  M=1,NSX
 110  F1(M)=ZERO
      DO 115  I=1,NSY
      V(I)=V(I)*CF1
      DO 115  M=1,NSX
 115  F1(M)=F1(M)+V(I)*A(I,M)
 120  IF(NDER.LE.1) GO TO 1000
      NNC=NSY*NSY
      NM=NMM
      NR=NNC/NMM
      NL=NNC-NMM*NR
      READ(INP2,6) M,N
      IF(M.NE.NA) THEN
           IFLAG=2
           RETURN
      END IF
      DO 125  N=1,NSX
      DO 125  I=1,NSY
 125  XS(I,N)=ZERO
      KK=0
      DO 128  II=1,NR+1
      IF(II.EQ.NR+1) NM=NL
      READ(INP2,3) (V(IK),IK=1,NM)
           DO 128  IK=1,NM
           V(IK)=V(IK)*CF2
           KK=KK+1
           J=(KK-1)/NSY+1
           I=KK-NSY*(J-1)
           DO 130  N=1,NSX
 130       XS(I,N)=XS(I,N)+A(J,N)*V(IK)
 128  CONTINUE
      DO 131  N=1,NSX
      DO 131  M=1,NSX
      XX=ZERO
      DO 132  I=1,NSY
 132  XX=XX+A(I,M)*XS(I,N)
 131  F2(M,N)=XX
      DO 134  M=2,NSX
      DO 134  N=1,M-1
      F2(M,N)=(F2(M,N)+F2(N,M))/TWO
 134  F2(N,M)=F2(M,N)
      IF(NDER.LE.2) GO TO 1000
      NM=NMM
      NNC=(NSY*(NSY+1)*(NSY+2))/6
      NR=NNC/NMM
      NL=NNC-NMM*NR
      READ(INP3,6) M,N
      IF(M.NE.NA.OR.N.NE.NNC) THEN
           IFLAG=3
           RETURN
      END IF
      DO 135  P=1,NSX
      DO 135  I=1,NSY
      DO 135  J=1,I
 135  F3(I,J,P)=ZERO
      I=0
      J=0
      K=0
      DO 138  II=1,NR+1
      IF(II.EQ.NR+1) NM=NL
      READ(INP3,3) (V(IK),IK=1,NM)
          DO 138  IK=1,NM
          V(IK)=V(IK)*CF3
          IF(K.LT.J) THEN
             K=K+1
             GO TO 139
          END IF
          IF(J.LT.I) THEN
             J=J+1
             K=1
             GO TO 139
          END IF
          I=I+1
          J=1
          K=1
 139  CONTINUE
          IF(I.NE.J) THEN
             IF(J.NE.K) THEN
                DO 140  P=1,NSX
                F3(I,J,P)=V(IK)*A(K,P)+F3(I,J,P)
                F3(I,K,P)=V(IK)*A(J,P)+F3(I,K,P)
 140            F3(J,K,P)=V(IK)*A(I,P)+F3(J,K,P)
             ELSE
                DO 141  P=1,NSX
                F3(I,J,P)=V(IK)*A(J,P)+F3(I,J,P)
 141            F3(J,J,P)=V(IK)*A(I,P)+F3(J,J,P)
             END IF
          ELSE
             IF(J.NE.K) THEN
                DO 142  P=1,NSX
                F3(I,I,P)=V(IK)*A(K,P)+F3(I,I,P)
 142            F3(I,K,P)=V(IK)*A(I,P)+F3(I,K,P)
             ELSE
                DO 143  P=1,NSX
 143            F3(I,I,P)=V(IK)*A(I,P)+F3(I,I,P)
             END IF
          END IF
 138  CONTINUE
      REWIND ISCR5
      WRITE(ISCR5,3)(((F3(I,J,P),P=1,NSX),J=1,I),I=1,NSY)
      DO 144  P=1,NSX
      DO 144  N=1,P
      DO 144  I=1,NSY
 144  F3(I,N,P)=ZERO
      NND=NSY*(NSY+1)*NSX/2
      NR=NND/NMM
      NL=NND-NMM*NR
      NM=NMM
      I=1
      J=1
      P=0
      REWIND ISCR5
      DO 146  II=1,NR+1
      IF(II.EQ.NR+1) NM=NL
      READ(ISCR5,3) (V(IK),IK=1,NM)
          DO 146  IK=1,NM
          IF(P.LT.NSX) THEN
             P=P+1
             GO TO 147
          END IF
          IF(J.LT.I) THEN
             J=J+1
             P=1
             GO TO 147
          END IF
          I=I+1
          J=1
          P=1
 147  CONTINUE
          IF(I.NE.J) THEN
                DO 148  N=1,P
                F3(I,N,P)=V(IK)*A(J,N)+F3(I,N,P)
 148            F3(J,N,P)=V(IK)*A(I,N)+F3(J,N,P)
          ELSE
                DO 149  N=1,P
 149            F3(I,N,P)=V(IK)*A(I,N)+F3(I,N,P)
          END IF
 146  CONTINUE
      REWIND ISCR5
      WRITE(ISCR5,3)(((F3(I,N,P),I=1,NSY),N=1,P),P=1,NSX)
      DO 150  P=1,NSX
      DO 150  N=1,P
      DO 150  M=1,N
 150  F3(M,N,P)=ZERO
      NND=NSY*NSX*(NSX+1)/2
      NR=NND/NMM
      NL=NND-NMM*NR
      NM=NMM
      I=0
      N=1
      P=1
      REWIND ISCR5
      DO 152  II=1,NR+1
      IF(II.EQ.NR+1) NM=NL
      READ(ISCR5,3) (V(IK),IK=1,NM)
          DO 152  IK=1,NM
          IF(I.LT.NSY) THEN
             I=I+1
             GO TO 153
          END IF
          IF(N.LT.P) THEN
             N=N+1
             I=1
             GO TO 153
          END IF
          P=P+1
          N=1
          I=1
 153  CONTINUE
          DO 154  M=1,N
 154      F3(M,N,P)=V(IK)*A(I,M)+F3(M,N,P)
 152  CONTINUE
      REWIND ISCR5
      CALL FILL3A(NC,NSX,F3)
      IF(NDER.LE.3) GO TO 1000
      NM=NMM
      NNC=(NSY*(NSY+1)*(NSY+2)*(NSY+3))/24
      NR=NNC/NMM
      NL=NNC-NMM*NR
      READ(INP4,6) M,N
      IF(M.NE.NA.OR.N.NE.NNC) THEN
           IFLAG=4
           RETURN
      END IF
      DO 170  Q=1,NSX
      DO 170  I=1,NSY
      DO 170  J=1,I
      DO 170  K=1,J
 170  F4(I,J,K,Q)=ZERO
      I=0
      J=0
      K=0
      L=0
      DO 178  II=1,NR+1
      IF(II.EQ.NR+1) NM=NL
      READ(INP4,3) (V(IK),IK=1,NM)
          DO 178  IK=1,NM
          V(IK)=V(IK)*CF4
          IF(L.LT.K) THEN
             L=L+1
             GO TO 180
          END IF
          IF(K.LT.J) THEN
             K=K+1
             L=1
             GO TO 180
          END IF
          IF(J.LT.I) THEN
             J=J+1
             K=1
             L=1
             GO TO 180
          END IF
             I=I+1
             J=1
             K=1
             L=1
 180      CONTINUE
          IF(I.NE.J) THEN
             IF(J.NE.K) THEN
                IF(K.NE.L) THEN
                DO 182  Q=1,NSX
                F4(I,J,K,Q)=V(IK)*A(L,Q)+F4(I,J,K,Q)
                F4(I,J,L,Q)=V(IK)*A(K,Q)+F4(I,J,L,Q)
                F4(I,K,L,Q)=V(IK)*A(J,Q)+F4(I,K,L,Q)
 182            F4(J,K,L,Q)=V(IK)*A(I,Q)+F4(J,K,L,Q)
                ELSE
                DO 184  Q=1,NSX
                F4(I,J,K,Q)=V(IK)*A(K,Q)+F4(I,J,K,Q)
                F4(I,K,K,Q)=V(IK)*A(J,Q)+F4(I,K,K,Q)
 184            F4(J,K,K,Q)=V(IK)*A(I,Q)+F4(J,K,K,Q)
                END IF
             ELSE
                IF(K.NE.L) THEN
                DO 186  Q=1,NSX
                F4(I,J,J,Q)=V(IK)*A(L,Q)+F4(I,J,J,Q)
                F4(I,J,L,Q)=V(IK)*A(J,Q)+F4(I,J,L,Q)
 186            F4(J,J,L,Q)=V(IK)*A(I,Q)+F4(J,J,L,Q)
                ELSE
                DO 188  Q=1,NSX
                F4(I,J,J,Q)=V(IK)*A(J,Q)+F4(I,J,J,Q)
 188            F4(J,J,J,Q)=V(IK)*A(I,Q)+F4(J,J,J,Q)
                END IF
             END IF
          ELSE
             IF(J.NE.K) THEN
                IF(K.NE.L) THEN
                DO 190  Q=1,NSX
                F4(I,I,K,Q)=V(IK)*A(L,Q)+F4(I,I,K,Q)
                F4(I,I,L,Q)=V(IK)*A(K,Q)+F4(I,I,L,Q)
 190            F4(I,K,L,Q)=V(IK)*A(I,Q)+F4(I,K,L,Q)
                ELSE
                DO 192  Q=1,NSX
                F4(I,I,K,Q)=V(IK)*A(K,Q)+F4(I,I,K,Q)
 192            F4(I,K,K,Q)=V(IK)*A(I,Q)+F4(I,K,K,Q)
                END IF
             ELSE
                IF(K.NE.L) THEN
                DO 194  Q=1,NSX
                F4(I,I,I,Q)=V(IK)*A(L,Q)+F4(I,I,I,Q)
 194            F4(I,I,L,Q)=V(IK)*A(I,Q)+F4(I,I,L,Q)
                ELSE
                DO 196  Q=1,NSX
 196            F4(I,I,I,Q)=V(IK)*A(I,Q)+F4(I,I,I,Q)
                END IF
             END IF
          END IF
 178  CONTINUE
      REWIND ISCR5
      WRITE(ISCR5,3) ((((F4(I,J,K,Q),Q=1,NSX),K=1,J),J=1,I),I=1,NSY)
      DO 198  Q=1,NSX
      DO 198  P=1,Q
      DO 198  I=1,NSY
      DO 198  J=1,I
 198  F4(I,J,P,Q)=ZERO
      NND=NSY*(NSY+1)*(NSY+2)*NSX/6
      NR=NND/NMM
      NL=NND-NMM*NR
      NM=NMM
      I=1
      J=1
      K=1
      Q=0
      REWIND ISCR5
      DO 200  II=1,NR+1
      IF(II.EQ.NR+1) NM=NL
      READ(ISCR5,3) (V(IK),IK=1,NM)
          DO 200  IK=1,NM
          IF(Q.LT.NSX) THEN
            Q=Q+1
            GO TO 202
          END IF
          IF(K.LT.J) THEN
            K=K+1
            Q=1
            GO TO 202
          END IF
          IF(J.LT.I) THEN
            J=J+1
            K=1
            Q=1
            GO TO 202
          END IF
            I=I+1
            J=1
            K=1
            Q=1
 202      CONTINUE
          IF(I.NE.J) THEN
             IF(J.NE.K) THEN
                DO 204  P=1,Q
                F4(I,J,P,Q)=V(IK)*A(K,P)+F4(I,J,P,Q)
                F4(I,K,P,Q)=V(IK)*A(J,P)+F4(I,K,P,Q)
 204            F4(J,K,P,Q)=V(IK)*A(I,P)+F4(J,K,P,Q)
             ELSE
                DO 206  P=1,Q
                F4(I,J,P,Q)=V(IK)*A(J,P)+F4(I,J,P,Q)
 206            F4(J,J,P,Q)=V(IK)*A(I,P)+F4(J,J,P,Q)
             END IF
          ELSE
             IF(J.NE.K) THEN
                DO 208  P=1,Q
                F4(I,I,P,Q)=V(IK)*A(K,P)+F4(I,I,P,Q)
 208            F4(I,K,P,Q)=V(IK)*A(I,P)+F4(I,K,P,Q)
             ELSE
                DO 210  P=1,Q
 210            F4(I,I,P,Q)=V(IK)*A(I,P)+F4(I,I,P,Q)
             END IF
          END IF
 200  CONTINUE
      REWIND ISCR5
      WRITE(ISCR5,3) ((((F4(I,J,P,Q),P=1,Q),Q=1,NSX),J=1,I),I=1,NSY)
      DO 212  Q=1,NSX
      DO 212  P=1,Q
      DO 212  N=1,P
      DO 212  I=1,NSY
 212  F4(I,N,P,Q)=ZERO
      NND=NSX*(NSX+1)*NSY*(NSY+1)/4
      NR=NND/NMM
      NL=NND-NMM*NR
      NM=NMM
      P=0
      Q=0
      I=1
      J=1
      REWIND ISCR5
      DO 214  II=1,NR+1
      IF(II.EQ.NR+1) NM=NL
      READ(ISCR5,3) (V(IK),IK=1,NM)
          DO 214  IK=1,NM
          IF(P.LT.Q) THEN
            P=P+1
            GO TO 216
          END IF
          IF(Q.LT.NSX) THEN
            Q=Q+1
            P=1
            GO TO 216
          END IF
          IF(J.LT.I) THEN
            J=J+1
            P=1
            Q=1
            GO TO 216
          END IF
            I=I+1
            J=1
            P=1
            Q=1
 216      CONTINUE
          IF(I.NE.J) THEN
                DO 218  N=1,P
                F4(I,N,P,Q)=V(IK)*A(J,N)+F4(I,N,P,Q)
 218            F4(J,N,P,Q)=V(IK)*A(I,N)+F4(J,N,P,Q)
          ELSE
                DO 220  N=1,P
 220            F4(I,N,P,Q)=V(IK)*A(I,N)+F4(I,N,P,Q)
          END IF
 214  CONTINUE
      REWIND ISCR5
      WRITE(ISCR5,3) ((((F4(I,N,P,Q),I=1,NSY),N=1,P),P=1,Q),Q=1,NSX)
      DO 222  Q=1,NSX
      DO 222  P=1,Q
      DO 222  N=1,P
      DO 222  M=1,N
 222  F4(M,N,P,Q)=ZERO
      NND=NSX*(NSX+1)*(NSX+2)*NSY/6
      NR=NND/NMM
      NL=NND-NMM*NR
      NM=NMM
      N=1
      P=1
      Q=1
      I=0
      REWIND ISCR5
      DO 224  II=1,NR+1
      IF(II.EQ.NR+1) NM=NL
      READ(ISCR5,3) (V(IK),IK=1,NM)
          DO 224  IK=1,NM
          IF(I.LT.NSY) THEN
            I=I+1
            GO TO 226
          END IF
          IF(N.LT.P) THEN
            N=N+1
            I=1
            GO TO 226
          END IF
          IF(P.LT.Q) THEN
            P=P+1
            N=1
            I=1
            GO TO 226
          END IF
            Q=Q+1
            P=1
            N=1
            I=1
 226      CONTINUE
          DO 228  M=1,N
 228      F4(M,N,P,Q)=V(IK)*A(I,M)+F4(M,N,P,Q)
 224  CONTINUE
      REWIND ISCR5
      CALL FILL4A(NC,NSX,F4)
C     ****************************************
 1000 IF(LPRT(1,NPRT).LT.5) RETURN
      IF(NDER.LE.1) RETURN
      IF(NINV.LE.0) THEN
         WRITE(IOUT,7)
      ELSE
         WRITE(IOUT,8)
      END IF
      WRITE(IOUT,9)
      DO 1015  M=1,NSX
 1015 WRITE(IOUT,5) (F2(M,N),N=1,NSX)
      IF(NDER.LE.2) RETURN
      WRITE(IOUT,10)
      DO 1020  M=1,NSX
      WRITE(IOUT,*) 'M=',M
      DO 1020  N=1,NSX
 1020 WRITE(IOUT,5) (F3(M,N,P),P=1,NSX)
      IF(NDER.LE.3) RETURN
      WRITE(IOUT,12)
      DO 1030  M=1,NSX
      DO 1030  N=1,NSX
      WRITE(IOUT,*) 'M=',M,'N=',N
      DO 1030  P=1,NSX
 1030 WRITE(IOUT,5) (F4(M,N,P,Q),Q=1,NSX)
      RETURN
      END
