C//////////////////////////////////////////////////////////////////////
      integer function getmas (na,m,a,i,xmass)
      IMPLICIT REAL*8 (A-H,O-Z)
#include <libparse/error.h>
      character a*12
      integer na,m,i
      DIMENSION XMASS(NA)
      PARAMETER(ZERO=0.0D0)
      PARAMETER(W1=1.007825D0, W2=2.014D0, W3A=3.01605D0)
      PARAMETER(W3B=3.01603D0, W4=4.0026D0, W5=5.0125D0,W6=6.01512D0)
      PARAMETER(W7=7.01600D0, W9=9.01218D0,W10=10.0129D0)
      PARAMETER(W11=11.00931D0, W12=12.0000D0,W13=13.00335D0)
      PARAMETER(W14=14.00307D0, W15=15.00011D0,W16=15.99491D0)
      PARAMETER(W17=16.999131D0,W18=17.999160D0)
      PARAMETER(W19=18.99840D0,W20=19.992435D0,W21=20.993843D0)
      PARAMETER(W22=21.991383D0,W23=22.989767D0,W24=23.985042D0)
      PARAMETER(W25=24.985837D0,W26=25.982593D0,W27=26.98153D0)
      PARAMETER(W28=27.976927D0,W29=28.976495D0,W30=29.973770D0)
      PARAMETER(W31=30.973762D0,W32=31.972070D0,W33=32.971456D0)
      PARAMETER(W34=33.967866D0,W35=34.968852D0,W36A=35.967080D0)
      PARAMETER(W36B=35.967545D0,W37=36.965903D0,W38=37.962732D0)
      PARAMETER(W40A=39.962384D0)

      XMASS(M)=-1.0D0

C HYDROGEN
      IF(A(I:I).EQ.'H') THEN
         IF(I.EQ.12) THEN
            XMASS(M)=W1
         ELSE
            J=I+1
            IF(A(J:J).EQ.'1'.OR.A(J:J).EQ.' ') XMASS(M)=W1
            IF(A(J:J).EQ.'2') XMASS(M)=W2
            IF(A(J:J).EQ.'3') XMASS(M)=W3A
         END IF
      END IF
      IF(A(I:I).EQ.'D') THEN
         IF(I.EQ.12) THEN
            XMASS(M)=W2
         ELSE
            J=I+1
            IF(A(J:J).EQ.' ') XMASS(M)=W2
         END IF
      END IF
      IF(A(I:I).EQ.'T') THEN
         IF(I.EQ.12) THEN
            XMASS(M)=W3A
         ELSE
            J=I+1
            IF(A(J:J).EQ.' ') XMASS(M)=W3A
         END IF
      END IF
C HELIUM
      IF(A(I:I+1).EQ.'HE') THEN
         IF(I.EQ.11) THEN
            XMASS(M)=W4
         ELSE
            J=I+2
            IF(A(J:J).EQ.'4'.OR.A(J:J).EQ.' ') XMASS(M)=W4
            IF(A(J:J).EQ.'3') XMASS(M)=W3B
         END IF
      END IF
C LITHIUM
      IF(A(I:I+1).EQ.'LI') THEN
         IF(I.EQ.11) THEN
            XMASS(M)=W7
         ELSE
            J=I+2
            IF(A(J:J).EQ.'7'.OR.A(J:J).EQ.' ') XMASS(M)=W7
            IF(A(J:J).EQ.'6') XMASS(M)=W6
         END IF
      END IF
C BERYLLIUM
      IF(A(I:I+1).EQ.'BE') THEN
         IF(I.EQ.11) THEN
            XMASS(M)=W9
         ELSE
            J=I+2
            IF(A(J:J).EQ.'9'.OR.A(J:J).EQ.' ') XMASS(M)=W9
         END IF
      END IF
C BORON
      IF(A(I:I).EQ.'B') THEN
         IF(I.EQ.12) THEN
            XMASS(M)=W11
         ELSE
            J=I+1
            K=I+2
            IF(A(J:K).EQ.'11'.OR.A(J:J).EQ.' ') XMASS(M)=W11
            IF(A(J:K).EQ.'10') XMASS(M)=W10
         END IF
      END IF
C CARBON
      IF(A(I:I).EQ.'C') THEN
         IF(I.EQ.12) THEN
            XMASS(M)=W12
         ELSE
            J=I+1
            K=I+2
            IF(A(J:K).EQ.'12'.OR.A(J:J).EQ.' ') XMASS(M)=W12
            IF(A(J:K).EQ.'13') XMASS(M)=W13
         END IF
      END IF
C NITROGEN
      IF(A(I:I).EQ.'N') THEN
         IF(I.EQ.12) THEN
            XMASS(M)=W14
         ELSE
            J=I+1
            K=I+2
            IF(A(J:K).EQ.'14'.OR.A(J:J).EQ.' ') XMASS(M)=W14
            IF(A(J:K).EQ.'15') XMASS(M)=W15
         END IF
      END IF
C OXYGEN
      IF(A(I:I).EQ.'O') THEN
         IF(I.EQ.12) THEN
            XMASS(M)=W16
         ELSE
            J=I+1
            K=I+2
            IF(A(J:K).EQ.'16'.OR.A(J:J).EQ.' ') XMASS(M)=W16
            IF(A(J:K).EQ.'17') XMASS(M)=W17
            IF(A(J:K).EQ.'18') XMASS(M)=W18
         END IF
      END IF
C FLUORINE
      IF(A(I:I).EQ.'F') THEN
         IF(I.EQ.12) THEN
            XMASS(M)=W19
         ELSE
            J=I+1
            K=I+2
            IF(A(J:K).EQ.'19'.OR.A(J:J).EQ.' ') XMASS(M)=W19
         END IF
      END IF
C NEON
      IF(A(I:I+1).EQ.'NE') THEN
         IF(I.EQ.11) THEN
            XMASS(M)=W20
         ELSE
            J=I+2
            K=I+3
            IF(A(J:K).EQ.'20'.OR.A(J:J).EQ.' ') XMASS(M)=W20
            IF(A(J:K).EQ.'21') XMASS(M)=W21
            IF(A(J:K).EQ.'22') XMASS(M)=W22
         END IF
      END IF
C SODIUM
      IF(A(I:I+1).EQ.'NA') THEN
         IF(I.EQ.11) THEN
            XMASS(M)=W23
         ELSE
            J=I+2
            K=I+3
            IF(A(J:K).EQ.'23'.OR.A(J:J).EQ.' ') XMASS(M)=W23
         END IF
      END IF
C MAGNESIUM
      IF(A(I:I+1).EQ.'MG') THEN
         IF(I.EQ.11) THEN
            XMASS(M)=W24
         ELSE
            J=I+2
            K=I+3
            IF(A(J:K).EQ.'24'.OR.A(J:J).EQ.' ') XMASS(M)=W24
            IF(A(J:K).EQ.'25') XMASS(M)=W25
            IF(A(J:K).EQ.'26') XMASS(M)=W26
         END IF
      END IF
C ALUMINUM
      IF(A(I:I+1).EQ.'AL') THEN
         IF(I.EQ.11) THEN
            XMASS(M)=W27
         ELSE
            J=I+2
            K=I+3
            IF(A(J:K).EQ.'27'.OR.A(J:J).EQ.' ') XMASS(M)=W27
         END IF
      END IF
C SILICON
      IF(A(I:I+1).EQ.'SI') THEN
         IF(I.EQ.11) THEN
            XMASS(M)=W28
         ELSE
            J=I+2
            K=I+3
            IF(A(J:K).EQ.'28'.OR.A(J:J).EQ.' ') XMASS(M)=W28
            IF(A(J:K).EQ.'29') XMASS(M)=W29
            IF(A(J:K).EQ.'30') XMASS(M)=W30
         END IF
      END IF
C PHOSPHORUS
      IF(A(I:I).EQ.'P') THEN
         IF(I.EQ.12) THEN
            XMASS(M)=W31
         ELSE
            J=I+1
            K=I+2
            IF(A(J:K).EQ.'31'.OR.A(J:J).EQ.' ') XMASS(M)=W31
         END IF
      END IF
C SULFUR
      IF(A(I:I).EQ.'S') THEN
         IF(I.EQ.12) THEN
            XMASS(M)=W32
         ELSE
            J=I+1
            K=I+2
            IF(A(J:K).EQ.'32'.OR.A(J:J).EQ.' ') XMASS(M)=W32
            IF(A(J:K).EQ.'33') XMASS(M)=W33
            IF(A(J:K).EQ.'34') XMASS(M)=W34
            IF(A(J:K).EQ.'36') XMASS(M)=W36A
         END IF
      END IF
C CHLORINE
      IF(A(I:I+1).EQ.'CL') THEN
         IF(I.EQ.11) THEN
            XMASS(M)=W35
         ELSE
            J=I+2
            K=I+3
            IF(A(J:K).EQ.'35'.OR.A(J:J).EQ.' ') XMASS(M)=W35
            IF(A(J:K).EQ.'37') XMASS(M)=W37
         END IF
      END IF
C ARGON
      IF(A(I:I+1).EQ.'AR') THEN
         IF(I.EQ.11) THEN
            XMASS(M)=W40A
         ELSE
            J=I+2
            K=I+3
            IF(A(J:K).EQ.'36') XMASS(M)=W36B
            IF(A(J:K).EQ.'38') XMASS(M)=W38
            IF(A(J:K).EQ.'40'.OR.A(J:J).EQ.' ') XMASS(M)=W40A
         END IF
      END IF

      IF(XMASS(M).LE.ZERO) IFLAG=1
      getmas = iflag
      END 

