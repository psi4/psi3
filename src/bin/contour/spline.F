      SUBROUTINE SPLINE(N,X,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N+1),Y(N+1),H(100),A(100),B(100),C(100),D(100),
     & E(100),F(100),P(100),Q(100)
   10 FORMAT(/,2X,10F10.5)
C
      x(n+1)=x(1)
      y(n+1)=y(1)
      DO 50 I=1,N
50    H(I) = DSQRT((X(I)-X(I+1))**2+(Y(I)-Y(I+1))**2)
      H(N+1)=H(1)
      DO 80 I=1,N
      W=H(I)
      T=H(I+1)
      P(I)=T/(W+T)
      Q(I)=T/(W+T)
      B(I)=1-P(I)
      E(I)=1-P(I)
      C(I)=2.0D0
      F(I)=2.0D0
      W=P(I)/W
      T=B(I)/T
      IF(I.LT.N) GO TO 60
      J=2
      GO TO 70
60    J=I+2
70    A(I)=3.0D0*(W*(X(I+1)-X(I))+T*(X(J)-X(I+1)))
      D(I)=3.0D0*(W*(Y(I+1)-Y(I))+T*(Y(J)-Y(I+1)))
80    CONTINUE
      CALL CYCTRI(N,P,C,B,A)
      CALL CYCTRI(N,Q,F,E,D)
      DO 120 I=1,N
      J=I-1
      IF (J.EQ.0) J=N
      T=1.0D0/H(I)
      W=2.0D0*T
      B(I)=W*(A(J)+2*A(I)-3.0D0*(X(I+1)-X(I))*T)
      E(I)=W*(D(J)+2*D(I)-3.0D0*(Y(I+1)-Y(I))*T)
      IF (J.EQ.N) GO TO 120
      C(J)=T*(B(I)-B(J))
      F(J)=T*(E(I)-E(J))
120   CONTINUE
      T=1.0D0/H(1)
      C(N)=T*(B(1)-B(N))
      F(N)=T*(E(1)-E(N))
      CALL CURVE(N,X,Y,A,B,C,D,E,F,H)
      RETURN
C
C
      END
      SUBROUTINE CYCTRI(N,G,AK,AL,AM)
C  to solve the cycling-tri-diagonal simultaneous equations.
      implicit real*8 (a-h,o-z)
      DIMENSION G(N),AK(N),AL(N),AM(N)
      ONE=1.0D0 
      V=ONE/AK(1)
      G(1)=G(1)*V
      AL(1)=AL(1)*V
      AM(1)=AM(1)*V
      DO 50 I=2,N-1
      II=I-1
      V=1/(AK(I)-G(I)*AL(II))
      AL(I)=AL(I)*V
      AM(I)=(AM(I)-G(I)*AM(II))*V
      G(I)=-G(I)*G(II)*V
      AM(N)=AM(N)-AM(II)*AL(N)
      AK(N)=AK(N)-AL(N)*G(II)
      AL(N)=-AL(II)*AL(N)
50    CONTINUE
      AK(N)=AK(N)-(AL(N)+G(N))*(AL(N-1)+G(N-1))
      AM(N)=(AM(N)-(AL(N)+G(N))*AM(N-1))/AK(N)
      DO 60 I=N-1,1,-1
60    AM(I)=AM(I)-AL(I)*AM(I+1)-G(I)*AM(N)
      RETURN
C
      END
      SUBROUTINE CURVE(N,X,Y,A,B,C,D,E,F,H)
      implicit real*8 (a-h,o-z)
      COMMON/FIG/AA1,AA2,BB1,BB2,S,OC
      character*4 oc
      DIMENSION X(N+1),Y(N+1),H(N+1),A(N),B(N),C(N),D(N),E(N),F(N),
     & XAY(2,100)
   10 FORMAT(2X,10F10.5)
      II=0
      NN=N
      if (oc.eq.'open') nn=n-1
      DO 50 I=1,NN
      NP=5
      IF (N.GT.20) NP=3
      ST=H(i)/NP
      L=I-1
      if (i.eq.1) l=n
      DO 40 J=1,NP
      STJ=ST*J-ST
      II=II+1
      XAY(1,II)=X(i)+STJ*A(l)+STJ*STJ*B(l)*5.0D-1+STJ*STJ*STJ*C(l)/6
      XAY(2,II)=Y(i)+STJ*D(l)+STJ*STJ*E(l)*5.0D-1+STJ*STJ*STJ*F(l)/6
40    CONTINUE
      write (6,*) ('spline results')
      write (6,10) (xay(1,m),m=ii-np+1,ii)
      write (6,10) (xay(2,m),m=ii-np+1,ii)
50    CONTINUE
      if (oc.eq.'open') then
      call points(xay,ii)
      else
      XAY(1,II+1)=XAY(1,1)
      XAY(2,II+1)=XAY(2,1)
      call points(xay,ii+1)
      endif
      RETURN
      END
