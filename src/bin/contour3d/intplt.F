      SUBROUTINE INTPLT(K,L,XY,FG,FH,IXY,ierr)
C    interpolation
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/FIG/XX1,XX2,YY1,YY2,S
c
C    linear interpolation
      IF (IXY.EQ.1) GO TO 10
      XYG=XX1+S*K-S
      GO TO 20
10    XYG=YY1+S*K-S
20    A=FG*S
      B=FH-FG
      IF (DABS(B).LT.1.0D-10) THEN
        WRITE (6,*) 'WARNING: B=0   fh=',fh,'fg=',fg
        WRITE (6,*) 'WARNING: B=0    k=', k,' l=', l,'xyg=', xyg
        call fioflu(6)
        ierr = 1
        return
        ENDIF
      XY=XYG+(K-L)*A/B
c
c   quadratic interpolation will be added when necessary.
      ierr = 0
      RETURN
C
      END
