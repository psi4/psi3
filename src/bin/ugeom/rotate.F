C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C This subroutine rotates the molecule to a certain direction.
C It reads 3 numbers from input.dat (finding the keyword: geomro).    
C This subroutine is called by ugeom. 
C       Created by Seung-Joon Kim and modefied by Yaoming Xie
c
c     installed to psi area by croag on sep 7 1993.  no credit taken
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine rotate(natom,coord)
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (matom=200)
      integer errcod,fcnt,frdia1
      DIMENSION ORG(3),COORD(3,MATOM),noxy(3)         
      DATA ZERO / 0.0D+00 /
    2 FORMAT(//,2X,'ROTATION INDICES'/                         
     2 2X,' N0  = ',I5/                                
     3 2X,' NX  = ',I5/                                 
     4 2X,' NXY = ',I5)                             
    3 FORMAT(3F20.10)                                   
    4 FORMAT(10X,3I5)                                
    5 FORMAT(/'  NEW COORDINATES BEFORE ROTATION ' /)
    6 FORMAT(/'  NEW COORDINATES AFTER ROTATION ' /)
Cyx
        if(natom.gt.matom) then
          write(6,*) ' '
          write(6,*) 'ERROR: Ugeom cant rotate over 200 atoms.  Sorry.'
          call qabort
        endif
c
      errcod = fcnt('GEOMRO',iro)
        if(errcod.eq.0.and.iro.ne.3) then
         write(6,*) ' '
         write(6,*) 'GEOMRO is a 3 element vector you dimwit!  RTFM.'
         call pabort
         call qabort
        endif
C
        errcod = frdia1('GEOMRO',noxy,3)
        if(errcod.ne.0) goto 109
C
      n0   = noxy(1)
      nx   = noxy(2)
      nxy  = noxy(3)
      WRITE(6,2) N0, NX, NXY
        if (n0.gt.natom.or.nx.gt.natom.or.nxy.gt.natom) goto 109
        if (n0.le.0.or.nx.le.0.or.nxy.le.0) goto 109
        if (n0.eq.nx.or.n0.eq.nxy.or.nx.eq.nxy) goto 109
C
      WRITE(6,5)
      DO 101 I=1, NATOM
      WRITE(6,3) COORD(1,I),COORD(2,I),COORD(3,I)
  101 CONTINUE
C
      DO 103 I=1,3           
      ORG(I)=COORD(I,N0)  
  103 CONTINUE             
C
      DO 104 I=1,NATOM    
      DO 104 J=1,3       
      COORD(J,I)=COORD(J,I)-ORG(J) 
  104 CONTINUE                    
C
      NABC=1
      NA=NX
  200 XNX=COORD(1,NA)
      YNX=COORD(2,NA)
      IF (NABC .EQ. 2) YNX=COORD(3, NA)
      RR=DSQRT(XNX**2 + YNX**2)
      COS=XNX/RR
      SIN=YNX/RR
C
      DO 105 I=1, NATOM
      XX=COORD(1,I)
      YY=COORD(2,I)
      IF (NABC .EQ. 2) GO TO 300
      COORD(1,I)=XX*COS + YY*SIN
      COORD(2,I)=XX*(-SIN) + YY*COS
      GO TO 105
  300 YY=COORD(3,I)
      COORD(1,I)=XX*COS + YY*SIN
      COORD(3,I)=XX*(-SIN) + YY*COS
  105 CONTINUE
      IF (NABC .GE. 2) GO TO 400
      NABC=NABC+1
      GO TO 200
  400 CONTINUE
C
      NA=NXY
      XNX=COORD(2,NA)
      YNX=COORD(3,NA)
      RR=DSQRT(XNX**2 + YNX**2)
      COS=XNX/RR
      SIN=YNX/RR
      DO 106 I=1, NATOM
      YY=COORD(2,I)
      ZZ=COORD(3,I)
      COORD(2,I)=YY*COS + ZZ*SIN
      COORD(3,I)=YY*(-SIN) + ZZ*COS
  106 CONTINUE
C
      WRITE(6,6)
      DO 108 I=1, NATOM
      WRITE(6,3) COORD(1,I),COORD(2,I),COORD(3,I)
  108 CONTINUE
      WRITE(6,*)
c
        goto 110
  109   write(6,*) ' '
        write(6,*)'NO geometry rotation was performed.'
  110 return                                  
      END                                    
