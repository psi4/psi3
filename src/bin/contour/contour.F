      PROGRAM CONTOUR
C    This is a contour map program for a molecular orbital, created by
c Yaoming Xie in 03/91.
c    Unlike CONTOUR3, overlap integrals are not checked.
c    Unlike CONTOUR4, boundary is read from INPUT file, and the plane
c other than XY can be treated so long as it is perpendicular to any
c axis (x, y, or z).
c**********************************************************************
c  Since CONTOURMAX and CONTOURMIN didn't work well, this version makes the
c    defaults of Psi values are -1.0 to -0.1 (with steps 0.05), 0.0, and 0.1
c    to 1.0 (with steps 0.05).          
c                                                   07/11/91   Y. Xie
c**********************************************************************
c  Atoms are indicated by big dots, other than cross's.  (see drawps.f)
c                                                   12/19/91   Y. Xie
c**********************************************************************
c  Fix the bug when the contours intersect the boundary.         
c  1. Do such a contour line twice (set IBD=1 or 2), and not print it at
c    first time.  At the second time, it would be from border to border, 
c    completely and non-redundantly.
c  2. Use array XAY to store data of all points.  Print them only when
c    a perfect line has been completed.   It will also be useful for the
c    adaption of interpolation later.        
c  3. Use subroutine POINTS, to replace subroutine POINT.  (see drawps.f)
c                                                08/26/92   Y. Xie
c**********************************************************************
c  Spline interpolation is used to make the contour lines smoother. 
c                                         09/02/92   Y. Xie
c  Spline interpolation for open lines, as well as for closed ones. 
c                                         09/20/92   Y. Xie
c**********************************************************************
c  f, and g functions added to psi.f .                               
c                                         12/01/93   Y. Xie
c  5d, 7f, 9g can be drawn, since program 'local' changed too.
c                                         12/01/93   Y. Xie
c**********************************************************************
c           **********************************************
C  It reads basis functions from FILE30 (matrix BF)
C  It reads geometry from FILE30  (matrix H)
C  It reads MO coefficients from INPUT, which is in terms of AO's.
c
c Current input format (by C Janssen)
c  string ORBTYPE is AO, CMO, or LMO (for CMO or LMO local must be run)
c  integer ORBNUM is the orbital for which output is generated
c  string PLANE reads XYZ (as described below)
c  boolean COLOR generates postscript output that requires a color
c                output device if TRUE (default: NO)
c  real GRIDSIZE is S (as below) (default: .05)
c  string LABEL is a label that appears in the plot
c  boolean ROTATE tells whether or not to rotate the drawing by
c                 90 degrees.  If it is not given, then the
c                 program decides.
c  real X0 is AA1 (as below) (default: smallest x - 1.5 Bohr)
c  real X1 is AA2 (as below) (default: largest  x + 1.5 Bohr)
c  real Y0 is BB1 (as below) (default: smallest y - 1.5 Bohr)
c  real Y1 is BB2 (as below) (default: largest  y + 1.5 Bohr)
c  real Z is C (as below)    (default: 0)
C
C  OLD Input format:
C  # CNTRMP #
C  1. FORMAT (A2)
C     XYZ     =  XY :  X is the abscissa, and Y is the ordinate.
C             =  YX :  Y is the abscissa, and X is the ordinate.
C             =  XZ :  X is the abscissa, and Z is the ordinate.
C             =  ZX :  Z is the abscissa, and X is the ordinate.
C             =  YZ :  Y is the abscissa, and Z is the ordinate.
C             =  ZY :  Z is the abscissa, and Y is the ordinate.
 
C  2. FORMAT (6F10.6)
C     S   : The step length (in bohr) for the grid of the map.
C     AA1 : left boundary of the abscissa.
C     AA2 : right boundary of the abscissa.
C     BB1 : bottom boundary of the ordinate.
C     BB2 : top boundary of the ordinate.
C     C   : the value of the axis perpendicular to the considered plane.
C
C  3. FORMAT (10F10.6)
C     (C(I),I=1,NBFAO)    coefficients of a MO. (NBFAO is the number
C                  of basis functions.)
c
C  4. FORMAT (F10.6)
c     FI      The psi value for a contour line.
c
c ... Repeat line 4 as many times as you want. End this with FI>90.0.
c
c**********************************************************************
c
      implicit real*8 (a-h,o-z)
c grddim is the dimension of the grid used to store the
c       values of the mo
c maxpnt is the maximum number of points in one contour
      integer grddim,maxpnt
      parameter (grddim=1000,maxpnt=400)
c
      DIMENSION F(grddim,grddim),IIJJ(maxpnt),XAY(2,maxpnt),
     & xay1(maxpnt),xay2(maxpnt)
      CHARACTER*2 XYZ
c
c  max number of atoms 112
c  number of shells is 200
c
      common/coord/chgs(112),GEOM(3,112),ict(112),ns(112),ks(112)
      common/funcs/natom,nbfao
      common/gauss/zs(200),CSPD(5,200)
      common/nshel/knuc(200),ktype(200),kprim(200),kloc(200),
     1 kmin(200),kmax(200),kng(200),nshell,nprim
c
      COMMON/MO/C(200),DC(50,200)
      COMMON/FIG/AA1,AA2,BB1,BB2,S,OC
clj
      integer errcod, prcntr, prdup, frdc, frdboo, color, frdf, frdi
      integer chrend
      character*80 label
      character*4 orbtyp,oc
      integer orbnum,rotate,splineyn
      real*8 cstep,cmin,cmax
C
    1 FORMAT(/,2X,' NPRIM IS ',I5)
    2 FORMAT(2X,I5,4F20.10)
    3 FORMAT(/,2X,' NSHELL IS ',I5)
    4 FORMAT(2X,8I5)
    5 FORMAT(/,2X,'The Number of Atoms is',I5)
    6 FORMAT(2X,I5,3F16.10,F8.4,5X,2I5)
20    FORMAT (A2)
30    FORMAT (10F10.6)
40    FORMAT ('This contour map is on the plane of ',A2)
50    FORMAT ('abscissa boundary:              from',f6.2,' to',f6.2,
     1        /,'ordinate boundary:              from',f6.2,' to',f6.2,
     2        /,'value of the third coordinate:',f6.2,
     3        /,'increment of the grid        :',f6.2)
55    FORMAT (/,'Coefficients of the MO')
60    FORMAT (22F6.4)
70    FORMAT (2F10.5,4I5)
80    FORMAT ('CONTOUR LINE VALUE',F10.4)
555   FORMAT(7F10.6)
 500  format(//1x,'not enough memory in routine overlap for s matrix'/)
c
      call psinit('APPEND')
      CALL TSTART(6)

      errcod = prcntr('ON BRIEF')
      errcod = prdup()

      label = ' '
      errcod = frdc('LABEL',label)

      color = 0
      errcod = frdboo('COLOR',color)

      rotate = -1
      errcod = frdboo('ROTATE',rotate)

      splineyn = 0
      errcod = frdboo('SPLINE',splineyn)

      errcod = frdi('ORBNUM',ORBNUM)
      if (errcod.ne.0) then
        write(6,*) 'ERROR: problem reading ''ORBNUM'''
        call pabort
        call qabort
        endif

      errcod = frdc('ORBTYPE',ORBTYP)
      if (errcod.ne.0) then
        write(6,*) 'ERROR: problem reading ''ORBTYPE'''
        call pabort
        call qabort
        endif

      errcod = frdc('PLANE',xyz)
      if (errcod.ne.0) then
        write(6,*) 'ERROR: problem reading ''PLANE'''
        call pabort
        call qabort
        endif

      cmax = 1.0
      cmin = -1.0
      cstep = .05
      errcod = frdf('CONTOURMAX',cmax)
      errcod = frdf('CONTOURMIN',cmin)
      errcod = frdf('CONTOURSTEP',cstep)

      fi = cmin - cstep
c
c
c   .........................  read in basis set information
      call basset

c Set up the default aa, bb

      aa1 = 1.d20
      aa2 = -1.d20
      bb1 = 1.d20
      bb2 = -1.d20
      do 1103 i = 1, natom
        if (xyz(1:1).eq.'X') then
          if (geom(1,i).lt.aa1) aa1 = geom(1,i)
          if (geom(1,i).gt.aa2) aa2 = geom(1,i)
        else if (xyz(1:1).eq.'Y') then
          if (geom(2,i).lt.aa1) aa1 = geom(2,i)
          if (geom(2,i).gt.aa2) aa2 = geom(2,i)
        else
          if (geom(3,i).lt.aa1) aa1 = geom(3,i)
          if (geom(3,i).gt.aa2) aa2 = geom(3,i)
          endif
 1103   continue
      do 1104 i = 1, natom
        if (xyz(2:2).eq.'X') then
          if (geom(1,i).lt.bb1) bb1 = geom(1,i)
          if (geom(1,i).gt.bb2) bb2 = geom(1,i)
        else if (xyz(2:2).eq.'Y') then
          if (geom(2,i).lt.bb1) bb1 = geom(2,i)
          if (geom(2,i).gt.bb2) bb2 = geom(2,i)
        else
          if (geom(3,i).lt.bb1) bb1 = geom(3,i)
          if (geom(3,i).gt.bb2) bb2 = geom(3,i)
          endif
 1104   continue

      aa1 = aa1 - 1.5d0
      aa2 = aa2 + 1.5d0
      bb1 = bb1 - 1.5d0
      bb2 = bb2 + 1.5d0
      c0 = 0.d0
      errcod = frdf('X0',aa1)
      errcod = frdf('X1',aa2)
      errcod = frdf('Y0',bb1)
      errcod = frdf('Y1',bb2)
      errcod = frdf('Z',c0)

c Use x0, x1, y0, y1 to determine the default gridsize
      if (aa2-aa1 .gt. bb2-bb1) then
        s = (aa2-aa1)/120
      else
        s = (bb2-bb1)/120
        endif
      errcod = frdf('GRIDSIZE',s)

c Write out the input parameters
 801  format(
     &  /,
     &  '  The Contour Plotting Program for Molecular Orbitals',/,
     &  '                    by Yaoming Xie',/,
     &  '           Postscript output by Curtis Janssen',/,
     &  /,
     &  'Input Parameters:',/,
     &  '  LABEL      = ',A,/,
     &  '  CONTOURMAX = ',F12.6,/,
     &  '  CONTOURMIN = ',F12.6,/,
     &  '  CONTOURSTEP= ',F12.6,/,
     &  '  ORBTYPE    = ',A,/,
     &  '  ORBNUM     = ',I5,/,
     &  '  PLANE      = ',A,/,
     &  '  X0         = ',F12.6,/,
     &  '  X1         = ',F12.6,/,
     &  '  Y0         = ',F12.6,/,
     &  '  Y1         = ',F12.6,/,
     &  '  Z          = ',F12.6,/,
     &  '  GRIDSIZE   = ',F12.6,/,
     &  '  COLOR      = ',I5,/,
     &  '  SPLINE     = ',I5
     &  )
 802  format(
     &  '  ROTATE     = ',I5
     &  )
      write(6,801) label(1:chrend(label)),cmax,cmin,cstep,
     &        orbtyp,orbnum,xyz,aa1,aa2,bb1,bb2,c0,s,color,splineyn
      if (rotate.ne.-1) write(6,802) rotate


      if (orbnum.gt.nbfao) then 
        write(6,*) 'ERROR: ''ORBNUM'' is greater than the number'
        write(6,*) '      of AO basis functions.'
        call qabort
        endif

      do 101 i = 1, nbfao
        c(i) = 0.d0
 101    continue
c Get the coefficients of the atomic orbitals
      if (orbtyp.eq.'AO') then
        c(orbnum) = 1.d0
      else
        if (orbtyp.eq.'CMO'.or.orbtyp.eq.'DENS') then
          call fioopn(8,'cmo.dat','USER','READ')
        else if (orbtyp.eq.'LMO') then
          call fioopn(8,'lmo.dat','USER','READ')
        else
          write(6,*) 'ERROR: ''ORBTYPE'' is invalid'
          write(6,*) '      it must be one of AO, CMO, LMO, or DENS'
          call qabort
          endif
        do 102 i = 1, orbnum
          read(8,'(I5)',end=901) ijunk
          if (ijunk.ne.i) then
            write(6,*) 'WARNING: orbitals in cmo.dat or lmo.dat were'
            write(6,*) '        not in the expected order'
            endif
          read(8,555,end=901) (c(j),j=1,nbfao)
            do 106 j=1,nbfao
 106        dc(i,j) = c(j)
 102      continue
        endif

c
clj contour no longer uses old input
c     CALL LOCATE(5,'# CNTRMP #',IERR)
c     READ (5,20) XYZ
c     READ (5,30) S,AA1,AA2,BB1,BB2,C0
c     READ (5,30) (C(I),I=1,NBFAO)
c     WRITE (6,40) XYZ
c     write (6,50) AA1,AA2,BB1,BB2,C0,S
c     write (6,55)
clj End of old stype input

      if (prcntr('IS_ON BRIEF').eq.0) then
        if (orbtyp.eq.'DENS') then
           do 107 i=1,orbnum
  107      write (6,30) (dc(i,j),j=1,nbfao)
        else
           write (6,30) (C(I),I=1,NBFAO)
        endif
      endif
c Initialize the postscript output file
      call initps(aa1,bb1,aa2,bb2,label,color,cmin,cmax,cstep,
     &            orbtyp,orbnum,rotate)
      call drawat(geom,chgs,natom,aa1,bb1,aa2,bb2,c0,xyz)
c
      if (prcntr('IS_ON BRIEF').eq.0) WRITE(6,1) NPRIM
      DO 103 I=1,NPRIM
        if (prcntr('IS_ON BRIEF').eq.0)
     &    WRITE(6,2) I,ZS(I),CSPD(1,I),CSPD(2,I),CSPD(3,I)
  103 continue
c     WRITE(6,3) NSHELL
c     DO 104 I=1,NSHELL
c 104 WRITE(6,4) I,KNUC(I),KTYPE(I),KPRIM(I),KLOC(I),KMIN(I),KMAX(I),
c    1 KNG(I)
      WRITE(6,5) NATOM
      DO 105 I=1,NATOM
  105 WRITE(6,6) I,GEOM(1,I),geom(2,I),geom(3,I),CHGS(I)
c
c  calculate psi values of all grid points within the boundary
      AA=AA2-AA1
      BB=BB2-BB1
      LL=AA/S+0.1
      MM=BB/S+0.1
      DO 120 L=1,LL+1
      A=AA1+S*L-S
      DO 120 M=1,MM+1
      B=BB1+S*M-S
      IF (XYZ.NE.'XY'.AND.XYZ.NE.'YX'.AND.XYZ.NE.'XZ'.AND.XYZ.NE.'ZX'
     1 .AND.XYZ.NE.'YZ'.AND.XYZ.NE.'ZY') GO TO 999
      IF (XYZ.EQ.'XY') THEN
      X=A
      Y=B
      Z=C0
      ENDIF
      IF (XYZ.EQ.'YX') THEN
      Y=A
      X=B
      Z=C0
      ENDIF
      IF (XYZ.EQ.'XZ') THEN
      X=A
      Z=B
      Y=C0
      ENDIF
      IF (XYZ.EQ.'ZX') THEN
      Z=A
      X=B
      Y=C0
      ENDIF
      IF (XYZ.EQ.'YZ') THEN
      Y=A
      Z=B
      X=C0
      ENDIF
      IF (XYZ.EQ.'ZY') THEN
      Z=A
      Y=B
      X=C0
      ENDIF
      if (l.gt.grddim.or.m.gt.grddim) then
        write(6,*) 'ERROR: maxgrd has been exceed'
        write(6,*) '     make X1-X0 or Y1-Y0 smaller'
        write(6,*) '     or GRIDSIZE larger'
        call qabort
        endif
      F(L,M)=PSI(X,Y,Z,orbtyp,orbnum)
120   CONTINUE
      if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,*)
      if (prcntr('IS_ON BRIEF').eq.0)
     &  WRITE (6,*) 'PSI VALUES OF THE GRID'
      M1=1
      M3=MM+1
125   M2=MM+1
      IF(M3.GE.22) M2=M1+21
      if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,*)
      DO 130 L=1,LL+1
130   if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,60) (F(L,M),M=M1,M2)
      M1=M2+1
      M3=MM+1-M2
      IF (M3.NE.0) GO TO 125
c
c   search for and record all possible pairs (in x-direction) psi value
c   is in between.
140   continue
clj   READ (5,30,END=1000)  FI
      fi = fi + cstep
cyx  07/11/91  the following 2 statements added.
      if (fi.eq.0.05) goto 140
      if (fi.eq.-0.05) goto 140
      if (fi.gt.cmax) goto 1000
clj 150   IF (FI.GT.90.0) GO TO 1000
      if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,*)
      if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,80) FI
clj   WRITE (7,80) FI
160   MA=0
      DO 170 J=1,MM+1
      DO 170 I=2,LL+1
      FG=F(I,J)-FI
      FH=F(I-1,J)-FI
      IF (FH*FG.GT.0.0D0) GO TO 170
      IJ=I*1000000+J*10000+(I-1)*100+J
      MA=MA+1
      if (ma.gt.maxpnt/2) then
        write(6,*) 'ERROR: too many points in a contour'
        write(6,*) '      increase ''GRIDSIZE'''
        call qabort
        endif
      IIJJ(MA)=IJ
170   CONTINUE
      IF (MA.EQ.0) GO TO 140
C
C    circle beginning
      if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,*) 'A NEW LINE BEGINNING'
      call startc(fi)
      I0=IIJJ(1)/1000000+0.1
      J0=IIJJ(1)/10000-I0*100+0.1
      M0=I0-1
      N0=J0
C     WRITE (6,*) I0,J0,M0,N0
      KQ=1
      IIJJ(1)=IIJJ(MA)
      MB=MA-1
180   continue
      oc='clos'
      IPNT=0
      IBD=1
      FG=F(I0,J0)-FI
      FH=F(M0,N0)-FI
      X=AA1+I0*S-S
      Y=BB1+J0*S-S
      IXY=2
      CALL INTPLT(I0,M0,XY,FG,FH,IXY,ierr)
      if (ierr.ne.0) then
        write(6,*) 'fg,fh,fi=',fg,fh,fi
        write(6,*) 'i0,m0=',i0,m0
        write(6,*) 'x,y = ', x, y
        call qabort
        endif
      X=XY
cyx   if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,70) X,Y, I0,J0,M0,N0
clj   WRITE (7,70) X,Y
cyx   call point(x,y)
      IPNT=IPNT+1
      XAY(1,IPNT)=X
      XAY(2,IPNT)=Y
      IJ1=I0*1000000+J0*10000+M0*100+N0
200   IM=M0-I0+(N0-J0)*KQ
      JN=N0-J0-(M0-I0)*KQ
      IF (IM.GT.1) IM=1
      IF (IM.LT.-1) IM=-1
      IF (JN.GT.1) JN=1
      IF (JN.LT.-1) JN=-1
      M1=I0+IM
      N1=J0+JN
      IF (M1.EQ.0.OR.N1.EQ.0) GO TO 310
      IF (M1.EQ.LL+2.OR.N1.EQ.MM+2) GO TO 310
210   FH=F(M1,N1)-FI
      FG=F(I0,J0)-FI
      IF (FH*FG.GT.0.0D0) GO TO 260
      IF (I0.NE.M1.AND.J0.NE.N1) GO TO 250
      X=AA1+I0*S-S
      Y=BB1+J0*S-S
      IF (I0.EQ.M1) THEN
      IXY=1
      K=J0
      L=N1
      ENDIF
      IF (J0.EQ.N1) THEN
      IXY=2
      K=I0
      L=M1
      ENDIF
      CALL INTPLT(K,L,XY,FG,FH,IXY,ierr)
      if (ierr.ne.0) then
        write(6,*) 'fg,fh,fi=',fg,fh,fi
        write(6,*) 'k,l=',k,l
        write(6,*) 'x,y = ', x, y
        call qabort
        endif
      IF (IXY.EQ.1) Y=XY
      IF (IXY.EQ.2) X=XY
cyx   if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,70) X,Y, I0,J0,M1,N1
clj   WRITE (7,70) X,Y
cyx   call point(x,y)
      IPNT=IPNT+1
      XAY(1,IPNT)=X
      XAY(2,IPNT)=Y
      IJ2=I0*1000000+J0*10000+M1*100+N1
      IJ3=M1*1000000+N1*10000+I0*100+J0
      IF (IBD.EQ.2) GO TO 215
      IF (IJ2.EQ.IJ1) GO TO 320
      IF (IJ3.EQ.IJ1) GO TO 320
215   DO 230 I=1,MB
      IF (IABS(IIJJ(I)).EQ.IJ2) GO TO 220
      IF (IABS(IIJJ(I)).EQ.IJ3) GO TO 220
      GO TO 230
220   IIJJ(I)=IIJJ(MB)
      MB=MB-1
230   CONTINUE
250   M0=M1
      N0=N1
      GO TO 200
260   IF (X.NE.AA1+S*M0-S) GO TO 270
      IF (Y.NE.BB1+S*N0-S) GO TO 270
      I0=M1
      J0=N1
      M1=M0
      N1=N0
      WRITE (6,*) 'WARNING: CONTOUR LINE MAY PASS THOUGH A GRID POINT.'
      GO TO 210
270   CONTINUE
C    The following paragrah was to solve the trouble that two contour
C  lines may pass through the same side.  But not in use currently.
C 270   IM=I0-M0+(N0-J0)*KQ
C     JN=J0-N0-(M0-I0)*KQ
C     IF (IM.GT.1) IM=1
C     IF (IM.LT.-1) IM=-1
C     IF (JN.GT.1) JN=1
C     IF (JN.LT.-1) JN=-1
C     M2=I0+IM
C     N2=J0+JN
C     IF (M2.EQ.0.OR.N2.EQ.0) GO TO 310
C     IF (M2.EQ.LL+2.OR.N2.EQ.MM+2) GO TO 310
C     FH=F(M2,N2)-FI
C     FG=F(M0,N0)-FI
C     IF (FH*FG.LT.0.0D0) GO TO 280
C     IF (DABS(F(I0,J0)).LT.DABS(F(M0,N0))) GO TO 290
280   I0=M0
      J0=N0
      KQ=KQ*(-1)
      GO TO 210
C290   M1=M2
C     N1=N2
C     GO TO 210
300   IF (MB.EQ.0) GO TO 140
      I0=IIJJ(1)/1000000+0.1
      J0=IIJJ(1)/10000-I0*100+0.1
      M0=I0-1
      N0=J0
      IIJJ(1)=IIJJ(MB)
      MB=MB-1
      if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,*) 'A NEW LINE BEGINNING'
      call startc(fi)
      GO TO 180
310   continue
      oc='open'
      if (prcntr('IS_ON BRIEF').eq.0)
     &  WRITE (6,*) 'THIS IS AN OPEN CONTOUR LINE'
clj   WRITE (7,*) 'END OF A LINE: OUT OF BOUNDARY'
clj   WRITE (7,*)
cyx   call bound()
cyx   go to 300
      IF (IBD.EQ.2) GO TO 320
      KQ=KQ*(-1)
      IPNT=0
      IBD=2
      GO TO 200
320   continue
      DO 330 I=1,IPNT
      WRITE (6,70) (XAY(J,I),J=1,2)
330   continue
c
cyx 09/02/92   added for spline interpolation
      if (ipnt.lt.50.and.splineyn.eq.1) go to 340
c
      call points(xay,ipnt)
      if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,*) 'END OF A LINE'
clj   WRITE (7,*)
      go to 360
c
cyx 09/02/92   added for spline interpolation
340   do 350 i=1,ipnt
      xay1(i)=xay(1,i)
      xay2(i)=xay(2,i)
350   continue
      ipntn=ipnt-1
      if (oc.eq.'open') ipntn=ipnt
      call spline (ipntn,xay1,xay2)
c
360   call closec()
      GO TO 300
1000  CONTINUE
      call doneps()
      call tstop(6)
      call psdone()
      STOP
 999  WRITE (6,*) 'ERROR: invalid ''PLANE'''
      call qabort
 901  continue
      write(6,*) 'ERROR: unexpected end of file in either cmo.dat or'
      write(6,*) '       lmo.dat'
      call qabort
C
      END
