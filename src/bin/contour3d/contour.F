c      This copy is good for use before two-way option in implemented.
c************************************************************************
      PROGRAM CONTOUR3D
C    This is a contour3d program was based on the 2d program (~xie/src
C   /contour2), which was copied and modified on 04/24/94 by Yaoming Xie.
c**********************************************************************
c
      implicit real*8 (a-h,o-z)
c grddim is the dimension of the grid used to store the
c       values of the mo
c maxpnt is the maximum number of points in one contour
      integer grddim,maxpnt
      parameter (grddim=1000,maxpnt=800)
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
      common/gauss/zs(500),CSPD(500,5)
      common/nshel/knuc(500),ktype(500),kprim(500),kloc(500),
     1 kmin(500),kmax(500),kng(500),nshell,nprim
c
      COMMON/MO/C(500),DC(200,500)
      COMMON/FIG/AA1,AA2,BB1,BB2,S,OC
clj
      integer errcod, prcntr, prdup, frdc, frdboo, color, frdf, frdi
      integer chrend
      character*80 label
      character*4 orbtyp,oc
      integer orbnum,rotate,splineyn
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
90    FORMAT ('Z coordinate',F10.4)
555   FORMAT(7F10.6)
 500  format(//1x,'not enough memory in routine overlap for s matrix'/)
c
      call psinit('APPEND')
      CALL TSTART(6)

C  PSI 2.0 Input...begin
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
        call fioflu(6)
        call pabort
        call qabort
        endif

      errcod = frdc('ORBTYPE',ORBTYP)
      if (errcod.ne.0) then
        write(6,*) 'ERROR: problem reading ''ORBTYPE'''
        call fioflu(6)
        call pabort
        call qabort
        endif

      errcod = frdc('PLANE',xyz)
      if (errcod.ne.0) then
        write(6,*) 'ERROR: problem reading ''PLANE'''
        call fioflu(6)
        call pabort
        call qabort
        endif

      fi = 1.d0
      errcod = frdf('CONTOURVAL',fi)
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
 1103   continue

      aa1 = aa1 - 1.5d0
      aa2 = aa2 + 1.5d0
      bb1 = bb1 - 1.5d0
      bb2 = bb2 + 1.5d0
      c0 = 0.d0
      errcod = frdf('X0',aa1)
      errcod = frdf('X1',aa2)
      errcod = frdf('Y0',bb1)
      errcod = frdf('Y1',bb2)
      errcod = frdf('Z0',cc1)
      errcod = frdf('Z1',cc2)
      ccd=0.20
      errcod = frdf('ZD',ccd)
      phi=10
      errcod = frdf('PHI',phi)
      theta=20
      errcod = frdf('THETA',theta)
      if (phi.lt.-90.0.or.phi.gt.90.0.or.theta.lt.-90.0.
     &  or.theta.gt.90.0)  then
        write(6,*) 'ERROR: ''PHI'' and ''THETA'' are required here'
        write(6,*) '  to be the values between -90.0 and 90.0 degrees.'
        call fioflu(6)
        call qabort
        endif
      c0=cc1-ccd
c Use x0, x1, y0, y1 to determine the default gridsize
      if (aa2-aa1 .gt. bb2-bb1) then
        s = (aa2-aa1)/120
      else
        s = (bb2-bb1)/120
        endif
      errcod = frdf('GRIDSIZE',s)

C  PSI 2.0 Input...end

c Write out the input parameters
 801  format(
     &  /,
     &  '       The 3-Dimension Plotting Program for',/,
     &  '      Molecular Orbital and Electron Density',/,
     &  '                    04/28/94',/,
     &  '                by Yaoming Xie',/,
     &  /,
     &  'Input Parameters:',/,
     &  '  LABEL      = ',A,/,
     &  '  CONTOURVAL = ',F12.6,/,
     &  '  ORBTYPE    = ',A,/,
     &  '  ORBNUM     = ',I5,/,
     &  '  PLANE      = ',A,/,
     &  '  PHI        = ',F12.6,/,
     &  '  THETA      = ',F12.6,/,
     &  '  X0         = ',F12.6,/,
     &  '  X1         = ',F12.6,/,
     &  '  Y0         = ',F12.6,/,
     &  '  Y1         = ',F12.6,/,
     &  '  Z0         = ',F12.6,/,
     &  '  Z1         = ',F12.6,/,
     &  '  ZD         = ',F12.6,/,
     &  '  GRIDSIZE   = ',F12.6,/,
     &  '  COLOR      = ',I5,/,
     &  '  SPLINE     = ',I5
     &  )
 802  format(
     &  '  ROTATE     = ',I5
     &  )

      write(6,801) label(1:chrend(label)),fi,orbtyp,orbnum,xyz,phi,
     &      theta,aa1,aa2,bb1,bb2,cc1,cc2,ccd,s,color,splineyn

      if (rotate.ne.-1) write(6,802) rotate
      call fioflu(6)

      IF (XYZ.NE.'XY'.AND.XYZ.NE.'YX'.AND.XYZ.NE.'XZ'.AND.XYZ.NE.'ZX'
     1 .AND.XYZ.NE.'YZ'.AND.XYZ.NE.'ZY') then
        WRITE (6,*) 'ERROR: invalid ''PLANE'''
        call fioflu(6)
        call qabort
        endif

      if (orbnum.gt.nbfao) then 
        write(6,*) 'ERROR: ''ORBNUM'' is greater than the number'
        write(6,*) '      of AO basis functions.'
        call fioflu(6)
        call qabort
        endif

      conv=180.0D0/3.1415926535D0
      theta=theta/conv
      phi=phi/conv
      ct=dcos(theta)
      cp=dcos(phi)
      sp=dsqrt(1-cp*cp)
      if(phi.lt.0.0) sp=-sp
      st=dsqrt(1-ct*ct)
      if(theta.lt.0.0) st=-st
      do 101 i = 1, nbfao
 101    c(i) = 0.d0

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
          call fioflu(6)
          call qabort
          endif
        do 102 i = 1, orbnum
          read(8,'(I5)',end=901) ijunk
          if (ijunk.ne.i) then
            write(6,*) 'WARNING: orbitals in cmo.dat or lmo.dat were'
            write(6,*) '        not in the expected order'
            call fioflu(6)
            endif
          read(8,555,end=901) (c(j),j=1,nbfao)
          do 106 j=1,nbfao
 106        dc(i,j) = c(j)
 102      continue
        endif

      goto 905
 901  write(6,*) 'ERROR: unexpected end of file in either cmo.dat or'
      write(6,*) '       lmo.dat'
      call fioflu(6)
      call qabort
      stop

clj contour no longer uses old input
c     CALL LOCATE(5,'# CNTRMP #',IERR)
c     READ (5,20) XYZ
c     READ (5,30) S,AA1,AA2,BB1,BB2,C0
c     READ (5,30) (C(I),I=1,NBFAO)
c     WRITE (6,40) XYZ
c     write (6,50) AA1,AA2,BB1,BB2,C0,S
c     write (6,55)
clj End of old stype input

 905  continue
      if (prcntr('IS_ON BRIEF').eq.0) then
        if (orbtyp.eq.'DENS') then
          do 107 i=1,orbnum
  107       write (6,30) (dc(i,j),j=1,nbfao)
          call fioflu(6)
        else
          write (6,30) (C(I),I=1,NBFAO)
          call fioflu(6)
          endif
        endif

c Initialize the postscript output file
      call initps(aa1,bb1,aa2,bb2,label,color,
     &            orbtyp,orbnum,rotate)

cyx   the next line has been moved to after lable 1000.
CYX   call drawat(geom,chgs,natom,aa1,bb1,aa2,bb2,c0,xyz,ct,cp)
 
      if (prcntr('IS_ON BRIEF').eq.0) then
        WRITE(6,1) NPRIM
        DO 103 I=1,NPRIM
  103     WRITE(6,2) I,ZS(I),CSPD(1,I),CSPD(2,I),CSPD(3,I)
        endif

c     WRITE(6,3) NSHELL
c     DO 104 I=1,NSHELL
c 104 WRITE(6,4) I,KNUC(I),KTYPE(I),KPRIM(I),KLOC(I),KMIN(I),KMAX(I),
c    1 KNG(I)

      WRITE(6,5) NATOM
      DO 105 I=1,NATOM
  105   WRITE(6,6) I,GEOM(1,I),geom(2,I),geom(3,I),CHGS(I)
      call fioflu(6)
 
c  calculate psi values of all grid points within the boundary
110   c0=c0+ccd
      if (c0.gt.cc2) goto 1000
      write (6,90) c0
      call fioflu(6)
      fi=dabs(fi)
      AA=AA2-AA1
      BB=BB2-BB1
      LL=AA/S+0.1
      MM=BB/S+0.1
      LL1=LL+1
      MM1=MM+1


      IF (XYZ.EQ.'XY') THEN
        DO 120 L=1,LL1
          A=AA1+S*L-S
          DO 120 M=1,MM1
            B=BB1+S*M-S
            if (l.gt.grddim.or.m.gt.grddim) goto 902
            F(L,M)=PSI(A,B,C0,orbtyp,orbnum)
120       CONTINUE
      else IF (XYZ.EQ.'YX') THEN
        DO 121 L=1,LL1
          A=AA1+S*L-S
          DO 121 M=1,MM1
            B=BB1+S*M-S
            if (l.gt.grddim.or.m.gt.grddim) goto 902
            F(L,M)=PSI(B,A,C0,orbtyp,orbnum)
121       CONTINUE
      else IF (XYZ.EQ.'XZ') THEN
        DO 122 L=1,LL1
          A=AA1+S*L-S
          DO 122 M=1,MM1
            B=BB1+S*M-S
            if (l.gt.grddim.or.m.gt.grddim) goto 902
            F(L,M)=PSI(A,C0,B,orbtyp,orbnum)
122       CONTINUE
      else IF (XYZ.EQ.'ZX') THEN
        DO 123 L=1,LL1
          A=AA1+S*L-S
          DO 123 M=1,MM1
            B=BB1+S*M-S
            if (l.gt.grddim.or.m.gt.grddim) goto 902
            F(L,M)=PSI(B,C0,A,orbtyp,orbnum)
123       CONTINUE
      else IF (XYZ.EQ.'YZ') THEN
        DO 124 L=1,LL1
          A=AA1+S*L-S
          DO 124 M=1,MM1
            B=BB1+S*M-S
            if (l.gt.grddim.or.m.gt.grddim) goto 902
            F(L,M)=PSI(C0,A,B,orbtyp,orbnum)
124       CONTINUE
      else IF (XYZ.EQ.'ZY') THEN
        DO 115 L=1,LL1
          A=AA1+S*L-S
          DO 115 M=1,MM1
            B=BB1+S*M-S
            if (l.gt.grddim.or.m.gt.grddim) goto 902
            F(L,M)=PSI(C0,B,A,orbtyp,orbnum)
115       CONTINUE
        ENDIF

      goto 904
 902  write(6,*) 'ERROR: maxgrd has been exceed'
      write(6,*) '  make X1-X0 or Y1-Y0 smaller or GRIDSIZE larger'
      call fioflu(6)
      call qabort
      stop

 904  continue
      if (prcntr('IS_ON BRIEF').eq.0) THEN
        WRITE (6,*)
        WRITE (6,*) 'PSI VALUES OF THE GRID'
        endif

      M1=1
      M3=MM+1
125   M2=MM+1
      IF(M3.GE.22) M2=M1+21
      if (prcntr('IS_ON BRIEF').eq.0) then
        WRITE (6,*)
        DO 130 L=1,LL+1
130       WRITE (6,60) (F(L,M),M=M1,M2)
        endif
      M1=M2+1
      M3=MM+1-M2
      IF (M3.NE.0) GO TO 125

c   search for and record all possible pairs (in x-direction) psi value
c   is in between.
140   continue

      if (prcntr('IS_ON BRIEF').eq.0) then
        WRITE (6,*)
        WRITE (6,80) FI
        endif

clj   WRITE (7,80) FI
160   MA=0
      DO 170 J=1,MM+1
        DO 170 I=2,LL+1
          FG=F(I,J)-FI
          FH=F(I-1,J)-FI
          IF (FH*FG.GT.0.0D0) GO TO 170
          IJ=I*10000+J
Ceeb-JX   IJ=I*1000000+J*10000+(I-1)*100+J
          MA=MA+1
          if (ma.gt.maxpnt/2) then
            write(6,*) 'ERROR: too many points in a contour'
            write(6,*) '      increase ''GRIDSIZE'''
            call fioflu(6)
            call qabort
            endif
          IIJJ(MA)=IJ
170     CONTINUE
      IF (MA.EQ.0) then
        if (fi.le.0.0) goto 110
        fi = -fi
        goto 140
      endif
 
C    circle beginning

      if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,*) 'A NEW LINE BEGINNING'

cYX   call startc(fi)
      call startz(fi,c0)
Ceeb-YX  I0=IIJJ(1)/1000000+0.1
Ceeb-YX  J0=IIJJ(1)/10000-I0*100+0.1
      I0=IIJJ(1)/10000+0.1
      J0=IIJJ(1)-I0*10000+0.1
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
        call fioflu(6)
        call qabort
        endif

      X=XY

cyx   if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,70) X,Y, I0,J0,M0,N0
clj   WRITE (7,70) X,Y
cyx   call point(x,y)

      xr=ct*x-st*c0
      yr=cp*y-sp*ct*c0-sp*st*x

c     xr=ct*x-st*sp*y-st*cp*c0
c     yr=cp*y-sp*c0

      IPNT=IPNT+1
      XAY(1,IPNT)=Xr
      XAY(2,IPNT)=Yr
Ceeb-YX  IJ1=I0*1000000+J0*10000+M0*100+N0
      IJ1=I0*10000+J0
200   IM=M0-I0+(N0-J0)*KQ
Ceeb  For possible endless loop
      ieeb_err = 0
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
      IF (FH*FG.GE.0.0D0) GO TO 260
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
        call fioflu(6)
        call qabort
        endif

      IF (IXY.EQ.1) Y=XY
      IF (IXY.EQ.2) X=XY

cyx   if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,70) X,Y, I0,J0,M1,N1
clj   WRITE (7,70) X,Y
cyx   call point(x,y)

      xr=ct*x-st*c0
      yr=cp*y-sp*ct*c0-sp*st*x

c     xr=ct*x-st*sp*y-st*cp*c0
c     yr=cp*y-sp*c0

      IPNT=IPNT+1
      XAY(1,IPNT)=Xr
      XAY(2,IPNT)=Yr
Ceeb-YX  IJ2=I0*1000000+J0*10000+M1*100+N1
Ceeb-YX  IJ3=M1*1000000+N1*10000+I0*100+J0
      IJ2=I0*10000+J0
      IJ3=M1*10000+N1
      IF (IBD.EQ.2) GO TO 215
Ceeb-YX  IF (IJ2.EQ.IJ1) GO TO 320
Ceeb-YX  IF (IJ3.EQ.IJ1) GO TO 320
      IF (I0.GT.M1.AND.IJ2.EQ.IJ1) GO TO 320
      IF (M1.GT.I0.AND.IJ3.EQ.IJ1) GO TO 320
215   DO 230 I=1,MB
        IF (IABS(IIJJ(I)).EQ.IJ2) GO TO 220
        IF (IABS(IIJJ(I)).EQ.IJ3) GO TO 220
        GO TO 230
220     IIJJ(I)=IIJJ(MB)
        MB=MB-1
230     CONTINUE
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
C270   IM=I0-M0+(N0-J0)*KQ
C      JN=J0-N0-(M0-I0)*KQ
C      IF (IM.GT.1) IM=1
C      IF (IM.LT.-1) IM=-1
C      IF (JN.GT.1) JN=1
C      IF (JN.LT.-1) JN=-1
C      M2=I0+IM
C      N2=J0+JN
C      IF (M2.EQ.0.OR.N2.EQ.0) GO TO 310
C      IF (M2.EQ.LL+2.OR.N2.EQ.MM+2) GO TO 310
C      FH=F(M2,N2)-FI
C      FG=F(M0,N0)-FI
C      IF (FH*FG.LT.0.0D0) GO TO 280
C      IF (DABS(F(I0,J0)).LT.DABS(F(M0,N0))) GO TO 290
280   if (i0.eq.m0.and.j0.eq.n0) then
        write(6,*) 'Possible endless loop...'
        write(6,*) 'Perhaps the contour line passes through',
     &    ' a grid point...'
        write(6,*) 'Will attempt a recovery...'
        ieeb_err = ieeb_err + 1
        if (ieeb_err.eq.6) then
          write(6,*) ' '
          write(6,*) 'Endless loop error!'
          write(6,*) ' '
          write(6,*) 'Aborting Program!'
          call fioflu(6)
          call qabort
          endif
        I0=M1
        J0=N1
        M1=M0
        N1=N0
        goto 210
        endif

      I0=M0
      J0=N0
      KQ=KQ*(-1)
      GO TO 210
C290   M1=M2
C      N1=N2
C      GO TO 210
300   IF (MB.EQ.0) then
        if (fi.le.0.0) goto 110
        fi = -fi
        GO TO 140
      endif
Ceeb-JX  I0=IIJJ(1)/1000000+0.1
Ceeb-JX  J0=IIJJ(1)/10000-I0*100+0.1
      I0=IIJJ(1)/10000+0.1
      J0=IIJJ(1)-I0*10000
      M0=I0-1
      N0=J0
      IIJJ(1)=IIJJ(MB)
      MB=MB-1

      if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,*) 'A NEW LINE BEGINNING'

cYX   call startc(fi)
      call startz(fi,c0)
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
330     continue
      call fioflu(6)
c
cyx 09/02/92   added for spline interpolation
      if (ipnt.lt.100.and.splineyn.eq.1) go to 340
c
      if (ipnt.gt.800) then
        write(6,*) 'ipnt = ', ipnt
        write(6,*) ' '
        write(6,*) ' '
        write(6,*) 'Error!  ipnt is greater than 800'
        write(6,*) ' '
        write(6,*) 'Please increase the maximum number ',
     1    'of points allowed'
        write(6,*) 'or decrease the grid size or turn ',
     1    'spline off (if on)'
        write(6,*) ' '
        write(6,*) ' '
        call fioflu(6)
        call qabort
      endif

      call points(xay,ipnt)

      if (prcntr('IS_ON BRIEF').eq.0) WRITE (6,*) 'END OF A LINE'

clj   WRITE (7,*)
      go to 360
c
cyx 09/02/92   added for spline interpolation
340   do 350 i=1,ipnt
        xay1(i)=xay(1,i)
        xay2(i)=xay(2,i)
350     continue
      ipntn=ipnt-1
      if (oc.eq.'open') ipntn=ipnt
      call spline (ipntn,xay1,xay2)
c
360   call closec()
      GO TO 300
1000  CONTINUE
      call drawat(geom,chgs,natom,aa1,bb1,aa2,bb2,c0,xyz,ct,cp,st,sp)
      call doneps()
      call tstop(6)
      call psdone()
      STOP
C
      END
