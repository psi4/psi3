C*** March 1995 version ***)
c
c  psi 2.0 style input parsing added 4/25/95
c  by justin t. fermann & c. david sherrill
c  everything we added is in lower case 
c
      PROGRAM INTDER95
      IMPLICIT REAL*8 (A-H,O-Z)
      character*80 lbls(100)
      character junk
      integer dum1, dum2
      integer iin1,iin2,iout,gmfile
      logical duminp
      DIMENSION AZ(1000000),TYPE(200),U(200,20),IA(200,6),IU(200,0:20)
      COMMON/MDIM/M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13
      common/input/duminp
 1    FORMAT(A10)
 2    FORMAT(2I5,25X,I5)
 3    FORMAT(1X,'NO INPUT FOUND.')
      MAXCOR=1000000
      IIN1=5
      IIN2=7
      IOUT=6
      gmfile=8
      call psinit('APPEND')
      CALL TSTART(IOUT)

      call ffile(iin1,'input',0)
      call ffile(iin2,'intco',0)

      CALL LOCATE(IIN1,'# INTDER #',IERR)
      if (ierr.eq.0) then
         duminp=.true.
         read(iin1,2) na,ns,ndum
      else
        duminp=.false.

         call ffile(gmfile,'geom',0)

         errcod = fwkclr()
         errcod = fwkadd(gmfile,'DUMMY')

         na = 0
         errcod = fcnt(':GEOMETRY',na)
         if(errcod.ne.EOK.or.na.eq.0) then
c          if geom.dat dosn't exist, get na from file11.dat
           call ffile(11,' ',0)
           rewind(11)
           read(11,770) junk
           read(11,'(i5)') na
 770       format(a)
           close(11)
           if(na.eq.0) then
             errcod = frdi('NATOMS',na)
             if(errcod.ne.EOK) then
               write(6,*) 'need NATOMS'
               call pabort
               call qabort
             endif
           endif
         endif
         close(gmfile)

         errcod = fwkclr()
         errcod = fwkadd(iin1,'DUMMY')
         errcod = fwkadd(iin1,'DEFAULT')
         errcod = fwkadd(iin2,'DUMMY')

         ndum =0 
         errcod = fcnt('DUMMY_ATOMS',ndum)

         dum2 = 0
         dum1 = 0
         errcod = fcnt(':INTCO:STRE',dum1)
         dum2 = dum2 + dum1
         dum1 = 0
         errcod = fcnt(':INTCO:BEND',dum1)
         dum2 = dum2 + dum1
         dum1 = 0
         errcod = fcnt(':INTCO:TORS',dum1)
         dum2 = dum2 + dum1
         dum1 = 0
         errcod = fcnt(':INTCO:OUT',dum1)
         dum2 = dum2 + dum1
         dum1 = 0
         errcod = fcnt(':INTCO:LIN1',dum1)
         dum2 = dum2 + dum1
         dum1 = 0
         errcod = fcnt(':INTCO:SPF',dum1)
         dum2 = dum2 + dum1
         dum1 = 0
         errcod = fcnt(':INTCO:LINX',dum1)
         dum2 = dum2 + dum1
         dum1 = 0
         errcod = fcnt(':INTCO:LINY',dum1)
         dum2 = dum2 + dum1
         dum1 = 0
         errcod = fcnt(':INTCO:RCOM',dum1)

         ns = dum2 + dum1
      endif
      
      NAD=NA+NDUM
      NC=3*NA
      NCD=3*NAD
      M1=1
      M2=M1+NCD
      M3=M2+NA
      M4=M3+NS
      M5=M4+NS
      M6=M5+NS*NC
      M7=M6+NC*NC
      M8=M7+NC*NC
      M9=M8+NC*NC
      M10=M9+NC*NC
      M11=M10+NC*NC
      M12=M11+3*NC
      M13=M12+3*NA
      NCORE=MAXCOR-M13+1
      CALL INTDER(NA,NAD,NC,NS,
     $ AZ(M1),AZ(M2),AZ(M3),AZ(M4),AZ(M5),AZ(M6),
     $ AZ(M7),AZ(M8),AZ(M9),AZ(M10),AZ(M11),AZ(M12),AZ(M13),
     $ TYPE,U,IA,IU,NCORE,IFLAG,lbls)
C      CALL TSTOP(3)
      CALL TSTOP(IOUT)
      STOP
 1000 WRITE(IOUT,3)
      STOP
      END
