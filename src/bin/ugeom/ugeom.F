      program ugeom
C----------------------------------------------------------------
c Updated by Seung-Joon Kim and Yaoming Xie to rotate the molecule     
c to a specific direction, so that the MO will be in a neat form.
C                                       (08/07/93)
c
c     installed to psi area by croag on sep 7 1993.  no credit taken
C----------------------------------------------------------------
c Read in the geometry from the geom.dat file and put in FILE30.
c matom = maximum number of atoms
C----------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (matom=200)
      DIMENSION I30(200),CHG(50),COORD(3,matom)
      integer prcntr,errcod,frdfa2,isizof,fwkadd
    1 FORMAT(/,2X,' THE GEOMETRY UPDATING PROGRAM')
    2 FORMAT(/,2X,' PARAMETERS'/
     1 2X,' MPOINT = ',I5/
     2 2X,' MCONST = ',I5/
     3 2X,' NCALCS = ',I5/
     4 2X,' NATOM  = ',I5/
     5 2X,' N3N    = ',I5/
     6 2X,' IGEOP  = ',I5)
    3 FORMAT(/,2X,' OLD COORDINATES IN FILE30'/
     1 39X,' OLD X IN A.U.'/
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z)
    4 FORMAT(3X,I5,1X,4F16.10)
    5 FORMAT(/,2X,' NEW COORDINATES IN FILE30'/
     1 39X,' NEW X IN A.U.'/
     2 4X,4H NO.,4X,7H CHARGE,9X,2H X,14X,2H Y,14X,2H Z)
    6 FORMAT(3F20.10)
    7 FORMAT(/,2X,' GEOMETRY IN FILE30 IS UPDATED')
C
      call psinit('APPEND')
      CALL TSTART(6)

      errcod = prcntr('ON BRIEF')
      call prdup()
C
      ITAPE6=6
      ITAP30=30
c
CRAWDAD (that's Tommy Dan) so that Russ won't fill /home/xerxes anymore
c     get rid of the old
c     call fioopn(7,'geom.dat','USER','READ')
c     add the new
      open(7,file='geom.dat',form='FORMATTED',status='OLD',iostat=ios)
            if (ios .ne. 0) then
              write(6,*) 'ERROR: you need a geom.dat'
              call pabort
              call qabort
            endif
c     end of changes - installation by jef 032593
C
      WRITE(6,1)
C
C************************************
C***UPDATE TAPE30 FOR NEW GEOMETRY***
C************************************
C
C   READ IN GEOMETRY FROM TAPE30
      CALL RFILE(ITAP30)
      CALL WREADW(ITAP30,I30,200,101,JUNK)
      MPOINT=I30(2)
      MCONST=I30(3)
      MCALCS=I30(4)
      NCALCS=I30(5)
      NATOM=I30(19)
      N3N=NATOM*3
      IGEOP=100+MCONST+MPOINT+NCALCS
      ichg=101+MCONST+MPOINT+MCALCS

c
      if (natom .gt. matom) then
         write(6,*) 'ERROR: You have too many atoms. Increase matom.'
         call qabort
         endif
C
      if (prcntr('IS_ON BRIEF').eq.0) then
        WRITE(6,2) MPOINT,MCONST,NCALCS,NATOM,N3N,IGEOP
        endif
      call wreadw(itap30,chg,natom*isizof('REAL'),ichg,junk)
      CALL WREADW(ITAP30,LOCCAL,1,IGEOP,JUNK)
      IGEOP=LOCCAL+80
      CALL WREADW(ITAP30,COORD,N3N*isizof('REAL'),IGEOP,JUNK)
C
      WRITE(6,3)
      DO 101 I=1,NATOM
      WRITE(6,4) I,CHG(I),COORD(1,I),COORD(2,I),COORD(3,I)
  101 CONTINUE
C
c
c Read in the next geometry.
      call fwkclr()
      errcod = fwkadd(7,'DUMMY')
      errcod = frdfa2(':GEOMETRY',coord,3,natom)
      if (errcod.ne.0) then
        write(6,*) 'ERROR: problem reading geometry from geom.dat'
        call pabort
        call qabort
        endif
      call fwkdef()
CWA
Cyx(08/06/93)    rotate the molecule to a certain direction
      if(natom.le.2) goto 103
      call rotate(natom,coord)
  103 WRITE(6,5)
      DO 102 I=1,NATOM
      WRITE(6,4) I,CHG(I),COORD(1,I),COORD(2,I),COORD(3,I)
  102 CONTINUE
      CALL WWRITW(ITAP30,COORD,N3N*2,IGEOP,JUNK)
      WRITE(6,7)
C
      CALL TSTOP(6)
      CALL RCLOSE(ITAP30,3)
      call psdone()
      STOP
      END

