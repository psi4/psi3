      subroutine fentry(ijnk1,ijnk2,ijnk3)
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
C Added calculation of spectroscopic constants, We, Be, De, Ae, WeXe
C  for diatomics.  Cleaned up some output.
C Daniel Crawford, David Sherrill, and Justin Fermann  7.24.95  
C search: cdc7.24.95
C**********************************************************************
c Program now enters through FENTRY, 3/3/89 - clj.
C**********************************************************************
c Moved to PSI distribution disk on 020389 - clj.
C**********************************************************************
C*  BY:  RICHARD REMINGTON                         search:  c4-10-88  *
C*  DATE:  April  10,  1988                                           *
C*  REASON: Migration to run in 7Mb on 9370                           *
C**********************************************************************
      IMPLICIT INTEGER (A-Z)
#include <error.h>
      integer frdi,frdboo,fcnt
      integer errcod
      integer ijnk1(1),ijnk2(2),ijnk3
      logical oldinp
      REAL*8 Z(10000000)
      INTEGER IZ(50000)
      CHARACTER LABEL*10
      COMMON/IO/IIN,IOUT,I12,I16,I21,IDER,ISCR1,ISCR2
    1 FORMAT(A10)
    2 FORMAT(1X,'NO INTDIF LABEL FOUND.')
    3 FORMAT(16I5)
    4 FORMAT(1X,'NOT ENOUGH CORE'/I10,' WORDS NEEDED BUT ONLY',
     $    I10,' WORDS AVAILABLE.')
    5 FORMAT(16I5)
    6 FORMAT(/1X,'CALCULATION TYPE IS ',A10)
    7 FORMAT(/1X,'    NS    NDER    NPERM   NPAR    NPRT   NMORSE')
    8 FORMAT(1X,8(I6,2X))
    9 FORMAT(/1X,'NPOINT    NDER   NGUESS')
cets041389
cets031191      call drum
      call psinit('APPEND')
      call tstart(6)
cets
      IIN=5
      IOUT=6
      IDER=8
      intco=9
      ISCR1=91
      ISCR2=92
      I12=12
      I16=16
      I21=21
cets041389
cets082891      call ffile(iin,'intder1',0)
cets082891      call ffile(iout,'intdifo',0)
c
      call ffile(ider,'ider',0)
      call ffile(iscr1,'file91a',0)
      call ffile(iscr2,'file92a',0)
cets041389
      NMAX=300000
      NIMAX=5000
      DO 50  I=1,NIMAX
  50  IZ(I)=0
      DO 52  I=1,NMAX
  52  Z(I)=0
c
cets072891 100  READ(IIN,1,END=105) LABEL
c      IF(LABEL.EQ.'# INTDIF #') GO TO 110
c      GO TO 100
c 105  WRITE(IOUT,2)
c      STOP
      oldinp=.false.
      call locate(iin,'# INTDIF #',ierr)
      if(ierr.eq.0) then
        oldinp=.true.

 110    READ(IIN,1) LABEL
        WRITE(IOUT,6) LABEL
 200    IF(LABEL.NE.' FCONSTDIF'.AND.LABEL.NE.'FCONSTDIF ') GO TO 220
        READ(IIN,3) NS,NDER,NPERM,NPAR,NPRT,NMORSE
      else
        idiatm=0
        errcod = frdboo('DIATOMIC',idiatm)

        if(idiatm.eq.0) then
          write(iout,6) 'FCONSTDIF'

          call ffile(intco,'intco',0)
          errcod = fwkadd(iin,'INTCO')
          errcod = fwkadd(intco,'INTCO')
       
          ns=0
          errcod = fcnt('SYMM',ns)
          nas=0
          errcod = fcnt('ASYMM',nas)

          ns=ns+nas

          errcod = fwkclr()
          errcod = fwkadd(iin,'INTDIF')
          errcod = fwkadd(iin,'DEFAULT')

          errcod = frdi('N_INTCO',ns)

          nder=1
          errcod = frdi('DERLVL',nder)

          nperm=0
          errcod = fcnt('SYMOPS',nperm)

          npar=0
          errcod = fcnt('POS_DISP',npar)

          nprt=0
          errcod = frdi('IPRINT',nprt)

          nmorse=2
          errcod = frdboo('MORSE',nmorse)
        else
          errcod = fwkclr()
          errcod = fwkadd(iin,'INTDIF')
          errcod = fwkadd(iin,'DEFAULT')

          npoint=0
          errcod = fcnt('POINTS',npoint)

          nguess=0
          errcod = frdboo('GUESS_F',nguess)

          nder=0
          errcod = frdi('DERLVL',nder)

          label = '    DIATOM'
          write(iout,6) label
          goto 220
        endif
      endif

      WRITE(IOUT,7)
      WRITE(IOUT,8) NS,NDER,NPERM,NPAR,NPRT,NMORSE
      if(nder.ge.1) call ffile(i12,'file12a',0)
      if(nder.ge.2) call ffile(i16,'file16a',0)
      if(nder.ge.3) call ffile(i21,'file21a',0)

      NX=NS
      F1=1
      F1A=F1+NX
      F1B=F1A+NX
      F1C=F1B+NX
      F1D=F1C+NX
      F1E=F1D+NX
      F2=F1E+NX
      NX=NX*NS
      F2A=F2+NX
      F2B=F2A+NX
      F3=F2B+NX
      NX=NX*NS
      NTOP=F3+NX
      IF(NDER.GE.2) THEN
          F3A=F3+NX
          F3B=F3A+NX
          F4=F3B+NX
          NX=NX*NS
          F4A=F4+NX
          F4B=F4A+NX
          NTOP=F4B+NX
      ELSE
          F3A=1
          F3B=1
          F4=1
          F4A=1
          F4B=1
      END IF
      IF(NTOP.GT.NMAX) THEN
          WRITE(IOUT,4) NTOP,NMAX
          STOP
      END IF
      IPAR=1
      IPP=IPAR+NS*NS
      IPERM=IPP+NS*NS
      IMORSE=IPERM+NS*(NS+1)
      ITOP=IMORSE+NS
      CALL DIFA(NS,NDER,NPERM,NPAR,NPRT,NMORSE,Z(F1),Z(F1A),Z(F1B),
     $  Z(F1C),Z(F1D),Z(F1E),Z(F2),Z(F2A),Z(F2B),Z(F3),Z(F3A),Z(F3B),
     $  Z(F4),Z(F4A),Z(F4B),IZ(IPAR),IZ(IPP),IZ(IPERM),IZ(IMORSE),
     +  oldinp)
      GO TO 1000

 220  IF(LABEL.NE.'    DIATOM') GO TO 240
      if(oldinp) then
        READ(IIN,5) NPOINT,NDER,NGUESS
      endif

      WRITE(IOUT,9)
      WRITE(IOUT,8) NPOINT,NDER,NGUESS
      S=1
      E=S+NPOINT

      CALL DIFB(NPOINT,NDER,NGUESS,Z(S),Z(E),oldinp)

 240  CONTINUE
1000  CONTINUE
cets072891
      close(iscr1,status='delete')
      close(iscr2,status='delete')
      call tstop(6)
      stop
      END
