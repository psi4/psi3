C     ////////////////////////////////////////////////////////////
      SUBROUTINE INTDER(NA,NAD,NC,NS,
     $        XA,XMASS,S,SS,B,BS,A,UGF,XS,XT,XKI,XAR,
     $        Z,TYPE,U,IA,IU,NCORE,IFLAG,lbls)
      IMPLICIT REAL*8 (A-H,O-Z)
#include <error.h>
      integer errcod, fwkclr, fwkadd, fcnt, frdi, dum1, dum2
      integer frdfa2
      integer nsyms, multi
      character junk*80
      CHARACTER LABEL*10,TYPE*5,NUMTST*5,ASTAT*6,LREF*4
      character lbls*80
      character*8 invrt,asdk,sfhess
      character*11 fqanal
      LOGICAL FLAG,RFLAG
      logical duminp, dumsqm, gmred
      integer gmfile
      real*8 gtmp(1500)
      DIMENSION IOPT(30),TYPE(NS),IA(NS,6),IU(NS,0:1),S(NS),SS(NS)
      dimension lbls(ns)
      DIMENSION Z(1),U(NS,1),IR(4),XR(4),XA(NAD,3),B(NS,NC)
      DIMENSION A(NC,NC),UGF(NC,1),BS(NC,NC),XS(NC,NC)
      DIMENSION XT(NC,NC),XMASS(1),DIP(3),XKI(NC,3),XAR(NA,3),CM(3)
      PARAMETER(RAD=57.29577951308232D0)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TOLINV=1.0D-10)
      common /input/ duminp,dumsqm
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /DIPOLE/DIP,QQ
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
 1    FORMAT(A10)
 3    FORMAT(17I5)
 4    FORMAT(A5,5I5,A5)
 5    FORMAT(I5,4(I4,F14.10))
 6    FORMAT(4F20.10)
 7    FORMAT(1X,'B MATRIX DETERMINANT = 0 .')
 8    FORMAT(1X,'AT EQUILIBRIUM ALL INTERNAL GRADIENTS = 0 .')
 9    FORMAT(1X,'NOT ENOUGH CORE'/1X,I10,' WORDS NEEDED AND ONLY',I10,
     $          ' WORDS AVAILABLE')
 10   FORMAT(3X,'PROGRAM OPTIONS'//
     $'   NA',5X,' NSYM',5X,'  NEQ',5X,' NINV',5X,'NTEST',5X,'NFREQ',
     $5X,' NVEC',5X,'NDISP',5x,'MULTI'/
     $5X,'   NS',5X,' NDER',5X,' NPRT',5X,' NDUM',5X,'NGEOM',5X,'IRINT',
     $5X,'NSTOP',5X,'NMODE')
 11   FORMAT(//1X,'NUCLEAR CARTESIAN COORDINATES (BOHR)'/)
 12   FORMAT(1X,3F20.10)
 13   FORMAT(//,1X,'DEFINITION OF INTERNAL COORDINATES')
 14   FORMAT(/,1X,'SIMPLE INTERNAL COORDINATES'/)
 15   FORMAT(/,1X,'SYMMETRY INTERNAL COORDINATES'/)
 16   FORMAT('S(',I2,')=',2X,4(F13.10,' L(',I2,')',3X),
     $     5(/,8X,4(F13.10,' L(',I2,')',3X)))
 17   FORMAT('L(',I2,')=',2X,A5,4I5,A5)
 18   FORMAT(///,1X,'VALUES OF SIMPLE INTERNAL COORDINATES ',
     $      '(ANG. OR DEG.)'/)
 19   FORMAT(4(I4,F16.10))
 20   FORMAT(/,1X,'B MATRIX FOR (SYMMETRY) INTERNAL COORDINATES'/)
 21   FORMAT(6F12.6)
 22   FORMAT(//' THIS VERSION OF THE PROGRAM USES NUMERICAL '/
     $  ' BR(I,J,K,L) MATRICES TO COMPUTE CONTRIBUTIONS TO THE'/
     $  ' QUARTIC FORCE CONSTANTS ARISING FROM NONZERO GRADIENTS.')
 23   FORMAT(/' B MATRIX INVERSION UNSUCCESSFUL FOR GIVEN TOLERANCE.')
 24   FORMAT(/' B*A MATRIX'/)
 25   FORMAT(/' NA FROM INPUT IS INCONSISTENT WITH THE HEADER IN A'/
     $    ' CARTESIAN DERIVATIVE FILE.   IFLAG=',I5)
 26   FORMAT(F20.10)
 27   FORMAT(///,1X,'VALUES OF SYMMETRY INTERNAL COORDINATES ',
     $      '(ANG. OR RAD.)'/)
 28   FORMAT(3F20.10)
 29   FORMAT(/,1X,'EXECUTION STOPPED DUE TO MASS INPUT ERROR IN ',
     $      'MASSIN.')
 30   FORMAT(//,1X,'NORMAL MODE ANALYSIS IN INTERNAL COORDINATES')
 31   FORMAT(//,1X,'NORMAL MODE ANALYSIS IN CARTESIAN COORDINATES')
 32   FORMAT(/,1X,'INTDER95'//
     $ ' A GENERAL PROGRAM DEVELOPED BY WESLEY D. ALLEN AND CO-WORKERS'/
     $ ' WHICH PERFORMS VARIOUS VIBRATIONAL ANALYSES AND HIGHER-ORDER'/
     $ ' NONLINEAR TRANSFORMATIONS AMONG FORCE FIELD REPRESENTATIONS'/)
 33   FORMAT(//'NUCLEAR MASSES'/)
 34   FORMAT(4(I4,F12.6,4X))
 35   FORMAT(I5,F20.10)
 36   FORMAT(/,1X,'EXECUTION STOPPED DUE TO ERROR IN GFMAT.',
     $      '  IFLAG=',I5)
 37   FORMAT(/,1X,'EXECUTION STOPPED DUE TO ERROR IN NORMCO.',
     $      '  IFLAG=',I5)
 38   FORMAT(1X,  'REFERENCE VALUE IS ',F20.10,' A.')
 39   FORMAT(//,1X,   'NDER=0   NO TRANSFORMATION REQUESTED.')
 40   FORMAT(//,1X,'SCALED QUANTUM MECHANICAL FORCE CONSTANT ANALYSIS')
 41   FORMAT(/,1X,'EXECUTION STOPPED DUE TO ERROR IN SQMFC',
     $      '  IFLAG=',I5)
 42   FORMAT(2I5,5X,I5)
 43   FORMAT(/1X,'DUMMY ATOM VECTORS (BOHR)'/)
 44   FORMAT(/1X,'EXECUTION STOPPED, NOT ENOUGH CORE TO PERFORM',
     $       ' SQM ANALYSIS.'/)
 45   FORMAT(//,'EXECUTION STOPPED.  REFERENCE GEOMETRY NOT FOUND.')
 46   FORMAT(/1X,'REFERENCE GEOMETRY (BOHR)'/)
 47   FORMAT(///' The input options specified are not compatible with ',
     $       'multiple'/' transformations, i.e. IOPT(17) > 0 (in the ',
     $       'current version of this'/' program).'/,
     $       ' Your values : NINV =',I7,',  NGEOM =',I3,
     $       ',  NFREQ =',I3,',  IRINT =',I3/
     $       ' should be   : NINV = 0 or 1,  NGEOM =  0, |NFREQ| < 4, ',
     $       '|IRINT| < 2'/)

c     start the new input parsing stuff here
c      call psinit('APPEND')

      IFLAG=0
      IIN1=5
      IOUT=6
      IIN2=7
      gmfile = 8
      ICHECK=3
      I11=11
      I12=12
      I15=15
      I16=16
      I17=17
      I18=18
      I20=20
      I21=21
      I24=24
      I25=25
CSA
      I31=31
      I32=32
      I33=33
      I35=35
      I36=36
      I37=37
CSA
      ISCR1=91
      ISCR2=92
      ISCR3=93
      ISCR4=94
      ISCR5=95
      ISCR6=96
      ISCR7=97
      ISCR8=98
      ISCR9=80
      ISCR10=81
      ISCR11=85
      ISCR12=86
      ISCR13=87
      ISCR14=88
CSA
      BOHR=0.529177249D0
      DEBYE=2.54176548D0
      HART=4.3597482D0
      WAVE0=1302.7910D0
      CINT=42.25472D0
CSA
      NCHUNK=1000
      ENERGY=ZERO
      CALL RFILE(ISCR1)
      CALL RFILE(ISCR2)
      CALL RFILE(ISCR3)
      CALL RFILE(ISCR4)
      CALL RFILE(ISCR6)
      CALL RFILE(ISCR7)
      CALL RFILE(ISCR8)
      CALL RFILE(ISCR9)
      CALL RFILE(ISCR10)

c     more new input parsing stuff

      errcod = fwkclr()
      errcod = fwkadd(iin1,'DUMMY')
      errcod = fwkadd(iin2,'DUMMY')
      errcod = fwkadd(iin1,'INTDER')
      errcod = fwkadd(iin1,'DEFAULT')


      CALL ffile(I11,'file11',0)
      CALL ffile(I12,'file12',0)
      CALL ffile(I15,'file15',0)
      CALL ffile(I16,'file16',0)
      CALL ffile(I17,'file17',0)
      CALL ffile(I18,'file18',0)
      CALL ffile(I20,'file20',0)
      CALL ffile(I21,'file21',0)
      CALL ffile(I24,'file24',0)
      CALL ffile(I25,'file25',0)
      CALL ffile(ISCR5,'file95',0)
      CALL ffile(ICHECK,'intderAUX',0)

      CALL LOCATE(IIN1,'# INTDER #',IERR)
      IF (IERR.EQ.0) THEN
         duminp=.true.
COMMENT111111...................
 105     READ(IIN1,3) (IOPT(K),K=1,16)
         iopt(17) = 0
      ELSE
         duminp=.false.

         iopt(1) = na
         iopt(2) = ns

         IOPT(3) = 0
         dum1=0
         errcod = fcnt(':INTCO:SYMM',IOPT(3))
         errcod = fcnt(':INTCO:ASYMM',dum1)
         iopt(3) = iopt(3) + dum1

         IOPT(4)=2
         errcod = frdi('DERLVL',IOPT(4))

         iopt(5)=0
         errcod = frdboo('TRANS_1DER',iopt(5))

         iopt(6)=3000
         errcod = frdi('IPRINT',iopt(6))

         iopt(7)=0
         invrt=' '
         errcod = frdc('TRANS_TYPE',invrt)
         if(invrt.eq.'C_TO_I  ') iopt(7)=0
         if(invrt.eq.'I_TO_C  ') iopt(7)=1
         if(invrt.eq.'I_TO_C_R') iopt(7)=2
         if(invrt.eq.'PROJECT') iopt(7)=3

         dum1 = 0
         errcod = frdboo('MASS_TRANS',dum1)
         if(dum1.gt.0) then
           iopt(7)=-iopt(7)
         endif

         iopt(8)=0
         errcod = fcnt('DUMMY_ATOMS',iopt(8))

         iopt(9)=0
         errcod = frdi('NTEST',iopt(9))
         dum1=0
         errcod = frdboo('CHECK_ORTHO',dum1)
         if(errcod.eq.EOK.and.dum1.eq.1) then
            if(iopt(7).gt.0) then
                iopt(7) = 3
                iopt(9) = 3
            else
                iopt(7) = -3
                iopt(9) = -3
            endif
         endif
         dum1=0
         errcod = frdboo('CHECK_INV',dum1)
         if(dum1.eq.1) iopt(9) = 4
         dum1=0
         errcod = frdboo('REDUN_TEST',dum1)
         if (dum1.eq.1) iopt(9)=5
          
         iopt(10)=0
         errcod = frdboo('READ_GEOM',iopt(10))


         iopt(11)=0
         errcod = frdc('FREQ_ANAL',fqanal)
         if(fqanal.eq.'NONE       ') iopt(11)=0
         if(fqanal.eq.'INT_CO     ') iopt(11)=1
         if(fqanal.eq.'CART_CO    ') iopt(11)=2
         if(fqanal.eq.'BOTH       ') iopt(11)=3
         if(fqanal.eq.'INT_CO_R   ') iopt(11)=4
         if(fqanal.eq.'SQMFC      ') iopt(11)=5
         if(fqanal.eq.'RXN_PATH   ') iopt(11)=10
         if(fqanal.eq.'INT_NO_T   ') iopt(11)=-1
         if(fqanal.eq.'CART_NO_T  ') iopt(11)=-2
         if(fqanal.eq.'BOTH_NO_T  ') iopt(11)=-3
         if(fqanal.eq.'INTCOR_NO_T') iopt(11)=-4
         if(fqanal.eq.'SQMFC_NO_T ') iopt(11)=-5

         if (iopt(11).eq.10) then
            dum1=0
            errcod = frdi('RXN_COORD',dum1)
            if(errcod.ne.EOK) then
               write(iout,*) 'ERROR: must specify RXN_COORD'
               call pabort
               call qabort
            endif
            iopt(11)=iopt(11)+dum1
         endif

         iopt(12)=0
         errcod = frdc('IRINT_T',invrt)
         if(invrt.eq.'NONE    ') iopt(12)=0
         if(invrt.eq.'DO      ') iopt(12)=1
         if(invrt.eq.'DO_R    ') iopt(12)=2

         iopt(13)=0
         errcod = frdboo('TRANS_DIPDER',iopt(13))

         iopt(14)=0
         errcod = frdi('STOP',iopt(14))

         iopt(15)=0
         errcod = fcnt('DISPLACEMENTS',dum1)
         if (errcod.eq.EOK) then
            iopt(4)=0
            iopt(15)=1
            dum1 = 0
            errcod = frdboo('HIGH_ORDER',dum1)
            if (dum1.ne.0) iopt(15)=2
         endif

         iopt(16)=0
         errcod = frdboo('PED',iopt(16))
 
         iopt(17)=0
         errcod = frdi('MULTI',iopt(17))


         dum1 = 0
         errcod = fcnt('MASSES',dum1)
         if(dum1.gt.0) then
            iopt(15)=-iopt(15)
         endif

      END IF

c     end of psi2.0 input parsing

      NA=IOPT(1)
      NS=IOPT(2)
      NSYM=IOPT(3)
      NDER=IOPT(4)
      NEQ=IOPT(5)
      NPRT=IOPT(6)
      NINV=IOPT(7)
      NDUM=IOPT(8)
      NTEST=IOPT(9)
      NGEOM=IOPT(10)
      NFREQ=IOPT(11)
      IRINT=IOPT(12)
      NVEC=IOPT(13)
      NSTOP=IOPT(14)
      NDISP=IOPT(15)
      NMODE=IOPT(16)
      MULTI=IOPT(17)

      IF ((MULTI.GT.1).AND.((NINV.GE.2).OR.(NGEOM.EQ.1).OR.
     $      (IABS(NFREQ).EQ.4).OR.(IABS(IRINT).GT.1))) THEN
            WRITE(IOUT,47) NINV,NGEOM,NFREQ,IRINT
            STOP
            END IF

      NREF=0
      NCOM=0
      IF(NVEC.NE.0) THEN
         IOPT(5)=1
         NEQ=1
      END IF
      NC=3*NA
      NAD=NA+NDUM
CSA
      IF(NTEST.EQ.5) THEN
       ITST=2
       NTEST=0
      ELSE IF (ABS(NTEST).EQ.3) THEN
       NSTOP=1
          NINV=ISIGN(1,NTEST)
      ELSE
          ITST=1
      END IF
CSA

c     continue old style input

CSA
      IF (ABS(NINV).EQ.3) THEN
         NEQ=0
      NVEC=0
      GO TO 117
      END IF

      if (duminp) then
         DO 100  J=1,NS
            IA(J,6)=0
CSA         IF(QTEST) IA(J,6)=NTEST
            IF(ABS(NTEST).LE.2.AND.NTEST.NE.0) IA(J,6)=NTEST
 100     CONTINUE

         DO 110  J=1,NS
          READ(IIN1,4) TYPE(J),(IA(J,K),K=1,5),NUMTST
          IF(TYPE(J).EQ.'STRE ') TYPE(J)=' STRE'
          IF(TYPE(J).EQ.'BEND ') TYPE(J)=' BEND'
          IF(TYPE(J).EQ.'LIN1 ') TYPE(J)=' LIN1'
          IF(TYPE(J).EQ.'LINX ') TYPE(J)=' LINX'
          IF(TYPE(J).EQ.'LINY ') TYPE(J)=' LINY'
          IF(TYPE(J).EQ.'TORS ') TYPE(J)=' TORS'
          IF(TYPE(J).EQ.' OUT '.OR.TYPE(J).EQ.'OUT  ')
     $      TYPE(J)='  OUT'
          IF(TYPE(J).EQ.' SPF '.OR.TYPE(J).EQ.'SPF  ')
     $      TYPE(J)='  SPF'
          IF(TYPE(J).EQ.'RCOM ') TYPE(J)=' RCOM'
      IF(NUMTST.EQ.'  ST '.OR.NUMTST.EQ.' ST  '.OR.NUMTST.EQ.'ST   ')
     +   NUMTST='   ST'
          IF(NTEST.NE.0.AND.NUMTST.EQ.'   ST') IA(J,6)=0
          IF(TYPE(J).EQ.'  SPF') THEN
                IF(IA(J,3).NE.0)  READ(IIN1,26) S(J)
          END IF
          IF(TYPE(J).EQ.' RCOM') NCOM=1
 110     CONTINUE

         IF(NSYM.EQ.0) GO TO 125
         DO 115  K=1,NSYM
 115     IU(K,0)=0
 116     READ(IIN1,5) L,(IR(K),XR(K),K=1,4)
         IF(L.EQ.0) GO TO 125
         IPOINT=IU(L,0)
         DO 120  K=1,4
           M=IR(K)
           IF(M.GT.0) THEN
             IPOINT=IPOINT+1
             IU(L,IPOINT)=M
             U(L,IPOINT)=XR(K)
           END IF
 120     CONTINUE
         IU(L,0)=IPOINT
         GO TO 116
 125  DO 121  K=1,NSYM
         M=IU(K,0)
         UNORM=ZERO
         DO 122  L=1,M
 122  UNORM=UNORM+U(K,L)*U(K,L)
         UNORM=ONE/DSQRT(UNORM)
         DO 123  L=1,M
 123  U(K,L)=U(K,L)*UNORM
 121  CONTINUE
CSA
c     possibly the end of old style input -

      ELSE
c     do all psi2.0 input
        nsyms = nsym
        call rdsimp(type,ia,s,ns)
        do 213  j=1,ns
            ia(j,6)=0
            if(abs(ntest).le.2.and.ntest.ne.0) ia(j,6)=ntest
            lbls(j) = '                     '
 213     continue

        call rdsymm(lbls,iu,u,ns,nsyms)
        
      endif
c  possibly end of psi2.0 input changes

 117  continue

      DO 800 IMULT=1,MAX(MULTI,1)
        if(ngeom.eq.0) then
          if(imult.eq.1.and.ndum.gt.0) then
            if (duminp) then
               do 132 j=na+1,nad
 132             read(iin1,28) (xa(j,k), k=1,3)
            else 
               ij=0
               errcod = frdfa2('DUMMY_ATOMS',gtmp,3,ndum)
               do 135 j=na+1,nad
                 do 135 k=1,3
                   ij=ij+1
                   xa(j,k)=gtmp(ij)
  135          continue
            endif
               
          endif
          IF (MULTI.EQ.0) THEN
            call bottom(i11)
            ik=2*na+1
            DO 128  K=1,IK
 128          BACKSPACE I11
          ELSE
            READ(I11,1) LABEL
          END IF
c
          READ(I11,35) NA11,ENERGY
          DO 130  J=1,NA
 130        READ(I11,6) xmass(j),(XA(J,K),K=1,3)
c
          IF ((MULTI.NE.0).AND.(NEQ.EQ.0)) THEN
c          Skip over 1st derivatives in File11
             DO 178 K=1,NA
 178           READ(I11,1) LABEL
          END IF
        end if
c
       IF (NGEOM.EQ.1) THEN
         if(duminp) then
           DO 210  J=1,NAD
 210       READ(IIN1,28) (XA(J,K),K=1,3)

           IF (NREF.NE.0) THEN
             READ(IIN1,FMT='(A4)') LREF
             IF (LREF.NE.'REFG') THEN
               WRITE (IOUT,45)
               GO TO 750
             END IF
             DO 211 J=1,NA
 211         READ (IIN1,28) (XAR(J,K),K=1,3)
           END IF
         else
c     read geom from geom.dat
           call ffile(gmfile, "geom", 0)
           errcod = fwkadd(gmfile,'DUMMY')
           errcod = frdfa2(':GEOMETRY',gtmp,3,na)
           close(gmfile)
           if (errcod.ne.EOK) then
             write(6,*) 'ERROR: problem reading geometry from geom.dat'
             call pabort
             call qabort
           endif
           ij=0
           do 780 j=1,na
             do 780 k=1,3
               ij=ij+1
               xa(j,k)=gtmp(ij)
  780      continue
           if(ndum.ne.0) then
             ij=0
             errcod = frdfa2('DUMMY_ATOMS',gtmp,3,ndum)
             do 781 j=na+1,nad
               do 781 k=1,3
                 ij=ij+1
                 xa(j,k)=gtmp(ij)
  781        continue
           endif
         endif
       END IF
CSA
      IF(ABS(NINV).EQ.3) GO TO 215
      DO 214  J=1,NS
      IF(TYPE(J).EQ.'  SPF'.AND.IA(J,3).EQ.0) THEN
         CALL VECT1(NAD,IA(J,1),IA(J,2),XT(1,1),XA,S(J))
         S(J)=S(J)*BOHR
      END IF
 214  CONTINUE
C
 215  CONTINUE
CSA
      IF((NFREQ.NE.0.AND.ABS(NFREQ).NE.5).OR.NVEC.GT.0.OR.NDISP.LT.0
     $    .OR.NCOM.GT.0.OR.NINV.LT.0) THEN
          CALL MASSIN (XMASS,NA,IFLAG,0,'MASSES  ',duminp)
          IF (IFLAG.NE.0) GO TO 1018
      END IF
      DO 212  I=1,NA
      IF (NVEC.GT.0.OR.NINV.LT.0) THEN
         XT(I,1) = XMASS(I)
      ELSE
         XT(I,1) = ONE
      END IF
      IF (NINV.EQ.3)    XMASS(I)= ONE
 212  CONTINUE
C
      IF ((ABS(NINV).EQ.3.OR.ABS(NTEST).EQ.3).AND.NREF.EQ.0) THEN
         CM(1)=ZERO
         CM(2)=ZERO
         CM(3)=ZERO
C
         TM=ZERO
         DO 240 J=1,NA
            TM=TM+XMASS(J)
         DO 240 K=1,3
 240        CM(K)=CM(K)+XA(J,K)*XMASS(J)
C
         DO 241 K=1,3
 241        CM(K)=CM(K)/TM
C
         DO 242 J=1,NA
         DO 242 K=1,3
 242        XAR(J,K)=XA(J,K)-CM(K)
C
      END IF
C
C     *********************************
      WRITE(IOUT,32)
      WRITE(IOUT,10)
      WRITE(IOUT,3)(IOPT(K),K=1,17)
      WRITE(IOUT,11)
      DO 133  J=1,NA
 133  WRITE(IOUT,12) (XA(J,K),K=1,3)
      IF(NDUM.GT.0) THEN
      WRITE(IOUT,43)
      DO 143  J=NA+1,NAD
 143  WRITE(IOUT,12) (XA(J,K),K=1,3)
      END IF
      DO 131  I=1,NA
      DO 131  J=1,3
 131  XA(I,J)=XA(I,J)*BOHR
CSA
      IF (ABS(NINV).EQ.3.OR.ABS(NTEST).EQ.3) THEN
         WRITE(IOUT,46)
         DO 144 J=1,NA
            WRITE (IOUT,12) (XAR(J,K),K=1,3)
         DO 144 K=1,3
 144     XAR(J,K)=XAR(J,K)*BOHR
      END IF
CSA
      IF(ABS(NINV).EQ.3) GO TO 700
CSA
      WRITE(IOUT,13)
      WRITE(IOUT,14)
      DO 134  J=1,NS
      NUMTST='     '
      IF(NTEST.NE.0.AND.ABS(NTEST).LE.2.AND.IA(J,6).EQ.0) NUMTST='   ST'
      WRITE(IOUT,17) J,TYPE(J),(IA(J,K),K=1,4),NUMTST
      IF(TYPE(J).EQ.'  SPF') WRITE(IOUT,38) S(J)
 134  CONTINUE
      DO 136  J=1,NS
      DO 136  I=1,NC
 136  UGF(I,J)=ZERO
      IF(NSYM.EQ.0) GO TO 145
      WRITE(IOUT,15)
      DO 137  I=1,NSYM
      LL=IU(I,0)
      WRITE(IOUT,16) I,(U(I,L),IU(I,L),L=1,LL)
      DO 137  L=1,LL
      K=IU(I,L)
 137  UGF(I,K)=U(I,L)
 145  IF(LPRT(4,NPRT).LT.4) GO TO 150
      DO 141  J=1,NAD
 141  WRITE(ICHECK,142) (XA(J,K),K=1,3)
 142  FORMAT(3F10.7)
      DO 138  I=1,NSYM
 138  WRITE(ICHECK,139) (UGF(I,K),K=1,NS)
 139  FORMAT(8F10.7)
C     **********************************
COMMENT2222222222.......................
 150  CONTINUE
      IF(NDISP.NE.0) THEN
         N1=1
         N2=N1+3*NAD
         N3=N2+3*NA
         NT=N3+NS
         CALL DISP(NAD,NA,NC,NS,NDISP,IOPT,XA,XMASS,S,SS,B,BS,A,
     $     TYPE,U,IA,IU,Z(N1),Z(N2),Z(N3),UGF,XS,XT,IFLAG)
         GO TO 750
      END IF
      CALL MACHB(NAD,NC,NS,XA,XMASS,TYPE,IA,B,S)
      IF(NSYM.LE.0) THEN
         DO 154  I=1,NS
         SS(I)=S(I)
         DO 154  J=1,NC
 154     BS(I,J)=B(I,J)
      ELSE
         DO 156  J=1,NC
         DO 156  I=1,NSYM
           BS(I,J)=ZERO
           DO 158  L=1,IU(I,0)
           K=IU(I,L)
 158       BS(I,J)=BS(I,J)+U(I,L)*B(K,J)
 156     CONTINUE
         DO 160  I=1,NSYM
           SS(I)=ZERO
           DO 160  L=1,IU(I,0)
           K=IU(I,L)
 160       SS(I)=SS(I)+U(I,L)*S(K)
      END IF
      IF(NVEC.GT.0.OR.NINV.LT.0) THEN
           WRITE(IOUT,33)
           WRITE(IOUT,34) (I,XMASS(I),I=1,NA)
      END IF
      NSX=NS
      IF(NSYM.GT.0) NSX=NSYM
      NSY=NSX
      DO 162  I=1,NS
      UGF(I,1)=S(I)
      IF(TYPE(I).EQ.' BEND'.OR.TYPE(I).EQ.' LIN1'.OR.TYPE(I).EQ.
     $   ' TORS'.OR.TYPE(I).EQ.'  OUT') THEN
      UGF(I,1)=UGF(I,1)*RAD
      END IF
 162  CONTINUE
      WRITE(IOUT,18)
      WRITE(IOUT,19) (I,UGF(I,1),I=1,NS)
      IF(NSYM.GT.0) THEN
      WRITE(IOUT,27)
      WRITE(IOUT,19) (I,SS(I),I=1,NSYM)
      END IF
      IF(LPRT(1,NPRT).GE.2) THEN
         WRITE(IOUT,20)
         DO 170  I=1,NSX
         WRITE(IOUT,*) 'I= ',I
 170     WRITE(IOUT,21) (BS(I,J),J=1,NC)
      END IF
      IF(NFREQ.LT.0) GO TO 600
      IF(NINV.GT.0) THEN
           NSY=NC
CSA        IF(NTEST.EQ.0) GO TO 175
      END IF
      CALL BINVRT(NSX,NC,XT(1,1),BS,UGF,A,XS,XT(1,2),IFLAG,TOLINV,ITST)
      IF(IFLAG.NE.0) GO TO 1002
COMMENT44444444444........................
 175  N1=1
      IF(NDER.EQ.0) GO TO 1014
      IF(NDER.EQ.1.AND.NEQ.EQ.0) GO TO 1004
      IF(NDER.EQ.1) THEN
         N2=1
         N3=1
         N4=1
         N5=N1+NC
         GO TO 180
      END IF
      IF(NEQ.EQ.0) THEN
         N2=1
      ELSE
         N2=NC+1
      END IF
      IF(NDER.EQ.2) THEN
         N3=1
         N4=1
         N5=N2+NC**2
         GO TO 180
      END IF
      N3=N2+NC**2
      IF(NDER.EQ.3) THEN
         N4=1
         N5=N3+NC**3
         GO TO 180
      END IF
      IF(NDER.EQ.4) THEN
         N4=N3+NC**3
         N5=N4+NC**4
         GO TO 180
      END IF
 180  NMM=NCHUNK*3
      NT=N5+NMM
      MDER=NDER+NEQ
      IF(MDER.GE.4) THEN
           NNN=NC**3
           N6=N5+NNN
           NMM=MAX0(NMM,2*NNN)
           NT=N5+NMM
      END IF
      IF(NT.GT.NCORE) GO TO 1006
      IF(NINV.EQ.2) THEN
        CALL INPFKM(NC,NDER,NEQ,NSX,Z(N1),Z(N2),Z(N3),Z(N4))
      END IF
      IF(MDER.GE.3) THEN
         IF(ABS(NTEST).EQ.1) THEN
            CALL NUMX(NAD,NC,NS,NSX,XA,XMASS,TYPE,IA,XT(1,1),A,S,UGF,XS)
         END IF
         CALL MACHX(NAD,NC,NS,NSX,IOPT,XA,XMASS,TYPE,IA,A,S,U,IU,UGF,XS)
         IF(NTEST.EQ.1) THEN
            CALL SRTST1(NC,NS,NSX,NSYM,IA,XS,XT)
         END IF
      END IF
      IF (MDER.GE.4) THEN
         IF(ABS(NTEST).EQ.2) THEN
            CALL NUMY(NAD,NC,NS,NSX,XA,TYPE,IA,A,S,Z(N5),Z(N6),UGF)
         END IF
       CALL MACHY(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,U,IU,Z(N5),Z(N6))
         IF(NTEST.EQ.2) THEN
            CALL SRTST2(NC,NS,NSX,NSYM,IA,Z(N5),Z(N6))
         END IF
      END IF
      IF (MDER.GE.5) THEN
       CALL NUMZ(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,Z(N5),Z(N6))
       CALL MACHZ(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,U,IU,Z(N5),Z(N6))
      END IF
      IF(NSTOP.EQ.1) GO TO 700
COMMENT5555555555...................................
      NRUN=0
      IF(NVEC.GT.0) THEN
         CALL VCXKI(NA,NAD,NC,XMASS,XA,XKI)
      END IF
 200  NRUN=NRUN+1
 220  IF(NINV.LE.0) THEN
      CALL LINTR(NC,NSX,NC,IOPT,A,XS,Z(N1),Z(N2),Z(N3),Z(N4),Z(N5),
     $         NMM,NRUN,IFLAG)
      ELSE
      CALL LINTR(NC,NC,NSX,IOPT,BS,XS,Z(N1),Z(N2),Z(N3),Z(N4),Z(N5),
     $         NMM,NRUN,IFLAG)
      END IF
 230  IF(IFLAG.NE.0) GO TO 1008
      IF(NVEC.GT.0) THEN
      CALL VCDER1(NA,NC,NSX,NINV,XMASS,XKI,A,Z(N1),Z(N2),
     $        Z(N3),Z(N4),Z(N5),NRUN)
      END IF
      IF(NDER.LE.1) GO TO 500
      IF(NDER.EQ.2.AND.NEQ.EQ.0) GO TO 500
      IF(NEQ.NE.0) THEN
         CALL MACHF2(NC,NSX,NINV,Z(N1),Z(N2),UGF,Z(N5))
      END IF
      IF(NDER.LE.2) GO TO 500
C      ADD X-V(M,N) CONTRIBUTIONS TO 3RD AND PERHAPS 4TH DERIVATIVES.
      CALL XF2(NDER,NC,NSX,NINV,BS,Z(N2),Z(N3),Z(N4),Z(N5),UGF,XS,XT)
      IF(NDER.EQ.3.AND.NEQ.EQ.0) GO TO 500
C      FORM Y MATRICES AND WRITE TO DISK.
C      COMPLETE V(M,N,P) COMPUTATION.
      IF(NEQ.NE.0) THEN
           CALL MACHF3(NC,NSX,NINV,Z(N1),Z(N3),Z(N5),Z(N6))
      END IF
      IF(NDER.LE.3) GO TO 500
C      ADD X-V(M,N,P) CONTRIBUTIONS TO 4TH DERIVATIVES.
      CALL XF3(NC,NSX,NINV,BS,UGF,XS,XT,Z(N3),Z(N4),Z(N5))
C      ADD Y-V(M,N) CONTRIBUTIONS TO 4TH DERIVATIVES.
      CALL YF2(NC,NSX,NINV,BS,UGF,XS,Z(N2),Z(N4),Z(N5))
      IF(NDER.EQ.4.AND.NEQ.EQ.0) GO TO 500
      WRITE(IOUT,22)
C      FORM Z MATRIX ELEMENTS NUMERICALLY AND
C        COMPLETE THE V(M,N,P,Q) COMPUTATION.
      CALL MACHF4(NC,NSX,NINV,Z(N1),Z(N4),Z(N5),Z(N6))
 500  CALL FCOUT(NC,NS,NSY,NEQ,NDER,NINV,NVEC,NRUN,SS,ENERGY,
     $           Z(N1),Z(N2),Z(N3),Z(N4))
      IF(NVEC.EQ.1.AND.NRUN.LT.3) GO TO 200
 600  NFREQ=IABS(NFREQ)
      IF(NFREQ.NE.0.AND.NFREQ.NE.2.AND.NFREQ.NE.5) THEN
         WRITE(IOUT,30)
         CALL GFMAT(NA,NAD,NC,NSX,NFREQ,IRINT,NMODE,IA,
     $           XMASS,XA,S,BS,UGF,A,XS,XT,XS,IFLAG,lbls)
         IF(IFLAG.NE.0) GO TO 1010
      END IF
      IF (NFREQ.EQ.2.OR.NFREQ.EQ.3) THEN
         WRITE(IOUT,31)
         CALL NORMCO(NA,NAD,NC,NFREQ,IRINT,
     $            XMASS,XA,BS,UGF,A,XS,XT,XS,IFLAG)
         IF(IFLAG.NE.0) GO TO 1012
      END IF
      IF (NFREQ.EQ.5) THEN
        WRITE(IOUT,40)
        CALL LOCATE(IIN1,'# SQMFC ##',IERR)
        if(ierr.eq.0) then
          dumsqm = .true.
          READ(IIN1,42) NSF,NISO,NH
        else
          dumsqm = .false.
          errcod = fwkclr()
          errcod = fwkadd(iin1,'SQMFC')
          errcod = fwkadd(iin1,'INTDER:SQMFC')
          errcod = fwkadd(iin1,'INTDER')
          errcod = fwkadd(iin1,'DEFAULT')
          errcod = fcnt('SCALE_FACTORS',nsf)
          errcod = fcnt('ISOTOPES',niso)
          nh = 0
          errcod = frdc('SF_HESS',sfhess)
          if(sfhess.eq.'GUESS   ') nh = 0
          if(sfhess.eq.'INPUT   ') nh = 1
          if(sfhess.eq.'ANAL    ') nh = 2
          if(sfhess.eq.'NUM_FD  ') nh = -1
          if(sfhess.eq.'LLS     ') nh = -2
          errcod = fwkclr()
        endif
c end of changes

         N1=1
         N2=N1+NA*NISO
         N3=N2+NSX*NISO
         N4=N3+NSX*NISO
         N5=N4+NSX*NSX*NISO
         N6=N5+NSX*NSX
         N7=N6+NSF
         N8=N7+NSF
         N9=N8+NSF
         N10=N9+NSF
         N11=N10+NSF*NSF
         N12=N11+NSF*NSF
         IF(ABS(NH).GE.2) THEN
           N13=N12+NSX*NSX*NSF
           NT=N13+NSX*NSF*NISO
         ELSE
           N13=N12
           NT=N12
         END IF
         IF(NT.GT.NCORE) GO TO 1020
         CALL SQMFC(NA,NAD,NC,NSX,NSF,NISO,NMODE,XA,BS,
     $   Z(N1),Z(N2),Z(N3),Z(N4),Z(N5),Z(N6),Z(N7),Z(N8),Z(N9),Z(N10),
     $   S,SS,XKI(1,1),XKI(1,2),XKI(1,3),A,XS,UGF,XT,UGF,Z(N11),
     $   Z(N12),Z(N13),IFLAG)
         IF(IFLAG.NE.0) GO TO 1016
      END IF
CSA
      GO TO 750
 700  CONTINUE
      IF (ABS(NTEST).EQ.3) THEN
         N2=NC*NC+1
         N3=N2+NC**3
         IF (N3.GT.NCORE) GO TO 1006
         CALL ORTHOG(NA,NC,NSX,IOPT,XAR,BS,Z(1),Z(N2),FLAG)
         IF (FLAG) GOTO 750
      END IF
      IF (NSTOP.EQ.1) GO TO 750
C
C IF ABS(NINV).EQ.3 CALCULATE THE PROJECTION
C MATRIX TO SEPARATE THE ROTATION
C AND THE TRANSLATION FROM THE FORCE CONSTANT MATRICES CALCULATED AT
C A NON-STATIONARY REFERENCE GEOMETRY
C
C ALLOCATE MEMORY FOR PROJECTION ROUTINES
C NDER=2
C       DK1     6*NC            N1
C       P2      NC,NC           N2
C       F2P     NC,NC           N3
C NDER=3
C       DK2     3*NC*NC         N4
C       P3      NC*NC*NC        N5
C       F3P     NC*NC*NC        N6
C       SCRTCH  NC**3           N6A
C       SCRTCH  NC**2           N6B
C NDER=4
C       DK3      3*NC*NC*NC     N7
C       P4      NC*NC*NC*NC     N8
C       F4P     NC*NC*NC*NC     N9
C       SCRATCH NC**4           N10
      NMM=NCHUNK*3
      IF(NDER.EQ.2) THEN
        N1  = 1
        N2  = N1  + 6*NC
        N3  = N2  + NC*NC
        N4  = N3  + NC*NC
        N5  = N4
        N6  = N5
        N6A = N6
        N6B = N6A
        N7  = N6B
        N8  = N7
        N9  = N8
        N10 = N9
        N11 = N10
      END IF
      IF(NDER.EQ.3) THEN
        N1  = 1
        N2  = N1  + 6*NC
        N3  = N2  + NC*NC
        N4  = N3  + NC*NC
        N5  = N4  + 3*NC*NC
        N6  = N5  + NC*NC*NC
        N6A = N6  + NC*NC*NC
        N6B = N6A + NC*NC*NC
        N7  = N6B + NC**2
        N8  = N7
        N9  = N8
        N10 = N9
        N11 = N10
      END IF
      IF(NDER.EQ.4) THEN
        N1  = 1
        N2  = N1  + 6*NC
        N3  = N2  + NC*NC
        N4  = N3  + NC*NC
        N5  = N4  + 3*NC*NC
        N6  = N5  + NC*NC*NC
        N6A = N6  + NC*NC*NC
        N6B = N6A + NC*NC*NC
        N7  = N6B + NC**2
        N8  = N7  + 3*NC**3
        N9  = N8  + NC**4
        N10 = N9  + NC**4
        N11 = N10 + NC**4
      END IF
CSA
      NT=N11+NMM
      IF(NT.GT.NCORE) GO TO 1006
CSA
CSA     CALL ROTIFLIN (NA,NC,XAR,XMASS,IOPT,
CSA  $                  Z(N1),Z(N3),Z(N6),Z(N9),RFLAG)
CSA     IF (NSTOP.EQ.4) GOTO 750
CSA
        RFLAG=.FALSE.
        CALL FORMDK0 (NA,NC,XAR,XMASS,RFLAG,IOPT,Z(N5),
     $                     Z(N6),Z(N8),Z(N9),
     $                     FLAG,Z(N1),Z(N4),Z(N7))
        IF (FLAG.OR.NSTOP.EQ.2) GOTO 750
        CALL FORMP (IOPT,XMASS,NA,NC,XAR,Z(N1),Z(N4),Z(N7),Z(N6),
     $              Z(N6A),Z(N6B),Z(N9),Z(N2),Z(N5),Z(N8))
        IF (NSTOP.EQ.3) GOTO 750
CSA     IF (ABS(NDISP).EQ.3) THEN
CSA     CALL PRJDISP(NC,IOPT,RFLAG,Z(N2),Z(N5),Z(N8),
CSA  $                  Z(N1), Z(N3), FLAG, Z(N4))
CSA     GOTO 750
CSA     END IF
        CALL LINTR (NC,NC,NC,IOPT,Z(N2),Z(N4),Z(N1),Z(N3),Z(N6),
     $              Z(N9),Z(N10),3000,1,IFLAG)
        CALL PROJV (NC,IOPT,Z(N2),Z(N5),Z(N8),Z(N1),Z(N4),Z(N6A),
     $              Z(N7),Z(N10),Z(N3),Z(N6),Z(N9))
CSA     CALL UNROTATE (NC,IOPT,Z(N3),Z(N6),Z(N9))
        CALL FCOUT (NC,NC,NC,0,NDER,3,0,1,Z(N1),ENERGY,Z(N1),
     $              Z(N3),Z(N6),Z(N9))
 720  CONTINUE
 750  CONTINUE
 800  continue
      ISTAT=4
      IF (NSTOP.EQ.1) ISTAT=3
      CALL RCLOSE(ISCR1, ISTAT)
      CALL RCLOSE(ISCR2, ISTAT)
      CALL RCLOSE(ISCR3, ISTAT)
      CALL RCLOSE(ISCR4, ISTAT)
      CALL RCLOSE(ISCR6, ISTAT)
      CALL RCLOSE(ISCR7, ISTAT)
      CALL RCLOSE(ISCR8, ISTAT)
      CALL RCLOSE(ISCR9, ISTAT)
      CALL RCLOSE(ISCR10,ISTAT)
C      CALL RCLOSE(ISCR11,ISTAT)
C      CALL RCLOSE(ISCR12,ISTAT)
C      CALL RCLOSE(ISCR13,ISTAT)
C      CALL RCLOSE(ISCR14,ISTAT)
CWA 4/28/99
      CALL NCLOSE(ICHECK)
      CALL NCLOSE(I11)
      CALL NCLOSE(I12)
      CALL NCLOSE(I15)
      CALL NCLOSE(I16)
      CALL NCLOSE(I17)
      CALL NCLOSE(I18)
      CALL NCLOSE(I20)
      CALL NCLOSE(I21)
      CALL NCLOSE(I24)
      CALL NCLOSE(I25)
CSA
      ASTAT='DELETE'
      CLOSE(ISCR5,STATUS=ASTAT)
      IF (NSTOP.EQ.2) ASTAT='KEEP'
         CLOSE(I31,STATUS=ASTAT)
         CLOSE(I32,STATUS=ASTAT)
         CLOSE(I33,STATUS=ASTAT)
      ASTAT='DELETE'
      IF (NSTOP.EQ.3) ASTAT='KEEP'
         CLOSE(I35,STATUS=ASTAT)
         CLOSE(I36,STATUS=ASTAT)
         CLOSE(I37,STATUS=ASTAT)
      ASTAT='DELETE'
C      IF (NSTOP.EQ.4) ASTAT='KEEP'
C         CLOSE(ISCR11,STATUS=ASTAT)
C         CLOSE(ISCR12,STATUS=ASTAT)
C         CLOSE(ISCR13,STATUS=ASTAT)
C         CLOSE(ISCR14,STATUS=ASTAT)
      RETURN
CSA
 1002 IF(IFLAG.EQ.1) THEN
         WRITE(IOUT,7)
      ELSE IF(IFLAG.EQ.2) THEN
         WRITE(IOUT,23)
         WRITE(IOUT,24)
         DO 1003  I=1,NSX
 1003    WRITE(IOUT,21) (UGF(I,J),J=1,NSX)
      END IF
      RETURN
 1004 WRITE(IOUT,8)
      RETURN
 1006 WRITE(IOUT,9)  NT,NCORE
      RETURN
 1008 WRITE(IOUT,25) IFLAG
      RETURN
 1010 WRITE(IOUT,36) IFLAG
      RETURN
 1012 WRITE(IOUT,37) IFLAG
      RETURN
 1014 WRITE(IOUT,39)
      RETURN
 1016 WRITE(IOUT,41) IFLAG
      RETURN
 1018 WRITE(IOUT,29)
      RETURN
 1020 WRITE(IOUT,44)
      RETURN
 1022 WRITE(IOUT,46)
      RETURN
      END
