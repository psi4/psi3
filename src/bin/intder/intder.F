      SUBROUTINE INTDER(NAD,NC,NS,XA,XMASS,S,SS,B,BS,A,UGF,XS,XT,XKI,
     $ TYPE,U,IA,IU,IFLAG,nst,nbnd,ntor,nout,nlin,nspf,nsym,
     2 na,lbls,z,ncore)
      IMPLICIT REAL*8 (A-H,O-Z)
#include <libparse/error.h>
      integer errcod, frdi, frdfa1, frdfa2, fwkadd, frdboo
      integer frdc, fcnt
      CHARACTER LABEL*10,TYPE*5,NUMTST*5
      character*80 lbls(nsym)
      character*8 invrt
      character*11 fqanal
      real*8 gtmp(1500)
      DIMENSION IOPT(16),TYPE(NS),IA(NS,6),IU(NS,0:1),S(NS),SS(NS)
      DIMENSION Z(ncore),U(NS,1),IR(4),XR(4),XA(NAD,3),B(NS,NC)
      DIMENSION A(NC,NC),UGF(NC,1),BS(NC,NC),XS(NC,NC)
      DIMENSION XT(NC,NC),XMASS(1),DIP(3),XKI(NC,3)
      PARAMETER(BOHR=0.52917706D0,RAD=57.29577951D0,TOLINV=1.0D-10)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
croag
      common/zpve/nvib, izvlim
croag
      COMMON /DIPOLE/DIP,QQ
Crawdad.7.18.93 Added masses through U(92)
      real*8 atm(92)
      DATA ATM / 1.007825D+00 , 4.00260D+00 , 7.01600D+00 ,
     *           9.01218D+00  , 11.00931D+00 , 12.00000D+00 ,
     *           14.00307D+00 , 15.99491D+00 , 18.99840D+00 ,
     *           19.99244D+00 , 22.98977D+00 , 23.98504D+00 ,
     *           26.98154D+00 , 27.97693D+00 , 30.97376D+00 ,
     *           31.97207D+00 , 34.96885D+00 , 39.96238D+00 ,
     *           38.96371D+00 , 39.96259D+00 , 44.95591D+00 ,
     *           47.94795D+00 , 50.94396D+00 , 51.94051D+00 ,
     *           54.93805D+00 , 55.93494D+00 , 58.93320D+00 ,
     *           57.93535D+00 , 62.92960D+00 , 63.92915D+00 ,
     *           68.92558D+00 , 73.92118D+00 , 74.92160D+00 ,
     *           79.91652D+00 , 78.91834D+00 , 83.91151D+00 ,
     + 84.9117996d+0, 87.9056249d+0, 88.9058560d+0, 89.9047080d+0,
     + 92.9063780d+0, 97.9054050d+0, 0.0d+0, 101.9043475d+0,
     + 102.905503d+0, 105.903475d+0, 106.905095d+0, 113.9033607d+0,
     + 114.903875d+0, 119.9021990d+0, 120.9038237d+0, 129.906229d+0,
     + 126.904477d+0, 131.904148d+0, 132.905433d+0, 137.905236d+0,
     + 138.906355d+0, 139.905442d+0, 140.907657d+0, 141.907731d+0,
     + 0.0d+0, 151.919741d+0, 152.921243d+0, 157.924111d+0,
     + 158.925350d+0, 163.929183d+0, 164.930332d+0, 165.930305d+0,
     + 168.934225d+0, 173.938873d+0, 174.940785d+0, 179.946561d+0,
     + 180.948014d+0, 183.950953d+0, 186.955765d+0, 191.961487d+0,
     + 192.962942d+0, 194.964785d+0, 196.966560d+0, 201.970632d+0,
     + 204.974410d+0, 207.976641d+0, 208.980388d+0, 0.0d+0,
     + 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 0.0d+0, 232.038053805d+0,
     + 0.0d+0, 238.050785782d+0/
 1    FORMAT(A10)
 3    FORMAT(16I5)
 4    FORMAT(A5,5I5,A5)
 5    FORMAT(I5,4(I4,F14.10))
 6    FORMAT(4F20.10)
 7    FORMAT(1X,'B MATRIX DETERMINANT = 0 .')
 8    FORMAT(1X,'AT EQUILIBRIUM ALL INTERNAL GRADIENTS = 0 .')
 9    FORMAT(1X,'NOT ENOUGH CORE'/1X,I10,' WORDS NEEDED AND ONLY',I10,
     $          ' WORDS AVAILABLE')
 10   FORMAT(3X,'PROGRAM OPTIONS'//
     $'   NA',5X,' NSYM',5X,'  NEQ',5X,' NINV',5X,'NTEST',5X,'NFREQ',
     $5X,' NDIP',5X,'MULTI'/
     $5X,'   NS',5X,' NDER',5X,' NPRT',5X,' NDUM',5X,'NGEOM',5X,'IRINT',
     $5X,'NSTOP')
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
     $     'ZR(M,N,P,Q) AND SR(I,J,K,L) MATRICES TO COMPUTE '/
     $     'CONTRIBUTIONS TO THE QUARTIC FORCE CONSTANTS ARISING '/
     $     'FROM NONZERO GRADIENTS.')
 23   FORMAT(/' B MATRIX INVERSION UNSUCCESSFUL FOR GIVEN TOLERANCE.')
 24   FORMAT(/' B*A MATRIX'/)
 25   FORMAT(/' NA FROM INPUT IS INCONSISTENT WITH THE HEADER IN A'/
     $    ' CARTESIAN DERIVATIVE FILE.   IFLAG=',I5)
 26   FORMAT(F20.10)
 27   FORMAT(///,1X,'VALUES OF SYMMETRY INTERNAL COORDINATES ',
     $      '(ANG. OR RAD.)'/)
 28   FORMAT(3F20.10)
 29   FORMAT(F12.6)
 30   FORMAT(//,1X,'NORMAL MODE ANALYSIS IN INTERNAL COORDINATES')
 31   FORMAT(//,1X,'NORMAL MODE ANALYSIS IN CARTESIAN COORDINATES')
 32   FORMAT(/,1X,'INTDER',/1X,'A PROGRAM TO PERFORM GENERAL ',
     $ 'CURVILINEAR TRANSFORMATIONS'/' AMONG HIGHER-ORDER DERIVATIVES'/
     $ 1X,'WRITTEN BY WESLEY D. ALLEN'/)
 33   FORMAT(//'NUCLEAR MASSES'/)
 34   FORMAT(4(I4,F12.6,4X))
 35   FORMAT(I5,F20.10)
 36   FORMAT(/,1X,'EXECUTION STOPPED DUE TO ERROR IN GFMAT.',
     $      '  IFLAG=',I5)
 37   FORMAT(/,1X,'EXECUTION STOPPED DUE TO ERROR IN NORMCO.',
     $      '  IFLAG=',I5)
 38   FORMAT(1X,  'REFERENCE VALUE IS ',F20.10,' A.')
cbfy
 39   FORMAT(///' The input options specified are not compatible with ',
     $       'multiple'/' transformations, i.e. IOPT(15) > 0 (in the ',
     $       'current version of this'/' program).'/,
     $       ' Your values : NINV =',I7,',  NGEOM =',I3,
     $       ',  NFREQ =',I3,',  IRINT =',I3/
     $       ' should be   : NINV = 0 or 1,  NGEOM =  0, |NFREQ| < 4, ',
     $       '|IRINT| < 2'/)
cbfy
      IFLAG=0
      iin1=7
      IOUT=6
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
c
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
cets050491
      iopt(1)=na
      iopt(2)=ns
      iopt(3)=nsym

      iopt(4)=2
      errcod = frdi('DERLVL',iopt(4))
      nder=iopt(4)

      iopt(5)=0
      errcod = frdboo('TRANS_1DER',iopt(5))
      neq=iopt(5)

      iopt(6)=3000
      errcod = frdi('IPRINT',iopt(6))
      nprt=iopt(6)

      iopt(7)=0
      invrt=' '
      errcod = frdc('TRANS_TYPE',invrt)
      if(invrt.eq.'C_TO_I  ') iopt(7)=0
      if(invrt.eq.'I_TO_C  ') iopt(7)=1
      if(invrt.eq.'I_TO_C_R') iopt(7)=2
      ninv=iopt(7)

      iopt(8)=0
      errcod = fcnt('DUMMY_ATOMS',iopt(8))
      ndum=iopt(8)

      iopt(9)=0
      errcod = frdi('NTEST',iopt(9))
      ntest=iopt(9)

      iopt(10)=0
      errcod = frdboo('READ_GEOM',iopt(10))
      ngeom=iopt(10)

      iopt(11)=1
      errcod = frdc('FREQ_ANAL',fqanal)
      if(fqanal.eq.'NONE       ') iopt(11)=0
      if(fqanal.eq.'INT_CO     ') iopt(11)=1
      if(fqanal.eq.'CART_CO    ') iopt(11)=2
      if(fqanal.eq.'BOTH       ') iopt(11)=3
      if(fqanal.eq.'NONE_R     ') iopt(11)=4
      if(fqanal.eq.'INT_NO_T   ') iopt(11)=-1
      if(fqanal.eq.'CART_NO_T  ') iopt(11)=-2
      if(fqanal.eq.'BOTH_NO_T  ') iopt(11)=-3
      if(fqanal.eq.'NONE_R_NO_T') iopt(11)=-4
      nfreq=iopt(11)

      iopt(12)=0
      errcod = frdc('IRINT_T',invrt)
      if(invrt.eq.'NONE    ') iopt(12)=0
      if(invrt.eq.'DO      ') iopt(12)=1
      if(invrt.eq.'DO_R    ') iopt(12)=2
      irint=iopt(12)

      iopt(13)=0
      errcod = frdboo('TRANS_DIPDER',iopt(13))
      nvec=iopt(13)

      iopt(14)=0
      errcod = frdboo('STOP',iopt(14))
      nstop=iopt(14)

croag
      nvib=3*na-5
      errcod = frdi('NVIB',nvib)

      izvlim=20
      errcod = frdi('ZPVE_LIM',izvlim)
croag

      call ffile(i11,' ',0)
      if(neq.eq.1) call ffile(i12,' ',0)

      if(nder.ge.2) then
        call ffile(i15,' ',0)
        call ffile(i16,' ',0)
      endif
      if(nder.ge.3) then
        call ffile(i20,' ',0)
        call ffile(i21,' ',0)
      endif
      if(nder.ge.4) then
        call ffile(i25,' ',0)
        call ffile(i24,' ',0)
      endif

      if(nvec.eq.1.or.irint.eq.1) then
        call ffile(i18,' ',0)
        call ffile(i17,' ',0)
      endif

cbfy
      iopt(15)=0
      errcod = frdi('MULTI',iopt(15))
      multi=iopt(15)
      IF ((MULTI.GT.1).AND.((NINV.GE.2).OR.(NGEOM.EQ.1).OR.
     $      (IABS(NFREQ).EQ.4).OR.(IABS(IRINT).GT.1))) THEN
         WRITE(IOUT,39) NINV,NGEOM,NFREQ,IRINT
         STOP
         END IF
cbfy
      if(nvec.eq.1) then
         iopt(5)=1
         neq=1
      end iF
      nc=3*na
      nad=na+ndum
      do 100  j=1,ns
        ia(j,6)=0
        if(ntest.ne.0) ia(j,6)=ntest
 100  continue

cets050591 now done from fentry
cets050491  read in simple internal coordinates
c    
c      call rdsimp(type,ia,s,ns,nst,nbnd,ntor,nlin,nout,nspf)
c
c      if(nsym.eq.0) go to 125
c
c      do 115  k=1,nsym
c 115  iu(k,0)=0
c
c now read in symmetrized internal coordinates
c      call rdsymm(lbls,iu,u,nsym,ns)
c
c 125  continue

      do 121 k=1,nsym
        m=iu(k,0)
        unorm=zero
        do 122 l=1,m
          unorm=unorm+u(k,l)*u(k,l)
 122    continue
        unorm=one/dsqrt(unorm)
        do 123 l=1,m
 123      u(k,l)=u(k,l)*unorm
 121  continue

cets
cbfy
      DO 800 IMULT=1,MAX(MULTI,1)
        IF (MULTI.EQ.0) THEN
          call bottom(i11)
          ik=2*na+1
          DO 128  K=1,IK
 128        BACKSPACE I11
        ELSE
          READ(I11,1) LABEL
        END IF
c
        READ(I11,35) NA11,ENERGY
        DO 130  J=1,NA
 130      READ(I11,6) xmass(j),(XA(J,K),K=1,3)
c
        IF ((MULTI.NE.0).AND.(NEQ.EQ.0)) THEN
c          Skip over 1st derivatives in File11
           DO 220 K=1,NA
 220          READ(I11,1) LABEL
        END IF
c
cets050491 get geometry from geom.dat file
        IF (NGEOM.EQ.1) THEN
          errcod = fwkadd(8,'DUMMY')
          errcod = frdfa2(':GEOMETRY',gtmp,3,na)
          if (errcod.ne.EOK) then
            write(6,*) 'ERROR: problem reading geometry from geom.dat'
            call pabort
            call qabort
          endif
          ij=0
          do 210 j=1,na
            do 210 k=1,3
              ij=ij+1
              xa(j,k)=gtmp(ij)
  210     continue
cets050491          DO 210  J=1,NAD
cets050491 210      READ(IIN1,28) (XA(J,K),K=1,3)
        END IF
        if(ndum.ne.0) then
          ij=0
          errcod = frdfa2('DUMMY_ATOMS',gtmp,3,ndum)
          do 211 j=na+1,nad
            do 211 k=1,3
              ij=ij+1
              xa(j,k)=gtmp(ij)
  211     continue
        endif
c
        DO 214  J=1,NS
          IF(TYPE(J).EQ.'  SPF'.AND.IA(J,3).EQ.0) THEN
            CALL VECT1(NAD,IA(J,1),IA(J,2),XT(1,1),XA,S(J))
            S(J)=S(J)*BOHR
          END IF
 214    CONTINUE
c
        if(nfreq.ne.0.or.nvec.gt.0) then
          errcod = frdfa1('MASSES',xmass,na)
          if(errcod.ne.EOK) then
            do 213 i=1,na
              j=dint(xmass(i))
              xmass(i)=atm(j)
 213        continue
          endif
        endif
        
          
        DO 212  I=1,NA
          IF(NVEC.GT.0) THEN
            XT(I,1)=XMASS(I)
          ELSE
            XT(I,1)=ONE
          END IF
 212    CONTINUE
c     
        WRITE(IOUT,32)
        WRITE(IOUT,10)
        WRITE(IOUT,3)(IOPT(K),K=1,16)
        WRITE(IOUT,11)
        DO 133  J=1,NAD
 133    WRITE(IOUT,12) (XA(J,K),K=1,3)
        DO 131  I=1,NAD
          DO 131  J=1,3
 131        XA(I,J)=XA(I,J)*BOHR
        WRITE(IOUT,13)
        WRITE(IOUT,14)
        DO 134  J=1,NS
          NUMTST='     '
          IF(NTEST.NE.0.AND.IA(J,6).EQ.0) NUMTST='   ST'
          WRITE(IOUT,17) J,TYPE(J),(IA(J,K),K=1,4),NUMTST
          IF(TYPE(J).EQ.'  SPF') WRITE(IOUT,38) S(J)
 134    CONTINUE
        DO 136  J=1,NS
          DO 136  I=1,NC
 136        UGF(I,J)=ZERO
        IF(NSYM.EQ.0) GO TO 145
          WRITE(IOUT,15)
          DO 137  I=1,NSYM
            LL=IU(I,0)
            WRITE(IOUT,16) I,(U(I,L),IU(I,L),L=1,LL)
            DO 137  L=1,LL
              K=IU(I,L)
 137          UGF(I,K)=U(I,L)
c
 145    IF(LPRT(4,NPRT).LT.4) GO TO 150
          DO 141  J=1,NAD
 141      WRITE(ICHECK,142) (XA(J,K),K=1,3)
 142      FORMAT(3F10.7)
          DO 138  I=1,NSYM
 138      WRITE(ICHECK,139) (UGF(I,K),K=1,NS)
 139      FORMAT(8F10.7)
C
 150    CALL MACHB(NAD,NC,NS,XA,TYPE,IA,B,S)
        IF(NSYM.LE.0) THEN
          DO 154  I=1,NS
            SS(I)=S(I)
            DO 154  J=1,NC
              BS(I,J)=B(I,J)
 154      continue
        ELSE
          DO 156  J=1,NC
            DO 156  I=1,NSYM
              BS(I,J)=ZERO
              DO 158  L=1,IU(I,0)
                K=IU(I,L)
                BS(I,J)=BS(I,J)+U(I,L)*B(K,J)
 158          CONTINUE
 156      CONTINUE
          DO 160  I=1,NSYM
            SS(I)=ZERO
            DO 160  L=1,IU(I,0)
            K=IU(I,L)
            SS(I)=SS(I)+U(I,L)*S(K)
 160      continue
        END IF
        IF(NVEC.GT.0) THEN
           WRITE(IOUT,33)
           WRITE(IOUT,34) (I,XMASS(I),I=1,NA)
        END IF
        NSX=NS
        IF(NSYM.GT.0) NSX=NSYM
        NSY=NSX
        DO 162  I=1,NS
          UGF(I,1)=S(I)
          IF(TYPE(I).EQ.' BEND'.OR.TYPE(I).EQ.' LIN1'.OR.TYPE(I).EQ.
     $                             ' TORS'.OR.TYPE(I).EQ.'  OUT') THEN
            UGF(I,1)=UGF(I,1)*RAD
          END IF
 162    CONTINUE
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
            WRITE(IOUT,21) (BS(I,J),J=1,NC)
 170      continue
        END IF
c
        IF(NFREQ.LT.0) GO TO 600
c
        IF(NINV.GT.0) THEN
          NSY=NC
          IF(NTEST.EQ.0) GO TO 175
        END IF
        CALL BINVRT(NSX,NC,XT(1,1),BS,UGF,A,IFLAG,TOLINV,1)
c
        IF(IFLAG.NE.0) GO TO 1002
C
 175    N1=1
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
 180    NMM=NCHUNK*3
        NT=N5+NMM
        MDER=NDER+NEQ
        IF(MDER.GE.4) THEN
          NNN=NC**3
          N6=N5+NNN
          NMM=MAX0(NMM,2*NNN)
          NT=N5+NMM
        END IF
        IF(MDER.GE.5) THEN
          NNN=NC**4
          N7=NT
          N8=N7+NNN
          NT=N8+NNN
        END IF

        IF(NT.GT.NCORE) GO TO 1006

        IF(NINV.EQ.2) THEN
          CALL INPFKM(NC,NDER,NEQ,NSX,Z(N1),Z(N2),Z(N3),Z(N4))
        END IF

        IF(MDER.GE.3) THEN
          IF(ABS(NTEST).EQ.1) THEN
            CALL NUMX(NAD,NC,NS,NSX,XA,TYPE,IA,XT(1,1),A,S,UGF,XS)
          END IF
          CALL MACHX(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,U,IU,UGF,XS)
          IF(NTEST.EQ.1) THEN
          CALL SRTST1(NC,NS,NSX,NSYM,IA,XS,XT)
          END IF
        END IF

        IF (MDER.GE.4) THEN
          IF(ABS(NTEST).EQ.2) THEN
            CALL NUMY(NAD,NC,NS,NSX,XA,TYPE,IA,A,S,Z(N5),Z(N6),UGF)
          END IF
          CALL MACHY(
     +             NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,U,IU,Z(N5),Z(N6))
          IF(NTEST.EQ.2) THEN
            CALL SRTST2(NC,NS,NSX,NSYM,IA,Z(N5),Z(N6))
          END IF
        END IF
c
        IF (MDER.GE.5) THEN
          CALL NUMZ(NAD,NC,NS,NSX,XA,TYPE,IA,A,S,Z(N5),Z(N7),Z(N8))
          CALL MACHZ(
     +              NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,U,IU,Z(N7),Z(N8))
        END IF
        IF(NSTOP.EQ.1) GO TO 700
C
        NRUN=0
        IF(NVEC.GT.0) THEN
          CALL VCXKI(NA,NAD,NC,XMASS,XA,XKI)
        END IF
 200    NRUN=NRUN+1
        IF(NINV.LE.0) THEN
            CALL LINTR(NC,NSX,NC,IOPT,A,XS,Z(N1),Z(N2),Z(N3),Z(N4),
     $         z(n5),NMM,NRUN,IFLAG)
        ELSE
        CALL LINTR(NC,NC,NSX,IOPT,BS,XS,Z(N1),Z(N2),Z(N3),Z(N4),
     $         z(n5),NMM,NRUN,IFLAG)
        END IF
c
        IF(IFLAG.NE.0) GO TO 1008
        IF(NVEC.GT.0) THEN
          CALL VCDER1(NA,NC,NSX,NINV,XMASS,XKI,A,Z(N1),Z(N2),
     $        Z(N3),Z(N4),Z(N5),NRUN)
        END IF
c
        IF(NDER.LE.1) GO TO 500
        IF(NDER.EQ.2.AND.NEQ.EQ.0) GO TO 500
c
        IF(NEQ.NE.0) THEN
           CALL MACHF2(NC,NSX,NINV,Z(N1),Z(N2),UGF,Z(N5))
        END IF
c
        IF(NDER.LE.2) GO TO 500
c
C add x-v(m,n) contributions to 3rd and perhaps 4th derivatives.
c
        CALL XF2(NDER,NC,NSX,NINV,BS,Z(N2),Z(N3),Z(N4),Z(N5),UGF,XS,XT)
        IF(NDER.EQ.3.AND.NEQ.EQ.0) GO TO 500
c
C FORM Y MATRICES AND WRITE TO DISK.
C COMPLETE V(M,N,P) COMPUTATION.
c
        IF(NEQ.NE.0) THEN
           CALL MACHF3(NC,NSX,NINV,Z(N1),Z(N3),Z(N5),Z(N6))
        END IF
        IF(NDER.LE.3) GO TO 500
c
C      ADD X-V(M,N,P) CONTRIBUTIONS TO 4TH DERIVATIVES.
c
        CALL XF3(NC,NSX,NINV,BS,UGF,XS,XT,Z(N3),Z(N4),Z(N5))
c
C      ADD Y-V(M,N) CONTRIBUTIONS TO 4TH DERIVATIVES.
c
        CALL YF2(NC,NSX,NINV,BS,UGF,XS,Z(N2),Z(N4),Z(N5))
        IF(NDER.EQ.4.AND.NEQ.EQ.0) GO TO 500
c
        WRITE(IOUT,22)
C   FORM Z MATRIX ELEMENTS NUMERICALLY AND
C   COMPLETE THE V(M,N,P,Q) COMPUTATION.
c
        CALL MACHF4(NC,NSX,NINV,Z(N1),Z(N4),Z(N7),Z(N5))
 
 500    CALL FCOUT(NC,NS,NSY,NEQ,NDER,NINV,NVEC,NRUN,SS,ENERGY,
     $           Z(N1),Z(N2),Z(N3),Z(N4))
 
        IF(NVEC.EQ.1.AND.NRUN.LT.3) GO TO 200
 
 600    NFREQ=IABS(NFREQ)

        IF(NFREQ.NE.0.AND.NFREQ.NE.2) THEN
           WRITE(IOUT,30)
cbfy
           IF (IOPT(11).GT.0) THEN
c           Reposition File16
             NLINES = (NSX*NSX +2)/3 + 1
             DO 610 I=1,NLINES
 610           BACKSPACE I16
           END IF
cbfy
           CALL GFMAT(NA,NAD,NC,NSX,NFREQ,IRINT,IA,
     $           XMASS,XA,S,BS,UGF,A,XS,XT,XS,IFLAG,lbls)
           IF(IFLAG.NE.0) GO TO 1010
        END IF
        IF (NFREQ.EQ.2.OR.NFREQ.EQ.3) THEN
          WRITE(IOUT,31)
cbfy
          IF (IOPT(11).GT.0) THEN
c           Reposition File15
            NLINES = (NC*NC +2)/3 + 1
            DO 620 I=1,NLINES
 620           BACKSPACE I15
          END IF
cbfy
          CALL NORMCO(NA,NAD,NC,NSX,NFREQ,IRINT,
     $            XMASS,XA,BS,UGF,A,XS,XT,XS,IFLAG)
          IF(IFLAG.NE.0) GO TO 1012
        END IF
c
 700    CONTINUE

        call srew(iscr1)
        call srew(iscr2)
        call srew(iscr3)
        call srew(iscr4)
        call srew(iscr6)
        call srew(iscr7)
        call srew(iscr8)
        call srew(iscr9)
        call srew(iscr10)
 800  CONTINUE

      CALL RCLOSE(ISCR1,4)
      CALL RCLOSE(ISCR2,4)
      CALL RCLOSE(ISCR3,4)
      CALL RCLOSE(ISCR4,4)
      CALL RCLOSE(ISCR6,4)
      CALL RCLOSE(ISCR7,4)
      CALL RCLOSE(ISCR8,4)
      CALL RCLOSE(ISCR9,4)
      CALL RCLOSE(ISCR10,4)
      RETURN
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
      END
