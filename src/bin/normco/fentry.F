cets The following is called by a C program which handles dynamic
cets allocation of memory.
      subroutine fentry(cc,ia,maxcor)
cets  PROGRAM NORMCO
C     The NORmal COordinate analysis program
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
C*  By: Evan Bolton						      *
C*  Date: Jan 20 1995                                                 *
C*  Reason: Put back in Tracy Hamilton's addition of GAPT charges     *
C*          from AUGUST 29, 1989.  SEE J. CIOSLOWSKI PAPERS.          *
C*          Uses the trace of the Dipole Moment Derivatives.          *
C*          See comments Ceeb.                                        *
C**********************************************************************
C*  BY:  George Vacek (croag)
C*  DATE: Aug 1 1994
C*  REASON: intder and normco read/write files with different format
C*********************************************************************
C*  BY:  T. Daniel Crawford (Crawdad.9.19.93)                         *
C*  DATE: September 19, 1993                                          *
C*  REASON: Force Raman Depolarization Ratios to zero when necessary. *
C**********************************************************************
C*  BY:  T. Daniel Crawford (Crawdad.7.18.93)                         *
C*  DATE: July 18, 1993                                               *
C*  REASON: Added elements 37-92 to mass list.  Some elements were    *
C*          excluded due to inavailability of accurate data.          *
C*          Same reference as BJDheen's from 1988.                    *
C**********************************************************************
C*  BY:  T. Daniel Crawford (Crawdad.6.10.93)			      *
C*  DATE: June 10, 1993						      *
C*  REASON: Correct inaccurate label for principal moments of         *
C*		inertia and rotational constants.  Mainly for         *
C*		beautification.	                                      *
C**********************************************************************
C*  BY:  Edward Seidl (cets030991)                                    *
C*  DATE:  March 9 1991						      *
C*  REASON: use new input and print control routines                  *
C**********************************************************************
C*  BY:  Edward Seidl (cets111690)                                    *
C*  DATE:  November  16,  1990                                        *
C*  REASON: change nvib to n3n-5 in calculation of zpve so others     *
c*          won't make the same stupid mistake I did                  *
c*          also changed code to reset idipol and ipolar if           *
c*          file17 or file18 are empty                                *
C**********************************************************************
C* BY:  BRADLEY DELEEUW                                               *
C* DATE:  JUNE 27, 1988                                               *
C* REASON:  UPDATE ATOMIC MASSES OF PRINCIPAL ISOTOPES.  ELEMENTS     *
C*        19-36 ADDED.  ELEMENT 11 CHANGED FROM 22.9898 TO 22.98977.  *
C*        ELEMENT 13 CHANGED FROM 26.98153 TO 26.98154.  ELEMENT 17   *
C*        CHANGED FROM 34.95885 TO 34.96885.  ELEMENT 18 CHANGED      *
C*        FROM 39.948 TO 39.96238.                                    *
C* DATA SOURCE:  J. CHEM. PHYS. REF. DATA, VOL.13, NO.3, 1984, p.309  *
C**********************************************************************
C*  BY:  RICHARD REMINGTON                         search:  c3-18-88  *
C*  DATE:  MARCH  18,  1988                                           *
C*  REASON: DECREASE CORE TO RUN IN 7MB ON 9370                       *
C**********************************************************************
C*  LAST UPDATED ON DECEMBER 28, 1987 BY YUKIO YAMAGUCHI              *
C**********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
cets030991
#include <libparse/error.h>
      integer frdc, frdboo, frdi, frdfa2, fcnt, fcnt1
      integer prdup, prcntr
      integer ip, errcod, niso, nat
      integer ifx, igf
      logical oldinp
c3-18-88  DIMENSION CC(360000),IA(1)
cets  DIMENSION CC(360000),IA(1)
cets  dimension cc(maxcor),ia(intowp(maxcor))
      dimension cc(maxcor),ia(1)
      DIMENSION I30(200),COORD(3,50),GRAD(3,50),WISO(20,50)
      real*8 wtemp(1000),gtemp(1500)
C   THE NORMAL COORDINATE ANALYSIS
C   FX---CARTESIAN COORDINATE SYSTEM
C   GF---INTERNAL COORDINATE SYSTEM
C   THIS PROGRAM IS able TO TREAT MOLECULES WITH '50' ATOMS
      LOGICAL CALEND,GEOINP,ISOMOL
Crawdad.7.18.93 Increase mass array dimension through U(92)
      DIMENSION ATM(92)
Ceeb  For GAPT charges.
      real*8 GAPT(50)
      COMMON/VIB101/NATOM,NAT1,N3N,ILIN,NVIB,IDIPOL,IPOLAR
      COMMON/VIB102/LDIM,MDIM,NDIM
      COMMON/VIB103/NST,NBND,NLIB,NDEF,NTORS
      COMMON/VIB104/PARA,WAVE,CONST,CYCL,CONV
      COMMON/VIB105/AIXX,AIYY,AIZZ,SUMW,ROTCM(3)
      COMMON/VIB106/NSTORE,NPLUS
      COMMON/VIB107/NIC,NINT,ISYM,NSYM
      COMMON/VIB108/IZVLIM,IPLOT
      COMMON/VIB109/ITHERM,IQELEC
      COMMON/VIB201/CHG(50),X(50),Y(50),Z(50),W(50)
      COMMON/VIB202/R(1275)
      COMMON/VIB203/ENGX(50),ENGY(50),ENGZ(50)
      COMMON/VIB204/DIPDR(3,150),POLDR(6,150)
      COMMON/VIB205/IOFF(150)
      COMMON/VIB206/SQM(150),ROOT(150),FREQ(150)
      COMMON/VIB207/VZERO
      COMMON/VIB301/KR(150),LR(150)
      COMMON/VIB302/KA(150),LA(150),MA(150)
      COMMON/VIB303/KB(150),LB(150),MB(150)
      COMMON/VIB304/KD(150),LD(150),MD(150),ND(150)
      COMMON/VIB305/KT(150),LT(150),MT(150),NT(150)
      COMMON/VIB306/AGLB(150),AGLT(150),AGLO(150),AGLR(150)
      COMMON/VIB307/AGLL(150),AGLD(150),AGLE(150),AGLF(150)
      COMMON/VIB308/IIS(150),JJS(150),KKS(150),LLS(150)
cets  EQUIVALENCE (CC,IA)
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
Crawdad.7.18.93.  Adding Atomic masses through Uranium (92)
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
      DATA PH,BK,AVN / 6.626176D+00 , 1.380662D+00 , 6.022045D+00 /
      DATA A0,HE,EV / 0.52917706D+00 , 4.359813653D+00 , 1.6021892D+00 /
      DATA CL,TEMP / 2.99792458D+00 , 298.15D+00 /
      DATA PI,ALN / 3.1415926536D+00 , 2.3025850930D+00 /
Ceeb  For GAPT Charges
      DATA DANG3 / 1.4409726D+01 /
Ceeb  
      DATA NRED / 10 /
      DATA ITAP15,ITAP17,ITAP18 / 15 , 17 , 18 /
    1 FORMAT(//,2X,' THE NORMAL COORDINATE ANALYSIS PROGRAM'/)
    2 FORMAT(14I5)
    3 FORMAT(//,2X,' PARAMETERS FOR THE CALCULATION'//
     * 2X,' IFXGF  = ',I5/
     * 2X,' ISOTOP = ',I5/
     * 2X,' NVIB   = ',I5/
     * 2X,' IFORCE = ',I5/
     * 2X,' IDIPOL = ',I5/
     * 2X,' IPOLAR = ',I5/
     * 2X,' IGEOMT = ',I5/
     * 2X,' ITHERM = ',I5/
     * 2X,' IQELEC = ',I5/
     * 2X,' IZVLIM = ',I5/
     * 2X,' IPLOT  = ',I5/
     * 2X,' IPRNT  = ',I5/
     * 2X,' NATOM  = ',I5/
     * 2X,' N3N    = ',I5/
     * 2X,' NSTORE = ',I5/
     * 2X,' NPLUS  = ',I5)
    4 FORMAT(//,2X,' SCF ENERGY = ',F20.10/)
    5 FORMAT(//,2X,' MOLECULAR GEOMETRY (IN A.U.)'/
     1 5X,4H NO.,13X,2H X,18X,2H Y,18X,2H Z,18X,2H W/)
    6 FORMAT(4F20.10)
    7 FORMAT(2X,I7,5X,4F20.10)
    8 FORMAT(//,2X,' MOLECULAR GEOMETRY (IN A)'/
     1 5X,4H NO.,13X,2H X,18X,2H Y,18X,2H Z,18X,2H W/)
    9 FORMAT(3F20.10)
   10 FORMAT(//,2X,' FIRST DERIVATIVES OF DIPOLE MOMENTS (IN DEBYE/A)'/)
   11 FORMAT(//,2X,' FIRST DERIVATIVES OF POLARIZABILITIES (IN A+2)'/)
   12 FORMAT(2I5)
   13 FORMAT(3I5)
   14 FORMAT(4I5)
   15 FORMAT(F20.10)
   16 FORMAT(/1H1//,2X,' THE NORMAL COORDINATE ANALYSIS PROPGRAM'/)
   17 FORMAT(//,3X,22(1H*)/2X,' VIBRATIONAL AMPLITUDES'/3X,22(1H*))
   18 FORMAT(//,3X,20(1H*)/2X,' INFRARED INTENSITIES'/3X,20(1H*))
   19 FORMAT(//,3X,17(1H*)/2X,' RAMAN INTENSITIES'/3X,17(1H*))
   20 FORMAT(//,2X,' THERE ARE NOT ENOUGH GRADIENTS TO CALCULATE NORMAL
     1COORDINATES'/)
Ceeb  BEGIN Necessary Formats for Tracy Hamilton's addition of GAPT charges
   21 FORMAT(/,'   ATOM      NET CHARGES',/,
     1         '   ----      -----------')
   22 FORMAT(1X,I5,6X,F10.5)
   23 FORMAT(/,'   THE SUM OF THE GAPT CHARGES IS ',F10.5,/)
   24 FORMAT(//,'   THE GAPT CHARGES DERIVED FROM THE DIPOLE MOMENT DERI
     1VATIVES',/)
Ceeb  END of GAPT Charge Format statements
C
cets030991      call drum
      call psinit('APPEND')
      CALL TSTART(6)
cets  CALL NOUNFL
C
cets  open (unit=15,status='unknown',file='file15')
cets  open (unit=17,status='unknown',file='file17')
cets  open (unit=18,status='unknown',file='file18') 
      call ffile(itap15,' ',0)
C   PHYSICAL CONSTANTS
      PARA=1.0D+00/AVN
      WAVE=1.0D+04*DSQRT(AVN)/(2.0D+00*PI*CL)
      CONST=1.0D+02*(PH*AVN)/(8.0D+00*PI*PI*CL)
      CYCL=1.0D+06*(PH*AVN)/(8.0D+00*PI*PI)
      CONV=1.0D-01*(PH*CL)/(2.0D+00*BK*TEMP)
C
      LDIM=1275
      MDIM=150
      NDIM=50
      INPUT=5
      ITAPE6=6
      ITAP30=30
c3-18-88  MAXCOR=360000
cets  MAXCOR=360000
C
      CALL LOCATE(INPUT,'# NORMCO #',IERR)
C
      WRITE(6,1)
C
C***INPUT PARAMETERS FOR THE CALCULATION***
C   IFXGF IS A PARAMETER FOR THE METHOD
C     IFXGF =1   :FX MATRIX METHOD
C     IFXGF =2   :GF MATRIX METHOD
C   ISOTOP IS A PARAMETER FOR ISOTOPIC MOLECULES
C     ISOTOP =0 AND 1 :USE REGULAR ATOMIC MASS
C     ISOTOP =2 OR MORE :NUMBER OF ISOTOPIC MOLECULES
C   NVIB IS DEGREES OF VIBRATIONAL FREEDOM
C   IFORCE IS A PARAMETER FOR GRADIENT SETS
C     IFORCE =1 :READ IN GRADIENT SETS FROM '15' FILE
C     IFORCE =2 :READ IN GRADIENT SETS FROM INPUT FILE
C     IFORCE =3 :READ IN A HESSIAN FROM '15' FILE
C   IDIPOL IS A PARAMETER FOR DIPOLE DERIVATIVES
C     IDIPOL =0  :NO DIPOLE DERIVATIVES
C     IDIPOL =1  :READ IN DIPOLE DERIVATIVES
C   IPOLAR IS A PARAMETER FOR POLARIZABILITY DRIVATIVES
C     IPOLAR =0  :NO POLARIZABILITY DERIVATIVES
C     IPOLAR =1  :READ IN POLARIZABILITY DERIVATIVES
C   IGEOMT IS A PARAMETER FOR INPUT DATA
C     IGEOMT =0 :READ IN GEOMETRY FROM TAPE 30
C     IGEOMT =1 :READ IN GEOMETRY FROM TAPE 5
C   ITHERM IS A PARAMETER FOR THERMOCHEMISTRY
C     ITHERM =0 AND 1 :USE 1 FOR SYMMETRY PARAMETER
C     ITHERM =N :FOR A MOLECULE WITH SYMMETRY NUMBER N
C   IQELEC IS A PARAMETER FOR MULTIPLICITY OF ELECTRONIC STATE
C     IQELEC =0 AND 1 :USE 1 FOR SINGLE MULTIPLICITY
C     IQELEC =N :FOR MULTIPLICITY OF ELECTRONIC STATE
C   IZVLIM IS A PARAMETER FOR LOWER LIMIT TO THE ZPVE CORRECTION
C     IZVLIM =0 :USE DEFAULT VALUE OF 20 CM-1
C     IZVLIM =N :SPECIFY LOWER LIMIT OF ZPVE
C   IPLOT  IS A PARAMETER FOR THERMOCHEMISTRY
C     IPLOT  =0 : DO NOT FORM FILE14
C     IPLOT  =1 : FORM FILE14 FOR PLOTTING PURPOSE
cets030991
      if(ierr.eq.0) then
        oldinp=.true.
        READ(5,2) IFXGF,ISOTOP,NVIB,IFORCE,IDIPOL,IPOLAR,IGEOMT,ITHERM,
     1          IQELEC,IZVLIM,IPLOT,IPRNT
        IF(ISOTOP.EQ.0) ISOTOP=1
        IF(ITHERM.EQ.0) ITHERM=1
        IF(IQELEC.EQ.0) IQELEC=1
        IF(IZVLIM.EQ.0) IZVLIM=20
      else
        oldinp=.false.
        iforce=3
        iplot=0

        idipol=1
        errcod = frdboo('IR_INT',idipol)
        
        ipolar=0
        errcod = frdboo('RAMAN_INT',ipolar)

        igf = 0
        errcod = frdboo('GF',igf)

        ifx = 0
        errcod = frdboo('FX',ifx)

        if(iand(igf,ifx).eq.1) then
          write(6,*) ' cannot do both fx and gf methods'
          write(6,*) ' defaulting to fx method'
          ifxgf=1
        else if(igf.eq.1) then
          write(6,*) ' cannot do gf method with new input yet'
          write(6,*) ' defaulting to fx method'
          write(6,*) ' if you really want gf, use the old-style input'
          ifxgf=1
        else
          ifxgf=1
        endif

        call rfile(30)
        call wreadw(itap30,i30,200,101,junk)
        call rclose(30,3)
        natom=i30(19)
        nvib=3*natom-5
        errcod = frdi('NVIB',nvib)

        iprnt = 0
        errcod = frdi('IPRINT',iprnt)

        itherm = 1
        errcod = frdi('ITHERM',itherm)
 
        iqelec = 1
        errcod = frdi('IQELEC',iqelec)

        izvlim=20
        errcod = frdi('ZPVE_LIM',izvlim)

        isotop=1
        errcod = fcnt('ISOTOPES',niso)
        if(errcod.eq.EOK) isotop=niso+1
        if(isotop.gt.1) then
          errcod = fcnt1('ISOTOPES',nat,niso)
          if(errcod.ne.EOK.or.nat.ne.natom) then
            write(6,*) ' number of atoms ne natom'
            call pabort
            call qabort
          endif

          errcod = frdfa2('ISOTOPES',wtemp,nat,niso)
          if(errcod.ne.EOK) then
            write(6,*) 'problems reading in isotopic masses'
            call pabort
            call qabort
          endif
          
          ij=0
          do 3001 j=2,niso+1
            do 3002 i=1,nat
              ij=ij+1
              wiso(j,i) = wtemp(ij)
 3002       continue
 3001     continue
        endif

        igeomt=0
        errcod = frdboo('READ_GEOM',igeomt)
        if(igeomt.eq.1) then
          call ffile(9,'geom',0)
          errcod = fwkadd(9,'DUMMY')
          errcod = frdfa2(':GEOMETRY',gtemp,3,natom)
        endif
      endif
c
      errcod = prcntr('ON BRIEF')
      errcod = prcntr('OFF VERBOSE DEBUG')
      errcod = prdup()

      if(prcntr('IS_ON VERBOSE').eq.1.or.prcntr('IS_ON DEBUG').eq.1)
     +  then
          errcod = prcntr('OFF BRIEF')
      endif
C
      REWIND ITAP15
      READ(ITAP15,2) NATOM,NSTORE
      NAT1=NATOM-1
      N3N=NATOM*3
      CALEND=NSTORE.GE.N3N*2
      NPLUS=NSTORE
      GEOINP=.NOT.CALEND
      GEOINP=IGEOMT.NE.0
      IF(.NOT.CALEND) NPLUS=NSTORE+1
      CALEND=NPLUS.GE.N3N*2
      WRITE(6,3) IFXGF,ISOTOP,NVIB,IFORCE,IDIPOL,IPOLAR,IGEOMT,ITHERM,
     1 IQELEC,IZVLIM,IPLOT,IPRNT,NATOM,N3N,NSTORE,NPLUS
C
      IOFF(1)=0
      DO 101 I=1,149
      IOFF(I+1)=IOFF(I)+I
  101 CONTINUE
C
      CALL RFILE(ITAP30)
      CALL WREADW(ITAP30,I30,200,101,JUNK)
      IEND=I30(1)
      MPOINT=I30(2)
      MCONST=I30(3)
      MCALCS=I30(4)
      NCALCS=I30(5)
      NATOM=I30(19)
      NAT1=NATOM-1
      N3N=NATOM*3
      JUNK=101+MCONST
      CALL WREADW(ITAP30,IPT,1,JUNK,JUNK)
      CALL WREADW(ITAP30,CHG,NATOM*2,IPT,JUNK)
      IGEOP=100+MCONST+MPOINT+NCALCS
      CALL WREADW(ITAP30,LOCCAL,1,IGEOP,JUNK)
      IGEOP=LOCCAL+60
      CALL WREADW(ITAP30,I30,20,IGEOP,IGEOP)
      CALL WREADW(ITAP30,COORD,N3N*2,IGEOP,JUNK)
      JUNK=JUNK+INTOWP(1)
      CALL WREADW(ITAP30,ESCF,2,JUNK,JUNK)
      JUNK=I30(2)
      CALL WREADW(ITAP30,GRAD,N3N*2,JUNK,JUNK)
C
C   REWIND INPUT FILE
      if(oldinp) then
        CALL LOCATE(INPUT,'# NORMCO #',IERR)
        READ(5,2) IDUM1,IDUM2,IDUM3,IDUM4,IDUM5,IDUM6,IDUM7,IDUM8,IDUM9,
     1          IDUM10
      endif
C
      WRITE(6,4) ESCF
cets111690
      if(nvib.ne.(n3n-5).and.nvib.ne.(n3n-6)) then
        write(6,*) ' '
        write(6,*) '        !!!!!!!! Warning !!!!!!!!!!!'
        write(6,*) '      nvib is not equal to 3n-5 or 3n-6 '
        write(6,*) '   are you sure you know what you are doing?'
        write(6,*) ' '
      endif
      WRITE(6,5)
      ij=0
      DO 102 I=1,NATOM
        if(geoinp.and.oldinp) go to 201
        if(geoinp) then
          x(i)=gtemp(ij+1)
          y(i)=gtemp(ij+2)
          z(i)=gtemp(ij+3)
          ij=ij+3
        else
          X(I)=COORD(1,I)
          Y(I)=COORD(2,I)
          Z(I)=COORD(3,I)
        endif
        II=CHG(I)
        W(I)=ATM(II)
        GO TO 202
  201   CONTINUE
        READ(5,6) COORD(1,I),COORD(2,I),COORD(3,I),W(I)
        X(I)=COORD(1,I)
        Y(I)=COORD(2,I)
        Z(I)=COORD(3,I)
  202   CONTINUE
        WISO(1,I)=W(I)
        WRITE(6,7) I,X(I),Y(I),Z(I),W(I)
  102 CONTINUE
      CALL DIST(1)
      WRITE(6,8)
      DO 103 I=1,NATOM
      X(I)=X(I)*A0
      Y(I)=Y(I)*A0
      Z(I)=Z(I)*A0
      WRITE(6,7) I,X(I),Y(I),Z(I),W(I)
  103 CONTINUE
      CALL DIST(1)
      DO 104 I=1,NATOM
      ENGX(I)=GRAD(1,I)
      ENGY(I)=GRAD(2,I)
      ENGZ(I)=GRAD(3,I)
  104 CONTINUE
C
C   FORM FORCE CONSTANT MATRIX
      N3N2=N3N*N3N
      IC1=1
      IC2=IC1+N3N2
      IC3=IC2+N3N2
      IC4=IC3+N3N2
      IC5=IC4+N3N2
      ICMAX=IC5+N3N2
C.................FX      FP      FM      FF      AA.................
      CALL FCONST(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),N3N,IFORCE)
C
      IF(.NOT.CALEND) GO TO 210
C
C   THE CALCULATION OF CENTER OF MASS AND MOMENTS OF INTERTIA
      CALL MOMENT
C
C   READ IN DIPOLE MOMENT DERIVATIVES
      IF(IDIPOL.EQ.0) GO TO 203
      call ffile(itap17,' ',0)
      REWIND ITAP17
croag
croag Intder writes these files with ICHG instead of N3N, so it's best
croag just not to read it in here.
croag
      READ(ITAP17,2,end=298) NATOM
      n3n = natom*3
      READ(ITAP17,9) ((DIPDR(I,J),J=1,N3N),I=1,3)
      WRITE(6,10)
      CALL MATOUT(DIPDR,3,150,3,N3N,6)

Ceeb  BEGIN OF Tracy Hamilton's GAPT charge calculation
CTPH                                                                    NOR02910
C   CALCULATE GAPT CHARGES (DANG3 IS DEBYE*ANGSTROM*3)                  NOR02920
      WRITE(6,24)                                                       NOR02930
      JOFF = 0                                                          NOR02940
      DO 212 J=1,NATOM                                                  NOR02950
        GAPT(J) = 0.0D0                                                 NOR02960
        DO 211 I=1,3                                                    NOR02970
          GAPT(J) = GAPT(J) + DIPDR(I,JOFF+I)                           NOR02980
  211   CONTINUE                                                        NOR02990
          JOFF = JOFF + 3                                               NOR03000
        GAPT(J) = GAPT(J)/DANG3                                         NOR03010
  212 CONTINUE                                                          NOR03020
      WRITE(6,21)                                                       NOR03030
      SUM = 0.0D0                                                       NOR03040
      DO 213 I=1,NATOM                                                  NOR03050
        WRITE(6,22) I,GAPT(I)                                           NOR03060
        SUM = SUM + GAPT(I)                                             NOR03070
  213 CONTINUE                                                          NOR03080
      WRITE(6,23) SUM                                                   NOR03090
CTPH                                                                    NOR03100
Ceeb  END OF Tracy Hamilton's GAPT charge calculation

      REWIND ITAP17
      goto 299
 298  continue
      write(6,*) ' '
      write(6,*) ' Cannot find dipole moment derivatives'
      write(6,*) ' so I will not calculate intensities'
      write(6,*) ' '
      idipol=0
 299  continue
C
C   READ IN POLARIZABILITY DERIVATIVES
  203 CONTINUE
      IF(IPOLAR.EQ.0) GO TO 204
      call ffile(itap18,' ',0)
      REWIND ITAP18
croag
croag Intder writes these files with ICHG instead of N3N, so it's best
croag just not to read it in here.
croag
      READ(ITAP18,2,end=296) NATOM
      n3n = natom*3
      READ(ITAP18,9) ((POLDR(I,J),J=1,N3N),I=1,6)
      WRITE(6,11)
      CALL MATOUT(POLDR,6,150,6,N3N,6)
      REWIND ITAP18
      goto 297
 296  continue
      write(6,*) ' '
      write(6,*) ' Cannot find polarizability derivatives'
      write(6,*) ' so I will not calculate raman intensities'
      write(6,*) ' '
      ipolar=0
 297  continue
C
C   READ IN INTERNAL COORDINATES
  204 CONTINUE
      IF(IFXGF.NE.2) GO TO 205
      READ(5,2) NST,NBND,NLIB,NDEF,NTORS,ISYM
      NINT=NST+NBND+NLIB+NDEF+NTORS
C   INPUT DATA FOR STRETCHINGS
      IF(NST.EQ.0) GO TO 301
      READ(5,12) (KR(I),LR(I),I=1,NST)
C   INPUT DATA FOR BENDINGS
  301 IF(NBND.EQ.0) GO TO 302
      READ(5,13) (KA(I),LA(I),MA(I),I=1,NBND)
C   INPUT DAT FOR LINEAR BENDINGS
  302 IF(NLIB.EQ.0) GO TO 303
      NLHF=NLIB/2
      READ(5,13) (KB(I),LB(I),MB(I),I=1,NLHF)
C   INPUT DATA FOR DEFORMATIONS
  303 IF(NDEF.EQ.0) GO TO 304
      READ(5,14) (KD(I),LD(I),MD(I),ND(I),I=1,NDEF)
C   INPUT DATA FOR TORSIONS
  304 IF(NTORS.EQ.0) GO TO 305
      READ(5,14) (KT(I),LT(I),MT(I),NT(I),I=1,NTORS)
C
C   INPUT UNITARY MATRIX ELEMENTS
  305 CONTINUE
      IF(ISYM.EQ.0) GO TO 205
      NSYM=0
  500 READ(5,2) II,JJ,KK,LL
      IJKL=II*JJ*KK*LL
      IF(IJKL.EQ.0) GO TO 501
      NSYM=NSYM+1
      IIS(NSYM)=II
      JJS(NSYM)=JJ
      KKS(NSYM)=KK
      LLS(NSYM)=LL
      GO TO 500
  501 CONTINUE
C
C   READ IN ISOTOPIC ATOMIC MASS
  205 CONTINUE
      IF(ISOTOP.LE.1) GO TO 206
      if(.not.oldinp) goto 206
      DO 105 ISO=2,ISOTOP
      READ(5,15) (WISO(ISO,I),I=1,NATOM)
  105 CONTINUE
C
C   VIBRATIONAL ANALYSIS FOR ISOTOPIC MOLECULES
  206 CONTINUE
      ISOMOL=.FALSE.
      NISO=0
  400 NISO=NISO+1
      IF(NISO.EQ.1) GO TO 401
      IF(NISO.GT.ISOTOP) GO TO 220
      ISOMOL=.TRUE.
      WRITE(6,5)
      DO 106 I=1,NATOM
      X(I)=COORD(1,I)
      Y(I)=COORD(2,I)
      Z(I)=COORD(3,I)
      W(I)=WISO(NISO,I)
      WRITE(6,7) I,X(I),Y(I),Z(I),W(I)
  106 CONTINUE
      WRITE(6,8)
      DO 107 I=1,NATOM
      X(I)=X(I)*A0
      Y(I)=Y(I)*A0
      Z(I)=Z(I)*A0
      WRITE(6,7) I,X(I),Y(I),Z(I),W(I)
  107 CONTINUE
C
      CALL MOMENT
C
C   THE NORMAL COORDINATE ANALYSIS
  401 CONTINUE
      IF(IFXGF.EQ.2) GO TO 402
C
C   FX MATRIX METHOD
      N3N2=N3N*N3N
      N3TRI=IOFF(N3N+1)
      IC1=1
      IC2=IC1+N3N2
      IC3=IC2+N3N2
      IC4=IC3+N3N2
      IC5=IC4+N3N2
      IC6=IC5+N3N
      IC7=IC6+N3N
      IC8=IC7+N3N
      ICMAX=IC8+N3TRI
C................FX      ELX     ELXM    EE      VALU    FV1.....
      CALL VIBFX(CC(IC1),CC(IC2),CC(IC3),CC(IC4),CC(IC5),CC(IC6),
C................FV2     FXM...............
     1           CC(IC7),CC(IC8),N3N,N3TRI)
C
C   EVALUATE THERMODYNAMIC FUNCTIONS
      CALL THERMO
      GO TO 400
C
C   GF MATRIX METHOD
  402 CONTINUE
      N3N2=N3N*N3N
      N3R=N3N+NRED
      N3R2=N3R*N3R
      NATOM2=NATOM*NATOM
      N3TRI=IOFF(N3N+1)
      IC1=1
      IC2=IC1+N3N2
      IC3=IC2+N3R2
      IC4=IC3+N3N2
      IC5=IC4+N3TRI
      IC6=IC5+N3R2
      IC7=IC6+NATOM2
      IC8=IC7+NATOM2
      ICMAX=IC8+NATOM2
C................B       EX      EY      EZ.......................
      CALL VIBGF(CC(IC5),CC(IC6),CC(IC7),CC(IC8),N3R,NATOM,ISOMOL)
C
C   FORM SYMMETRY INTERNAL COORDINATE MATRIX IF NECESSARY
      IF(ISYM.EQ.0) GO TO 403
      IC1=1
      IC2=IC1+N3N2
      IC3=IC2+N3R2
      IC4=IC3+N3N2
      IC5=IC4+N3TRI
      IC6=IC5+N3R2
      IC7=IC6+N3R2
      ICMAX=IC7+N3R2
C................U       B       BS      EE.................
      CALL SYMCO(CC(IC2),CC(IC5),CC(IC6),CC(IC7),N3R,ISOMOL)
C
C   FORM G AND F MATRICES
  403 CONTINUE
      IC1=1
      IC2=IC1+N3N2
      IC3=IC2+N3R2
      IC4=IC3+N3N2
      IC5=IC4+N3TRI
      IC6=IC5+N3R2
      IC7=IC6+N3N2
      IC8=IC7+N3N2
      IC9=IC8+NATOM
      IA9=IC9+IC9-1
      ICMAX=IC9+NATOM
C................FX      F       G       B       BIN     EE......
      CALL GFMAT(CC(IC1),CC(IC3),CC(IC4),CC(IC5),CC(IC6),CC(IC7),
C................PIVOT   INDEX........................
     1           CC(IC8),IA(IC9),N3N,N3R,N3TRI,ISOMOL)
C
C   NORMAL COORDINATE ANALYSIS
      NVIB2=NVIB*NVIB
      IC1=1
      IC2=IC1+N3N2
      IC3=IC2+N3R2
      IC4=IC3+N3N2
      IC5=IC4+N3TRI
      IC6=IC5+N3R2
      IC7=IC6+N3N2
      IC8=IC7+NVIB2
      IC9=IC8+NVIB2
      IC10=IC9+NVIB2
      IC11=IC10+NVIB2
      IC12=IC11+NVIB2
      IC13=IC12+NVIB2
      IC14=IC13+NVIB2
      IC15=IC14+NVIB2
      IC16=IC15+NVIB
      IC17=IC16+NVIB
      IC18=IC17+NVIB
      ICMAX=IC18+NVIB
C.................F       G       S       WW      H       AA......
      CALL NORMGF(CC(IC3),CC(IC4),CC(IC4),CC(IC7),CC(IC8),CC(IC9),
C.................Y        ELIN     ELL      POT      EKED     EE.......
     1            CC(IC10),CC(IC11),CC(IC12),CC(IC13),CC(IC13),CC(IC13),
C.................GIN      VALU     ROOT     FV1      FV2......
     2            CC(IC14),CC(IC15),CC(IC16),CC(IC17),CC(IC18),
C................................
     3            N3N,N3TRI,NVIB)
C
C   THE CALCULATION OF VIBRATIONAL AMPLITUDES
      WRITE(6,17)
      IC1=1
      IC2=IC1+N3N2
      IC3=IC2+N3R2
      IC4=IC3+N3N2
      IC5=IC4+N3TRI
      IC6=IC5+N3R2
      IC7=IC6+N3N2
      IC8=IC7+NVIB2
      IC9=IC8+NVIB2
      IC10=IC9+NVIB2
      IC11=IC10+NVIB2
      IC12=IC11+NVIB2
      IC13=IC12+NVIB2
      IC14=IC13+N3N2
      ICMAX=IC14+N3N2
      CALL AMPGF(CC(IC5),CC(IC11),CC(IC13),CC(IC14),N3R)
C
C   THE CALCULATION OF INFRARED INTENSITIES
      IF(IDIPOL.EQ.0) GO TO 404
      WRITE(6,18)
      CALL IRINTR(CC(IC6),CC(IC12))
C
C   THE CALCULATION OF RAMAN INTENSITIES
  404 IF(IPOLAR.EQ.0) GO TO 405
      WRITE(6,19)
      CALL RAMANR(CC(IC6),CC(IC12))
C
C   END OF AN ISOTOPIC MOLECULE
  405 CONTINUE
      GO TO 400
C
  210 WRITE(6,20)
C
  220 CONTINUE
      CALL RCLOSE(ITAP30,3)
      CALL TSTOP(6)
      stop
      RETURN
      END
