      SUBROUTINE GIVENS (NX,NROOTX,NJX,A,B,ROOT,VECT)                   
C     DOUBLE PRECISION VERSION BY MEL LEVY 8/72                         
C 62.3  GIVENS  -EIGENVALUES AND EIGENVECTORS BY THE GIVENS METHOD.     
C     BY FRANKLIN PROSSER, INDIANA UNIVERSITY.                          
C     SEPTEMBER, 1967                                                   
C     CALCULATES EIGENVALUES AND EIGENVECTORS OF REAL SYMMETRIC MATRIX  
C     STORED IN PACKED UPPER TRIANGULAR FORM.                           
C                                                                       
C     THANKS ARE DUE TO F. E. HARRIS (STANFORD UNIVERSITY) AND H. H.    
C     MICHELS (UNITED AIRCRAFT RESEARCH LABORATORIES) FOR EXCELLENT     
C     WORK ON NUMERICAL DIFFICULTIES WITH EARLIER VERSIONS OF THIS      
C     PROGRAM.                                                          
C                                                                       
C     THE PARAMETERS FOR THE ROUTINE ARE...                             
C         NX     ORDER OF MATRIX                                        
C         NROOTX NUMBER OF ROOTS WANTED.  THE NROOTX SMALLEST (MOST     
C                 NEGATIVE) ROOTS WILL BE CALCULATED.  IF NO VECTORS    
C                 ARE WANTED, MAKE THIS NUMBER NEGATIVE.                
C         NJX    ROW DIMENSION OF VECT ARRAY.  SEE -VECT- BELOW.        
C                 NJX MUST BE NOT LESS THAN NX.                         
C         A      MATRIX STORED BY COLUMNS IN PACKED UPPER TRIANGULAR    
C                FORM, I.E. OCCUPYING NX*(NX+1)/2 CONSECUTIVE           
C                LOCATIONS.                                             
C         B      SCRATCH ARRAY USED BY GIVENS.  MUST BE AT LEAST        
C                 NX*5 CELLS.                                           
C         ROOT   ARRAY TO HOLD THE EIGENVALUES.  MUST BE AT LEAST       
C                NROOTX CELLS LONG.  THE NROOTX SMALLEST ROOTS ARE      
C                 ORDERED LARGEST FIRST IN THIS ARRAY.                  
C         VECT   EIGENVECTOR ARRAY.  EACH COLUMN WILL HOLD AN           
C                 EIGENVECTOR FOR THE CORRESPONDING ROOT.  MUST BE      
C                 DIMENSIONED WITH -NJX- ROWS AND AT LEAST -NROOTX-     
C                 COLUMNS, UNLESS NO VECTORS                            
C                 ARE REQUESTED (NEGATIVE NROOTX).  IN THIS LATTER      
C                 CASE, THE ARGUMENT VECT IS JUST A DUMMY, AND THE      
C                 STORAGE IS NOT USED.                                  
C                 THE EIGENVECTORS ARE NORMALIZED TO UNITY.             
C                                                                       
C     THE ARRAYS A AND B ARE DESTROYED BY THE COMPUTATION.  THE RESULTS 
C     APPEAR IN ROOT AND VECT.                                          
C     FOR PROPER FUNCTIONING OF THIS ROUTINE, THE RESULT OF A FLOATING  
C     POINT UNDERFLOW SHOULD BE A ZERO.                                 
C     TO CONVERT THIS ROUTINE TO DOUBLE PRECISION (E.G. ON IBM 360      
C     MACHINES), BE SURE THAT ALL REAL VARIABLES AND FUNCTION           
C     REFERENCES ARE PROPERLY MADE DOUBLE PRECISION.                    
C     THE VALUE OF -ETA- (SEE BELOW) SHOULD ALSO BE CHANGED, TO REFLECT 
C     THE INCREASED PRECISION.                                          
C                                                                       
C     THE ORIGINAL REFERENCE TO THE GIVENS TECHNIQUE IS IN OAK RIDGE    
C     REPORT NUMBER ORNL 1574 (PHYSICS), BY WALLACE GIVENS.             
C     THE METHOD AS PRESENTED IN THIS PROGRAM CONSISTS OF FOUR STEPS,   
C     ALL MODIFICATIONS OF THE ORIGINAL METHOD...                       
C     FIRST, THE INPUT MATRIX IS REDUCED TO TRIDIAGONAL FORM BY THE     
C     HOUSEHOLDER TECHNIQUE (J. H. WILKINSON, COMP. J. 3, 23 (1960)).   
C     THE ROOTS ARE THEN LOCATED BY THE STURM SEQUENCE METHOD (J. M.    
C     ORTEGA (SEE REFERENCE BELOW).  THE VECTORS OF THE TRIDIAGONAL     
C     FORM ARE THEN EVALUATED (J. H. WILKINSON, COMP. J. 1, 90 (1958)), 
C     AND LAST THE TRIDIAGONAL VECTORS ARE ROTATED TO VECTORS OF THE    
C     ORIGINAL ARRAY (FIRST REFERENCE).                                 
C     VECTORS FOR DEGENERATE (OR NEAR-DEGENERATE) ROOTS ARE FORCED      
C     TO BE ORTHOGONAL, USING A METHOD SUGGESTED BY B. GARBOW, ARGONNE  
C     NATIONAL LABS (PRIVATE COMMUNICATION, 1964).  THE GRAM-SCHMIDT    
C     PROCESS IS USED FOR THE ORTHOGONALIZATION.                        
C                                                                       
C     AN EXCELLENT PRESENTATION OF THE GIVENS TECHNIQUE IS FOUND IN     
C     J. M. ORTEGA-S ARTICLE IN -MATHEMATICS FOR DIGITAL COMPUTERS,-    
C     VOLUME 2, ED. BY RALSTON AND WILF, WILEY (1967), PAGE 94.         
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DOUBLE PRECISION B(NX,5),A(1),ROOT(NROOTX),VECT(NJX,NROOTX)       
      DOUBLE PRECISION ETA,THETA,DEL1,DELTA,SMALL,DELBIG,THETA1,TOLER   
      DOUBLE PRECISION RPOWER,RPOW1,RAND1,FACTOR,ANORM,ALIMIT,SUM,TEMP  
      DOUBLE PRECISION AK,ROOTL,ROOTX,TRIAL,F0,SAVE,AROOT               
      DOUBLE PRECISION ELIM1,ELIM2                                      
C                                                                       
C ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C **  USERS PLEASE NOTE...                                              
C **  THE FOLLOWING TWO PARAMETERS, ETA AND THETA, SHOULD BE ADJUSTED   
C **  BY THE USER FOR HIS PARTICULAR MACHINE.                           
C **  ETA IS AN INDICATION OF THE PRECISION OF THE FLOATING POINT       
C **  REPRESENTATION ON THE COMPUTER BEING USED (ROUGHLY 10**(-M),      
C **  WHERE M IS THE NUMBER OF DECIMALS OF PRECISION ).                 
C **  THETA IS AN INDICATION OF THE RANGE OF NUMBERS THAT CAN BE        
C **  EXPRESSED IN THE FLOATING POINT REPRESENTATION (ROUGHLY THE       
C **  LARGEST NUMBER).                                                  
C **  SOME RECOMMENDED VALUES FOLLOW.                                   
C **  FOR IBM 7094, UNIVAC 1108, ETC. (27-BIT BINARY FRACTION, 8-BIT    
C **  BINARY EXPONENT), ETA=1.E-8, THETA=1.E37.                         
C **  FOR CONTROL DATA 3600 (36-BIT BINARY FRACTION, 11-BIT BINARY      
C **  EXPONENT), ETA=1.E-11, THETA=1.E307.                              
C **  FOR CONTROL DATA 6600 (48-BIT BINARY FRACTION, 11-BIT BINARY      
C **  EXPONENT), ETA=1.E-14, THETA=1.E307.                              
C **  FOR IBM 360/50 AND 360/65 DOUBLE PRECISION (56-BIT HEXADECIMAL    
C **  FRACTION, 7-BIT HEXADECIMAL EXPONENT), ETA=1.E-16, THETA=1.E75.   
C **                                                                    
      THETA = 1.D75                                                     
      ETA = 1.D-16                                                      
C ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C                                                                       
      DEL1=ETA/100.0D0                                                  
      DELTA=ETA**2*100.0D0                                              
      SMALL=ETA**2/100.0D0                                              
      DELBIG=THETA*DELTA/1000.0D0                                       
      THETA1=1000.0D0/THETA                                             
C     TOLER  IS A FACTOR USED TO DETERMINE IF TWO ROOTS ARE CLOSE       
C     ENOUGH TO BE CONSIDERED DEGENERATE FOR PURPOSES OF ORTHOGONALI-   
C     ZING THEIR VECTORS.  FOR THE MATRIX NORMED TO UNITY, IF THE       
C     DIFFERENCE BETWEEN TWO ROOTS IS LESS THAN TOLER, THEN             
C     ORTHOGONALIZATION WILL OCCUR.                                     
      TOLER = ETA*100.                                                  
C                                                                       
C     INITIAL VALUE FOR PSEUDORANDOM NUMBER GENERATOR... (2**23)-3      
      RPOWER=8388608.0D0                                                
      RPOW1 = RPOWER/2.                                                 
      RAND1=RPOWER-3.0D0                                                
C                                                                       
      N = NX                                                            
      NROOT = IABS(NROOTX)                                              
      IF (NROOT.EQ.0) GO TO 1001                                        
      IF (N-1) 1001,1003,105                                            
1003  ROOT(1) = A(1)                                                    
      IF(NROOTX .GT. 0)VECT(1,1)=1.0D0                                  
      GO TO 1001                                                        
105   CONTINUE                                                          
C     NSIZE    NUMBER OF ELEMENTS IN THE PACKED ARRAY                   
      NSIZE = (N*(N+1))/2                                               
      NM1 = N-1                                                         
      NM2 = N-2                                                         
C                                                                       
C     SCALE MATRIX TO EUCLIDEAN NORM OF 1.  SCALE FACTOR IS ANORM.      
      FACTOR=0.0D-0                                                     
      DO 70 I=1,NSIZE                                                   
 70   FACTOR=DMAX1(FACTOR,DABS(A(I)))                                   
      IF(FACTOR .NE. 0.0D0)GO TO 72                                     
C     NULL MATRIX.  FIX UP ROOTS AND VECTORS, THEN EXIT.                
      DO 78 I=1,NROOT                                                   
           IF (NROOTX.LT.0) GO TO 78                                    
           DO 77 J=1,N                                                  
 77        VECT(J,I)=0.0D0                                              
           VECT(I,I)=1.0D0                                              
 78   ROOT(I)=0.0D0                                                     
      GO TO 1001                                                        
C                                                                       
 72   ANORM=0.0D0                                                       
      J = 1                                                             
      K = 1                                                             
      DO 80 I=1,NSIZE                                                   
           IF (I.NE.J) GO TO 81                                         
           ANORM=ANORM+(A(I)/FACTOR)**2/2.0D0                           
           K = K+1                                                      
           J = J+K                                                      
           GO TO 80                                                     
81         ANORM = ANORM + (A(I)/FACTOR)**2                             
80    CONTINUE                                                          
      ANORM=DSQRT(ANORM*2.0D0)*FACTOR                                   
      DO 91 I=1,NSIZE                                                   
91    A(I) = A(I)/ANORM                                                 
      ALIMIT=1.0D-0                                                     
C                                                                       
C     TRIDIA SECTION.                                                   
C     TRIDIAGONALIZATION OF SYMMETRIC MATRIX                            
      ID = 0                                                            
      IA = 1                                                            
      IF (NM2.EQ.0) GO TO 201                                           
      DO 200  J=1,NM2                                                   
C     J       COUNTS ROW  OF A-MATRIX TO BE DIAGONALIZED                
C     IA      START OF NON-CODIAGONAL ELEMENTS IN THE ROW               
C     ID      INDEX OF CODIAGONAL ELEMENT ON ROW BEING CODIAGONALIZED.  
           IA = IA+J+2                                                  
           ID = ID + J + 1                                              
           JP2 = J+2                                                    
C     SUM SQUARES OF NON-CODIAGONAL ELEMENTS IN ROW J                   
           II = IA                                                      
           SUM=0.0D-0                                                   
           DO 100 I=JP2,N                                               
                SUM = SUM + A(II)**2                                    
100        II = II + I                                                  
           TEMP = A(ID)                                                 
           IF (SUM.GT.SMALL) GO TO 110                                  
C     NO TRANSFORMATION NECESSARY IF ALL THE NON-CODIAGONAL             
C     ELEMENTS ARE TINY.                                                
           B(J,1) = TEMP                                                
           A(ID)=0.0D0                                                  
           GO TO 200                                                    
C     NOW COMPLETE THE SUM OF OFF-DIAGONAL SQUARES                      
 110       SUM=DSQRT(SUM+TEMP**2)                                       
C     NEW CODIAGONAL ELEMENT                                            
           B(J,1)=-DSIGN(SUM,TEMP)                                      
C     FIRST NON-ZERO ELEMENT OF THIS W-VECTOR                           
           B(J+1,2)=DSQRT((1.0D0+DABS(TEMP)/SUM)/2.0D0)                 
C     FORM REST OF THE W-VECTOR ELEMENTS                                
           TEMP=DSIGN(0.5D0/(B(J+1,2)*SUM),TEMP)                        
           II = IA                                                      
           DO 130 I=JP2,N                                               
                B(I,2) = A(II)*TEMP                                     
130        II = II + I                                                  
C     FORM P-VECTOR AND SCALAR.  P-VECTOR = A-MATRIX*W-VECTOR.          
C     SCALAR = W-VECTOR*P-VECTOR.                                       
           AK=0.0D0                                                     
C     IC      LOCATION OF NEXT DIAGONAL ELEMENT                         
           IC = ID + 1                                                  
           J1 = J + 1                                                   
           DO 190  I=J1,N                                               
                JJ = IC                                                 
                TEMP=0.0D0                                              
                DO 180  II=J1,N                                         
C     I       RUNS OVER THE NON-ZERO P-ELEMENTS                         
C     II      RUNS OVER ELEMENTS OF W-VECTOR                            
                     TEMP = TEMP + B(II,2)*A(JJ)                        
C     CHANGE INCREMENTING MODE AT THE DIAGONAL ELEMENTS.                
                     IF (II.LT.I) GO TO 210                             
                     JJ = JJ + II                                       
                     GO TO 180                                          
210                  JJ = JJ + 1                                        
180             CONTINUE                                                
C     BUILD UP THE K-SCALAR (AK)                                        
                AK = AK + TEMP*B(I,2)                                   
                B(I,1) = TEMP                                           
C     MOVE IC TO TOP OF NEXT A-MATRIX -ROW-                             
190        IC = IC + I                                                  
C     FORM THE Q-VECTOR                                                 
           DO 150  I=J1,N                                               
150        B(I,1) = B(I,1) - AK*B(I,2)                                  
C     TRANSFORM THE REST OF THE A-MATRIX                                
C     JJ      START-1 OF THE REST OF THE A-MATRIX                       
           JJ = ID                                                      
C     MOVE W-VECTOR INTO THE OLD A-MATRIX LOCATIONS TO SAVE SPACE       
C     I       RUNS OVER THE SIGNIFICANT ELEMENTS OF THE W-VECTOR        
           DO 160  I=J1,N                                               
                A(JJ) = B(I,2)                                          
                DO 170  II=J1,I                                         
                     JJ = JJ + 1                                        
 170            A(JJ)=A(JJ)-2.0D0*(B(I,1)*B(II,2)+B(I,2)*B(II,1))       
160        JJ = JJ + J                                                  
200   CONTINUE                                                          
C     MOVE LAST CODIAGONAL ELEMENT OUT INTO ITS PROPER PLACE            
201   CONTINUE                                                          
      B(NM1,1) = A(NSIZE-1)                                             
      A(NSIZE-1)=0.0D-0                                                 
C                                                                       
C     STURM SECTION.                                                    
C     STURM SEQUENCE ITERATION TO OBTAIN ROOTS OF TRIDIAGONAL FORM.     
C     MOVE DIAGONAL ELEMENTS INTO SECOND N ELEMENTS OF B-VECTOR.        
C     THIS IS A MORE CONVENIENT INDEXING POSITION.                      
C     ALSO, PUT SQUARE OF CODIAGONAL ELEMENTS IN THIRD N ELEMENTS.      
      JUMP=1                                                            
      DO 320 J=1,N                                                      
           B(J,2)=A(JUMP)                                               
           B(J,3) = B(J,1)**2                                           
320   JUMP = JUMP+J+1                                                   
      DO 310 I=1,NROOT                                                  
310   ROOT(I) = +ALIMIT                                                 
      ROOTL = -ALIMIT                                                   
C     ISOLATE THE ROOTS.  THE NROOT LOWEST ROOTS ARE FOUND, LOWEST FIRST
      DO 330 I=1,NROOT                                                  
C     FIND CURRENT -BEST- UPPER BOUND                                   
           ROOTX = +ALIMIT                                              
           DO 335 J=I,NROOT                                             
 335       ROOTX=DMIN1(ROOTX,ROOT(J))                                   
           ROOT(I) = ROOTX                                              
C     GET IMPROVED TRIAL ROOT                                           
 500       TRIAL=(ROOTL+ROOT(I))*0.5D0                                  
           IF (TRIAL.EQ.ROOTL.OR.TRIAL.EQ.ROOT(I)) GO TO 330            
C     FORM STURM SEQUENCE RATIOS, USING ORTEGA-S ALGORITHM (MODIFIED).  
C     NOMTCH IS THE NUMBER OF ROOTS LESS THAN THE TRIAL VALUE.          
           NOMTCH=N                                                     
           J=1                                                          
360        F0 = B(J,2) - TRIAL                                          
370        CONTINUE                                                     
           IF(DABS(F0) .LT. THETA1)GO TO 380                            
           IF(F0 .GE. 0.0D0)NOMTCH=NOMTCH-1                             
           J = J + 1                                                    
           IF (J.GT.N) GO TO 390                                        
C     SINCE MATRIX IS NORMED TO UNITY, MAGNITUDE OF B(J,3) IS LESS THAN 
C     ONE, AND SINCE F0 IS GREATER THAN THETA1, OVERFLOW IS NOT POSSIBLE
C     AT THE DIVISION STEP.                                             
           F0 = B(J,2) - TRIAL - B(J-1,3)/F0                            
           GO TO 370                                                    
380        J = J + 2                                                    
           NOMTCH = NOMTCH - 1                                          
           IF (J.LE.N) GO TO 360                                        
390        CONTINUE                                                     
C     FIX NEW BOUNDS ON ROOTS                                           
           IF (NOMTCH.GE.I) GO TO 540                                   
           ROOTL = TRIAL                                                
           GO TO 500                                                    
540        ROOT(I) = TRIAL                                              
           NOM = MIN0(NROOT,NOMTCH)                                     
           ROOT(NOM) = TRIAL                                            
           GO TO 500                                                    
330   CONTINUE                                                          
C     REVERSE THE ORDER OF THE EIGENVALUES, SINCE CUSTOM DICTATES       
C     -LARGEST FIRST-.  THIS SECTION MAY BE REMOVED IF DESIRED WITHOUT  
C     AFFECTING THE REMAINDER OF THE ROUTINE.                           
C     NRT = NROOT/2                                                     
C     DO 10 I=1,NRT                                                     
C     SAVE = ROOT(I)                                                    
C     NMIP1 = NROOT - I + 1                                             
CC    ROOT(I) = ROOT(NMIP1)                                             
C10   ROOT(NMIP1) = SAVE                                                
C                                                                       
C     TRIVEC SECTION.                                                   
C     EIGENVECTORS OF CODIAGONAL FORM                                   
C807  CONTINUE                                                          
C     QUIT NOW IF NO VECTORS WERE REQUESTED.                            
      IF (NROOTX.LT.0) GO TO 1002                                       
C     INITIALIZE VECTOR ARRAY.                                          
      DO 15 I=1,N                                                       
           DO 15 J=1,NROOT                                              
 15   VECT(I,J)=1.0D-0                                                  
      DO 700 I=1,NROOT                                                  
           AROOT = ROOT(I)                                              
C     ORTHOGONALIZE IF ROOTS ARE CLOSE.                                 
           IF (I.EQ.1) GO TO 710                                        
C     THE ABSOLUTE VALUE IN THE NEXT TEST IS TO ASSURE THAT THE TRIVEC  
C     SECTION IS INDEPENDENT OF THE ORDER OF THE EIGENVALUES.           
           IF(DABS(ROOT(I-1)-AROOT) .LT. TOLER)GO TO 720                
710        IA = -1                                                      
720        IA = IA + 1                                                  
           ELIM1 = A(1) - AROOT                                         
           ELIM2 = B(1,1)                                               
           JUMP = 1                                                     
           DO 750  J=1,NM1                                              
                JUMP = JUMP+J+1                                         
C     GET THE CORRECT PIVOT EQUATION FOR THIS STEP.                     
                IF(DABS(ELIM1) .LE. DABS(B(J,1)))GO TO 760              
C     FIRST (ELIM1) EQUATION IS THE PIVOT THIS TIME.  CASE 1.           
                B(J,2) = ELIM1                                          
                B(J,3) = ELIM2                                          
                B(J,4)=0.0D0                                            
                TEMP = B(J,1)/ELIM1                                     
                ELIM1 = A(JUMP) - AROOT - TEMP*ELIM2                    
                ELIM2 = B(J+1,1)                                        
                GO TO 755                                               
C     SECOND EQUATION IS THE PIVOT THIS TIME.  CASE 2.                  
760             B(J,2) = B(J,1)                                         
                B(J,3) = A(JUMP) - AROOT                                
                B(J,4) = B(J+1,1)                                       
                TEMP=1.0D0                                              
                IF(DABS(B(J,1)) .GT. THETA1)TEMP=ELIM1/B(J,1)           
                ELIM1 = ELIM2 - TEMP*B(J,3)                             
                ELIM2 = -TEMP*B(J+1,1)                                  
C     SAVE FACTOR FOR THE SECOND ITERATION.                             
755             B(J,5) = TEMP                                           
750        CONTINUE                                                     
           B(N,2) = ELIM1                                               
           B(N,3)=0.0D-0                                                
           B(N,4)=0.0D-0                                                
           B(NM1,4)=0.0D-0                                              
           ITER = 1                                                     
           IF (IA.NE.0) GO TO 801                                       
C     BACK SUBSTITUTE TO GET THIS VECTOR.                               
790        L = N + 1                                                    
           DO 780 J=1,N                                                 
                L = L - 1                                               
786             CONTINUE                                                
                ELIM1=VECT(L,I)-VECT(L+1,I)*B(L,3)-VECT(L+2,I)*B(L,4)   
C     IF OVERFLOW IS CONCEIVABLE, SCALE THE VECTOR DOWN.                
C     THIS APPROACH IS USED TO AVOID MACHINE-DEPENDENT AND SYSTEM-      
C     DEPENDENT CALLS TO OVERFLOW ROUTINES.                             
                IF(DABS(ELIM1) .GT. DELBIG)GO TO 782                    
                TEMP = B(L,2)                                           
                IF(DABS(B(L,2)) .LT. DELTA)TEMP=DELTA                   
                VECT(L,I) = ELIM1/TEMP                                  
                GO TO 780                                               
C     VECTOR IS TOO BIG.  SCALE IT DOWN.                                
782             DO 784 K=1,N                                            
784             VECT(K,I) = VECT(K,I)/DELBIG                            
                GO TO 786                                               
780        CONTINUE                                                     
           GO TO (820,800), ITER                                        
C     SECOND ITERATION.  (BOTH ITERATIONS FOR REPEATED-ROOT VECTORS).   
820        ITER = ITER + 1                                              
890        ELIM1 = VECT(1,I)                                            
           DO 830 J=1,NM1                                               
                IF (B(J,2).EQ.B(J,1)) GO TO 840                         
C     CASE ONE.                                                         
                VECT(J,I) = ELIM1                                       
                ELIM1 = VECT(J+1,I) - ELIM1*B(J,5)                      
                GO TO 830                                               
C     CASE TWO.                                                         
840             VECT(J,I) = VECT(J+1,I)                                 
                ELIM1 = ELIM1 - VECT(J+1,I)*TEMP                        
830        CONTINUE                                                     
           VECT(N,I) = ELIM1                                            
           GO TO 790                                                    
C     PRODUCE A RANDOM VECTOR                                           
801        CONTINUE                                                     
           DO 802 J=1,N                                                 
C     GENERATE PSEUDORANDOM NUMBERS WITH UNIFORM DISTRIBUTION IN (-1,1).
C     THIS RANDOM NUMBER SCHEME IS OF THE FORM...                       
C     RAND1 = AMOD((2**12+3)*RAND1,2**23)                               
C     IT HAS A PERIOD OF 2**21 NUMBERS.                                 
                RAND1=DMOD(4099.0D0*RAND1,RPOWER)                       
 802       VECT(J,I)=RAND1/RPOW1-1.0D0                                  
           GO TO 790                                                    
C                                                                       
C     ORTHOGONALIZE THIS REPEATED-ROOT VECTOR TO OTHERS WITH THIS ROOT. 
800        IF (IA.EQ.0) GO TO 885                                       
           DO 860 J1=1,IA                                               
                K = I - J1                                              
                TEMP=0.0D0                                              
                DO 870 J=1,N                                            
870             TEMP = TEMP + VECT(J,I)*VECT(J,K)                       
                DO 880 J=1,N                                            
880             VECT(J,I) = VECT(J,I) - TEMP*VECT(J,K)                  
860        CONTINUE                                                     
885        GO TO (890,900), ITER                                        
C     NORMALIZE THE VECTOR                                              
 900       ELIM1=0.0D0                                                  
           DO 904 J=1,N                                                 
 904       ELIM1=DMAX1(DABS(VECT(J,I)),ELIM1)                           
           TEMP=0.0D-0                                                  
           DO 910 J=1,N                                                 
                ELIM2=VECT(J,I)/ELIM1                                   
910        TEMP = TEMP + ELIM2**2                                       
           TEMP=1.0D0/(DSQRT(TEMP)*ELIM1)                               
           DO 920 J=1,N                                                 
                VECT(J,I) = VECT(J,I)*TEMP                              
                IF(DABS(VECT(J,I)) .LT. DEL1)VECT(J,I)=0.0D0            
920        CONTINUE                                                     
700   CONTINUE                                                          
C                                                                       
C     SIMVEC SECTION.                                                   
C     ROTATE CODIAGONAL VECTORS INTO VECTORS OF ORIGINAL ARRAY          
C     LOOP OVER ALL THE TRANSFORMATION VECTORS                          
      IF (NM2.EQ.0) GO TO 1002                                          
      JUMP = NSIZE - (N+1)                                              
      IM = NM1                                                          
      DO 950  I=1,NM2                                                   
           J1 = JUMP                                                    
C     MOVE A TRANSFORMATION VECTOR OUT INTO BETTER INDEXING POSITION.   
           DO 955  J=IM,N                                               
                B(J,2) = A(J1)                                          
955        J1 = J1 + J                                                  
C     MODIFY ALL REQUESTED VECTORS.                                     
           DO 960  K=1,NROOT                                            
                TEMP=0.0D0                                              
C     FORM SCALAR PRODUCT OF TRANSFORMATION VECTOR WITH EIGENVECTOR     
                DO 970  J=IM,N                                          
970             TEMP = TEMP + B(J,2)*VECT(J,K)                          
                TEMP = TEMP + TEMP                                      
                DO 980  J=IM,N                                          
980             VECT(J,K) = VECT(J,K) - TEMP*B(J,2)                     
960        CONTINUE                                                     
           JUMP = JUMP - IM                                             
950   IM = IM - 1                                                       
1002  CONTINUE                                                          
C     RESTORE ROOTS TO THEIR PROPER SIZE.                               
      DO 95 I=1,NROOT                                                   
95    ROOT(I) = ROOT(I)*ANORM                                           
1001  RETURN                                                            
      END                                                               
