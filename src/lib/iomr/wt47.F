      SUBROUTINE WT47(IA,IADD)                                          
      IMPLICIT INTEGER (A-Z)                                            
      DIMENSION IA(1)                                                   
      COMMON/LOC47/LOCS(192)                                            
      COMMON/T47SEC/ITAP47                                              
C                                                                       
C   NO.  :  CONTENTS                                                    
C    TRIA. = TRIANGLE                                                   
C    SQUA. = SQUARE MATRIX                                              
C    RECT. = RECTANGULAR MATRIX                                         
C    1   : PARAMETERS(1-100),LOCATIONS(101-142),LENGTHS(151-192)        
C    2   : LORD                                                         
C    3   : EIGENVALUES                                                  
C    4   : FOCC                                                         
C    5   : NSORB                                                        
C    6   : ALPHA (TRIA.)                                                
C    7   : BETA  (TRIA.)                                                
C    8   : SO-AO EIGENVECTORS IN PITZER'S SCF (SQUA.)                   
C    9   : SO-MO EIGENVECTORS IN PITZER'S SCF (SQUA.)                   
C   10   : AO-MO EIGENVECTORS IN PITZER'S SCF (SQUA.)                   
C   11   : WMOT (TRIA.)                                                 
C   12   : WAO (TRIA.)                                                  
C   13   : CM                                                           
C   14   : DC (TRIA.)                                                   
C   15   : DO (TRIA.)                                                   
C   16   : DO2 (TRIA.)                                                  
C   17   : DROL                                                         
C   18   : SO-MO EIGENVECTORS IN EIGENVALUE ORDER (SQUA.)               
C   19   : AO-MO EIGENVECTORS IN EIGENVALUE ORDER (SQUA.)               
C   20   : X (MO CI LAGRANGIAN) (SQUA.)                                 
C   21   : XAO (SYMMETRIZED AO CI LAGRANGIAN) (TRIA.)                   
C   22   : QAO (AO ONE PARTICLE DENSITY MATRIX) (TRIA.)                 
C   23   : QMO (MO ONE PARTICLE DENSITY MATRIX) (TRIA.)                 
C   24   : ZVECT (TRIA.)                                                
C   25   : U (FIRST-ORDER CPHF COEFFICIENTS)                            
C   26   : FORCE MATRIX                                                 
C   27   : IJORD (TRIA.)                                                
C   28   : HIJ (SQUARE)                                                 
C   29   : MO DERIVATIVE CI LAGRANGIAN (SQUA.)                          
C   30   : ELECTRONIC CONTRIBUTION TO CI GRADIENT (3*NAT)               
C                                                                       
      IADDR=LOCS(IADD+100)                                              
      NLEN=LOCS(IADD+150)                                               
      CALL WWRITW(ITAP47,IA,NLEN,IADDR,JUNK)                            
      RETURN                                                            
C                                                                       
      ENTRY RT47(IA,IADD)                                               
C                                                                       
CC    WRITE(*,*) ' IADD  = ',IADD                                       
      IADDR=LOCS(IADD+100)                                              
      NLEN=LOCS(IADD+150)                                               
CC    WRITE(*,*) ' IADDR= ',IADDR                                       
CC    WRITE(*,*) ' NLEN = ',NLEN                                        
      CALL WREADW(ITAP47,IA,NLEN,IADDR,JUNK)                            
      RETURN                                                            
      END                                                               
