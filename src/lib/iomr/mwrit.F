      SUBROUTINE MWRIT(IA,IADD)                                         
      IMPLICIT INTEGER (A-Z)                                            
      DIMENSION IA(1)                                                   
      COMMON/LOCAT/LOCS(1024)                                           
      COMMON/MFSEC/MFILE,NSECT                                          
C                                                                       
C  NO.  :  CONTENTS                                                     
C                                                                       
C   TRIA. = LOWER TRIANGULAR MATRIX                                     
C   SQUA. = SQUARE MATRIX                                               
C   RECT. = RECTANGULAR MATRIX                                          
C   1   : LOCATION                                                      
C   2   : PARAMETERS                                                    
C   3   : ENERGIES AND REAL CONSTANTS                                   
C   4   : MO ORDERING IN SCF                                            
C   5   : MO ORDERING IN DRT                                            
C   6   : MO CODE FOR CI                                                
C   7   : MO CODE FOR MCSCF                                             
C   8   : NUCLEAR CHARGE AND GEOMETRY                                   
C   9   : AO-SO TRANSFORMATION MATRIX (RECT.)                           
C  10   : EIGENVALUES IN PITZER'S SCF (SQUA.)                           
C  11   : OCCUPATION IN PITZER'S SCF                                    
C  12   : SO-MO EIGENVECTORS IN PITZER'S SCF (SQUA.)                    
C  13   : AO-MO EIGENVECTORS IN PITZER'S SCF (RECT.)                    
C  14   : SO OVERLAP INTEGRALS IN PITZER'S SCF (TRIA.)                  
C  15   : ONE ELECTRON SO INTEGRALS IN PITZER'S SCF (TRIA.)             
C  16   : ONE ELECTRON MO INTEGRALS IN PITZER'S SCF (TRIA.)             
C  17   : EIGENVALUES IN SORTED SCF                                     
C  18   : OCCUPATION IN SORTED SCF                                      
C  19   : SO-MO EIGENVECTORS IN SORTED SCF (SQUA.)                      
C  20   : AO-MO EIGENVECTORS IN SORTED SCF (RECT.)                      
C  21   : SO OVERLAP INTEGRALS IN SORTED SCF (TRIA.)                    
C  22   : ONE ELECTRON SO INTEGRALS IN SORTED SCF (TRIA.)               
C  23   : ONE ELECTRON MO INTEGRLAS IN SORTED SCF (TRIA.)               
C  24   : LAGRANGIAN MATRIX FOR SCF IN AO BASIS (TRIA.)                 
C  25   : LAGRANGIAN MATRIX FOR SCF IN MO BASIS (TRIA.)                 
C  26   : K MATRIX FOR HIGH SPIN OPEN-SHELL SCF (TRIA.)                 
C  27   : FIRST ZETA MATRIX  FOR GENERALIZED OPEN-SHELL SCF (TRIA.)     
C  28   : SECOND ZETA MATRIX FOR GENERALIZED OPEN-SHELL SCF (TRIA.)     
C  29   : THIRD ZETA MATRIX  FOR GENERALIZED OPEN-SHELL SCF (TRIA.)     
C  30   : FOURTH ZETA MATRIX FOR GENERALIZED OPEN-SHELL SCF (TRIA.)     
C  31   : FIFTH ZETA MATRIX  FOR GENERALIZED OPEN-SHELL SCF (TRIA.)     
C  32   : LAGRANGIAN MATRIX FOR CI IN AO BASIS (TRIA.)                  
C  33   : LAGRANGIAN MATRIX FOR CI IN MO BASIS (SQUA.)                  
C  34   : ONE PDM IN AO BASIS FOR CI (TRIA.)                            
C  35   : ONE PDM IN MO BASIS FOR CI (TRIA.)                            
C  36   : LAGRANGIAN MATRIX FOR MCSCF IN AO BASIS (TRIA.)               
C  37   : LAGRANGIAN MATRIX FOR MCSCF IN MO BASIS (SQUA.)               
C  38   : ONE PDM IN AO BASIS FOR MCSCF (TRIA.)                         
C  39   : ONE PDM IN MO BASIS FOR MCSCF (TRIA.)                         
C                                                                       
      IADDR=LOCS(IADD)                                                  
      NSECTH=NSECT/2                                                    
      NLEN=LOCS(IADD+NSECTH)                                            
      CALL RWRIT(MFILE,IA,NLEN,IADDR)                                   
      RETURN                                                            
C                                                                       
      end
