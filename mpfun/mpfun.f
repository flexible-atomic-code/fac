CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C THIS SIMPLE PACKAGE IMPLEMENTS THE QUAD-PRECISION
C ARITHMETICS FOR *, /, +, AND -, ACCORDING TO THE 
C ALGORITHMS OF LINNAINMMA, S. ACM TOMS. VOL. 7, PP. 272.
C M. F. GU
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MULTIPLICATION OF TWO DOUBLE PRECISION TO YEILD A QUAD 
C PRECISION NUMBER. THIS IS THE BASIC BUILDING BLOCK OF 
C THE MULTIPLICATION ALGORITHM.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE EXACTMUL(A, C, R)
      DOUBLE PRECISION A, C, R(2)
      
      DOUBLE PRECISION A1, A2, C1, C2, C21, C22, T

      INTEGER NBITS
      DOUBLE PRECISION CONST

C     IEEE FLOATING POINT HAS 53 BITS ACCURACY.
      PARAMETER (NBITS = 53)
      PARAMETER (CONST = 2.0**(NBITS - NBITS/2))

      T = CONST*A
      A1 = A - T
      A1 = A1 + T
      A2 = A - A1
      
      T = CONST*C
      C1 = C - T
      C1 = C1 + T
      C2 = C - C1

      T = C2*CONST
      C21 = C2 - T
      C21 = C21 + T
      C22 = C2 - C21

      R(1) = A*C
      T = A1*C1
      R(2) = T - R(1)
      T = A1*C2
      R(2) = R(2) + T
      T = C1*A2
      R(2) = R(2) + T
      T = C21*A2
      R(2) = R(2) + T
      T = C22*A2
      R(2) = R(2) + T
      
      END

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MULTIPLICATION OF TWO QUAD PRECISION NUMBERS.            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MPMUL(A, C, R)                            
      DOUBLE PRECISION A(2), C(2), R(2)                    
      DOUBLE PRECISION Z(2), T, D                          
                                                           
      CALL EXACTMUL(A(1), C(1), Z)                         
      T = A(1)*C(2)                                        
      D = A(2)*C(1)                                        
      T = T + D                                            
      T = T + Z(2)                                         
                                                           
      R(1) = Z(1) + T                                      
      D = Z(1) - R(1)
      R(2) = D + T                                         
                                                           
      END                                                  
                                                           
                                                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C DIVISION OF TWO QUAD PRECISION NUMBERS.                  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MPDIV(A, C, R)                            
      DOUBLE PRECISION A(2), C(2), R(2)                    
      DOUBLE PRECISION Q(2), Z(2), T, D                    
                                                           
      Z(1) = A(1)/C(1)                                     
      CALL EXACTMUL(C(1), Z(1), Q)                         
      T = A(1) - Q(1)                                      
      T = T - Q(2)                                         
      T = T + A(2)                                         
      D = Z(1)*C(2)                                        
      T = T - D                                            
      Z(2) = T/C(1)                                        
                                                           
      R(1) = Z(1) + Z(2)                                   
      D = Z(1) - R(1)                                      
      R(2) = D + Z(2)                                      
                                                           
      END                                                  
                                                           
                                                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ADDITION OF TWO QUAD PRECISION NUMBERS.                  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MPADD(A, C, R)                            
      DOUBLE PRECISION A(2), C(2), R(2)                    
      DOUBLE PRECISION Z(2), Q, T, D                       
                                                           
      Z(1) = A(1) + C(1)                                   
      Q = A(1) - Z(1)                                      
      T = Q + Z(1)                                         
      D = A(1) - T                                         
      T = Q + C(1)                                         
      D = T + D                                            
      D = D + A(2)                                         
      Z(2) = D + C(2)                                      
                                                           
      R(1) = Z(1) + Z(2)                                   
      T = Z(1) - R(1)                                      
      R(2) = T + Z(2)                                      
                                                           
      END                                                  
                                                           
                                                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SUBTRACTION OF TWO QUAD PRECISION NUMBERS.               
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MPSUB(A, C, R)                            
      DOUBLE PRECISION A(2), C(2), R(2)                    
      DOUBLE PRECISION D(2)                                
                                                           
      D(1) = -C(1)                                         
      D(2) = -C(2)                                         
      CALL MPADD(A, D, R)                                  
                                                           
      END                                                  
                                                           
                                                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ASIGN 1 QUAD PRECISION NUMBER TO ANOTHER                 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MPEQ(S, D)                                
      DOUBLE PRECISION D(2), S(2)                          
                                                           
      D(1) = S(1)                                          
      D(2) = S(2)                                          
                                                           
      END                                                  
                                                           
                                                           
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C CONVERT DOUBLE TO QUAD PRECISION                         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MPDMC(H, L, R)
      DOUBLE PRECISION H, L, R(2)

      R(1) = H
      R(2) = L

      END
                    
                                       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C ABS VALUE OF A QUAD PRECISION NUMBER
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MPABS(A, R)
      DOUBLE PRECISION A(2), R(2)
      
      IF (A(1) .LT. 0) THEN
         R(1) = -A(1)
         R(2) = -A(2)
      ELSE 
         R(1) = A(1)
         R(2) = A(2)
      ENDIF

      END
