        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 29 22:07:03 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INTEGRAL_1D__genmod
          INTERFACE 
            SUBROUTINE INTEGRAL_1D(FNX,A,B,DEL,P1,P2,P3,PT,IDP1,IDP2,   &
     &IDP3,CENTR,SUMF)
              REAL(KIND=8) :: FNX
              EXTERNAL FNX
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: DEL
              REAL(KIND=8) :: P1
              REAL(KIND=8) :: P2
              REAL(KIND=8) :: P3
              REAL(KIND=8) :: PT
              INTEGER(KIND=4) :: IDP1
              INTEGER(KIND=4) :: IDP2
              INTEGER(KIND=4) :: IDP3
              REAL(KIND=8) :: CENTR
              REAL(KIND=8) :: SUMF
            END SUBROUTINE INTEGRAL_1D
          END INTERFACE 
        END MODULE INTEGRAL_1D__genmod
