        !COMPILER-GENERATED INTERFACE MODULE: Wed Jul 13 10:46:34 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INTEGRAL__genmod
          INTERFACE 
            SUBROUTINE INTEGRAL(A,B,F,EPS,T,IDP)
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: F
              EXTERNAL F
              REAL(KIND=4) :: EPS
              REAL(KIND=8) :: T
              INTEGER(KIND=4) :: IDP
            END SUBROUTINE INTEGRAL
          END INTERFACE 
        END MODULE INTEGRAL__genmod
