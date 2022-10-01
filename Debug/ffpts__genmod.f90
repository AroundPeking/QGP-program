        !COMPILER-GENERATED INTERFACE MODULE: Wed Jul 13 10:46:34 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FFPTS__genmod
          INTERFACE 
            SUBROUTINE FFPTS(A,B,F,EPS,T,IDP,PT)
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: F
              EXTERNAL F
              REAL(KIND=4) :: EPS
              REAL(KIND=8) :: T
              INTEGER(KIND=4) :: IDP
              REAL(KIND=8) :: PT
            END SUBROUTINE FFPTS
          END INTERFACE 
        END MODULE FFPTS__genmod
