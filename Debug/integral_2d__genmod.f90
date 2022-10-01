        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 29 22:07:03 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INTEGRAL_2D__genmod
          INTERFACE 
            SUBROUTINE INTEGRAL_2D(FNX,PT,IDP1,IDP2,IDP3,CENTR,SUMF)
              REAL(KIND=8) :: FNX
              EXTERNAL FNX
              REAL(KIND=8) :: PT
              INTEGER(KIND=4) :: IDP1
              INTEGER(KIND=4) :: IDP2
              INTEGER(KIND=4) :: IDP3
              REAL(KIND=8) :: CENTR
              REAL(KIND=8) :: SUMF
            END SUBROUTINE INTEGRAL_2D
          END INTERFACE 
        END MODULE INTEGRAL_2D__genmod
