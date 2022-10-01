        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 29 22:27:54 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZKSCL__genmod
          INTERFACE 
            SUBROUTINE ZKSCL(ZRR,ZRI,FNU,N,YR,YI,NZ,RZR,RZI,ASCLE,TOL,  &
     &ELIM)
              REAL(KIND=8) :: ZRR
              REAL(KIND=8) :: ZRI
              REAL(KIND=8) :: FNU
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: YR(1)
              REAL(KIND=8) :: YI(1)
              INTEGER(KIND=4) :: NZ
              REAL(KIND=8) :: RZR
              REAL(KIND=8) :: RZI
              REAL(KIND=8) :: ASCLE
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: ELIM
            END SUBROUTINE ZKSCL
          END INTERFACE 
        END MODULE ZKSCL__genmod
