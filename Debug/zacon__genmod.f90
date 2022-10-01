        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 29 22:27:54 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZACON__genmod
          INTERFACE 
            SUBROUTINE ZACON(ZR,ZI,FNU,KODE,MR,N,YR,YI,NZ,RL,FNUL,TOL,  &
     &ELIM,ALIM)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: ZR
              REAL(KIND=8) :: ZI
              REAL(KIND=8) :: FNU
              INTEGER(KIND=4) :: KODE
              INTEGER(KIND=4) :: MR
              REAL(KIND=8) :: YR(N)
              REAL(KIND=8) :: YI(N)
              INTEGER(KIND=4) :: NZ
              REAL(KIND=8) :: RL
              REAL(KIND=8) :: FNUL
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: ELIM
              REAL(KIND=8) :: ALIM
            END SUBROUTINE ZACON
          END INTERFACE 
        END MODULE ZACON__genmod
