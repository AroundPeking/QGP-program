        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 29 22:27:54 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZBINU__genmod
          INTERFACE 
            SUBROUTINE ZBINU(ZR,ZI,FNU,KODE,N,CYR,CYI,NZ,RL,FNUL,TOL,   &
     &ELIM,ALIM)
              REAL(KIND=8) :: ZR
              REAL(KIND=8) :: ZI
              REAL(KIND=8) :: FNU
              INTEGER(KIND=4) :: KODE
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: CYR(1)
              REAL(KIND=8) :: CYI(1)
              INTEGER(KIND=4) :: NZ
              REAL(KIND=8) :: RL
              REAL(KIND=8) :: FNUL
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: ELIM
              REAL(KIND=8) :: ALIM
            END SUBROUTINE ZBINU
          END INTERFACE 
        END MODULE ZBINU__genmod
