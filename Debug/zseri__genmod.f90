        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 29 22:27:54 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZSERI__genmod
          INTERFACE 
            SUBROUTINE ZSERI(ZR,ZI,FNU,KODE,N,YR,YI,NZ,TOL,ELIM,ALIM)
              REAL(KIND=8) :: ZR
              REAL(KIND=8) :: ZI
              REAL(KIND=8) :: FNU
              INTEGER(KIND=4) :: KODE
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: YR(1)
              REAL(KIND=8) :: YI(1)
              INTEGER(KIND=4) :: NZ
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: ELIM
              REAL(KIND=8) :: ALIM
            END SUBROUTINE ZSERI
          END INTERFACE 
        END MODULE ZSERI__genmod
