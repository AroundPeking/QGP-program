        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 29 22:27:54 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZBUNI__genmod
          INTERFACE 
            SUBROUTINE ZBUNI(ZR,ZI,FNU,KODE,N,YR,YI,NZ,NUI,NLAST,FNUL,  &
     &TOL,ELIM,ALIM)
              REAL(KIND=8) :: ZR
              REAL(KIND=8) :: ZI
              REAL(KIND=8) :: FNU
              INTEGER(KIND=4) :: KODE
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: YR(1)
              REAL(KIND=8) :: YI(1)
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: NUI
              INTEGER(KIND=4) :: NLAST
              REAL(KIND=8) :: FNUL
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: ELIM
              REAL(KIND=8) :: ALIM
            END SUBROUTINE ZBUNI
          END INTERFACE 
        END MODULE ZBUNI__genmod
