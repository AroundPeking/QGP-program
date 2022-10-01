        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 29 22:27:54 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZUOIK__genmod
          INTERFACE 
            SUBROUTINE ZUOIK(ZR,ZI,FNU,KODE,IKFLG,N,YR,YI,NUF,TOL,ELIM, &
     &ALIM)
              REAL(KIND=8) :: ZR
              REAL(KIND=8) :: ZI
              REAL(KIND=8) :: FNU
              INTEGER(KIND=4) :: KODE
              INTEGER(KIND=4) :: IKFLG
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: YR(1)
              REAL(KIND=8) :: YI(1)
              INTEGER(KIND=4) :: NUF
              REAL(KIND=8) :: TOL
              REAL(KIND=8) :: ELIM
              REAL(KIND=8) :: ALIM
            END SUBROUTINE ZUOIK
          END INTERFACE 
        END MODULE ZUOIK__genmod
