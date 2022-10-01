        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 29 22:27:54 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZUNIK__genmod
          INTERFACE 
            SUBROUTINE ZUNIK(ZRR,ZRI,FNU,IKFLG,IPMTR,TOL,INIT,PHIR,PHII,&
     &ZETA1R,ZETA1I,ZETA2R,ZETA2I,SUMR,SUMI,CWRKR,CWRKI)
              REAL(KIND=8) :: ZRR
              REAL(KIND=8) :: ZRI
              REAL(KIND=8) :: FNU
              INTEGER(KIND=4) :: IKFLG
              INTEGER(KIND=4) :: IPMTR
              REAL(KIND=8) :: TOL
              INTEGER(KIND=4) :: INIT
              REAL(KIND=8) :: PHIR
              REAL(KIND=8) :: PHII
              REAL(KIND=8) :: ZETA1R
              REAL(KIND=8) :: ZETA1I
              REAL(KIND=8) :: ZETA2R
              REAL(KIND=8) :: ZETA2I
              REAL(KIND=8) :: SUMR
              REAL(KIND=8) :: SUMI
              REAL(KIND=8) :: CWRKR(16)
              REAL(KIND=8) :: CWRKI(16)
            END SUBROUTINE ZUNIK
          END INTERFACE 
        END MODULE ZUNIK__genmod
