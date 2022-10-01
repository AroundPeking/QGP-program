        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 29 22:27:54 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ZBESK__genmod
          INTERFACE 
            SUBROUTINE ZBESK(ZR,ZI,FNU,KODE,N,CYR,CYI,NZ,IERR)
              REAL(KIND=8) :: ZR
              REAL(KIND=8) :: ZI
              REAL(KIND=8) :: FNU
              INTEGER(KIND=4) :: KODE
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: CYR(1)
              REAL(KIND=8) :: CYI(1)
              INTEGER(KIND=4) :: NZ
              INTEGER(KIND=4) :: IERR
            END SUBROUTINE ZBESK
          END INTERFACE 
        END MODULE ZBESK__genmod
