        !COMPILER-GENERATED INTERFACE MODULE: Tue Mar 29 22:07:45 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VEGAS__genmod
          INTERFACE 
            SUBROUTINE VEGAS(REGION,NDIM,FXN,INIT,NCALL,ITMX,NPRN,TGRAL,&
     &SD,CHI2A,ID,Q,C)
              INTEGER(KIND=4) :: NDIM
              REAL(KIND=8) :: REGION(2*NDIM)
              REAL(KIND=8) :: FXN
              EXTERNAL FXN
              INTEGER(KIND=4) :: INIT
              INTEGER(KIND=4) :: NCALL
              INTEGER(KIND=4) :: ITMX
              INTEGER(KIND=4) :: NPRN
              REAL(KIND=8) :: TGRAL
              REAL(KIND=8) :: SD
              REAL(KIND=8) :: CHI2A
              INTEGER(KIND=4) :: ID
              REAL(KIND=8) :: Q
              REAL(KIND=8) :: C
            END SUBROUTINE VEGAS
          END INTERFACE 
        END MODULE VEGAS__genmod
