!usage example:
!    program main
!    external f
!    DOUBLE PRECISION f,s,a,b
!    a=0.
!    b=5.
!    e=0.1
!    call integral(a,b,f,e,s)
!    write(*,*)s
!    WRITE(*,*)125./3
!    end
!    
!function f(x)
!double precision f,x
!f=x**2
!return
!end
SUBROUTINE integral(A,B,F,EPS,T,IDP)
DIMENSION Y(10)
DOUBLE PRECISION A,B,F,T,Y,H,P,S,Q
integer IDP
    H=B-A
    Y(1)=H*(F(A,IDP)+F(B,IDP))/2.0
    M=1
    N=1
10  P=0.0
    DO 20 I=0,N-1
20      P=P+F(A+(I+0.5)*H,IDP)
    P=(Y(1)+H*P)/2.0
    S=1.0
    DO 30 K=1,M
        S=4*S
        Q=(S*P-Y(K))/(S-1)
        Y(K)=P
        P=Q
30  CONTINUE
    IF ((ABS(Q-Y(M)).GE.EPS).AND.(M.LE.9))THEN
        M=M+1
        Y(M)=Q
        N=N+N
        H=H/2.0
        GOTO 10
    ENDIF
    T=Q

    RETURN
END

SUBROUTINE integral_S(A,B,F,EPS,T,IDP,pt)
DIMENSION Y(10)
DOUBLE PRECISION A,B,F,T,Y,H,P,S,Q,pt
integer IDP
    H=B-A
    Y(1)=H*(F(pt,A,IDP)+F(pt,B,IDP))/2.0
    M=1
    N=1
10  P=0.0
    DO 20 I=0,N-1
20      P=P+F(pt,A+(I+0.5)*H,IDP)
    P=(Y(1)+H*P)/2.0
    S=1.0
    DO 30 K=1,M
        S=4*S
        Q=(S*P-Y(K))/(S-1)
        Y(K)=P
        P=Q
30  CONTINUE
    IF ((ABS(Q-Y(M)).GE.EPS).AND.(M.LE.9))THEN
        M=M+1
        Y(M)=Q
        N=N+N
        H=H/2.0
        GOTO 10
    ENDIF
    T=Q

    RETURN
END   