
!========================================
!    Note
!1.
!If you would like to change the cutoff(integral interval) the shower SS and TS(S_j) 
!to improve results,Please ctrl+F to find "the cutoff" to change them
!
!Suggested values:
!for Jpsi, we choose q=2.7-30 for SS integral, q=2.8-50 for Shower integral to reproduce.
!for D0, now choose q=3.2-30 for SS integral, q=2.8-50 for Shower integral to reproduce.
!for Ds, now choose q=3.2-30 for SS integral, q=3.0-50 for Shower integral to reproduce.
!
!2.
!for Jpsi bL=0.01 to charm, but bL=2.39 to charm and gluon for D meson
!
!3.
!before run this program, check the energy of system, i.e. /sqrt s_NN=200 GeV, 2.76 TeV, 5.02 TeV, 5.5 TeV
!which will lead to different parameters in fik.
!(parameters for 2.76 and 5.02 are given except for c/cbar, which are linearly interpolated by ln(fik) and sNN between 5.5 and 0.2)
!
!4.SS(2j) should be considered at LHC energy, which is negligible at RHIC energy.
!=====================================================
!Two free parameters to fit
!1.v_T=tanh(eta_T), we will use sinh(eta_T) and cosh(eta_T)
!2.gamac
!
!
    program charmed_meson_production
    implicit none
    !include 'Recom.h'
    
    ! functions to calculate different components for meson
	external meson_TT, meson_TS, meson_SS,int_SS, meson_SS2j
	real*8 meson_TT, meson_TS, meson_SS,int_SS, meson_SS2j
    
    real*8 Dc,Dg,x
    external S_ij, thermal,int_Sj,S_j
    double precision S_ij, thermal,int_Sj,S_j
    ! IDP1, parton 1 id
    ! The parton ID
    ! 1, 2, 3,    4,    5,   6,  7,     8, 9
    ! u, d, ubar, dbar, g,   s,  sbar,  c, cbar
    ! IDmeson, represent the particle
    ! particle 
    ! J/psi (c, cbar) (8, 9) --1
    ! D0  (c,ubar) (8, 3)    --2
    ! Ds   (c, sbar)  (8, 7) --3
	real*8 pt, q, dNptdpt(5)   !4 components for meson
	integer IDP1, IDP2, IDmeson
	integer i,j, N
    !-------------------------------------------
    !---------------------------------------------
! scan the parameter space
    REAL*8 DATAX(500), DATAY(500),plow,phigh
	INTEGER IPOINT
	REAL(8) CHISQ


!	open(unit=101, file='C:\Users\86177\Desktop\innovation\charmed meson program\charm\AuAu200Jpsi_data.dat', status='old')
!	!open(unit=102, file='p05_pb502_chisq.dat',status='unknown')
!
!	IPOINT = 1
!	
!200	READ(101, *, END=210)DATAX(IPOINT), DATAy(IPOINT)
!	PRINT *, IPOINT, DATAX(IPOINT), DATAY(IPOINT)*1.d-6
!	IPOINT = IPOINT + 1
!	GOTO 200
!210	CONTINUE
!
!	IPOINT = IPOINT-1
!    do i=1, IPOINT
!        pt = DATAX(i)
!        dNptdpt(1)  = meson_TT(pt, 8, 9, 1) 
!	    dNptdpt(2)  = meson_TS(pt, 8, 9, 1) 
!	    dNptdpt(3)  = meson_SS(pt, 8, 9, 1) 
!	    dNptdpt(4)  = meson_SS2j(pt, 8, 9, 1)
!	    dNptdpt(5)  = dNptdpt(1)+dNptdpt(2)+dNptdpt(3)+dNptdpt(4)
!
!	!CHISQ = CHISQ + (dNptdpt(5)/DATAY(I)-1.0)**2
!      ! write(*,*) "dNptdpt1"
!       ! write(*,*) dNptdpt(1)
!       ! write(*,*) "dNptdpt2"
!       ! write(*,*) dNptdpt(2)
!      !  write(*,*) "dNptdpt3"
!     !   write(*,*) dNptdpt(3)
!     !   write(*,*) "all:dNptdpt"
!        write(*,*) dNptdpt(5)
!    enddo
!    CHISQ = 0.0
    
    !test SS----------------------------
    !open(unit = 17, file = "SS_Jpsi_c020.dat", status = "unknown")
    !write(17,111)
    !testify thermal-------------------
    !open(unit = 14, file = "Tu_Jpsi_c020.dat", status = "unknown")
    !open(unit = 15, file = "Ts_Jpsi_c020.dat", status = "unknown")
    !open(unit = 16, file = "Tc_Jpsi_c020.dat", status = "unknown")
    !write(14, 98)
    !write(15, 98)
    !write(16, 98)
    !---------------------------------calculate Jpsi begin-------------------------------------------
    open(unit = 11, file = "dNptdpt_Jpsi_c020_55.dat", status = "unknown")
    write(11, 108)
    phigh=20.00
    plow=10.00
    N=20
    do i = 1, N
	    pt = plow+(i-1)*(phigh-plow)/N
        dNptdpt(1)  = meson_TT(pt, 8, 9, 1) 
	    dNptdpt(2)  = meson_TS(pt, 8, 9, 1) 
	    dNptdpt(3)  = meson_SS(pt, 8, 9, 1) 
	    dNptdpt(4)  = meson_SS2j(pt, 8, 9, 1)
	    dNptdpt(5)  = dNptdpt(1)+dNptdpt(2)+dNptdpt(3)+dNptdpt(4)
        
        write(11, 110) pt, dNptdpt(1:5)
    enddo
    !!---------------------------------calculate Jpsi end-------------------------------------------
    
    !---------------------------------calculate D0 begin---------------------------------------------
    !open(unit = 18, file = "dNptdpt_D0_c020_276.dat", status = "unknown")
    !write(18, 108)
    !phigh=20.00
    !plow=0.05
    !N=100
    !do i = 1, N
	   ! pt = plow+(i-1)*(phigh-plow)/N
    !    dNptdpt(1)  = meson_TT(pt, 8, 3, 2) 
	   ! dNptdpt(2)  = meson_TS(pt, 8, 3, 2) 
	   ! dNptdpt(3)  = meson_SS(pt, 8, 3, 2) 
	   ! dNptdpt(4)  = meson_SS2j(pt, 8, 3, 2)
	   ! dNptdpt(5)  = dNptdpt(1)+dNptdpt(2)+dNptdpt(3)+dNptdpt(4)
    !    
    !    write(18, 110) pt, dNptdpt(1:5)
    !    write(*, *) pt, dNptdpt(1:5)
    !enddo
    !!---------------------------------calculate D0 end---------------------------------------------
    
    !---------------------------------calculate Ds begin---------------------------------------------
    !open(unit = 22, file = "dNptdpt_Ds_c020_276.dat", status = "unknown")
    !write(22, 108)
    !phigh=20.00
    !plow=0.1
    !N=100
    !do i = 1, N
	   ! pt = plow+(i-1)*(phigh-plow)/N
    !    dNptdpt(1)  = meson_TT(pt, 8, 7, 3) 
	   ! dNptdpt(2)  = meson_TS(pt, 8, 7, 3) 
	   ! dNptdpt(3)  = meson_SS(pt, 8, 7, 3) 
	   ! dNptdpt(4)  = meson_SS2j(pt, 8, 7, 3)
	   ! dNptdpt(5)  = dNptdpt(1)+dNptdpt(2)+dNptdpt(3)+dNptdpt(4)
    !    
    !    write(22, 110) pt, dNptdpt(1:5)
    !enddo
    !---------------------------------calculate Ds end---------------------------------------------
    
    !write(14, 110) pt, Thermal(pt, 1)
    !write(15, 110) pt, Thermal(pt, 6)
    !write(16, 110) pt, Thermal(pt, 8)
    
        !test int_SS ,  integral by person
    !dNptdpt(2)=0
    !dNptdpt(5)=0
    !dNptdpt(1)=0
    !do j=1,1000
    !    if(pt.le.3.5)then
    !        q=4.-1.14286*pt +0.285714*pt**2+(j-1)*(30.0-pt)/1000!3.5+dexp(-pt)+(j-1)*(30.0-pt)/1000
    !    elseif(pt.gt.3.5)then
    !        q=pt+(j-1)*(30.0-pt)/1000
    !    endif
    !    dNptdpt(5) = int_SS(pt,q,1)
    !    dNptdpt(1)=dNptdpt(5)*(30.0-pt)/1000
    !    dNptdpt(2)=dNptdpt(2)+dNptdpt(1)
    !   ! write(*,*)pt, dNptdpt(5) 
    !enddo
    !   ! write(*,*)"all---------------------------"
    !    write(17,110)pt, dNptdpt(2)/pt*(1-dexp(-pt/2.0))
    !
    !enddo
    
    !open(unit = 25, file = "int_Sj_c020.dat", status = "unknown")
    !write(25,112)
    !do i=1,1000
    !    pt=2.3+0.048*i
    !    dNptdpt(1)=int_Sj(2.0d0,pt,8)
    !    write(25,110)2+0.048*i,dNptdpt(1)
    !enddo
    
    !open(unit = 26, file = "Sj_c020.dat", status = "unknown")
    !write(26,113)
    !do i=1,1000
    !    pt=0.05+i*6.5/1000
    !    dNptdpt(1)=S_j(pt,9)
    !    write(26,110)pt,dNptdpt(1)
    !enddo
    
    !having testified D and SPD right
    !open(unit = 12, file = "Dc_Jpsi_c020.dat", status = "unknown")
    !write(12, 100)
    !open(unit = 13, file = "Dg_Jpsi_c020.dat", status = "unknown")
    !write(13, 100)
    !open(unit = 55, file = "S_c_cbar_Jpsi_c020.dat", status = "unknown")
    !open(unit = 66, file = "S_g_ccbar_Jpsi_c020.dat", status = "unknown")
    !open(unit = 77, file = "S_c_l_Jpsi_c020.dat", status = "unknown")
    !open(unit = 88, file = "S_c_c_Jpsi_c020.dat", status = "unknown")
    !do i=1, 100
    !    x=0+0.01*i
    !    Dc=0.5/x*                                            &
    !    ( S_ij(x/2,8,8)*S_ij(x/(2-x),8,9)+S_ij(x/2,8,9)*S_ij(x/(2-x),8,8))
    !    Dg=0.5/x*                                                   &
    !    ( S_ij(x/2,5,8)*S_ij(x/(2-x),5,9)+S_ij(x/2,5,9)*S_ij(x/(2-x),5,8))
    !    write(12,110)x,Dc
    !    write(13,110)x,Dg
        !write(55,110)x,S_ij(x,8,9)
        !write(66,110)x,S_ij(x,5,8)
        !write(77,110)x,S_ij(x,8,1)
        !write(88,110)x,S_ij(x,8,8)
    !enddo
!99  format(10x, "%x", 12x, "S_c_cbar/S_g_ccbar/S_c_l/S_c_c")
!100 format(10x, "%x", 12x, "Dc/Dg")
!98  format(10x, "%x", 12x, "T")
!111  format(10x, "%pt", 12x, "SS(by person)")
!112 format(10x, "%pt", 12x, "int_Sj")
!113 format(10x, "%pt", 12x, "Sj")
108 format(10x, "%pt", 12x, "TT", 12x,  "TS", 12x,  "SS", 12x,  "SS2j", 12x, "Tot")
110	format((1x, F12.3), 5(1x, ES14.5E2)) 

    end!!!!!!!!!!!!main program end
    
    !=======================================================
	! meson spectrum from TT
	function meson_TT(pt, IDP1, IDP2, IDmeson)
	implicit none
	!include 'Recom.h'
	
	real*8 meson_TT
	real*8 pt
	integer IDP1, IDP2, IDmeson, i
    double precision :: pi=3.141592653589793D0
	real*8 gamc, gams, gaml
    real*8 :: C_M=(2*3)**2, T=0.175
    real*8 x, M_T, m_h, tau, A_T, mc, mu, ms, resu, q,tmp
    real*8 Ir, Ii, Kr, Ki, Kr2
    DOUBLE PRECISION I0r(1),I0i(1),K1r(1),K1i(1),K1r2(1)
    integer ierr, nz
    !200 GeV for J/psi vT=tanh(etaT)=0.3c,  gamc=0.26        eta=0.30952
    !200 GeV for D0    vT=tanh(etaT)=0.42c, gamc=0.26        eta=0.447692

    gaml=1.d0
    !gamc=0.26d0!200 GeV
    !gamc=1.2d0!test for 2.76 TeV
    gamc=0.26d0
    gams=0.8d0


!----------------------------------------------	
	if(IDmeson.eq.1)then ! for J/psi
    m_h=3.097
    mc=1.5
    M_T=Sqrt(m_h**2+pt**2)
    tau=25.34d0
    A_T=(45.61d0)**2*pi
    x=0.5d0!wavefunction**2 = delta(x-0.5) for J/psi
    
    !----------------free parameters for TT/TS (v_T=tanh eta_T)-----------------------------
    !Ir=sinh(0.100335)*pt/T                                                       !test for 2.76 TeV(vT=0.1)
    Ir=sinh(0.30952)*pt/T                                                  !for 200 GeV and 5.5 TeV(vT=0.3)
    Ii=0.d0
    !Kr=cosh(0.100335)*(Sqrt(mc**2+(x*pt)**2)+Sqrt(mc**2+((1-x)*pt)**2))/T        !test for 2.76 TeV(vT=0.1)
    Kr=cosh(0.30952)*(Sqrt(mc**2+(x*pt)**2)+Sqrt(mc**2+((1-x)*pt)**2))/T   !for 200 GeV(vT=0.3)
    Ki=0.d0
    !----------------free parameters for TT/TS (v_T=tanh eta_T)-----------------------------
    
    call ZBESI(Ir,Ii,0.d0,1,1,I0r,I0i,nz,ierr)
    call ZBESK(Kr,Ki,0.d0,1,1,K1r,K1i,nz,ierr)
	
	meson_TT = C_M*M_T*tau*A_T/(2*pi)**3*2*gamc*gamc*I0r(1)*K1r(1)!over
	
!-------------------------------------------	
    elseif(IDmeson.eq.2)then ! for D0
    m_h=1.8
    M_T=Sqrt(m_h**2+pt**2)
    tau=25.34d0
    A_T=(45.6d0)**2*pi
    mc=1.5
    mu=0.26

    resu=0.d0
    
		do i=1,1000
            x=i*(1.-0.)/1000
            tmp=x+(1.-0.)/1000
            Ir=sinh(0.447692)*pt/T
            Ii=0.d0
            Kr=cosh(0.447692)*(Sqrt(mu**2+(x*pt)**2)+Sqrt(mc**2+((1-x)*pt)**2))/T
            Kr2=cosh(0.447692)*(Sqrt(mu**2+(tmp*pt)**2)+Sqrt(mc**2+((1-tmp)*pt)**2))/T
            Ki=0.d0
            call ZBESI(Ir,Ii,0.d0,1,1,I0r,I0i,nz,ierr)
            call ZBESK(Kr,Ki,0.d0,1,1,K1r,K1i,nz,ierr)
            call ZBESK(Kr2,Ki,0.d0,1,1,K1r2,K1i,nz,ierr)
            resu=resu+5*( K1r(1)*(1-x)**4 + K1r2(1)*(1-(tmp))**4 )*(1.-0.)/1000/2
        enddo
	meson_TT =  C_M*M_T*tau*A_T/(2*pi)**3*2*gamc*gaml*I0r(1)*resu
	
	elseif(IDmeson.eq.3)then ! for Ds
	
	m_h=1.96
    M_T=Sqrt(m_h**2+pt**2)
    tau=25.34d0
    A_T=(45.6d0)**2*pi
    mc=1.5
    ms=0.46

    resu=0.d0
    
		do i=1,1000
            x=i*(1.-0.)/1000
            tmp=x+(1.-0.)/1000
            Ir=sinh(0.447692)*pt/T
            Ii=0.d0
            Kr=cosh(0.447692)*(Sqrt(ms**2+(x*pt)**2)+Sqrt(mc**2+((1-x)*pt)**2))/T
            Kr2=cosh(0.447692)*(Sqrt(ms**2+(tmp*pt)**2)+Sqrt(mc**2+((1-tmp)*pt)**2))/T
            Ki=0.d0
            call ZBESI(Ir,Ii,0.d0,1,1,I0r,I0i,nz,ierr)
            call ZBESK(Kr,Ki,0.d0,1,1,K1r,K1i,nz,ierr)
            call ZBESK(Kr2,Ki,0.d0,1,1,K1r2,K1i,nz,ierr)
            resu=resu+660*( K1r(1)*x**2 *(1-x)**9 + K1r2(1)*tmp**2 *(1-(tmp))**9 )*(1.-0.)/1000/2
        enddo
	meson_TT =  C_M*M_T*tau*A_T/(2*pi)**3*2*gamc*gams*I0r(1)*resu
	
	endif
	
	return

    end!meson_TT end
	
	function meson_TS(pt, IDP,IDP2, IDmeson)	
	implicit none
!	include 'Recom.h'
	
	real*8 meson_TS
	real*8 pt, C_M, gamc,gaml,gams,resu,tmp,x
	integer IDP, IDP2, IDmeson,g,i
    external Thermal, S_j
    double precision Thermal, S_j

    C_M=(2*3)**2
    g=6
    gamc=1.2!2.76 TeV
    !gamc=0.26!200 GeV
    gaml=1.0
    gams=0.8
    resu=0.d0
!------------------------------------------------------------	
	if(IDmeson .eq.1)then !for J/psi

	meson_TS = (1-dexp(-pt/2.0))*C_M/g/gamc*2*    &
    (S_j(pt/2.0,IDP)*Thermal(pt/2.0,IDP2)+S_j(pt/2.0,IDP2)*Thermal(pt/2.0,IDP))
	
!------------------------------------------------

    elseif(IDmeson.eq.2)then !for D0
        do i=1,998
            x=i*(1.-0.)/1000
            tmp=x+1./1000
            resu=resu+(( Thermal(pt*x,IDP2)*S_j(pt*(1-x),IDP)/gaml/x              &
            +S_j(pt*x,IDP2)*Thermal(pt*(1-x),IDP)/gamc/(1-x) )*(1-x)**4 +         &
                ( Thermal(pt*tmp,IDP2)*S_j(pt*(1-tmp),IDP)/gaml/tmp               &
            +S_j(pt*tmp,IDP2)*Thermal(pt*(1-tmp),IDP)/gamc/(1-tmp) )*(1-tmp)**4)*(1.-0.)/1000/2
        enddo
	
    meson_TS = (1-dexp(-pt/2.0))*C_M/g *5*resu
	
	
	elseif(IDmeson.eq.3)then ! for Ds
	
	do i=1,998
            x=i*(1.-0.)/1000
            tmp=x+1./1000
            resu=resu+(( Thermal(pt*x,IDP2)*S_j(pt*(1-x),IDP)/gams/x                    &
            +S_j(pt*x,IDP2)*Thermal(pt*(1-x),IDP)/gamc/(1-x) )*x**2 *(1-x)**9 +         &
                ( Thermal(pt*tmp,IDP2)*S_j(pt*(1-tmp),IDP)/gams/tmp                     &
            +S_j(pt*tmp,IDP2)*Thermal(pt*(1-tmp),IDP)/gamc/(1-tmp) )*tmp**2 *(1-tmp)**9)*(1.-0.)/1000/2
        enddo
	
    meson_TS = (1-dexp(-pt/2.0))*C_M/g *660*resu
	
	endif
	
	return
	
    end ! end meson_TS

!---------------------------------------------------
    function meson_SS(pt, IDP,IDP2, IDmeson)	
	implicit none
!	include 'Recom.h'
	external int_SS,S_j
	real*8 meson_SS
	real*8 pt,q,m_h,p0,tmp
	integer IDP, IDP2, IDmeson,i
    double precision fa,fb,fsum, int_SS,S_j
    
    if((IDmeson .eq. 1))then !jpsi
        m_h=3.097
    elseif(IDmeson.eq.2)then !D0
        m_h=1.8
    elseif(IDmeson.eq.3)then !D0
        m_h=1.96
    endif
    p0=Sqrt(m_h**2+pt**2)
    fsum=0.0d0

        fa=2.7!the cutoff
        fb=30.0
    if(pt.gt.fa) then
        tmp= pt
        fa = tmp
        fb = 30.0
    endif

 !   call integral_S(2.0d0,30.0d0,int_SS,0.000001,fsum,IDmeson,0.05d0)
 !   write(*,*)fsum
!------------------------------------------------------------	
	if(IDmeson .eq.1)then !for J/psi
	!call FFPTS(fa,fb,int_SS,0.000001,fsum,IDmeson,pt)
    do i=1,1000
        q=fa+i*(fb-fa)/1000.
        fsum=fsum+( int_SS(pt,q,IDmeson)+int_SS(pt,q+(fb-fa)/1000.,IDmeson) )*(fb-fa)/1000./2.
    enddo
	meson_SS = (1-dexp(-pt/2.0))*(fsum)/p0
	!(1-dexp(-pt/2.0))*   S_j(pt/2,8)*S_j(pt/2,9)
!------------------------------------------------

	elseif(IDmeson.eq.2)then !for D0
	    fa=3.2!the cutoff
        fb=30.0
    if(pt.gt.fa) then
        tmp= pt
        fa = tmp
        fb = 30.0
    endif
    
     do i=1,1000
        q=fa+i*(fb-fa)/1000
        fsum=fsum+( int_SS(pt,q,IDmeson)+int_SS(pt,q+(fb-fa)/1000,IDmeson) )*(fb-fa)/1000/2
    enddo
	meson_SS = (1-dexp(-pt/2.0))*(fsum)/p0
	
    elseif(IDmeson.eq.3)then ! for Ds
        
        fa=3.2!the cutoff
        fb=30.0
    if(pt.gt.fa) then
        tmp= pt
        fa = tmp
        fb = 30.0
    endif
    
     do i=1,1000
        q=fa+i*(fb-fa)/1000
        fsum=fsum+( int_SS(pt,q,IDmeson)+int_SS(pt,q+(fb-fa)/1000,IDmeson) )*(fb-fa)/1000/2
    enddo
	meson_SS = (1-dexp(-pt/2.0))*(fsum)/p0 
	endif
	
	return
    end ! end meson_SS
!---------------------------------------------------
!---------------------------------------------------
    function meson_SS2j(pt, IDP1,IDP2, IDmeson)	
	implicit none
!	include 'Recom.h'
	external int_SS2j
	real*8 meson_SS2j
	real*8 pt,q1, q2, k1, k2, m_h,p0,tmp
	integer IDP1, IDP2, IDmeson,i, j
    double precision fa,fb,fsum,fsum1,fsum2, int_SS2j
    
    if((IDmeson .eq. 1))then !jpsi
        m_h=3.097
    elseif(IDmeson.eq.2)then !D0
        m_h=1.8
    elseif(IDmeson.eq.3)then !Ds
        m_h=1.96
    endif
    p0=Sqrt(m_h**2+pt**2)
    fsum=0.0d0
    fsum1=0.0d0
    fsum2=0.0d0

	if(IDmeson .eq.1)then !for J/psi
        fa=2.7!the cutoff
        fb=30.0
        if(pt.gt.fa) then
            tmp= pt
            fa = tmp
        endif
    
    do i=1,200
        q1=fa+i*(fb-fa)/200.
        q2=q1+(fb-fa)/200.
        do j=1,200
            k1=fa+j*(fb-fa)/200.
            k2=k1+(fb-fa)/200.
            fsum1=fsum1+( int_SS2j(pt,q1,k1,IDmeson)+   &
                          int_SS2j(pt,q1,k2,IDmeson) )*(fb-fa)/200./2.
            fsum2=fsum2+( int_SS2j(pt,q2,k1,IDmeson)+   &
                          int_SS2j(pt,q2,k2,IDmeson) )*(fb-fa)/200./2.
        enddo
        fsum=fsum+(fsum1+fsum2)/2.*(fb-fa)/200.
        fsum1=0.0d0
        fsum2=0.0d0
    enddo
    
	meson_SS2j = (1-dexp(-pt/2.0))*(fsum)/p0/pt*10.**(-2)

	elseif(IDmeson.eq.2)then !for D0
	    fa=3.2!the cutoff
        fb=30.0
    if(pt.gt.fa) then
        tmp= pt
        fa = tmp
        fb = 30.0
    endif
    
    do i=1,100
        q1=fa+i*(fb-fa)/100
        q2=q1+(fb-fa)/100
        do j=1,100
            k1=fa+j*(fb-fa)/100
            k2=k1+(fb-fa)/100
            fsum=fsum+( int_SS2j(pt,q1,k1,IDmeson)+   &
                        int_SS2j(pt,q2,k2,IDmeson) )*(fb-fa)/100/2
        enddo
    enddo
	meson_SS2j = (1-dexp(-pt/2.0))*(fsum)/p0/pt**6*5*10.**(-2)
	
    elseif(IDmeson.eq.3)then ! for Ds
        
        fa=3.2!the cutoff
        fb=30.0
    if(pt.gt.fa) then
        tmp= pt
        fa = tmp
        fb = 30.0
    endif
    
    do i=1,100
        q1=fa+i*(fb-fa)/100
        q2=q1+(fb-fa)/100
        do j=1,100
            k1=fa+j*(fb-fa)/100
            k2=k1+(fb-fa)/100
            fsum=fsum+( int_SS2j(pt,q1,k1,IDmeson)+   &
                        int_SS2j(pt,q2,k2,IDmeson) )*(fb-fa)/100/2
        enddo
    enddo
	meson_SS2j = (1-dexp(-pt/2.0))*(fsum)/p0/pt**13*660*10.**(-2)
    endif
    
	return
    end ! end meson_SS(2j)
!---------------------------------------------------
!---------------------------------------------------
    function Thermal(pt, IDP)
    implicit none
    double precision thermal
    real*8  m_T, m_h, gama, T,pt
    integer :: g=6
    integer IDP
    double precision :: pi=3.141592653589793D0
    real*8 Ir, Ii, Kr, Ki
    DOUBLE PRECISION I0r(1),I0i(1),K1r(1),K1i(1)
    integer ierr, nz
    T=0.175
    !for J/psi vT=tanh(etaT)=0.3c,  gamc=0.26        eta=0.30952
    !for D0    vT=tanh(etaT)=0.42c, gamc=0.26        eta=0.447692
    if((IDP .eq. 8) .or. (IDP .eq. 9))then !c,cbar
        m_h=1.5
        m_T=Sqrt(m_h**2+pt**2)
        !Ir=sinh(0.1)*pt/T!test for 2.76 TeV
        Ir=sinh(0.30952)*pt/T!200 GeV
        Ii=0.d0
        !Kr=cosh(0.1)*m_T/T!test for 2.76 TeV
        Kr=cosh(0.30952)*m_T/T!200 GeV
        Ki=0.d0
        !gama=1.2!2.76 TeV
        gama=0.26!200 GeV and 5.5 TeV
    elseif((IDP .eq. 6) .or. (IDP .eq. 7))then !s,sbar
        m_h=0.46
        m_T=Sqrt(m_h**2+pt**2)
        Ir=sinh(0.30952)*pt/T
        Ii=0.d0
        Kr=cosh(0.30952)*m_T/T
        Ki=0.d0
        gama=0.8
    elseif((IDP .eq. 1) .or. (IDP .eq. 3))then !u,ubar
        m_h=0.26
        m_T=Sqrt(m_h**2+pt**2)
        Ir=sinh(0.447692)*pt/T
        Ii=0.d0
        Kr=cosh(0.447692)*m_T/T
        Ki=0.d0
        gama=1
    elseif((IDP .eq. 2) .or. (IDP .eq. 4))then !d,dbar
        m_h=0.26
        m_T=Sqrt(m_h**2+pt**2)
        Ir=sinh(0.447692)*pt/T
        Ii=0.d0
        Kr=cosh(0.447692)*m_T/T
        Ki=0.d0
        gama=1
    endif
    
    call ZBESI(Ir,Ii,0.d0,1,1,I0r,I0i,nz,ierr)
    call ZBESK(Kr,Ki,0.d0,1,1,K1r,K1i,nz,ierr)
    
    thermal = 2*g*gama*m_T/(2*pi)**3*I0r(1)*K1r(1)


    return
    end!function thermal(pt,IDP) end
    
    function S_ij(z, IDP1, IDP2)!!IDP1:inside i; IDP2:outside j
    implicit none
    integer IDP1, IDP2
    real*8 z, s1, s2, S_ij, A0, a, b, c, d
    s1=0.1d0
    s2=0.d0
    if(IDP1.eq.5)then!i=g
        if((IDP2.eq.8) .or. (IDP2.eq.9))then!j=c,cbar
            A0=5.020d-3
            a=0.702
            b=1.201
            c=1.952
            d=4.415
        elseif((IDP2.eq.1) .or. (IDP2.eq.2) .or. (IDP2.eq.3) .or. (IDP2.eq.4))then!j=u,d,ubar,dbar
            A0=0.919
            a=-0.0646
            b=2.0425
            c=-0.977
            d=1.528
        elseif((IDP2.eq.6) .or. (IDP2.eq.7))then!j=s,sbar
            A0=0.0151
            a=-0.8754
            b=2.128
            c=3.7674
            d=1.148
        endif
    elseif(IDP1.eq.8)then!i=c
        if((IDP2.eq.1) .or. (IDP2.eq.2) .or. (IDP2.eq.3) .or. (IDP2.eq.4))then!j=u,d,ubar,dbar
            A0=2.185
            a=0.102
            b=4.017
            c=-0.879
            d=1.066
        elseif(IDP2.eq.8)then!j=c
            A0=6.483
            a=2.526
            b=2.375
            c=22.819
            d=3.233
        elseif(IDP2.eq.9)then!j=cbar
            A0=5.879d-4
            a=0.315
            b=2.577
            c=23.505
            d=16.078
        elseif((IDP2.eq.6) .or. (IDP2.eq.7))then!j=s,sbar
            A0=0.118
            a=-0.138
            b=2.3
            c=0.9
            d=0.1
        endif
    elseif(IDP1.eq.1)then
        if(IDP2.eq.3)then!S_u_ubar
            A0=0.7316
            a=-0.0813
            b=3.22
            c=5.5386
            d=5.1416
        elseif(IDP2.eq.7)then!S_u_sbar
            A0=2.0842
            a=0.8722
            b=8.377
            c=94.412
            d=3.0606
        endif
    elseif(IDP1.eq.2)then
        if(IDP2.eq.3)then!S_d_ubar
            A0=0.7316
            a=-0.0813
            b=3.22
            c=5.5386
            d=5.1416
        elseif(IDP2.eq.7)then!S_d_sbar
            A0=2.0842
            a=0.8722
            b=8.377
            c=94.412
            d=3.0606
        endif
    elseif(IDP1.eq.3)then
        if(IDP2.eq.3)then!S_ubar_ubar=K_NS+L=K
            A0=0.0101
            a=0.0256
            b=1.3418
            c=174.471
            d=1.093
            s1=A0*z**a*(1-z)**b*(1+c*z**d)!K_NS
            A0=0.7316
            a=-0.0813
            b=3.22
            c=5.5386
            d=5.1416
            s2=A0*z**a*(1-z)**b*(1+c*z**d)!L
        elseif(IDP2.eq.7)then!S_ubar_sbar
            A0=2.0842
            a=0.8722
            b=8.377
            c=94.412
            d=3.0606
        endif
    elseif(IDP1.eq.4)then
        if(IDP2.eq.3)then!S_dbar_ubar
            A0=0.7316
            a=-0.0813
            b=3.22
            c=5.5386
            d=5.1416
        elseif(IDP2.eq.7)then!S_dbar_sbar
            A0=2.0842
            a=0.8722
            b=8.377
            c=94.412
            d=3.0606
        endif
    elseif(IDP1.eq.6)then
        if(IDP2.eq.7)then!S_s_sbar
            A0=2.0842
            a=0.8722
            b=8.377
            c=94.412
            d=3.0606
        endif
    elseif(IDP1.eq.7)then
        if(IDP2.eq.7)then!S_sbar_sbar
            A0=0.0101
            a=0.0256
            b=1.3418
            c=174.471
            d=1.093
            s1=A0*z**a*(1-z)**b*(1+c*z**d)!K_NS
            A0=2.0842
            a=0.8722
            b=8.377
            c=94.412
            d=3.0606
            s2=A0*z**a*(1-z)**b*(1+c*z**d)!L_s
        endif
    endif
    
    if(s1.eq.0.1d0)then
        S_ij = A0*z**a*(1-z)**b*(1+c*z**d)
    else
        S_ij = s1+s2
    endif
            
    return
    end!function Shower_ij(z, IDP1, IDP2) end
    
    function S_j(pt, IDP2)
    implicit none

    integer IDP2,i
    external int_Sj
    double precision pt, S_j, fa, fb, int_Sj
    real*8 plowlim, phighlim,q!plowlim=0.05,phighlim=6.5
    plowlim=3.0!for Ds  !2.8!for D0!  2.8!2.8 for Jpsi!the cutoff
    phighlim=50.0
    fa=plowlim
    fb=phighlim
    S_j=0.0d0
    if(pt.ge.plowlim) then
        fa = pt
    endif
    
    if(IDP2 .eq. 8)then!S_c
        !call FFPTS(fa,fb,int_Sj,0.000001,S_j,IDP2,pt)
        do i=1,1000
            q=i*(fb-fa)/1000
            S_j=S_j+( int_Sj(pt,fa+q,IDP2)+int_Sj(pt,fa+q+(fb-fa)/1000,IDP2) )*(fb-fa)/1000/2
        enddo
    elseif(IDP2 .eq. 9)then!S_cbar
        !call FFPTS(fa,fb,int_Sj,0.000001,S_j,IDP2,pt)
        do i=1,1000
            q=i*(fb-fa)/1000
            S_j=S_j+( int_Sj(pt,fa+q,IDP2)+int_Sj(pt,fa+q+(fb-fa)/1000,IDP2) )*(fb-fa)/1000/2
        enddo
    elseif(IDP2 .eq. 3)then!S_ubar
        do i=1,1000
            q=i*(fb-fa)/1000
            S_j=S_j+( int_Sj(pt,fa+q,IDP2)+int_Sj(pt,fa+q+(fb-fa)/1000,IDP2) )*(fb-fa)/1000/2
        enddo
    elseif(IDP2 .eq. 7)then!S_sbar
        do i=1,1000
            q=i*(fb-fa)/1000
            S_j=S_j+( int_Sj(pt,fa+q,IDP2)+int_Sj(pt,fa+q+(fb-fa)/1000,IDP2) )*(fb-fa)/1000/2
        enddo
    endif
    S_j=S_j
    
    !test--------------------------------
    !call FFPTS(3.0d0,50.0d0,int_Sj,0.000001,S_j,IDP2,3.0d0)
    !write(*,*)3.0,S_j
    !call FFPTS(3.2d0,50.0d0,int_Sj,0.000001,S_j,IDP2,3.0d0)
    !write(*,*)3.2,S_j
    !------------------------------------
    
    return
    end!function S_j(pt, IDP2) end
    
    function int_Sj(pt,q, IDP)!function for integral
    implicit none
    external Fiq, S_ij
    real*8 Fiq, S_ij, pt, q
    double precision int_Sj
    integer IDP
    
    if(IDP .eq. 8)then!S_c from g/c
        int_Sj = (Fiq(q, 8)*S_ij(pt/q, 8, 8)+ Fiq(q, 5)*S_ij(pt/q, 5,8))/q
    elseif(IDP .eq. 9)then!S_cbar from g/c
        int_Sj = (Fiq(q, 8)*S_ij(pt/q, 8, 9)+ Fiq(q, 5)*S_ij(pt/q, 5,9))/q
    elseif(IDP .eq. 3)then!S_ubar from light(bar)/g
        int_Sj =( Fiq(q, 1)*S_ij(pt/q, 1, 3)+ Fiq(q, 2)*S_ij(pt/q, 2, 3)  &
                + Fiq(q, 3)*S_ij(pt/q, 3, 3)+ Fiq(q, 4)*S_ij(pt/q, 4, 3)  &
                + Fiq(q, 5)*S_ij(pt/q, 5, 3) )/q
    elseif(IDP .eq. 7)then!S_sbar from light(bar)/g/s(bar)
        int_Sj =( Fiq(q, 1)*S_ij(pt/q, 1, 7)+ Fiq(q, 2)*S_ij(pt/q, 2, 7)  &
                + Fiq(q, 3)*S_ij(pt/q, 3, 7)+ Fiq(q, 4)*S_ij(pt/q, 4, 7)  &
                + Fiq(q, 5)*S_ij(pt/q, 5, 7)+ Fiq(q, 6)*S_ij(pt/q, 6, 7)  &
                + Fiq(q, 7)*S_ij(pt/q, 7, 7))/q
        
    endif

    return
    end!end int_Sj(pt,q, IDP)
    
    
    
!    function fik(k,IDP)!AuAu@200GeV at RHIC
!    implicit none
!    real*8 Ai, ki, ni
!    double precision k, fik
!    integer IDP
! 	select case(IDP)
!! u quark		
!	case(1)
!	
!	Ai = 9.113d2
!	ki = 1.459d0
!	ni = 7.679d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! d quark		
!	case(2)
!	
!	Ai = 9.596d2
!	ki = 1.467d0
!	ni = 7.662d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! ubar quark		
!	case(3)
!	
!	Ai = 2.031d2
!	ki = 1.767d0
!	ni = 8.546d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! dbar quark		
!	case(4)
!	
!	Ai = 2.013d2
!	ki = 1.759d0
!	ni = 8.566d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! gluon		
!	case(5)
!	
!	Ai = 4.455d3
!	ki = 1.7694
!	ni = 8.610d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! s, sbar quark			
!	case(6, 7)
!	
!	Ai = 1.038d2
!	ki = 1.868d0
!	ni = 8.642d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! c, cbar quark in fik[15]
!    case(8, 9)
!    
!    fik = 19.2*(1+(k/6)**2)/(1+k/3.7)**12/(1+dexp(0.9-2*k))
!
!	case default
!	
!	print *, "Need to check the parton IDP in fik"
!	stop
!
!    end select
!
!    return
!    end!function fik(k,IDP) end
    
!    function fik(k,IDP)!PbPb@2.76 TeV at LHC
!    implicit none
!    real*8 Ai, ki, ni
!    double precision k, fik
!    integer IDP
! 	select case(IDP)
!! u quark		
!	case(1)
!	
!	Ai = 1.138d4
!	ki = 0.687d0
!	ni = 5.67d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! d quark		
!	case(2)
!	
!	Ai = 1.266d4
!	ki = 0.677d0
!	ni = 5.66d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! ubar quark		
!	case(3)
!	
!	Ai = 0.24d4
!	ki = 0.87d0
!	ni = 5.97d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! dbar quark		
!	case(4)
!	
!	Ai = 0.23d4
!	ki = 0.88d0
!	ni = 5.99d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! gluon		
!	case(5)
!	
!	Ai = 6.2d4
!	ki = 0.98d0
!	ni = 6.22d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! s, sbar quark			
!	case(6, 7)
!	
!	Ai = 0.093d4
!	ki = 1.05d0
!	ni = 6.12d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! c, cbar quark by linear relation between fik and sNN, interpolating between 5.5 and 0.2 TeV.
!    case(8, 9)
!        
!    !fik = (2497*(1 +                                 &
!    !        39.37226/k**5.5))/                       &
!    !    (k*(1 + 16/(0.1 + k)**2)) +                  &
!    !   0.51698*                                      &
!    !    ((19.2*(1 + k**2/36.))/                      &
!    !       ((1 + Exp(0.9 - 2*k))*                    &
!    !         (1 + 0.27027*k)**12)                    &
!    !        - (2497*                                 &  
!    !         (1 + 39.37226/k**5.5))                  &
!    !        /(k*(1 + 16/(0.1 + k)**2)))!fik vs s_NN
!            
!    fik = (201.5857*(0.1 + k)**2*                   &
!          ((1 + k**2/36.)/                          &
!             ((1 + 2.4596/Exp(2*k))*                &
!               (1 + 0.27027*k)**12)                 &
!             )**0.51698*                            &
!          (39.37226 + k**5.5))/                     &
!        (k**6.5*((1 +                               &
!               39.37226/k**5.5)/                    &
!             (k + (16.*k)/(0.1 + k)**2))**          &
!           0.51698*                                 &
!          (16.01 + k*(0.2 + k)))! ln fik vs s_NN
!
!	case default
!	
!	print *, "Need to check the parton IDP in fik"
!	stop
!
!    end select
!
!    return
!    end!function fik(k,IDP) end
    
!    function fik(k,IDP)!PbPb@5.02 TeV at LHC
!    implicit none
!    real*8 Ai, ki, ni
!    double precision k, fik
!    integer IDP
! 	select case(IDP)
!! u quark		
!	case(1)
!	
!	Ai = 2.023d4
!	ki = 0.578d0
!	ni = 5.295d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! d quark		
!	case(2)
!	
!	Ai = 2.279d4
!	ki = 0.567d0
!	ni = 5.278d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! ubar quark		
!	case(3)
!	
!	Ai = 0.42d4
!	ki = 0.74d0
!	ni = 5.51d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! dbar quark		
!	case(4)
!	
!	Ai = 0.397d4
!	ki = 0.75d0
!	ni = 5.52d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! gluon		
!	case(5)
!	
!	Ai = 11.2d4
!	ki = 0.789d0
!	ni = 5.67d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! s, sbar quark			
!	case(6, 7)
!	
!	Ai = 0.154d4
!	ki = 0.925d0
!	ni = 5.616d0
!	
!	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
!! c, cbar quark by linear relation between fik and sNN, interpolating between 5.5 and 0.2 TeV.
!    case(8, 9)
!    
    !fik = (2497*(1 +                            &
    !        39.37226/k**5.5))/                  &
    !    (k*(1 + 16/(0.1 + k)**2)) +             &
    !   0.090566*                                &
    !    ((19.2*(1 + k**2/36.))/                 &
    !       ((1 + Exp(0.9 - 2*k))*               &
    !         (1 + 0.27027*k)**12)               &
    !        - (2497*                            &
    !         (1 + 39.37226/k**5.5))             &
    !        /(k*(1 + 16/(0.1 + k)**2)))!fik vs s_NN
!
!    fik = (1606.76*(0.1 + k)**2*               &
!         ((1 + k**2/36.)/                      &
!            ((1 + 2.4596/Exp(2*k))*            &
!              (1 + 0.27027*k)**12)             &
!            )**0.090566*                       &
!         (39.37226 + k**5.5))/                 &
!       (k**6.5*((1 +                           &
!              39.37226/k**5.5)/                &
!            (k + (16.*k)/(0.1 + k)**2))**      &
!          0.090566*                            &
!         (16.01 + k*(0.2 + k)))! ln fik vs s_NN
!
!	case default
!	
!	print *, "Need to check the parton IDP in fik"
!	stop
!
!    end select
!
!    return
!    end!function fik(k,IDP) end

    function fik(k,IDP)!PbPb@5.5 TeV at LHC
    implicit none
    real*8 Ai, ki, ni
    double precision k, fik
    integer IDP
 	select case(IDP)
! u quark		
	case(1)
	
	Ai = 2.209d4
	ki = 0.5635d0
	ni = 5.240d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! d quark		
	case(2)
	
	Ai = 2.493d4
	ki = 0.5522d0
	ni = 5.223d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! ubar quark		
	case(3)
	
	Ai = 4.581d3
	ki = 0.7248d0
	ni = 5.437d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! dbar quark		
	case(4)
	
	Ai = 4.317d3
	ki = 0.7343d0
	ni = 5.448d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! gluon		
	case(5)
	
	Ai = 1.229d5
	ki = 0.7717d0
	ni = 5.600d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! s, sbar quark			
	case(6, 7)
	
	Ai = 1.662d3
	ki = 0.9064d0
	ni = 5.548d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! c, cbar quark in Ref.meson_fik_charm_5.5
    case(8, 9)
    
    fik = 2497*(1+k/1.95)**(-5.5)*(   k*(1+(4/(1+dexp(0.1+k)))**2)     )**(-1)
	case default
	
	print *, "Need to check the parton IDP in fik"
	stop

    end select

    return
    end!function fik(k,IDP) end
    
    
    function f1ik(k,IDP)!f'ik
    implicit none
    external fik
    double precision k, fik, f1ik
    double precision :: pi=3.141592653589793D0
    real*8 E, m_h
    integer IDP
    if((IDP .eq. 8) .or. (IDP .eq. 9))then !c,cbar
        m_h=1.5
    elseif((IDP .eq. 6) .or. (IDP .eq. 7))then !s,sbar
        m_h=0.46
    elseif((IDP .eq. 1) .or. (IDP .eq. 3))then !u,ubar
        m_h=0.26
    elseif((IDP .eq. 2) .or. (IDP .eq. 4))then !d,dbar
        m_h=0.26
    elseif(IDP.eq.5)then!gluon
        m_h=0.1
    endif
    E=Sqrt(m_h**2+k**2)!k**2/2/m_h!Sqrt(m_h**2+k**2)!？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
    f1ik = (2*pi)**3/E*fik(k,IDP)/k

    return
    end!function f1ik(k,IDP)!f'ik end
    
    function Fiq(q,IDP)!with energy
    implicit none
    external f1ik!f'ik
    double precision q, fsum, f1ik, q1, tmp
    real*8  Fiq
    real*8  bL
    integer IDP,i

    if(IDP.eq.5)then
        bL=2.39
    elseif(IDP.eq.8 .or. IDP.eq.9)then
        bL=0.01!for Jpsi bL=0.01 while for D meson bL=2.39, the same as light quark 
        !bL=2.39!choose the same energy loss factor for charm            fik with energy
    elseif(IDP.eq.1 .or. IDP.eq.2 .or. IDP.eq.3 .or. IDP.eq.4 .or. IDP.eq.6.or. IDP.eq.7)then
        bL=2.39
    endif
    
    fsum=0.d0
    do i=1,500
            q1=q+i*(q*exp(bL)-q)/500.
            tmp=q1+(q*exp(bL)-q)/500.
            fsum=fsum+(f1ik(q1,IDP)+f1ik(tmp,IDP))/2*(q*exp(bL)-q)/500.
    enddo
    Fiq = 1/bL*fsum
    return
    end!Fiq(q,IDP) end
    
    function int_SS(p,q,IDmeson)!SS for integral
    implicit none 
    external S_ij, F1iq, int_DD0c, int_DD0g, int_DDsc, int_DDsg
    real*8 S_ij, int_DD0c,int_DD0g,x1,tmp, int_DDsc, int_DDsg
    double precision  F1iq, fsum, p, q,int_SS, Dc, Dg
    integer IDmeson,i
    
    Dc=0.d0
    Dg=0.d0
    if(IDmeson.eq.1)then
        Dc=0.5/(p/q)*                                            &
        ( S_ij(p/q/2,8,8)*S_ij(p/q/(2-p/q),8,9)+S_ij(p/q/2,8,9)*S_ij(p/q/(2-p/q),8,8))
        Dg=0.5/(p/q)*                                            &
        ( S_ij(p/q/2,5,8)*S_ij(p/q/(2-p/q),5,9)+S_ij(p/q/2,5,9)*S_ij(p/q/(2-p/q),5,8))
        
        int_SS =( (F1iq(q,8)*Dc)+(F1iq(q,5)*Dg)  )/q**2
    
    elseif(IDmeson.eq.2)then
        do i=1,498
            x1=i*(p/q-0.)/500
            tmp=x1+(p/q-0.)/500
            Dc=Dc+(int_DD0c(p/q,x1)+int_DD0c(p/q,tmp))/2*(p/q-0.)/500
        enddo
        
        do i=1,498
            x1=i*(p/q-0.)/500
            tmp=x1+(p/q-0.)/500
            Dg=Dg+(int_DD0g(p/q,x1)+int_DD0g(p/q,tmp))/2*(p/q-0.)/500
        enddo
        !write(*,*)p/q,Dc,Dg
        int_SS =( (F1iq(q,8)*Dc)+(F1iq(q,5)*Dg)  )/q**2
        
    elseif(IDmeson.eq.3)then
        do i=1,498
            x1=i*(p/q-0.)/500
            tmp=x1+(p/q-0.)/500
            Dc=Dc+(int_DDsc(p/q,x1)+int_DDsc(p/q,tmp))/2*(p/q-0.)/500
        enddo
        
        do i=1,498
            x1=i*(p/q-0.)/500
            tmp=x1+(p/q-0.)/500
            Dg=Dg+(int_DDsg(p/q,x1)+int_DDsg(p/q,tmp))/2*(p/q-0.)/500
        enddo
        !write(*,*)p/q,Dc,Dg
        int_SS =( (F1iq(q,8)*Dc)+(F1iq(q,5)*Dg)  )/q**2
    endif

    return
    end!int_SS(p,q,IDmeson) end
    
    
    function int_DD0c(x,x1)
    implicit none
    real*8 int_DD0c, x,x1
    external S_ij
    double precision S_ij
    
    int_DD0c = ( S_ij(x1,8,3)*S_ij((x-x1)/(1-x1),8,8)+       &
    S_ij(x1/(1-x+x1),8,3)*S_ij(x-x1,8,8) )/2*5*(x-x1)**4/x**6
   !write(*,*)x,x1,int_DD0c
    return
    end
    

    function int_DD0g(x,x1)
    implicit none
    real*8 int_DD0g, x,x1
    external S_ij
    double precision S_ij
    
    int_DD0g = ( S_ij(x1,5,3)*S_ij((x-x1)/(1-x1),5,8)+       &
    S_ij(x1/(1-x+x1),5,3)*S_ij(x-x1,5,8) )/2*5*(x-x1)**4/x**6
    !write(*,*)(x-x1)/(1-x1),S_ij((x-x1)/(1-x1),5,8)
    return
    end
    
    function int_DDsc(x,x1)
    implicit none
    real*8 int_DDsc, x,x1
    external S_ij
    double precision S_ij
    
    int_DDsc = ( S_ij(x1,8,7)*S_ij((x-x1)/(1-x1),8,8)+       &
    S_ij(x1/(1-x+x1),8,7)*S_ij(x-x1,8,8) )*330*x1**2 *(x-x1)**9/x**13
   !write(*,*)x,x1,int_DD0c
    return
    end
    

    function int_DDsg(x,x1)
    implicit none
    real*8 int_DDsg, x,x1
    external S_ij
    double precision S_ij
    
    int_DDsg = ( S_ij(x1,5,7)*S_ij((x-x1)/(1-x1),5,8)+       &
    S_ij(x1/(1-x+x1),5,7)*S_ij(x-x1,5,8) )*330*x1**2 *(x-x1)**9/x**13
    !write(*,*)(x-x1)/(1-x1),S_ij((x-x1)/(1-x1),5,8)
    return
    end
    
    
    function F1iq(q,IDP)!F'iq(without energy E)
    implicit none
    external int_fik
    double precision q, fsum, int_fik, q1, tmp
    real*8  F1iq
    real*8  bL
    integer IDP, i
    
    if(IDP.eq.5)then
        !bL=2.39
        bL=2.9 !for gluon
    elseif(IDP.eq.8 .or. IDP.eq.9)then
        bL=0.01!for Jpsi bL=0.01 while for D meson bL=2.39, the same as light quark 
        !bL=2.39
    elseif(IDP.eq.1 .or. IDP.eq.2 .or. IDP.eq.3 .or. IDP.eq.4 .or. IDP.eq.6.or. IDP.eq.7)then
        bL=2.39
    endif
    
    fsum=0.d0
    do i=1,500
            q1=q+i*(q*exp(bL)-q)/500.
            tmp=q1+(q*exp(bL)-q)/500.
            fsum=fsum+(int_fik(q1,IDP)+int_fik(tmp,IDP))/2.*(q*exp(bL)-q)/500.
    enddo
    F1iq = 1/bL*fsum

    return
    end!F1iq(q,IDP) end
    
    function int_fik(k,IDP)
    implicit none
    external fik
    double precision int_fik, fik
    real*8 k
    integer IDP
    int_fik=fik(k, IDP)*k
    return
    end!fik/k for integral end
    
    function int_SS2j(pt, q, k,IDmeson)!SS2j for integral, q=q and k=q'
    implicit none 
    external S_ij, F1iq, int_fss,Fiq
    real*8 S_ij, int_fss, p2, p2_next
    double precision  F1iq, fsum, pt, q, k, int_SS2j,Fiq
    real*8 fs1s5, fs2s5, fs3s5, fs4s5, fs5s5, fs6s5, fs7s5,  fs1s8, fs2s8, fs3s8, fs4s8, fs5s8, fs6s8,fs7s8
    integer IDmeson,i

    fs1s5=0.d0
    fs2s5=0.d0
    fs3s5=0.d0
    fs4s5=0.d0
    fs5s5=0.d0
    fs6s5=0.d0
    fs7s5=0.d0
    
    fs1s8=0.d0
    fs2s8=0.d0
    fs3s8=0.d0
    fs4s8=0.d0
    fs5s8=0.d0
    fs6s8=0.d0
    fs7s8=0.d0
    
    if(IDmeson.eq.1)then!Jpsi
        !sum over i=c/g, i'=c/g, totally 4 terms
        !int_SS2j = Fiq(q, 8)*Fiq(k, 8)*S_ij(pt/2/q,8,8)*S_ij(pt/2/k,8,9)+ &
        !           Fiq(q, 5)*Fiq(k, 8)*S_ij(pt/2/q,5,8)*S_ij(pt/2/k,8,9)+ &
        !           Fiq(q, 8)*Fiq(k, 5)*S_ij(pt/2/q,8,8)*S_ij(pt/2/k,5,9)+ &
        !           Fiq(q, 5)*Fiq(k, 5)*S_ij(pt/2/q,5,8)*S_ij(pt/2/k,5,9)
        int_SS2j = F1iq(q, 8)*F1iq(k, 8)*S_ij(pt/2/q,8,8)*S_ij(pt/2/k,8,9)+ &
                   F1iq(q, 5)*F1iq(k, 8)*S_ij(pt/2/q,5,8)*S_ij(pt/2/k,8,9)+ &
                   F1iq(q, 8)*F1iq(k, 5)*S_ij(pt/2/q,8,8)*S_ij(pt/2/k,5,9)+ &
                   F1iq(q, 5)*F1iq(k, 5)*S_ij(pt/2/q,5,8)*S_ij(pt/2/k,5,9)
        
    elseif(IDmeson.eq.2)then!D0
        do i=1,98!dp_2:0-pt
            p2=i*pt/100
            p2_next=p2+pt/100
            fs1s5=fs1s5+(int_fss(pt,q,k,p2,IDmeson,1,5)+int_fss(pt,q,k,p2_next,IDmeson,1,5))/2*pt/100
            fs2s5=fs2s5+(int_fss(pt,q,k,p2,IDmeson,2,5)+int_fss(pt,q,k,p2_next,IDmeson,2,5))/2*pt/100
            fs3s5=fs3s5+(int_fss(pt,q,k,p2,IDmeson,3,5)+int_fss(pt,q,k,p2_next,IDmeson,3,5))/2*pt/100
            fs4s5=fs4s5+(int_fss(pt,q,k,p2,IDmeson,4,5)+int_fss(pt,q,k,p2_next,IDmeson,4,5))/2*pt/100
            fs5s5=fs5s5+(int_fss(pt,q,k,p2,IDmeson,5,5)+int_fss(pt,q,k,p2_next,IDmeson,5,5))/2*pt/100
            
            fs1s8=fs1s8+(int_fss(pt,q,k,p2,IDmeson,1,8)+int_fss(pt,q,k,p2_next,IDmeson,1,8))/2*pt/100
            fs2s8=fs2s8+(int_fss(pt,q,k,p2,IDmeson,2,8)+int_fss(pt,q,k,p2_next,IDmeson,2,8))/2*pt/100
            fs3s8=fs3s8+(int_fss(pt,q,k,p2,IDmeson,3,8)+int_fss(pt,q,k,p2_next,IDmeson,3,8))/2*pt/100
            fs4s8=fs4s8+(int_fss(pt,q,k,p2,IDmeson,4,8)+int_fss(pt,q,k,p2_next,IDmeson,4,8))/2*pt/100
            fs5s8=fs5s8+(int_fss(pt,q,k,p2,IDmeson,5,8)+int_fss(pt,q,k,p2_next,IDmeson,5,8))/2*pt/100
            
        enddo
        !sum over i=1-5, i'=c/g, totally 10 terms
        int_SS2j = F1iq(q, 1)*F1iq(k, 5)*fs1s5+ F1iq(q, 2)*F1iq(k, 5)*fs2s5+    &
                   F1iq(q, 3)*F1iq(k, 5)*fs3s5+ F1iq(q, 4)*F1iq(k, 5)*fs4s5+    &
                   F1iq(q, 5)*F1iq(k, 5)*fs5s5+ F1iq(q, 1)*F1iq(k, 8)*fs1s8+    &
                   F1iq(q, 2)*F1iq(k, 8)*fs2s8+ F1iq(q, 3)*F1iq(k, 8)*fs3s8+    &
                   F1iq(q, 4)*F1iq(k, 8)*fs4s8+ F1iq(q, 5)*F1iq(k, 8)*fs5s8
        
        
    elseif(IDmeson.eq.3)then!Ds
        do i=1,98!dp_2:0-pt
            p2=i*pt/100
            p2_next=p2+pt/100
            fs1s5=fs1s5+(int_fss(pt,q,k,p2,IDmeson,1,5)+int_fss(pt,q,k,p2_next,IDmeson,1,5))/2*pt/100
            fs2s5=fs2s5+(int_fss(pt,q,k,p2,IDmeson,2,5)+int_fss(pt,q,k,p2_next,IDmeson,2,5))/2*pt/100
            fs3s5=fs3s5+(int_fss(pt,q,k,p2,IDmeson,3,5)+int_fss(pt,q,k,p2_next,IDmeson,3,5))/2*pt/100
            fs4s5=fs4s5+(int_fss(pt,q,k,p2,IDmeson,4,5)+int_fss(pt,q,k,p2_next,IDmeson,4,5))/2*pt/100
            fs5s5=fs5s5+(int_fss(pt,q,k,p2,IDmeson,5,5)+int_fss(pt,q,k,p2_next,IDmeson,5,5))/2*pt/100
            fs6s5=fs6s5+(int_fss(pt,q,k,p2,IDmeson,6,5)+int_fss(pt,q,k,p2_next,IDmeson,6,5))/2*pt/100
            fs7s5=fs7s5+(int_fss(pt,q,k,p2,IDmeson,7,5)+int_fss(pt,q,k,p2_next,IDmeson,7,5))/2*pt/100
            
            fs1s8=fs1s8+(int_fss(pt,q,k,p2,IDmeson,1,8)+int_fss(pt,q,k,p2_next,IDmeson,1,8))/2*pt/100
            fs2s8=fs2s8+(int_fss(pt,q,k,p2,IDmeson,2,8)+int_fss(pt,q,k,p2_next,IDmeson,2,8))/2*pt/100
            fs3s8=fs3s8+(int_fss(pt,q,k,p2,IDmeson,3,8)+int_fss(pt,q,k,p2_next,IDmeson,3,8))/2*pt/100
            fs4s8=fs4s8+(int_fss(pt,q,k,p2,IDmeson,4,8)+int_fss(pt,q,k,p2_next,IDmeson,4,8))/2*pt/100
            fs5s8=fs5s8+(int_fss(pt,q,k,p2,IDmeson,5,8)+int_fss(pt,q,k,p2_next,IDmeson,5,8))/2*pt/100
            fs6s8=fs6s8+(int_fss(pt,q,k,p2,IDmeson,6,8)+int_fss(pt,q,k,p2_next,IDmeson,6,8))/2*pt/100
            fs7s8=fs7s8+(int_fss(pt,q,k,p2,IDmeson,7,8)+int_fss(pt,q,k,p2_next,IDmeson,7,8))/2*pt/100
            
        enddo
        !sum over i=1-7, i'=c/g, totally 14 terms
        int_SS2j = F1iq(q, 1)*F1iq(k, 5)*fs1s5+ F1iq(q, 2)*F1iq(k, 5)*fs2s5+    &
                   F1iq(q, 3)*F1iq(k, 5)*fs3s5+ F1iq(q, 4)*F1iq(k, 5)*fs4s5+    &
                   F1iq(q, 5)*F1iq(k, 5)*fs5s5+ F1iq(q, 6)*F1iq(k, 5)*fs6s5+    &
                   F1iq(q, 7)*F1iq(k, 5)*fs7s5+ F1iq(q, 1)*F1iq(k, 8)*fs1s8+    &
                   F1iq(q, 2)*F1iq(k, 8)*fs2s8+ F1iq(q, 3)*F1iq(k, 8)*fs3s8+    &
                   F1iq(q, 4)*F1iq(k, 8)*fs4s8+ F1iq(q, 5)*F1iq(k, 8)*fs5s8+    &
                   F1iq(q, 6)*F1iq(k, 8)*fs6s8+ F1iq(q, 7)*F1iq(k, 8)*fs7s8
                   
    endif

    return
    end!int_SS2j(p,q,IDmeson) end
    
    
    function int_fss(pt, q, k ,p2,IDmeson, ID1, ID2)!function to integral to form SS2j
    implicit none
    real*8 int_fss, pt, q, k ,p2
    integer IDmeson, ID1, ID2!ID1:i, ID2:i'
    external S_ij
    double precision S_ij
    
    if(IDmeson .eq. 2)then ! D0
        int_fss = p2**4*S_ij((pt-p2)/q, ID1,3)*S_ij(p2/k, ID2,8)
    
    elseif(IDmeson .eq. 3)then ! Ds
        int_fss = p2**9*(pt-p2)**2 *S_ij((pt-p2)/q, ID1,7)*S_ij(p2/k, ID2,8)
        
    endif
    
    return
    end!int_fss end
    
!===================================================
! beta function	
        SUBROUTINE BETA(P,Q,BT)
!       ==========================================
!       Purpose: Compute the beta function B(p,q)
!       Input :  p  --- Parameter  ( p > 0 )
!                q  --- Parameter  ( q > 0 )
!       Output:  BT --- B(p,q)
!       Routine called: GAMMA for computing ‚(x)
!       ==========================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        CALL GAMMA(P,GP)
        CALL GAMMA(Q,GQ)
        PPQ=P+Q
        CALL GAMMA(PPQ,GPQ)
        BT=GP*GQ/GPQ
        RETURN
        END!beta function end

! gamma function	
        SUBROUTINE GAMMA(X,GA)
!
!       ==================================================
!       Purpose: Compute gamma function ‚(x)
!       Input :  x  --- Argument of ‚(x)
!                       ( x is not equal to 0,-1,-2,˙˙˙)
!       Output:  GA --- ‚(x)
!       ==================================================
!
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           
              G=(/1.0D0,0.5772156649015329D0,                       &
               -0.6558780715202538D0, -0.420026350340952D-1,        &  
               0.1665386113822915D0,-.421977345555443D-1,           &
               -.96219715278770D-2, .72189432466630D-2,             &
               -.11651675918591D-2, -.2152416741149D-3,             & 
               .1280502823882D-3, -.201348547807D-4,                &
               -.12504934821D-5, .11330272320D-5,                   &
               -.2056338417D-6, .61160950D-8,                       &
               .50020075D-8, -.11812746D-8,                         &
               .1043427D-9, .77823D-11,                             &
               -.36968D-11, .51D-12,                                &
               -.206D-13, -.54D-14, .14D-14, .1D-15/)               
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
    END!gamma function end
    