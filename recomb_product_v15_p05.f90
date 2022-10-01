! This program is written by hua to calcuate the yields of meson and baryon from Recombination model
! treat u and d separately Dec-26-2013
! We can calculate pion, kaon and proton results Dec-28-2013
! We decide to fix gamma_q = 0.45 and gamma_g = 0.9. Jan-7-2014
! \gamma_g=.8/(1+q/10) and \gamma_g=\gamma_q/2. Apr-8-14, lilin
! Add the contributions meson_SS2j, baryon_TSS2j, baryon_SSS2j, Apr-8-2014, hua
! Add subroutine Integral_2d(fnx, pt, IDP1, IDP2, IDP3, centr, sumf) to optimize the code, Apr-9-2014, hua
! Add the function FF_Wang, Apr-9-2014
! Add code for Lambda particle
! debug the NAN in function Intf_SS1J, Apr-10-2014
! Add the code for Xi and Phi, Apr-17-2014, hua
! Add the periodic condition when we calculate the xibar, Aug-10-2014
! Add flags in module Iflag, Feb-5-2015
! Add fik for pPb at 5.02 TeV from v12, Mar-11-2015 by Hua
! Add fik for PbPb at 5.02 TeV and XeXe at 5.44 TeV, May-15-2020 by Lilin
! Add plowlim and phighlim
! Before running, fix system (ISys), centrality (c), particle (ISu), T_q, T_s and p_T range
! Then adjust gamma factor, gamma_0 and q_0 to fit the experimental data
! Fail to use Minuit to search for optimized parameters because of time cosuming
! Change the strategy to scan for the optimized parameters, May-25-2020 by Hua
! Data should be d^2N/2piptdpt,  May-25-2020 by Hua
!===============================================================
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Check the several places before you run the code
! 1) gamai in function Intf_fik(intvar, wgt, IDH, q, c) ! make sure it is right
! 2) Tabc, TabcRatio, in function Intf_fik(intvar, wgt, IDH, q, c)
!    This is fitted from experimental data. It is different for different system
! 3) choose the right range for pt in main program
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!========================================
! The parton ID
! 1, 2, 3,    4,    5,   6,  7
! u, d, ubar, dbar, g,   s,  sbar

! particle 
! pi+ (u dbar) (1, 4)
! k+  (u,sbar) (1, 7)
! p   (u u d)  (1, 1, 2)
! lambda^0 (u, d, s) (1, 2, 6)
!========================================
! IYangZhu, =1 use Yang's equations; =2 use Zhu's equations
! ISys, =1 AuAu, =2 PbPb, =3 dAu
! IResPi, =0 no resonance contribution, =1 with resonance contribution

    module Iflag
    implicit none

    integer IYangZhu, ISys, IResPi

    end module Iflag

!---------------------------------------
! module to include the integral limits for hard jets

    module limitsHJ
    implicit none

    real*8 plowlim, phighlim

    end module limitsHJ

!========================================

	program recomb_product
    use Iflag
    use limitsHJ
	implicit none
	include 'Recom.h'

! functions to calculate different components for meson
	external meson_TT, meson_TS, meson_SS, meson_SS2j
	real*8 meson_TT, meson_TS, meson_SS, meson_SS2j

! IDP1, parton 1 id
! Isu, represent the particle
	real*8 pt, dNptdpt(7)
	integer IDP1, IDP2, IDP3, Isu
	
	integer i
!-------------------------------------------

! functions to calculate different components for baryon
	external baryon_TTT, baryon_TTS, baryon_TSS, baryon_TSS2j, baryon_SSS, baryon_SSS2j
	real*8 baryon_TTT, baryon_TTS,baryon_TSS, baryon_TSS2j, baryon_SSS, baryon_SSS2j
	
	real*8 betaB1, betaB2, betaB3, betaB4
	
!-------------------------------------------
	
	character(len=32) ctemp
    character(len=7) CollSys

!-------------------------------------------
! test the precision of integral

	external Fiqc
	real*8 Fiqc
	real*8 q

!---------------------------------------------
! scan the parameter space
      	REAL*8 DATAX(500), DATAY(500)
	INTEGER IPOINT, IPAR1, IPAR2
	REAL(8) CHISQ


	open(unit=101, file='C:\Users\86177\Desktop\AuAu39_fit.dat', status='old')
	open(unit=102, file='p05_pb502_chisq.dat',status='unknown')

	IPOINT = 1
	
200	READ(101, *, END=210)DATAX(IPOINT), DATAy(IPOINT)
	PRINT *, IPOINT, DATAX(IPOINT), DATAY(IPOINT)
	IPOINT = IPOINT + 1
	GOTO 200
210	CONTINUE

	IPOINT = IPOINT-1

!--------------------------------------------------------	
!	call Set_Fiqc_Matrix
	
!	do i=1, 20
	
!	q = 0.1d0*i
!	print *, q, Fiqc(6, q, 0.1d0)
	
!	enddo

!------------------------------------------------
!------------------------------------------------
! We only need to change the parameters here

	c = 0.025D0        ! centrality
	Npart = 382.7D0   ! N participant from experimental data

	Isu = 3 ! 1 for pion, 2 for Kaon, 3 for proton, 4 for lambda^0, 5 for Xi, 6 for Phi, 7 for Omega
	IDparticle = Isu

    IYangZhu = 2 ! =1 use Yang's equations,  =2 use Zhu's equations
    ISys = 3 ! =1 for AuAu, =2 for PbPb at 2.76, =3 for PbPb at 5.02, =4 for XeXe at 5.44, =5 for dAu, =6 for pPb
    IResPi = 0 ! =1 with resonance contribution for pi; =0, no in meson_TT
    if(IYangZhu.eq.1) IResPi = 0 ! no resonance contribution in Yang's case

    if(ISys.eq.1)then
    CollSys = 'AuAu'
    elseif(ISys.eq.2)then
    CollSys = 'PbPb276'
    elseif(ISys.eq.3)then
    CollSys = 'PbPb502'
    elseif(ISys.eq.4)then
    CollSys = 'XeXe'
    elseif(ISys.eq.5)then
    CollSys = 'dAu'
    elseif(ISys.eq.6)then
    CollSys = 'pPb'
    elseif(ISys.eq.7)then
    else
    print *, 'Which system do you investigate?'
    call sleep(20)
    endif

    print *, 'The study system is: ', CollSys

	write(ctemp, '(F5.3)')c	
!------------------------------------------------
! input data, they are the common variables
! the thermal distribution temperature T and coefficient cT for u, d	
! the thermal distribution temperature Ts and coefficient cS for s

! parameters for Au-Au collisions at 200 GeV Lilin, PRC 88,044919(2013)
!	T = 0.31D0
!	cT = 23.2D0    ! c=0.05

!	Ts = 0.283D0
!	cS = 23.2D0    ! c=0.05

! parameters for Pb-Pb collisions at 2.76 TeV ! lilin Apr-8-2014
! T and T_s for 2.76 and 5.02 TeV are from Eqs. (4)-(7) of HZ2018, PRC97, 054902(2018)
!	T = 0.39D0
!	cT=23.2D0    ! c=0.025

!	Ts = 0.51D0
!	cS = 11.0D0    ! c=0.025

 ! parameters for Pb-Pb collisions at 5.02 TeV ! lilin May-16-2020
         T = 0.415D0
!         cT=23.2D0    ! c=0.025   ! We search the parameters cT and cS in v15
 
         Ts = 0.545D0
!         cS = 11.0D0    ! c=0.025
	
	gamma2j = 0.1D0 ! factor for two jets
    gammaYang = 0.07D0 ! the average xi factor in Yang's paper
    if(ISys.eq.5)gammaYang = 1.0D0 ! no energy loss in dAu case

! set up the hard jet momentum limits
    if(IYangZhu.eq.1)then
    plowlim = 3.0d0
    phighlim = 20.0d0
    else
    plowlim = 2.0d0
    phighlim = 30.0d0
    endif
!-------------------------------------------
! parameters for baryon	
	
	if(IDparticle.eq.3) then
	balp = 1.75D0               !common variables
	bbet = 1.05D0
	
	call Beta(balp+1.0D0, bbet+1.0D0, betaB1)
	call Beta(balp+1.0D0, balp+bbet+2.0D0, betaB2)
	call Beta(balp+2.0D0, bbet+2.0D0, betaB3)
	call Beta(balp+2.0D0, balp+bbet+4.0D0, betaB4)
	
	Nb = 1.0D0/betaB1/betaB2    !common variables
	Nbprime = betaB3*betaB4
	gbst = 1.0D0/6.0D0
	
	elseif(IDparticle.eq.4) then
	
	balp = 1.0D0               !common variables
	bbet = 2.0D0
	
	call Beta(balp+1.0D0, bbet+1.0D0, betaB1)
	call Beta(balp+1.0D0, balp+bbet+2.0D0, betaB2)
	
	Nb = 1.0D0/betaB1/betaB2    !common variables
	gbst = 1.0D0/4.0D0

!  The statistical factors for the three particles are from HZ2020, JPG47, 055102(2020)
	elseif(IDparticle.eq.5) then
	
	gbst = 0.03D0
	
	elseif(IDParticle.eq.6) then
	
	gbst = 0.423D0

    elseif(IDParticle.eq.7) then

    gbst = 0.01D0

	endif

!------------------------------------------------------
! gamma0  is PAR(1), q0 is PAR(2)

	PAR(1) = 4.5D0
	PAR(2) = 7.0D0

!------------------------------------------------------
! search the parameters by computer

	DO IPAR1 = 0, 10

	cT = 20.0d0 + REAL(IPAR1)*0.3d0

	DO IPAR2 = 0, 0

	cS = 7.0d0 + REAL(IPAR2)*0.4d0

!-------------------------------------------
! call this subroutine after we assign the value for c

	call Set_Fiqc_Matrix ! initialize MatFiqc	
	
!--------------------------------------------------
! open the data file to save output data
	if(Isu.eq.1)then
!------------------------------------------------	
! pion
	IDP1 = 1
	IDP2 = 4	
	print*, 'The is for pion case'
	open(unit = 11, file = "dNptdpt_"//CollSys//"_pion_c025.txt", status = "unknown")
	
	elseif(Isu.eq.2)then
	
!------------------------------------------------	
! kaon
	IDP1 = 1
	IDP2 = 7	
	print*, 'The is for kaon case'
	open(unit = 11, file = "dNptdpt_"//CollSys//"_kaon_c025.txt", status = "unknown")
	
	elseif(Isu.eq.3)then
	
!------------------------------------
! proton
	IDP1 = 1  ! (IDP1!=IDP2)
	IDP2 = 2
	IDP3 = 1	
	print*, 'The is for proton case'
	open(unit = 11, file = "dNptdpt_"//CollSys//"_p_c025.txt", status = "unknown")
	
	elseif(Isu.eq.4)then
	
!------------------------------------
! lambda^0
	IDP1 = 1 
	IDP2 = 2
	IDP3 = 6	
	print*, 'The is for lambda0 case'
	open(unit = 11, file = "dNptdpt_"//CollSys//"_lambda_c025.txt", status = "unknown")
	
	elseif(Isu.eq.5)then
	
!------------------------------------
! Xi
	IDP1 = 1 
	IDP2 = 6
	IDP3 = 6	
	print*, 'The is for Xi case'
	open(unit = 11, file = "dNptdpt_"//CollSys//"_Xi_c025.txt", status = "unknown")

	elseif(Isu.eq.6)then
!------------------------------------
! phi
	IDP1 = 6 
	IDP2 = 7
		
	print*, 'The is for Phi case'
	open(unit = 11, file = "dNptdpt_"//CollSys//"_phi_c025.txt", status = "unknown")

    elseif(Isu.eq.7)then
!------------------------------------
! Omega
    IDP1 = 6
    IDP2 = 6
    IDP3 = 6

    print*, 'The is for Omega case'
    open(unit = 11, file = "dNptdpt_"//CollSys//"_Omega_c025.txt", status = "unknown")

	endif

!---------------------------------------------	
	if(Isu.eq.1.or.Isu.eq.2.or.Isu.eq.6)then
	write(11, 108)
	elseif(Isu.eq.3.or.Isu.eq.4.or.Isu.eq.5.or.Isu.eq.7)then
	write(11, 112)
	endif

!---------------------------------------------	
! we need to modify here to choose the pt range

	CHISQ = 0.0

    do i = 1, IPOINT

	pt = DATAX(i)
	
!	pt = 0.225D0+0.3D0*(i-1) ! n=10 for proton to adjust the pt range
!    pt = 3.5D0+0.5D0*(i-1) ! n=7 for proton to adjust the pt range
!	pt = 7.5D0+1.0D0*(i-1) !n=17 for proton
!	pt = 13.5D0+1.0D0*(i-1) !n=7 for proton
!    pt = 0.11D0+0.2D0*(i-1) !n=40 for pion
!	pt = 12.0D0+2.0D0*(i-1) !n=9 for pion
!	pt = 0.65D0+1.0D0*(i-1) ! n=27 for Lambda
!	pt = 0.7D0+0.3D0*(i-1)	! n=24 for Xi
!	pt = 0.65D0+0.3D0*(i-1) ! n=15 for phi
!	pt= 1.3D0+0.5D0*(i-1)  !n=12 for Omega

	if(Isu.eq.1.or.Isu.eq.2.or.Isu.eq.6)then
	
	dNptdpt(1)  = meson_TT(pt, IDP1, IDP2, Isu) 
	dNptdpt(2)  = meson_TS(pt, IDP1, IDP2, Isu) 
	dNptdpt(3)  = meson_SS(pt, IDP1, IDP2, Isu) 
	dNptdpt(4)  = meson_SS2j(pt, IDP1, IDP2, Isu)
	dNptdpt(5)  = dNptdpt(1)+dNptdpt(2)+dNptdpt(3)+dNptdpt(4)

	CHISQ = CHISQ + (dNptdpt(5)/DATAY(I)-1.0)**2

    if(IYangZhu.eq.1)then
    dNptdpt(1)  = dNptdpt(1)
    dNptdpt(2)  = dNptdpt(2)*gammaYang
    dNptdpt(3)  = dNptdpt(3)*gammaYang
    dNptdpt(4)  = dNptdpt(4)*gammaYang**2
    dNptdpt(5)  = dNptdpt(1)+dNptdpt(2)+dNptdpt(3)+dNptdpt(4)
    endif

!	print *, pt, (dNptdpt(2:3))
	
!	write(11, 110) pt, dNptdpt(1:5)
	
!----------------------------------------------------------	
! for baryon

	elseif(Isu.eq.3.or.Isu.eq.4.or.Isu.eq.5.or.Isu.eq.7)then
	
	dNptdpt(1) = baryon_TTT(pt, IDP1, IDP2, IDP3, Isu)
	dNptdpt(2) = baryon_TTS(pt, IDP1, IDP2, IDP3, Isu)
	dNptdpt(3) = baryon_TSS(pt, IDP1, IDP2, IDP3, Isu)
	dNptdpt(4) = baryon_TSS2j(pt, IDP1, IDP2, IDP3, Isu)
	dNptdpt(5) = baryon_SSS(pt, IDP1, IDP2, IDP3, Isu)
	dNptdpt(6) = baryon_SSS2j(pt, IDP1, IDP2, IDP3, Isu)
	dNptdpt(7) = dNptdpt(1)+dNptdpt(2)+dNptdpt(3)+dNptdpt(4)+dNptdpt(5)+dNptdpt(6)

    CHISQ = CHISQ + (dNptdpt(7)/DATAY(I)-1.0)**2

    if(IYangZhu.eq.1)then
    dNptdpt(1)  = dNptdpt(1)
    dNptdpt(2)  = dNptdpt(2)*gammaYang
    dNptdpt(3)  = dNptdpt(3)*gammaYang
    dNptdpt(4)  = dNptdpt(4)*gammaYang**2
    dNptdpt(5)  = dNptdpt(5)*gammaYang
    dNptdpt(6)  = dNptdpt(6)*gammaYang**2
    dNptdpt(7) = dNptdpt(1)+dNptdpt(2)+dNptdpt(3)+dNptdpt(4)+dNptdpt(5)+dNptdpt(6)
    endif
!	print *, pt, dNptdpt(2),  dNptdpt(3)
	
	write(11, 114) pt, dNptdpt(1:7)
	else
	
	print *, 'Check the Isu in main program'
	stop
	
	endif
	
	enddo  ! enddo of pt

	PRINT *, PAR(1), PAR(2), CHISQ
	WRITE(102, *)PAR(1), PAR(2), CHISQ

	enddo  ! enddo of IPAR2

	ENDDO ! ENDDO OF IPAR1

108	format(10x, "%pt", 12x, "TT", 12x,  "TS", 12x,  "SS", 12x,  "SS2j", 12x, "Tot")
!110	format(5(1x, F14.9))	
110	format((1x, F12.3), 5(1x, ES14.5E2)) 

112	format(10x, "%pt", 11x, "TTT", 11x,  "TTS", 11x,  "TSS" , 11x,  "TSS2j",  11x,  "SSS", 11x,  "SSS2j", 11x, "Tot")
!114	format(6(1x, F14.9))
114	format((1x, F12.3), 7(1x, ES14.5E2))

		
	end

!==============================================
! Contributions from hard partons for one-shower parton

	function HardInt(IDP, p, centr)
    use limitsHJ
	implicit none
	
	real*8 HardInt, p, centr
	integer IDP
	
	external Intf_S1J
	real*8 Intf_S1J

	real*8 fa, fb, fdel   ! f means these variables for function
	real*8 fp1, fp2, fp3, fpt
	integer fIDP1, fIDP2, fIDP3
	real*8 fcentr, fsumf
	
	fa = plowlim
	fb = phighlim
	fdel = 0.02D0
	fp1 = p
	fp2 = 0.0D0
	fp3 = 0.0D0
	fpt = 0.0D0
	fIDP1 = IDP
	fIDP2 = 0
	fIDP3 = 0
	fcentr = centr
		
	if(fp1.ge.plowlim) then
	fa = fp1
	endif
	
	call Integral_1d(Intf_S1J, fa, fb, fdel, fp1, fp2, fp3, fpt, fIDP1, fIDP2, fIDP3, fcentr, fsumf)	
	
	HardInt = fsumf
	
	return
	
	end 
	
!--------------------------------------------------------------
! integral function for one shower from 1 jet distribution
! inputs x, p1, IDP1, centr
	
	function Intf_S1J(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	implicit none
	real*8 Intf_S1J
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
!------------------------------------------------------
	
	external Fiqc, QuenchFiqc, SPD, QuenchSPD
	real*8 Fiqc, QuenchFiqc, SPD, QuenchSPD

!------------------------------------------------------	
	
	integer IDH
	
	Intf_S1J = 0.0D0
	
	do IDH = 1, 7
	
	Intf_S1J = Intf_S1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x)*QuenchSPD(p1)/x
	
	enddo
	
	return
	
	end	
	
!=======================================================
! Contributions from hard partons for two-shower partons

	function HardInt2(IDP1, IDP2, p1, p2, centr)
    use limitsHJ
	implicit none
	
	real*8 HardInt2, p1, p2, centr
	integer IDP1, IDP2
	
	external Intf_SS1J
	real*8 Intf_SS1J

	real*8 fa, fb, fdel   ! f means these variables for function
	real*8 fp1, fp2, fp3, fpt
	integer fIDP1, fIDP2, fIDP3
	real*8 fcentr, fsumf
	
	fa = plowlim
	fb = phighlim
	fdel = 0.1D0
	fp1 = p1
	fp2 = p2
	fp3 = 0.0D0
	fpt = 0.0D0
	fIDP1 = IDP1
	fIDP2 = IDP2
	fIDP3 = 0
	fcentr = centr
		
	if(fp1+fp2.ge.plowlim) then
	fa = fp1+fp2
	endif
	
	call Integral_1d(Intf_SS1J, fa, fb, fdel, fp1, fp2, fp3, fpt, fIDP1, fIDP2, fIDP3, fcentr, fsumf)	
	
	HardInt2 = fsumf
	
	return
	
	end 
	
	
	
!--------------------------------------------------------------
! integral function for two showers from 1 jet distribution
! inputs x, p1, p2, IDP1, IDP2, centr
	
	function Intf_SS1J(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	implicit none
	real*8 Intf_SS1J
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
!------------------------------------------------------
	
	external Fiqc, QuenchFiqc, SPD, QuenchSPD
	real*8 Fiqc, QuenchFiqc, SPD, QuenchSPD
	
	external ll
	real*8 ll, xll
	
	external squ
	real*8 squ, xsqu

!------------------------------------------------------	
	
	integer IDH
	
	Intf_SS1J = 0.0D0
	
	if(IDP1.eq.IDP2)then ! two shower partons have same flavors in one jet
	
!----------------------------------
	
	if(IDP1.lt.6)then ! one shower parton is from L
		
	xll = p2/(x-p1)
	if(xll.gt.0.9999D0) xll = 0.9999D0  ! need to enforce this condition
	
	do IDH = 1, 7
	
	Intf_SS1J = Intf_SS1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x) &
	*QuenchSPD(p1)*ll(xll)*QuenchSPD(p2)/x
	
	enddo
	
	else ! one shower parton is from Ls
	
	xsqu = p2/(x-p1)
	if(xsqu.gt.0.9999D0) xsqu = 0.9999D0  ! need to enforce this condition
	
	do IDH = 1, 7
	
	Intf_SS1J = Intf_SS1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x) &
	*QuenchSPD(p1)*squ(xsqu)*QuenchSPD(p2)/x
	
	enddo
	
	
	endif
!------------------------------------
	
	else
	
	do IDH = 1, 7
	
	Intf_SS1J = Intf_SS1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x)*QuenchSPD(p1) &
	*SPD(IDH, IDP2, p2/(x-p1))*QuenchSPD(p2)/x
		
	enddo
	
	endif
	
	return
	
	end		

!=======================================================
! Contributions from hard partons for three-shower partons

	function HardInt3(IDP1, IDP2, IDP3, p1, p2, p3, centr)
    use limitsHJ
	implicit none
	
	real*8 HardInt3, p1, p2, p3, centr
	integer IDP1, IDP2, IDP3
	
	external Intf_SSS1J
	real*8 Intf_SSS1J

	real*8 fa, fb, fdel   ! f means these variables for function
	real*8 fp1, fp2, fp3, fpt
	integer fIDP1, fIDP2, fIDP3
	real*8 fcentr, fsumf
	
	fa = plowlim
	fb = phighlim
	fdel = 0.1D0
	fp1 = p1
	fp2 = p2
	fp3 = p3
	fpt = 0.0D0
	fIDP1 = IDP1
	fIDP2 = IDP2
	fIDP3 = IDP3
	fcentr = centr
		
	if(fp1+fp2+fp3.ge.plowlim) then
	fa = fp1+fp2+fp3
	endif
	
	call Integral_1d(Intf_SSS1J, fa, fb, fdel, fp1, fp2, fp3, fpt, fIDP1, fIDP2, fIDP3, fcentr, fsumf)	
	
	HardInt3 = fsumf
	
	return
	
	end 


!--------------------------------------------------------------
! integral function for three showers from 1 jet distribution
! inputs x, p1, p2, p3, IDP1, IDP2, IDP3, centr
	
	function Intf_SSS1J(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	implicit none
	real*8 Intf_SSS1J
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
!------------------------------------------------------
	
	external Fiqc, QuenchFiqc, SPD, QuenchSPD
	real*8 Fiqc, QuenchFiqc, SPD, QuenchSPD
	
	external ll
	real*8 ll, xll1, xll2
	
	external squ
	real*8 squ, xsqu1, xsqu2

!------------------------------------------------------	
	
	integer IDH
	
	Intf_SSS1J = 0.0D0
	
!--------------------------------

	if((IDP1.eq.IDP2).and.(IDP1.eq.IDP3))then

!+++++++++++++++++++++++++++++++++	
	
	if(IDP1.lt.6)then ! two shower parton are from L
	
	xll1 = p2/(x-p1)
	if(xll1.gt.0.9999D0) xll1 = 0.9999D0  ! need to enforce this condition
	xll2 = p3/(x-p1-p2)
	if(xll2.gt.0.9999D0) xll2 = 0.9999D0  ! need to enforce this condition
	
	do IDH = 1, 7
	
	Intf_SSS1J = Intf_SSS1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x) &
	*QuenchSPD(p1)*ll(xll1)*QuenchSPD(p2)*ll(xll2)*QuenchSPD(p3)/x
	
	enddo
	
	else ! two shower partons are from Ls
	
	xsqu1 = p2/(x-p1)
	if(xsqu1.gt.0.9999D0) xsqu1 = 0.9999D0  ! need to enforce this condition
	xsqu2 = p3/(x-p1-p2)
	if(xsqu2.gt.0.9999D0) xsqu2 = 0.9999D0  ! need to enforce this condition
	
	do IDH = 1, 7
	
	Intf_SSS1J = Intf_SSS1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x) &
	*QuenchSPD(p1)*squ(xsqu1)*QuenchSPD(p2)*squ(xsqu2)*QuenchSPD(p3)/x
	
	enddo
	
	endif

	
!+++++++++++++++++++++++++++++++++

	elseif((IDP1.eq.IDP2).or.(IDP1.eq.IDP3).or.(IDP2.eq.IDP3))then

!===================================

	if(IDP1.eq.IDP2)then
!IDP1=IDP2----------------------------------
	
	if(IDP1.lt.6)then ! one shower parton is from L
	
	xll1 = p2/(x-p1)
	if(xll1.gt.0.9999D0) xll1 = 0.9999D0  ! need to enforce this condition
	xll2 = p3/(x-p1-p2)
	if(xll2.gt.0.9999D0) xll2 = 0.9999D0  ! need to enforce this condition
	
	do IDH = 1, 7
	
	Intf_SSS1J = Intf_SSS1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x) &
	*QuenchSPD(p1)*ll(xll1)*QuenchSPD(p2)*SPD(IDH, IDP3, xll2)*QuenchSPD(p3)/x
	
	enddo
	
	else ! one shower parton are from Ls
	
	xsqu1 = p2/(x-p1)
	if(xsqu1.gt.0.9999D0) xsqu1 = 0.9999D0  ! need to enforce this condition
	xsqu2 = p3/(x-p1-p2)
	if(xsqu2.gt.0.9999D0) xsqu2 = 0.9999D0  ! need to enforce this condition
	
	do IDH = 1, 7
	
	Intf_SSS1J = Intf_SSS1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x) &
	*QuenchSPD(p1)*squ(xsqu1)*QuenchSPD(p2)*SPD(IDH, IDP3, xsqu2)*QuenchSPD(p3)/x
	
	enddo
	
	endif

!IDP1=IDP2----------------------------------	
	
	elseif(IDP1.eq.IDP3)then
	
!IDP1=IDP3----------------------------------

	if(IDP1.lt.6)then ! one shower parton is from L
	
	xll1 = p2/(x-p1)
	if(xll1.gt.0.9999D0) xll1 = 0.9999D0  ! need to enforce this condition
	xll2 = p3/(x-p1-p2)
	if(xll2.gt.0.9999D0) xll2 = 0.9999D0  ! need to enforce this condition
	
	do IDH = 1, 7
	
	Intf_SSS1J = Intf_SSS1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x) &
	*QuenchSPD(p1)*SPD(IDH, IDP2, xll1)*QuenchSPD(p2)*ll(xll2)*QuenchSPD(p3)/x
	
	enddo
	
	else ! one shower parton are from Ls
	
	xsqu1 = p2/(x-p1)
	if(xsqu1.gt.0.9999D0) xsqu1 = 0.9999D0  ! need to enforce this condition
	xsqu2 = p3/(x-p1-p2)
	if(xsqu2.gt.0.9999D0) xsqu2 = 0.9999D0  ! need to enforce this condition
	
	do IDH = 1, 7
	
	Intf_SSS1J = Intf_SSS1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x) &
	*QuenchSPD(p1)*SPD(IDH, IDP2, xsqu1)*QuenchSPD(p2)*squ(xsqu2)*QuenchSPD(p3)/x
	
	enddo
	
	endif
	
!IDP1=IDP3----------------------------------	
	else

!IDP2=IDP3----------------------------------
! this case is exactly same as IDP1=IDP3		
	
	if(IDP2.lt.6)then ! one shower parton is from L
	
	xll1 = p2/(x-p1)
	if(xll1.gt.0.9999D0) xll1 = 0.9999D0  ! need to enforce this condition
	xll2 = p3/(x-p1-p2)
	if(xll2.gt.0.9999D0) xll2 = 0.9999D0  ! need to enforce this condition
	
	do IDH = 1, 7
	
	Intf_SSS1J = Intf_SSS1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x) &
	*QuenchSPD(p1)*SPD(IDH, IDP2, xll1)*QuenchSPD(p2)*ll(xll2)*QuenchSPD(p3)/x
	
	enddo
	
	else ! one shower parton are from Ls
	
	xsqu1 = p2/(x-p1)
	if(xsqu1.gt.0.9999D0) xsqu1 = 0.9999D0  ! need to enforce this condition
	xsqu2 = p3/(x-p1-p2)
	if(xsqu2.gt.0.9999D0) xsqu2 = 0.9999D0  ! need to enforce this condition
	
	do IDH = 1, 7
	
	Intf_SSS1J = Intf_SSS1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x) &
	*QuenchSPD(p1)*SPD(IDH, IDP2, xsqu1)*QuenchSPD(p2)*squ(xsqu2)*QuenchSPD(p3)/x
	
	enddo
	
	endif
	
!IDP2=IDP3----------------------------------		
	endif
	
!====================================	

	else
	
	do IDH = 1, 7
	
	Intf_SSS1J = Intf_SSS1J + Fiqc(IDH, x, centr)*QuenchFiqc(x)*SPD(IDH, IDP1, p1/x)*QuenchSPD(p1) &
	*SPD(IDH, IDP2, p2/(x-p1))*QuenchSPD(p2)*SPD(IDH, IDP2, p3/(x-p1-p2))*QuenchSPD(p3)/x
		
	enddo
			
!----------------------------------	
	endif

!--------------------------------	

	return
	
	end

!==============================================================================
! one dimensional integral subroutine using expotential integral algorithm
! integral function fnx
! integral range from a to b, and the integral step is del
! momentum p1, p2, p3, pt
! parton ID IDP1, IDP2, IDP3
! centrality centr
! integral result sumf
	            	
	subroutine Integral_1d(fnx, a, b, del, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr, sumf)
	implicit none
!-------------------------------------------------	
	external fnx
	real*8 fnx
	
	real*8 a, b, del
	real*8 p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr, sumf
!-------------------------------------------------

	integer NP
	real*8 fa, fb, xi, aa
	integer i

!-------------------------------------------------	
! calculate how many steps we will take

	sumf = 0.0D0
	
	if((b-a)/del-int((b-a)/del).le.0.5D0)then
	NP = int((b-a)/del)
	else
	NP = int((b-a)/del)+1
	endif	
	
!--------------------------------------------------
	xi = a
	fa = fnx(xi, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	if(fa.le.0.0D0) fa = 1.0D-18
		
	do i = 1, NP
	
	xi = xi+del
	if(xi.gt.b) xi = b
	
	fb = fnx(xi, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	if(fb.le.0.0D0) fb = 1.0D-18
	
	aa = log(fb/fa)/del
		
	if(abs(aa).lt.0.0001D0)then
	sumf = sumf+fa*del
	else 
	sumf = sumf+(fb-fa)/aa
	
	endif
	
	fa = fb
	
	enddo
	
	return
	
	end	

!==============================================================================
! two dimensional integral subroutine using expotential integral algorithm
! this subroutine is specifically written for recombnation model to calculate the spectra of particles 
! the first integral range from fa to fb, and the integral step is fdel
! the second integral up limit is not fixed, it depends on the value of first integral variable
! integral function fnx
! momentum pt
! parton ID IDP1, IDP2, IDP3
! centrality centr
! integral result sumf
	            	
	subroutine Integral_2d(fnx, pt, IDP1, IDP2, IDP3, centr, sumf)
	implicit none
!-------------------------------------------------	
	external fnx
	real*8 fnx
	
	real*8 p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr, sumf
!-------------------------------------------------

	integer NP
	real*8 resa, resb, x, aa
	integer i
	
	real*8 fa, fb, fdel  ! f means first
	real*8 sa, sb, sdel  ! s means second
	real*8 fp1, fp2, fp3
	real*8 ssumf

!-------------------------------------------------	
! calculate how many steps we will take

	sumf = 0.0D0
	fa = 0.01D0
	fb = pt-0.01D0
	fdel = 0.05D0
	
	if((fb-fa)/fdel-int((fb-fa)/fdel).le.0.5D0)then
	NP = int((fb-fa)/fdel)
	else
	NP = int((fb-fa)/fdel)+1
	endif	
	
!--------------------------------------------------
	
	sa = 0.01D0
	sdel = fdel
	
	fp2 = 0.0D0
	fp3 = 0.0D0
	
	x = fa  ! x is the second integral variable
	sb = pt-x-0.01D0 ! the up limit of second integral variable which depends on value of the first integral variable 
	fp1 = x  
	
	call Integral_1d(fnx, sa, sb, sdel, fp1, fp2, fp3, pt, IDP1, IDP2, IDP3, centr, ssumf)
	
	resa = ssumf
	
	if(resa.le.0.0D0) resa = 1.0D-18
	
	do i =1, NP
	
	x = x + fdel
	sb = pt-x
	fp1 = x
	
	call Integral_1d(fnx, sa, sb, sdel, fp1, fp2, fp3, pt, IDP1, IDP2, IDP3, centr, ssumf)
	
	resb = ssumf
	
	if(resb.le.0.0D0) resb = 1.0D-18
	
	aa = log(resb/resa)/fdel
		
	if(abs(aa).lt.0.0001D0)then
	sumf = sumf+resa*fdel
	else 
	sumf = sumf+(resb-resa)/aa
	endif
	
	resa = resb

	enddo

	return
	
	end	
	

!=======================================================
! Isu = 1 for pion
! Isu = 2 for kaon

	subroutine meson_coe(am, bm, p0, pt, Isu)
	implicit none
	
	real*8 am, bm, p0, pt
	integer Isu
	
	if(Isu.eq.1)then
	
	am = 0.0D0
	bm = 0.0D0
	
!	p0 = sqrt(p**2+0.139D0**2)
	p0 = pt
	
	elseif(Isu.eq.2)then
	
	am = 1.0D0
	bm = 2.0D0
	
	p0 = sqrt(pt**2+0.494D0**2)
	
	elseif(Isu.eq.6)then
	
	am = 0.0D0
	bm = 0.0D0
	
	p0 = sqrt(pt**2+1.0194D0**2)
	
	else
	
	print *, 'Check the flag in meson_coe'
	
	endif
	
	end	

!=======================================================
! meson spectrum from TT

	function meson_TT(pt, IDP1, IDP2, Isu)
    use Iflag
	implicit none
	include 'Recom.h'
	
	real*8 meson_TT
	real*8 pt
	integer IDP1, IDP2, Isu
	
	external Intf_kaonTT
	real*8 Intf_kaonTT
	
	real*8 am, bm, p0
	real*8 betaM1, betaM2
	
	real*8 fa, fb, fdel   ! f means these variables for function
	real*8 fp1, fp2, fp3, fpt
	integer fIDP1, fIDP2, fIDP3
	real*8 fcentr, fsumf
	
	fa = 0.0D0
	fb = pt
	fdel = 0.02D0
	fp1 = 0.0D0
	fp2 = 0.0D0
	fp3 = 0.0D0
	fpt = pt
	fIDP1 = IDP1
	fIDP2 = IDP2
	fIDP3 = 0
	fcentr = c

	call meson_coe(am, bm, p0, pt, Isu) ! calculate am, bm, p0
!----------------------------------------------	
	if(Isu.eq.1)then ! for pion
	
	call Beta(am+2.0D0, bm+2.0D0, betaM1) ! calculate betaM1
	call Beta(am+1.0D0, bm+1.0D0, betaM2) ! calculate betaM2
	
	meson_TT = betaM1/betaM2*pt/p0*cT**2*dexp(-pt/T)

! including resonances Eq. (46) in arXiv: 1307.3328

    if(IResPi.eq.1)then
	
	meson_TT = (1.0D0+(2.8D0+0.003D0*Npart)*dexp(-pt/0.65D0))*meson_TT

    endif
	
!-------------------------------------------	
	elseif(Isu.eq.2)then ! for kaon
	
	call  Integral_1d(Intf_kaonTT, fa, fb, fdel, fp1, fp2, fp3, fpt, fIDP1, fIDP2, fIDP3, fcentr, fsumf)
		
	
    = 12.0*cT*cS/p0/pt**5*fsumf
	
	elseif(Isu.eq.6)then ! for phi
	
	meson_TT = gbst*cS**2*pt/4.0D0/p0*dexp(-pt/Ts)
	
	endif
	
	return

	end
	

!---------------------------------------------------
!---------------------------------------------------
! integral function for Kaon from TT
! inputs, x, pt

	function Intf_kaonTT(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	implicit none
	include 'Recom.h'
	real*8 Intf_kaonTT
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
!--------------------------------------------
	
	Intf_kaonTT = x*(pt-x)**2*x*dexp(-x/T)*(pt-x)*dexp(-(pt-x)/Ts)	
	
	return
	
	end

!=======================================================
! meson spectrum from TS

	function meson_TS(pt, IDP1,IDP2, Isu)	
	implicit none
	include 'Recom.h'
	
	real*8 meson_TS
	real*8 pt
	integer IDP1, IDP2, Isu

	external Intf_pionTS, Intf_kaonTS
	real*8 Intf_pionTS, Intf_kaonTS
	
	external f_PhiTS
	real*8 f_PhiTS
	
	real*8 am, bm, p0
	real*8 betaM
	
	real*8 fa, fb, fdel   ! f means these variables for function
	real*8 fp1, fp2, fp3, fpt
	integer fIDP1, fIDP2, fIDP3
	real*8 fcentr, fsumf
	
	fa = 0.01D0
	fb = pt-0.01D0
	fdel = 0.02D0
	fp1 = 0.0D0
	fp2 = 0.0D0
	fp3 = 0.0D0
	fpt = pt
	fIDP1 = IDP1
	fIDP2 = IDP2
	fIDP3 = 0
	fcentr = c
	
	call meson_coe(am, bm, p0, pt, Isu) ! calculate am, bm, p0
!------------------------------------------------------------	
	if(Isu.eq.1)then
	
	call Beta(am+1.0D0, bm+1.0D0, betaM)
	call Integral_1d(Intf_pionTS, fa, fb, fdel, fp1, fp2, fp3, fpt, fIDP1, fIDP2, fIDP3, fcentr, fsumf)

	meson_TS = 1.0D0/betaM/pt**(am+bm+2.0D0)/p0*fsumf
	
!------------------------------------------------

	elseif(Isu.eq.2)then !for kaon
	
	call Integral_1d(Intf_kaonTS, fa, fb, fdel, fp1, fp2, fp3, fpt, fIDP1, fIDP2, fIDP3, fcentr, fsumf)
	
	meson_TS = 12.0D0/p0/pt**5*fsumf
	
	elseif(Isu.eq.6)then ! for phi
	
	meson_TS = gbst*cS/2.0D0/p0*f_PhiTS(fpt, fIDP1, fIDP2, fIDP3, fcentr)
	
	endif
	
	return
	
	end

!---------------------------------------------------
! integral function for pion from TS
! inputs, x, pt, IDP1, IDP2, centr
	
	function Intf_pionTS(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	implicit none
	include 'Recom.h'
	real*8 Intf_pionTS
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt
	real*8 HardInt
	
!--------------------------------------------
	
	Intf_pionTS = x*cT*dexp(-x/T)*HardInt(IDP2, pt-x, centr)
	Intf_pionTS = Intf_pionTS + HardInt(IDP1, x, centr)*cT*(pt-x)*dexp(-(pt-x)/T)	 
	
	return
	
	end

!---------------------------------------------------
! integral function for Kaon from TS
! inputs, x, pt, IDP1, IDP2, centr
	
	function Intf_kaonTS(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	implicit none
	include 'Recom.h'
	real*8 Intf_kaonTS
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt
	real*8 HardInt
!--------------------------------------------
	
	Intf_kaonTS = x**2*(pt-x)**2*cT*dexp(-x/T)*HardInt(IDP2, pt-x, centr)
	Intf_kaonTS = Intf_kaonTS + x**2*(pt-x)**2*HardInt(IDP1, x, centr)*cS*(pt-x)/x*dexp(-(pt-x)/Ts)	 
		
	return
	
	end

!---------------------------------------------------
! function for Phi from TS
! inputs, pt, IDP1, IDP2, centr
	
	function f_PhiTS(pt, IDP1, IDP2, IDP3, centr)
	implicit none
	include 'Recom.h'
	real*8 f_PhiTS
	real*8 pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt
	real*8 HardInt
	
	real*8 p1, p2
!--------------------------------------------

	p1 = pt/2.0d0
	p2 = p1
	
	f_PhiTS = dexp(-p2/Ts)*HardInt(IDP1, p1, centr)
		 
	return
	
	end



!=======================================================
! meson spectrum from SS

	function meson_SS(pt, IDP1,IDP2, Isu)
    use limitsHJ
	implicit none
	include 'Recom.h'
	
	real*8 meson_SS
	real*8 pt
	integer IDP1, IDP2, Isu
		
	external Intf_mesonSS
	real*8 Intf_mesonSS	
			
	real*8 am, bm, p0
	real*8 betaM
	
	real*8 fa, fb, fdel   ! f means these variables for function
	real*8 fp1, fp2, fp3, fpt
	integer fIDP1, fIDP2, fIDP3
	real*8 fcentr, fsumf
	
	external HardInt2
	real*8 HardInt2
	
	fa = plowlim
	fb = phighlim
	fdel = 0.1D0
	fp1 = 0.0D0
	fp2 = 0.0D0
	fp3 = 0.0D0
	fpt = pt
	fIDP1 = IDP1
	fIDP2 = IDP2
	fIDP3 = 0
	fcentr = c
	
	if(pt.ge.plowlim) then
	fa = pt
	endif
		
!---------------------------------------------

	call meson_coe(am, bm, p0, pt, Isu) ! calculate am, bm, p0

!---------------------------------------------	

	if(Isu.eq.1.or.Isu.eq.2) then
	
	call Integral_1d(Intf_mesonSS, fa, fb, fdel, fp1, fp2, fp3, fpt, fIDP1, fIDP2, fIDP3, fcentr, fsumf)
	
	meson_SS = 1.0D0/p0*fsumf
	
	elseif(Isu.eq.6) then
	
	meson_SS = gbst/pt/p0*HardInt2(IDP1, IDP2, pt/2.0D0, pt/2.0D0, fcentr)
	
	else
	
	print *, 'check the Isu in meson_SS'
	stop
	
	endif
	
	return
	
	end

!---------------------------------------------------
! integral function for meson from SS
! inputs, x, pt, IDP1, IDP2, centr
	
	function Intf_mesonSS(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	implicit none
	include 'Recom.h'
	real*8 Intf_mesonSS
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external Fiqc, QuenchFiqc, FF, QuenchFF
	real*8 Fiqc, QuenchFiqc, FF, QuenchFF	
	
	integer IDH
!--------------------------------------------
	
	
	Intf_mesonSS = 0.0D0
	
	do IDH = 1, 7
	
	Intf_mesonSS = Intf_mesonSS + Fiqc(IDH, x, centr)*QuenchFiqc(x)*FF(IDH, IDparticle, pt/x)*QuenchFF(pt)/x**2/2.0D0
	
	enddo	 
	
	return
	
	end
	
!=======================================================
! meson spectrum from SS where SS from 2 jets

	function meson_SS2j(pt, IDP1,IDP2, Isu)	
	implicit none
	include 'Recom.h'
	
	real*8 meson_SS2j
	real*8 pt
	integer IDP1, IDP2, Isu
		
	external Intf_pionSS2j, Intf_kaonSS2j
	real*8 Intf_pionSS2j, Intf_kaonSS2j	
			
	real*8 am, bm, p0
	real*8 betaM
	
	real*8 fa, fb, fdel   ! f means these variables for function
	real*8 fp1, fp2, fp3, fpt
	integer fIDP1, fIDP2, fIDP3
	real*8 fcentr, fsumf
	
	external HardInt
	real*8 HardInt
	
	fa = 0.01D0
	fb = pt-0.01d0
	fdel = 0.02D0
	fp1 = 0.0D0
	fp2 = 0.0D0
	fp3 = 0.0D0
	fpt = pt
	fIDP1 = IDP1
	fIDP2 = IDP2
	fIDP3 = 0
	fcentr = c
			
!---------------------------------------------

	call meson_coe(am, bm, p0, pt, Isu) ! calculate am, bm, p0

!---------------------------------------------	
	
	if(Isu.eq.1)then
	
	call Integral_1d(Intf_pionSS2j, fa, fb, fdel, fp1, fp2, fp3, fpt, fIDP1, fIDP2, fIDP3, fcentr, fsumf)
	meson_SS2j = 1.0D0*gamma2j/pt**3*fsumf
	
	elseif(Isu.eq.2)then
	
	call Integral_1d(Intf_kaonSS2j, fa, fb, fdel, fp1, fp2, fp3, fpt, fIDP1, fIDP2, fIDP3, fcentr, fsumf)
	meson_SS2j = 12.0d0*gamma2j/p0/pt**5*fsumf
	
	elseif(Isu.eq.6)then
	
	meson_SS2j = gbst*gamma2j/pt/p0*HardInt(IDP1, pt/2.0D0, fcentr)*HardInt(IDP2, pt/2.0D0, fcentr)
	
	else
	
	print *, 'need to check the code in meson_SS2j'
	
	endif
	
	return
	
	end

!---------------------------------------------------
! integral function for pion from SS2j
! inputs, x, pt, IDP1, IDP2, centr
	
	function Intf_pionSS2j(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	implicit none
	include 'Recom.h'
	real*8 Intf_pionSS2j
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt
	real*8 HardInt	
	
!--------------------------------------------
	
	Intf_pionSS2j = HardInt(IDP1, x, centr)*HardInt(IDP2, pt-x, centr)
	
	return
	
	end

!---------------------------------------------------
! integral function for kaon from SS2j
! inputs, x, pt, IDP1, IDP2, centr
	
	function Intf_kaonSS2j(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	implicit none
	include 'Recom.h'
	real*8 Intf_kaonSS2j
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt
	real*8 HardInt	
	
!--------------------------------------------
	
	Intf_kaonSS2j = x*(pt-x)**2*HardInt(IDP1, x, centr)*HardInt(IDP2, pt-x, centr)
	
	return
	
	end	

!=============================================================
!=============================================================
! for baryon

	subroutine baryon_coe(p0, pt, Isu)
	implicit none
	
	real*8 p0, pt
	integer Isu
	
	if(Isu.eq.3) then
	
	p0 = dsqrt(pt**2+0.938D0**2)
	
	elseif(Isu.eq.4) then
	
	p0 = dsqrt(pt**2+1.1157D0**2)
	
	elseif(Isu.eq.5) then
	
	p0 = dsqrt(pt**2+1.3148D0**2)

    elseif(Isu.eq.7) then

    p0 = dsqrt(pt**2+1.672D0**2)

	else
	
	print *, 'check the baryon particle ID'
	
	endif
	
	end

!=======================================================
! baryon spectrum from TTT

	function baryon_TTT(pt, IDP1, IDP2, IDP3, Isu)
	implicit none
	include 'Recom.h'
	
	real*8 baryon_TTT
	real*8 pt
	integer IDP1, IDP2, IDP3
	integer Isu
	
	external Intf_lambdaTTT
	real*8 Intf_lambdaTTT
	
	real*8 p0, centr, sumf
	
	
	centr = c ! common variable
	call baryon_coe(p0, pt, Isu) ! calculate p0
	
!-------------------------------------------	
	if(Isu.eq.3) then
	
	baryon_TTT = gbst*Nb*Nbprime*cT**3*pt**2/p0*dexp(-pt/T)
	
	elseif(Isu.eq.4) then 
	
	call Integral_2d(Intf_lambdaTTT, pt, IDP1, IDP2, IDP3, centr, sumf)
	
	baryon_TTT = gbst*Nb*cT**2*cS/p0/pt**(2.0d0*balp+bbet+3.0d0)*sumf 
	
	elseif(Isu.eq.5) then
	
	baryon_TTT = gbst*cT*cS**2*pt**2/27.0d0/p0*dexp(-pt/3.0d0/T)*dexp(-2.0d0*pt/3.0d0/Ts)

!For Omega by Lilin
    elseif(Isu.eq.7) then

    baryon_TTT = gbst*cS**3*pt**2/27.0d0/p0*dexp(-pt/Ts)

	endif
!----------------------------------------------	
	
	return

	end
	
!----------------------------------------------------------
! inputs x, p1, pt

	function Intf_lambdaTTT(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 Intf_lambdaTTT
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
		
	real*8 delta_ptp1x
!----------------------------------------------------------	
	delta_ptp1x = pt -p1-x
	if(delta_ptp1x.le.0.0D0) delta_ptp1x = 1.0D-10
	
	Intf_lambdaTTT = (p1*x)**(balp+1.0d0)*dexp(-(p1+x)/T)*delta_ptp1x**(bbet+1.0d0)*dexp(-delta_ptp1x/Ts)
		
	return
	
	end			
	
!===========================================================
! baryon spectrum from TTS

	function baryon_TTS(pt, IDP1, IDP2, IDP3, Isu)
	implicit none
	include 'Recom.h'
	
	real*8 baryon_TTS
	real*8 pt
	integer IDP1, IDP2, IDP3
	integer Isu
	
	external Intf_protonTTS, Intf_lambdaTTS
	real*8 Intf_protonTTS, Intf_lambdaTTS
	
	external  f_XiTTS, f_OmegaTTS
	real*8  f_XiTTS,f_OmegaTTS
	
	real*8 p0, centr, sumf

!--------------------------------------------------------

	centr = c	
	call baryon_coe(p0, pt, Isu) ! calculate p0
	
	if(Isu.eq.3)then
	
	call Integral_2d(Intf_protonTTS, pt, IDP1, IDP2, IDP3, centr, sumf)
	
	baryon_TTS = gbst*Nb*cT**2/p0/pt**(2.0D0*balp+bbet+3.0D0)*sumf
	
	elseif(Isu.eq.4) then
		
	call Integral_2d(Intf_lambdaTTS, pt, IDP1, IDP2, IDP3, centr, sumf)
	
	baryon_TTS = gbst*Nb/p0/pt**(2.0d0*balp+bbet+3.0d0)*sumf 
	
	elseif(Isu.eq.5) then
	
	baryon_TTS = gbst/9.0d0/p0*pt*f_XiTTS(pt, IDP1, IDP2, IDP3, centr)

! For Omega by Lilin
    elseif(Isu.eq.7) then

    baryon_TTS = gbst/9.0d0/p0*pt*f_OmegaTTS(pt, IDP1, IDP2, IDP3, centr)

	endif
	
	end
!-----------------------------------------------------------
! input x, p1, pt, IDP1, IDP2

	function Intf_protonTTS(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 Intf_protonTTS
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt
	real*8 HardInt
	
	real*8 delta_ptp1x
!----------------------------------------------------------
	
	delta_ptp1x = pt-p1-x
	if(delta_ptp1x.le.0.0D0) delta_ptp1x = 1.0D-10
	
	Intf_protonTTS = dexp(-(p1+x)/T)*(p1*x)**(balp+1.0D0)*delta_ptp1x**bbet*HardInt(IDP2, delta_ptp1x, centr)
	Intf_protonTTS = Intf_protonTTS + dexp(-(p1+x)/T)*p1**(balp+1.0D0)*x**(bbet+1.0D0)*delta_ptp1x**balp &
	*HardInt(IDP1, delta_ptp1x, centr)
	
	return
	
	end

!-----------------------------------------------------------
! input x, p1, pt, IDP1, IDP3

	function Intf_lambdaTTS(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 Intf_lambdaTTS
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt
	real*8 HardInt
	
	real*8 delta_ptp1x
!----------------------------------------------------------
	
	delta_ptp1x = pt-p1-x
	if(delta_ptp1x.le.0.0D0) delta_ptp1x = 1.0D-10
	
	Intf_lambdaTTS = cT**2*p1*x*dexp(-(p1+x)/T)*(p1*x)**balp*delta_ptp1x**bbet*HardInt(IDP3, delta_ptp1x, centr)
	Intf_lambdaTTS = Intf_lambdaTTS + cT*cS*p1*dexp(-p1/T)*x*dexp(-x/Ts)*p1**balp*x**bbet*delta_ptp1x**balp &
	*HardInt(IDP1, delta_ptp1x, centr)
	
	return
	
	end	
	
!-----------------------------------------------------------
! input pt, IDP1, IDP2, centr

	function f_XiTTS(pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 f_XiTTS
	real*8 pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt
	real*8 HardInt
	
	real*8 p1, p2, p3
	
!----------------------------------------------------------

	p1 = pt/3.0d0
	p2 = p1
	
	f_XiTTS = cS**2*dexp(-2.0d0*p2/Ts)*HardInt(IDP1, p1, centr)
	f_XiTTS = f_XiTTS + cT*cS*dexp(-pt/3.0d0*(1.0d0/T+1.0d0/Ts))*HardInt(IDP2, p2, centr)
	
	return
	
	end	

!For Omega by Lilin -----------------------------------------------------------
! input pt, IDP1, centr

    function f_OmegaTTS(pt, IDP1, IDP2, IDP3, centr)
    implicit none
    include 'Recom.h'

    real*8 f_OmegaTTS
    real*8 pt
    integer IDP1, IDP2, IDP3
    real*8 centr

    external HardInt
    real*8 HardInt

    real*8 p1, p2, p3

!----------------------------------------------------------

    p1 = pt/3.0d0
    p2 = p1

    f_OmegaTTS = cS**2*dexp(-2.0d0*p2/Ts)*HardInt(IDP1, p1, centr)


return

end


!===========================================================
! baryon spectrum from TSS

	function baryon_TSS(pt, IDP1, IDP2, IDP3, Isu)
	implicit none
	include 'Recom.h'
	
	real*8 baryon_TSS
	real*8 pt
	integer IDP1, IDP2, IDP3
	integer Isu
	
	external Intf_protonTSS, Intf_lambdaTSS
	real*8 Intf_protonTSS, Intf_lambdaTSS
	
	external f_XiTSS, f_OmegaTSS
	real*8 f_XiTSS, f_OmegaTSS
	
	real*8 p0, centr, sumf

!------------------------------------------------	
	
	centr = c
	call baryon_coe(p0, pt, Isu) ! calculate p0
	
	if(Isu.eq.3)then
	
	call Integral_2d(Intf_protonTSS, pt, IDP1, IDP2, IDP3, centr, sumf)
	
	baryon_TSS = gbst*Nb*cT/p0/pt**(2.0D0*balp+bbet+3.0D0)*sumf
	
	elseif(Isu.eq.4) then
		
	call Integral_2d(Intf_lambdaTSS, pt, IDP1, IDP2, IDP3, centr, sumf)
	
	baryon_TSS = gbst*Nb/p0/pt**(2.0d0*balp+bbet+3.0d0)*sumf 


	elseif(Isu.eq.5) then
	
	baryon_TSS = gbst/3.0d0/p0*f_XiTSS(pt, IDP1, IDP2, IDP3, centr)

! For Omega by Lilin
    elseif(Isu.eq.7) then

    baryon_TSS = gbst/3.0d0/p0*f_OmegaTSS(pt, IDP1, IDP2, IDP3, centr)

	endif
	
	end
!-----------------------------------------------------------
! input x, p1, pt, IDP1, IDP2, IDP3

	function Intf_protonTSS(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 Intf_protonTSS
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt2
	real*8 HardInt2
	
	real*8 delta_ptp1x
!----------------------------------------------------------
	
	delta_ptp1x = pt-p1-x
	if(delta_ptp1x.le.0.0D0) delta_ptp1x = 1.0D-10
	
	Intf_protonTSS = dexp(-p1/T)*p1**(bbet+1.0)*x**balp*delta_ptp1x**balp*HardInt2(IDP1, IDP3, x, delta_ptp1x, centr)
	Intf_protonTSS = Intf_protonTSS+dexp(-p1/T)*p1*(p1*x)**(balp)*delta_ptp1x**bbet &
	*HardInt2(IDP1, IDP2, x, delta_ptp1x, centr)
	
	return
	
	end

!-----------------------------------------------------------
! input x, p1, pt, IDP1, IDP2, IDP3

	function Intf_lambdaTSS(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 Intf_lambdaTSS
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt2
	real*8 HardInt2
	
	real*8 delta_ptp1x
!----------------------------------------------------------
	
	delta_ptp1x = pt-p1-x
	if(delta_ptp1x.le.0.0D0) delta_ptp1x = 1.0D-10
	
	Intf_lambdaTSS = cS*p1*dexp(-p1/Ts)*p1**bbet*x**balp*delta_ptp1x**balp*HardInt2(IDP1, IDP2, x, delta_ptp1x, centr)
	Intf_lambdaTSS = Intf_lambdaTSS+cT*p1*dexp(-p1/T)*(p1*x)**balp*delta_ptp1x**bbet &
	*HardInt2(IDP2, IDP3, x, delta_ptp1x, centr)
	
	return
	
	end	
	
!-----------------------------------------------------------
! input pt, IDP1, IDP2, IDP3

	function f_XiTSS(pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 f_XiTSS
	real*8 pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt2
	real*8 HardInt2
	
	real*8 p1, p2, p3
!----------------------------------------------------------
	
	p1 = pt/3.0d0
	p2 = p1
	p3 = p1
	
	f_XiTSS = cT*dexp(-p1/T)*HardInt2(IDP2, IDP3, p2, p3, centr)+cS*dexp(-p2/Ts)*HardInt2(IDP1, IDP3, p1, p3, centr)
	
	return
	
	end

!For Omega by Lilin-----------------------------------------------------------
! input pt, IDP1, IDP2, IDP3

    function f_OmegaTSS(pt, IDP1, IDP2, IDP3, centr)
    implicit none
    include 'Recom.h'

    real*8 f_OmegaTSS
    real*8 pt
    integer IDP1, IDP2, IDP3
    real*8 centr

    external HardInt2
    real*8 HardInt2

    real*8 p1, p2, p3
!----------------------------------------------------------

    p1 = pt/3.0d0
    p2 = p1
    p3 = p1

    f_OmegaTSS = cS*dexp(-p1/Ts)*HardInt2(IDP2, IDP3, p2, p3, centr)

    return

    end
!===========================================================
! baryon spectrum from TSS which SS is from 2 jets

	function baryon_TSS2j(pt, IDP1, IDP2, IDP3, Isu)
	implicit none
	include 'Recom.h'
	
	real*8 baryon_TSS2j
	real*8 pt
	integer IDP1, IDP2, IDP3
	integer Isu
	
	external Intf_protonTSS2j, Intf_lambdaTSS2j
	real*8 Intf_protonTSS2j, Intf_lambdaTSS2j
	
	external f_XiTSS2j, f_OmegaTSS2j
	real*8 f_XiTSS2j, f_OmegaTSS2j
	
	real*8 p0, centr, sumf

!------------------------------------------------	
	
	centr = c
	call baryon_coe(p0, pt, Isu) ! calculate p0
	
	if(Isu.eq.3)then
	
	call Integral_2d(Intf_protonTSS2j, pt, IDP1, IDP2, IDP3, centr, sumf)
	
	baryon_TSS2j = gbst*Nb*gamma2j*cT/p0/pt**(2.0D0*balp+bbet+3.0D0)*sumf
	
	elseif(Isu.eq.4) then
		
	call Integral_2d(Intf_lambdaTSS2j, pt, IDP1, IDP2, IDP3, centr, sumf)
	
	baryon_TSS2j = gbst*Nb*gamma2j/p0/pt**(2.0d0*balp+bbet+3.0d0)*sumf
	
	elseif(Isu.eq.5) then
	
	baryon_TSS2j = gbst*gamma2j/3.0d0/p0*f_XiTSS2j(pt, IDP1, IDP2, IDP3, centr)

! For Omega by Lilin
    elseif(Isu.eq.7) then

    baryon_TSS2j = gbst*gamma2j/3.0d0/p0*f_OmegaTSS2j(pt, IDP1, IDP2, IDP3, centr)

	endif
	
	end
!-----------------------------------------------------------
! input x, p1, pt, IDP1, IDP2, IDP3

	function Intf_protonTSS2j(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 Intf_protonTSS2j
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt
	real*8 HardInt
	
	real*8 delta_ptp1x
!----------------------------------------------------------
	
	delta_ptp1x = pt-p1-x
	if(delta_ptp1x.le.0.0D0) delta_ptp1x = 1.0D-10
	
	Intf_protonTSS2j = dexp(-p1/T)*p1**(bbet+1.0)*x**balp*delta_ptp1x**balp &
	*HardInt(IDP1, x,centr)*HardInt(IDP3, delta_ptp1x,centr)
	Intf_protonTSS2j = Intf_protonTSS2j+dexp(-p1/T)*p1*(p1*x)**(balp)*delta_ptp1x**bbet &
	*HardInt(IDP1, x, centr)*HardInt(IDP2, delta_ptp1x,centr)
	
	return
	
	end	

!-----------------------------------------------------------
! input x, p1, pt, IDP1, IDP2, IDP3

	function Intf_lambdaTSS2j(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 Intf_lambdaTSS2j
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt
	real*8 HardInt
	
	real*8 delta_ptp1x
!----------------------------------------------------------
	
	delta_ptp1x = pt-p1-x
	if(delta_ptp1x.le.0.0D0) delta_ptp1x = 1.0D-10
	
	Intf_lambdaTSS2j = cS*p1*dexp(-p1/Ts)*p1**bbet*x**balp*delta_ptp1x**balp &
	*HardInt(IDP1, x,centr)*HardInt(IDP2, delta_ptp1x,centr)
	Intf_lambdaTSS2j = Intf_lambdaTSS2j+cT*p1*dexp(-p1/T)*(p1*x)**(balp)*delta_ptp1x**bbet &
	*HardInt(IDP2, x, centr)*HardInt(IDP3, delta_ptp1x,centr)
	
	return
	
	end	

!-----------------------------------------------------------
! input pt, IDP1, IDP2, IDP3

	function f_XiTSS2j(pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 f_XiTSS2j
	real*8  pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt
	real*8 HardInt
	
	real*8 p1, p2, p3
!----------------------------------------------------------
	
	p1 = pt/3.0d0
	p2 = p1
	p3 = p1
	
	f_XiTSS2j = cT*dexp(-p1/T)*HardInt(IDP2, p2,centr)*HardInt(IDP3, p3,centr)
	f_XiTSS2j = f_XiTSS2j + cS*dexp(-p2/Ts)*HardInt(IDP1, p1,centr)*HardInt(IDP3, p3,centr)
	
	return
	
	end	


!-----------------------------------------------------------
! input pt, IDP1, IDP2, IDP3

    function f_OmegaTSS2j(pt, IDP1, IDP2, IDP3, centr)
    implicit none
    include 'Recom.h'

    real*8 f_OmegaTSS2j
    real*8  pt
    integer IDP1, IDP2, IDP3
    real*8 centr

    external HardInt
    real*8 HardInt

    real*8 p1, p2, p3
!----------------------------------------------------------

    p1 = pt/3.0d0
    p2 = p1
    p3 = p1


    f_OmegaTSS2j = cS*dexp(-p1/Ts)*HardInt(IDP2, p2,centr)*HardInt(IDP3, p3,centr)

    return

    end

!=======================================================
! baryon spectrum from SSS

	function baryon_SSS(pt, IDP1,IDP2,IDP3, Isu)
    use limitsHJ
	implicit none
	include 'Recom.h'
	
	real*8 baryon_SSS
	real*8 pt
	integer IDP1, IDP2, IDP3
	integer Isu
		
	external Intf_protonSSS, Intf_lambdaSSS 
	real*8 Intf_protonSSS, Intf_lambdaSSS	
	
	external f_XiSSS, f_OmegaSSS
	real*8 f_XiSSS, f_OmegaSSS
			
	real*8 p0
	
	real*8 fa, fb, fdel   ! f means these variables for function
	real*8 fp1, fp2, fp3, fpt
	integer fIDP1, fIDP2, fIDP3
	real*8 fcentr, fsumf
	
	fa = plowlim
	fb = phighlim
	fdel = 0.1D0
	fp1 = 0.0D0
	fp2 = 0.0D0
	fp3 = 0.0D0
	fpt = pt
	fIDP1 = IDP1
	fIDP2 = IDP2
	fIDP3 = IDP3
	fcentr = c
	
	if(pt.ge.plowlim) then
	fa = pt
	endif
		
!---------------------------------------------
	
	call baryon_coe(p0, pt, Isu) ! calculate p0
	
	If(Isu.eq.3)then
	
	call Integral_1d(Intf_protonSSS, fa, fb, fdel, fp1, fp2, fp3, fpt, fIDP1, fIDP2, fIDP3, fcentr, fsumf)
	
	elseif(Isu.eq.4)then
	
	call Integral_1d(Intf_lambdaSSS, fa, fb, fdel, fp1, fp2, fp3, fpt, fIDP1, fIDP2, fIDP3, fcentr, fsumf)
	
	elseif(Isu.eq.5)then
	
	fsumf = f_XiSSS(pt, IDP1, IDP2, IDP3, fcentr)

! For Omega by Lilin
    elseif(Isu.eq.7)then

    fsumf = f_OmegaSSS(pt, IDP1, IDP2, IDP3, fcentr)

	endif
	
	baryon_SSS = 1.0D0/p0*fsumf
	
	return
	
	end

!---------------------------------------------------
! integral function for proton from SSS
! inputs, x, pt, centr
	
	function Intf_protonSSS(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	implicit none
	include 'Recom.h'
	real*8 Intf_protonSSS
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external Fiqc, QuenchFiqc, FF, QuenchFF
	real*8 Fiqc, QuenchFiqc, FF, QuenchFF	
	
	integer IDH
!--------------------------------------------
	
	
	Intf_protonSSS = 0.0D0
	
	do IDH = 1, 7
	
	Intf_protonSSS = Intf_protonSSS + Fiqc(IDH, x, centr)*QuenchFiqc(x)*FF(IDH, IDparticle, pt/x)*QuenchFF(pt)/x**2
	
	enddo	 
	
	return
	
	end
	
!---------------------------------------------------
! integral function for lambda from SSS
! inputs, x, pt, centr
	
	function Intf_lambdaSSS(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)
	implicit none
	include 'Recom.h'
	real*8 Intf_lambdaSSS
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external Fiqc, QuenchFiqc, FF_Wang, QuenchFF
	real*8 Fiqc, QuenchFiqc, FF_Wang, QuenchFF	
	
	integer IDH
!--------------------------------------------
	
	
	Intf_lambdaSSS = 0.0D0
	
	do IDH = 1, 7
	
	Intf_lambdaSSS = Intf_lambdaSSS + Fiqc(IDH, x, centr)*QuenchFiqc(x)*FF_Wang(IDH, IDparticle, pt/x)*QuenchFF(pt)/x**2
	
	enddo	 
	
	return
	
	end

!-----------------------------------------------------------
! input pt, IDP1, IDP2, IDP3

	function f_XiSSS(pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 f_XiSSS
	real*8  pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt3
	real*8 HardInt3
	
	real*8 p1, p2, p3
!----------------------------------------------------------
	
	p1 = pt/3.0d0
	p2 = p1
	p3 = p1
	
	f_XiSSS= gbst/pt*HardInt3(IDP1, IDP2, IDP3, p1, p2, p3, centr)
	
	return
	
	end

! For Omega by Lilin-----------------------------------------------------------
! input pt, IDP1, IDP2, IDP3

    function f_OmegaSSS(pt, IDP1, IDP2, IDP3, centr)
    implicit none
    include 'Recom.h'

    real*8 f_OmegaSSS
    real*8  pt
    integer IDP1, IDP2, IDP3
    real*8 centr

    external HardInt3
    real*8 HardInt3

    real*8 p1, p2, p3
!----------------------------------------------------------

    p1 = pt/3.0d0
    p2 = p1
    p3 = p1

    f_OmegaSSS= gbst/pt*HardInt3(IDP1, IDP2, IDP3, p1, p2, p3, centr)

    return

    end


!===========================================================
! baryon spectrum from SSS which S, SS are from 2 jets

	function baryon_SSS2j(pt, IDP1, IDP2, IDP3, Isu)
	implicit none
	include 'Recom.h'
	
	real*8 baryon_SSS2j
	real*8 pt
	integer IDP1, IDP2, IDP3
	integer Isu
	
	external Intf_protonSSS2j, Intf_lambdaSSS2j
	real*8 Intf_protonSSS2j, Intf_lambdaSSS2j
	
	external f_XiSSS2j, f_OmegaSSS2j
	real*8 f_XiSSS2j, f_OmegaSSS2j
	
	real*8 p0, centr, sumf
!------------------------------------------------	
	
	centr = c	
	call baryon_coe(p0, pt, Isu) ! calculate p0
	
	if(Isu.eq.3)then
	
	call Integral_2d(Intf_protonSSS2j, pt, IDP1, IDP2, IDP3, centr, sumf)
	
	baryon_SSS2j = gbst*Nb*gamma2j/p0/pt**(2.0D0*balp+bbet+3.0D0)*sumf
	
	elseif(Isu.eq.4) then
		
	call Integral_2d(Intf_lambdaSSS2j, pt, IDP1, IDP2, IDP3, centr, sumf)
	
	baryon_SSS2j = gbst*Nb*gamma2j/p0/pt**(2.0d0*balp+bbet+3.0d0)*sumf
	
	elseif(Isu.eq.5) then
	
	baryon_SSS2j = gbst*gamma2j/pt/p0*f_XiSSS2j(pt, IDP1, IDP2, IDP3, centr)

! For Omega by Lilin
    elseif(Isu.eq.7) then

    baryon_SSS2j = gbst*gamma2j/pt/p0*f_OmegaSSS2j(pt, IDP1, IDP2, IDP3, centr)

	endif	
	
	end
!-----------------------------------------------------------
! input x, p1, pt, IDP1, IDP2, IDP3

	function Intf_protonSSS2j(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 Intf_protonSSS2j
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt, HardInt2
	real*8 HardInt, HardInt2
	
	real*8 delta_ptp1x
!----------------------------------------------------------
	
	delta_ptp1x = pt-p1-x
	if(delta_ptp1x.le.0.0D0) delta_ptp1x = 1.0D-10
	
	Intf_protonSSS2j = p1**bbet*x**balp*delta_ptp1x**balp &
	*HardInt(IDP2, p1,centr)*HardInt2(IDP1, IDP3, x, delta_ptp1x, centr)
	Intf_protonSSS2j = Intf_protonSSS2j+(p1*x)**(balp)*delta_ptp1x**bbet &
	*HardInt(IDP1, p1, centr)*HardInt2(IDP2, IDP3, x, delta_ptp1x, centr)
	
	return
	
	end	

!-----------------------------------------------------------
! input x, p1, pt, IDP1, IDP2, IDP3

	function Intf_lambdaSSS2j(x, p1, p2, p3, pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 Intf_lambdaSSS2j
	real*8 x, p1, p2, p3, pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt, HardInt2
	real*8 HardInt, HardInt2
	
	real*8 delta_ptp1x
!----------------------------------------------------------
	
	delta_ptp1x = pt-p1-x
	if(delta_ptp1x.le.0.0D0) delta_ptp1x = 1.0D-10
	
	Intf_lambdaSSS2j = p1**bbet*x**balp*delta_ptp1x**balp &
	*HardInt(IDP3, p1,centr)*HardInt2(IDP1, IDP2, x, delta_ptp1x, centr)
	Intf_lambdaSSS2j = Intf_lambdaSSS2j+(p1*x)**(balp)*delta_ptp1x**bbet &
	*HardInt(IDP1, p1, centr)*HardInt2(IDP2, IDP3, x, delta_ptp1x, centr)
	
	return
	
	end			

!-----------------------------------------------------------
! input pt, IDP1, IDP2, IDP3

	function f_XiSSS2j(pt, IDP1, IDP2, IDP3, centr)	
	implicit none
	include 'Recom.h'
	
	real*8 f_XiSSS2j
	real*8 pt
	integer IDP1, IDP2, IDP3
	real*8 centr
	
	external HardInt, HardInt2
	real*8 HardInt, HardInt2
	
	real*8 p1, p2, p3
!----------------------------------------------------------
	
	p1 = pt/3.0d0
	p2 = p1
	p3 = p1
	
	f_XiSSS2j = HardInt(IDP1, p1, centr)*HardInt2(IDP2, IDP3, p2, p3, centr)
	f_XiSSS2j = f_XiSSS2j+ HardInt(IDP2, p2, centr)*HardInt2(IDP1, IDP3, p1, p3, centr)
	
	return
	
	end			

!For Omega by Lilin-----------------------------------------------------------
! input pt, IDP1, IDP2, IDP3

    function f_OmegaSSS2j(pt, IDP1, IDP2, IDP3, centr)
    implicit none
    include 'Recom.h'

    real*8 f_OmegaSSS2j
    real*8 pt
    integer IDP1, IDP2, IDP3
    real*8 centr

    external HardInt, HardInt2
    real*8 HardInt, HardInt2

    real*8 p1, p2, p3
!----------------------------------------------------------

    p1 = pt/3.0d0
    p2 = p1
    p3 = p1

    f_OmegaSSS2j = HardInt(IDP1, p1, centr)*HardInt2(IDP2, IDP3, p2, p3, centr)

    return

    end

!=========================================================================
! The quenching functions
!=========================================================================
! Eq. (C3) in ZH, PhysRevC.88.044919
! A cut-off for shower parton distribution

	function QuenchSPD(p1)
    use Iflag
	implicit none
    include 'Recom.h'
	real*8 QuenchSPD, p1

    if(IYangZhu.eq.2.and.ISys.eq.1)then
	QuenchSPD = 1.0D0 - dexp(-(p1/0.3D0)**2)  ! for AuAu
    elseif(IYangZhu.eq.2.and.Isys.eq.2)then
	QuenchSPD = 1.0D0 - dexp(-(p1/0.5D0)**2) ! for PbPb at 2.76 TeV
    elseif(IYangZhu.eq.2.and.Isys.eq.3)then
    QuenchSPD = 1.0D0 - dexp(-(p1/0.5D0)**2) ! for PbPb at 5.02 TeV Should be modified or not?
    elseif(IYangZhu.eq.1)then
    QuenchSPD = 1.0D0 ! In Yang's equations, there is no SPD quenching
    else
    print *, 'Check your values of IYangZhu and ISys'
    print *, 'In QuenchSPD'
    stop
    endif
	
	return
	
	end
	
!----------------------------------------------------------------------------
! Eq. (C4) in ZH, PhysRevC.88.044919
! This factor \gamma_1 (a cut-off on FF) in PhysRevC.88.044919 for SSS^1j and SSS^1j
	function QuenchFF(pt)
    use Iflag
	implicit none
    include 'Recom.h'
	real*8 QuenchFF, pt

    if(IYangZhu.eq.2)then
	QuenchFF = 1.0d0-dexp(-pt**2)
    else
    QuenchFF = 1.0D0 !In Yang's equations, there is no FF quenching
    endif
	
	return
	
	end
		
!------------------------------------------------------------------------
! Eq. (44) in PRC 84 064914 (2011)
! The parton distribution \bar F is not reliable for q<3, the cut-off is for minijet distribution
	function QuenchFiqc(q)
	implicit none
	real*8 QuenchFiqc, q
	
! we don't need this quench factor Apr-7-2014	
	QuenchFiqc = 1.0D0/(1.0D0 + dexp((3.5D0-q)/0.5D0))
	
! put quench factor equal 1	
	
!	QuenchFiqc = 1.0D0
	
	return
	
	end		

!===================================================
!===================================================
! This part we don't need to change	
!===================================================
!===================================================
! shower parton distribution functions which are 
! derived from fragmentation functions
! From HY2004, PRC 70, 024904
	function kns(x)  ! valence quark Kns
	implicit none
	
	real*8 kns, x
	
	real*8 nv, av, bv, cv, dv
	
	nv = 0.3328D0
	av = 0.45D0
	bv = 2.1D0
	cv = 5.0D0
	dv = 0.5D0
	
	kns = nv*x**av*(1.0D0-x)**bv*(1.0D0+cv*x**dv)
	
	return
	
	end 	
!---------------------------------------------------
! From HY2004, PRC 70, 024904

	function ll(x)  ! L
	implicit none
	
	real*8 ll, x
	
	real*8 ns, as, bs, cs, ds
	
	ns = 1.881D0
	as = 0.13323D0
	bs = 3.3842D0
	cs = -0.99098D0
	ds = 0.30986D0
	
	ll = ns*x**as*(1.0D0-x)**bs*(1.0D0+cs*x**ds)
	
	return
	
	end
	
!---------------------------------------------------
		
	function gg(x)  ! G
	implicit none
	
	real*8 gg, x
	
	real*8 ns, as, bs, cs, ds
! From HY2004, PRC70, 024904(2004)

!	ns = 0.8109D0
!	as = -0.0561D0
!	bs = 2.5470D0
!	cs = -0.1757D0
!	ds = 1.2D0

! From ZH2011, PRC84, 064914(2011)
! Z-5-18-SPD-ZH2011
!	ns = 0.739D0
!	as = -0.28D0
!	bs = 4.387D0
!	cs = 4.502D0
!	ds = 10.469D0

! From ZH2011, PRC84, 064914(2011)
! Z-5-26a-SPD-ZH2011
! The scale-breaking effect due to a cut-off at low p_1 is included.
	ns = 0.0739D0
	as = -0.6882D0
	bs = 3.5245D0
	cs = 16.4849D0
	ds = 0.7561D0
	
	gg = ns*x**as*(1.0D0-x)**bs*(1.0D0+cs*x**ds)
	
	return
	
	end
	
!-------------------------------------------------
! From HY2004, PRC 70, 024904

	function squ(x)  ! Ls
	implicit none
	
	real*8 squ, x
	
	real*8 ns, as, bs, cs, ds
	
	ns = 0.1182D0
	as = -0.1378D0
	bs = 2.3D0
	cs = 0.9D0
	ds = 0.1D0
	
	squ = ns*x**as*(1.0D0-x)**bs*(1.0D0+cs*x**ds)
	
	return
	
	end
	
!-------------------------------------------------
! From HY2004, PRC 70, 024904

	function sqg(x)  ! Gs
	implicit none
	
	real*8 sqg, x
	
	real*8 ns, as, bs, cs, ds
	
	ns = 0.0687D0
	as = -0.425D0
	bs = 2.489D0
	cs = -0.5D0
	ds = 1.1D0
	
	sqg = ns*x**as*(1.0D0-x)**bs*(1.0D0+cs*x**ds)
	
	return
	
	end
	
!-------------------------------------------------
! From HY2004, PRC 70, 024904

	function kk(x)  ! K
	implicit none
	
	real*8 kk, x
	
	external kns, ll
	real*8 kns, ll
	
	kk = kns(x)+ll(x)
	
	return
	
	end	
	
!===================================================
! shower parton distribution of different types
! Sij from IDH to IDP
	
	function SPD(IDH, IDP, x)
	implicit none
	
	real*8 SPD, x
	integer  IDH, IDP
	
	external kns, ll, gg, squ, sqg, kk
	real*8  kns, ll, gg, squ, sqg, kk

!-----------------------------------------
	
	if(x.ge.0.9999D0) x = 0.9999D0
	if(x.le.1.0D-6) x = 1.0D-6
	
	if(IDH.eq.IDP)then
	
	if(IDP<6)then
	SPD = kk(x)
	else
	SPD = kns(x)+squ(x)
	endif
	
!------------------------------------------	
	else
	
	if(IDH.lt.5)then
	
	if(IDP.lt.5) then
	SPD = ll(x)
	else
	SPD = squ(x)
	endif
	
	elseif(IDH.eq.5)then
	
	if(IDP.lt.5)then 
	SPD = gg(x)
	else
	SPD = sqg(x)
	endif
	
	else
	
	if(IDP.lt.5)then
	SPD = ll(x)
	else
	SPD = squ(x)
	endif

	endif
!--------------------------------------------	
	endif
	
	return
	
	end	
	
!=================================================================
! NLO parameterizations of fragmentation functions
! NPB 597, 337 (2001), KKP
! IDH
! 1 u, 2 d, 3 ubar, 4 dbar, 5 g, 6 s, 7 sbar
! IDparticle
! 1 pi+-, 2 k+-, 3 p+-

	function FF(IDH, IDparticle, x)
	implicit none
	
	real*8 FF
	integer IDH, IDparticle
	real*8 x
	
	real*8 mu2, mu02, Lambda5MS 
	real*8 ss
	
	real*8 N, alp, bet, gam
	
	mu2 = 100.0D0  ! GeV^2
	mu02 = 2.0D0   ! GeV^2
	Lambda5MS = 0.213D0**2  !GeV^2
	
	ss = log(log(mu2/Lambda5MS)/log(mu02/Lambda5MS))
	
	if(IDparticle.eq.1)then
!-------------------------------------------------------------------------	
	if(IDH.eq.1.or.IDH.eq.2.or.IDH.eq.3.or.IDH.eq.4)then   ! u, d->pi
	
	N   = 0.44809D0 - 0.13828D0*ss - 0.06951D0*ss**2 + 0.01354D0*ss**3
	alp = -1.47598D0 - 0.30498D0*ss - 0.01863D0*ss**2 - 0.12529D0*ss**3
	bet = 0.91338D0 + 0.64145D0*ss + 0.07270D0*ss**2 - 0.16989D0*ss**3
	gam =             0.07396D0*ss - 0.07757D0*ss**2
	
	elseif(IDH.eq.6.or.IDH.eq.7)then  ! s->pi
	
	N   = 16.5987D0 - 18.3856D0*ss + 2.44225D0*ss**2 + 2.13225D0*ss**3
	alp = 0.13345D0 + 0.22712D0*ss - 0.83625D0*ss**2 + 0.38526D0*ss**3
	bet = 5.89903D0 - 0.16911D0*ss + 0.59886D0*ss**2 - 0.25630D0*ss**3
	gam =           - 0.18619D0*ss + 0.87362D0*ss**2
	
	elseif(IDH.eq.5)then ! g->pi
	
	N   = 3.73331D0 - 3.16946D0*ss - 0.47683D0*ss**2 + 0.70270D0*ss**3
	alp = -0.74159D0 - 0.51377D0*ss - 0.19705D0*ss**2 - 0.17917D0*ss**3
	bet = 2.33092D0 + 2.03394D0*ss - 0.50764D0*ss**2 - 0.08565D0*ss**3
	gam =             0.09466D0*ss - 0.10222D0*ss**2
	
	else
	
	print *, 'Check the Hard jet ID in FF'
	
	endif
!--------------------------------------------------------------------------
	
	elseif(IDparticle.eq.2)then
	
	if(IDH.eq.1.or.IDH.eq.3.or.IDH.eq.6.or.IDH.eq.7)then   ! u, s->K
	
	N   = 0.17806D0 - 0.10988D0*ss - 0.02524D0*ss**2 + 0.03142D0*ss**3
	alp = -0.53733D0 - 0.60058D0*ss + 0.07863D0*ss**2 + 0.13276D0*ss**3
	bet = 0.75940D0 + 0.61356D0*ss - 0.43886D0*ss**2 + 0.23942D0*ss**3
	gam =             0.10742D0*ss + 0.12800D0*ss**2
	
	elseif(IDH.eq.2.or.IDH.eq.4)then  ! d->K
	
	N   = 4.96269D0 + 1.54098D0*ss - 9.06376D0*ss**2 + 4.94791D0*ss**3
	alp = 0.05562D0 + 1.88660D0*ss - 2.94350D0*ss**2 + 1.04227D0*ss**3
	bet = 2.79926D0 + 3.02991D0*ss - 4.14807D0*ss**2 + 1.91494D0*ss**3
	gam =             0.85450D0*ss - 0.61016D0*ss**2
	
	elseif(IDH.eq.5)then ! g->K
	
	N   = 0.23140D0 - 0.33644D0*ss + 0.16204D0*ss**2 - 0.02598D0*ss**3
	alp = -1.3640D0 + 0.97182D0*ss - 0.02908D0*ss**2 - 0.43195D0*ss**3
	bet = 1.79761D0 + 1.57116D0*ss + 0.71847D0*ss**2 - 0.68331D0*ss**3
	gam =             0.36906D0*ss + 2.39060D0*ss**2
	else
	
	print *, 'Check the Hard jet ID in FF'
	
	endif
!--------------------------------------------------------------------------
	
	elseif(IDparticle.eq.3)then	
	
	if(IDH.eq.1.or.IDH.eq.2)then   ! u, d-> p
	
	N   = 1.25946D0 - 1.17505D0*ss + 0.37550D0*ss**2 - 0.01416D0*ss**3
	alp = 0.07124D0 - 0.29533D0*ss - 0.24540D0*ss**2 + 0.16543D0*ss**3
	bet = 4.12795D0 + 0.98867D0*ss - 0.46846D0*ss**2 + 0.20750D0*ss**3
	gam =             0.18957D0*ss - 0.01116D0*ss**2
	
	elseif(IDH.eq.3.or.IDH.eq.4.or.IDH.eq.6.or.IDH.eq.7)then !  ubar, dbar, s->p
	
	N   = 4.01135D0 + 8.67124D0*ss - 22.7888D0*ss**2 + 11.4720D0*ss**3
	alp = 0.17258D0 + 4.57608D0*ss - 9.64835D0*ss**2 + 4.61792D0*ss**3
	bet = 5.20766D0 + 7.25144D0*ss - 12.6313D0*ss**2 + 6.07314D0*ss**3
	gam =             0.16931D0*ss - 0.09541D0*ss**2
	
	elseif(IDH.eq.5)then ! g->p
	
	N   = 1.56255D0 - 1.48158D0*ss - 0.39439D0*ss**2 + 0.51249D0*ss**3
	alp = 0.01567D0 - 2.16232D0*ss + 2.47127D0*ss**2 - 0.93259D0*ss**3
	bet = 3.57583D0 + 3.33958D0*ss - 3.05265D0*ss**2 + 1.21042D0*ss**3
	gam =           - 0.84816D0*ss + 1.23583D0*ss**2

	else
	
	print *, 'Check the Hard jet ID in FF'
	stop
	
	endif
	
!--------------------------------------------------------------------------
	endif
	
	if(x.eq.0.0D0) x = 1.0D-9
	
	FF = N*x**alp*(1.0D0-x)**bet*(1+gam/x)
	
	if(IDparticle.eq.3.and.(IDH.eq.2)) FF = FF/2.0D0
	
	return	
	
	end	
	
	

!=================================================================
!fragmentation function for Lambda
! PRC 58, 2321 (1998)
! IDH
! 1 u, 2 d, 3 ubar, 4 dbar, 5 g, 6 s, 7 sbar
! IDparticle
! 4 Lambda, lambdabar

	function FF_Wang(IDH, IDparticle, x)
	implicit none
	
	real*8 FF_Wang
	integer IDH, IDparticle
	real*8 x
	
	real*8 nhQ, ss
	real*8 Q2, Q02, LambdaQCD2
	
	real*8 N, alp, bet, gam
	real*8 N1, alp1, bet1, gam1
	real*8 N2, alp2, bet2
	real*8 coea, coeb, coec
	
	Q2 = 100.0D0  ! GeV^2 
	Q02 = 2.0D0   ! GeV^2
	LambdaQCD2 = 0.213D0**2  !GeV^2
	
	ss = log(log(Q2/LambdaQCD2)/log(Q02/LambdaQCD2))
	
	if(IDparticle.eq.4)then
!-------------------------------------------------------------------------	
	if(IDH.eq.1.or.IDH.eq.2.or.IDH.eq.6)then   ! u, d, s->Lambda
	
	coea = 0.0033D0
	coeb = 0.0232D0
	coec = 0.0265D0
	nhQ = coea+coeb*ss+coec*ss**2
	
	N = 0.318D0
	alp = -0.989D0
	bet = 4.956D0
	gam = 5.186D0
	
	N1 = 0.0D0
	alp1 = 0.0D0
	bet1 = 0.0D0
	gam1 = 0.0D0
	
	N2 = 0.0D0
	alp2 = 0.0D0
	bet2 = 0.0D0 
	
	elseif(IDH.eq.3.or.IDH.eq.4.or.IDH.eq.7)then  ! ubar, dbar, sbar->Lambda
	
	coea = 0.0033D0
	coeb = 0.0232D0
	coec = 0.0265D0
	nhQ = coea+coeb*ss+coec*ss**2
	
	N = 0.318D0
	alp = -0.989D0
	bet = 4.956D0
	gam = 5.186D0
	
	N1 = 0.0D0
	alp1 = 0.0D0
	bet1 = 0.0D0
	gam1 = 0.0D0
	
	N2 = 0.0D0
	alp2 = 0.0D0
	bet2 = 0.0D0 
		
	elseif(IDH.eq.5)then ! g->Lambda, Lambdabar
	
	coea = 0.0215D0
	coeb = 0.0454D0
	coec = 0.0568D0
	nhQ = coea+coeb*ss+coec*ss**2
	
	N = 3.378D0
	alp = -0.1660D0
	bet = 4.394D0
	gam = 0.105D0
	
	N1 = 0.0D0
	alp1 = 0.0D0
	bet1 = 0.0D0
	gam1 = 0.0D0
	
	N2 = 0.0D0
	alp2 = 0.0D0
	bet2 = 0.0D0 
		
	
	else
	
	print *, 'Check the Hard jet ID in FF_Wang'
	stop
	
	endif
	
!--------------------------------------------------------------------------
	endif
	
	if(x.eq.0.0D0) x = 1.0D-9
	
	FF_Wang = nhQ*N*x**alp*(1.0D0-x)**bet*(1.0D0+x)**gam
	FF_Wang = FF_Wang+nhQ*(N1*x**alp1*(1.0D0-x)**bet1*(1.0D0+x)**gam1+N2*x**alp2*(1.0D0-x)**bet2)
		
	return	
	
	end
			
	
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
        END


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
        END
        
!------------------------------------------------------------------------
	function Fiqc(ID, q, centr)
	implicit none
	include 'Recom.h'
	
	real*8 Fiqc
	integer ID
	real*8 q, centr      

	integer IQ
	
	IQ = floor((q+0.001d0)/0.01d0)
	
	if(IQ.eq.0) then
	print *, 'IQ = 0 in function Fiqc'
	stop
	endif
		
	centr = c
	
	Fiqc = MatFiqc(IQ, ID)
	
	return

	end
!------------------------------------------------------------------------
	subroutine Set_Fiqc_Matrix
	implicit none
	include 'Recom.h'
	integer ID
	integer IQ
	real*8 q
	
	external Cal_Fiqc
	real*8 Cal_Fiqc
	
	print *, 'Calculating MatFiqc'
	print *, 'centrality c is a common variable: ', c
	
	MatFiqc = 0.0d0
	
	do ID = 1, 7
	do IQ = 1, 3100  ! the maximum of pt 31 GeV
	
	q = 0.01d0*IQ
	
	MatFiqc(IQ, ID) = Cal_Fiqc(ID, q, c)
	
	enddo
	enddo
	
	print *, 'Finish calculating MatFiqc'

	end
!------------------------------------------------------------------------
	include 'vegas.f90' ! include the vegas integral subroutine
	
!------------------------------------------------------------------------

	function Cal_Fiqc(ID, q, centr)
    use Iflag
	implicit none
	
	real*8 Cal_Fiqc
	integer ID
	real*8 q, centr
	
	real*8 reg(4), pi, result, sdr, chi
     	integer ndim,ncall,itmx,nprn,idum
     	common /ranno/ idum  ! the common in vegas.f90
     	external Intf_fik     ! you must declare the fcn to be integrated as external
     	real*8 Intf_fik
        external fik
        real*8 fik

        if(IYangZhu.eq.1)then  ! In Yang's case, Fiqc = q**2 * fik
        Cal_Fiqc = q**2*fik(ID, q, centr)
        return
        endif

	pi = 3.14159265d0
	
	idum=9812371     ! a random seed for vegas
      	ndim=2           ! the number of dimensions of the integral
      	itmx=8           ! block the evaluations into itmx groups (for error estimation)
      	nprn=-1          ! something to do with information echoing

      	reg(1) = 0.d0    ! the integration region is defined in an
      	reg(2) = 0.d0    ! array of dim 2*d. In this case variable 1 i
       	reg(3) = pi/2.0d0   ! 0-pi, and variable 3 from 0-2pi.
      	reg(4) = 2.4d0

     	ncall = 2000  ! this is the total number of fcn calls
     	
      	call vegas(reg,ndim,Intf_fik,0,ncall,itmx,nprn,result,sdr,chi, ID, q, centr)
     	
     	Cal_Fiqc = result
     	
	return

	end 
	
!-----------------------------------------------------------------
	function Intf_fik(intvar, wgt, IDH, q, centr)
    use Iflag
	implicit none
	include 'Recom.h'
	real*8 Intf_fik
	real*8 intvar(2), wgt, q, centr
	integer IDH
	
	external Psi, fik
	real*8 Psi, fik 
	
	real*8 phi, z, pi, k
	real*8 gama, gamai
	real*8 xibari
	real*8 Tabc, TabcRatio
	
	pi = 3.14159265d0
	phi = intvar(1)
	z = intvar(2)
	
	gama = 0.11d0
	if(IDH.eq.5)then
!   For Au+Au, 200 GeV, Phys. Rev. C 88 044919
!	gamai = 0.14d0
!   For Pb+Pb, 2.76 TeV JPG47, 055102(2020)
!   gamai = 2.8d0/(1.0d0+q**2/49.0d0) ! consider the gammai depends on q
!	gamai = 3.2d0/(1.0d0+q**2/81.0d0) ! consider the gammai depends on q
	gamai = PAR(1)/(1.0d0+(q/PAR(2))**2) ! consider the gammai depends on q, search parameters by computer
	else
!   For Au+Au, 200 GeV, Phys. Rev. C 88 044919
!	gamai = 0.07d0
!   For Pb+Pb, 2.76 TeV JPG47, 055102(2020)
!   gamai = 2.8d0/(1.0d0+q**2/49.0d0) ! consider the gammai depends on q
!   For Pb+Pb, 5.02 TeV
!	gamai = 1.6d0/(1.0d0+q**2/81.0d0)
	gamai = PAR(1)/2.0d0/(1.0d0+(q/PAR(2))**2) ! consider the gammai depends on q, search parameters by computer
	endif
	
	if(.True.)then
	
	if(phi.gt.pi/2.0.and.phi.le.pi)then  !using the periodic condition
	phi = pi - phi
	elseif(phi.gt.pi.and.phi.le.1.5*pi)then
	phi = phi - pi
	elseif(phi.gt.1.5*pi.and.phi.le.2.0*pi)then
	phi = 2.0*pi - phi
	endif
	
	endif

!   ZH2013,Phys. Rev. C 88 044919
	xibari = gamai/gama*0.655d0*(1.0d0-centr-0.32*cos(phi)*sin(centr**(0.71d0)*pi)) ! depends on c
    if(IYangZhu.eq.1)then
    xibari = 0.0
    endif

	k = q*dexp(z*xibari)

!============================================
! T_AA(c) for different systems.
! Add by Lilin on May-20-2020
!============================================
    if(ISys.eq.5)then
!   For dAu, obtained by fitting T_AA vs centrality of dAu at 200 GeV
!---------------------------------------------------------------------------------------
!   From Hua's note: Recombination_model.pdf
!	Tabc = 0.45761d0-0.9933d0*centr+1.2564d0*centr**2-0.84044d0*centr**3 !need to modify for different systems
!   From LiLin by fit T_dAu from TBALE II PRC75,024909(2007)
!    Tabc = 0.4353d0-0.8382d0*centr+0.9157d0*centr**2-0.6008d0*centr**3
!    TabcRatio = Tabc/0.36d0
!---------------------------------------------------------------------------------------
    elseif(ISys.eq.1)then
!   For AuAu, obtained by fitting T_AA vs centrality of AuAu at 200 GeV
!   The data for T_AA vs centrality are from TABLE I, PRC69, 034909(2004)
!---------------------------------------------------------------------------------------
!    Tabc = 26.6381d0-97.3355d0*centr+122.0414d0*centr**2-51.9970d0*centr**3
!    TabcRatio = Tabc/24.28d0
!---------------------------------------------------------------------------------------
    elseif(ISys.eq.2.or.ISys.eq.3.)then
!   For PbPb, obtained by fitting T_AA vs centrality of PbPb
!   The data for T_AA vs centrality are from TABLE I, PRC88, 044909(2013)
!---------------------------------------------------------------------------------------
    Tabc = 27.3731d0-100.8426d0*centr+126.5849d0*centr**2-53.6066d0*centr**3
    TabcRatio = Tabc/24.93d0
!---------------------------------------------------------------------------------------
    endif
!	TabcRatio = 1.0D0 ! for central collision

	Intf_fik = Psi(z)*q**2*dexp(2.0d0*z*xibari)*fik(IDH, k, centr)*TabcRatio*2.0d0/pi
	
	end	
	
!-----------------------------------------------------------------	
! Psi(z) is very insensitive to phi and c, so we use this simplist expression for all systems and centralities.
! Appendix B, Phys. Rev. C 88 044919(2013), ZH2013
    function Psi(z)
	implicit none
	real*8 Psi, z

	Psi = (z/2.4d0)**0.5d0*(1.0d0-z/2.4d0)/0.64d0
	
	return
	
	end 
	
!-----------------------------------------------------------------	
! for different systems, we need to choose the right fik
! the minijet distribution 1/(2piptdptdy)

	function fik(IDH, k, centr)
    use Iflag
	implicit none
	real*8 fik
	integer IDH
	real*8 k, centr
	
	real*8 Ai, ki, ni

!============================================	
! for AuAu system 0-10% centrality
! D.K. Strivastava, C. Gale and R.J. Fries, PRC 67, 034903(2003)
!=============================================
	
	if(ISys.eq.1)then
	
!	print *, 'system: AuAu'
!-----------------------------	
	select case(IDH)
! u quark		
	case(1)
	
	Ai = 9.113d2
	ki = 1.459d0
	ni = 7.679d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! d quark		
	case(2)
	
	Ai = 9.596d2
	ki = 1.467d0
	ni = 7.662d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! ubar quark		
	case(3)
	
	Ai = 2.031d2
	ki = 1.767d0
	ni = 8.546d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! dbar quark		
	case(4)
	
	Ai = 2.013d2
	ki = 1.759d0
	ni = 8.566d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! gluon		
	case(5)
	
	Ai = 4.455d3
	ki = 1.7694
	ni = 8.610d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! s, sbar quark			
	case(6, 7)
	
	Ai = 1.038d2
	ki = 1.868d0
	ni = 8.642d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
	
	case default
	
	print *, "Need to check the parton IDH in fik"
	stop
	
	end select
!----------------------------	
	endif
	
!========================	
! for PbPb at 2.76 TeV system
! Logarithmic interpolation of the parameters between 200 and 5.5 TeV in PRC 67, 034903(2003).
!========================

	if(ISys.eq.2)then
	
!	print *, 'system: PbPb276'
!----------------------------
	select case(IDH)
! u quark	
	case(1)
	
	Ai = 1.138d4
	ki = 0.6868d0
	ni = 5.6736d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! d quark	
	case(2)
	
	Ai = 1.266d4
	ki = 0.6767d0
	ni = 5.6565d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! ubar quark
	case(3)
	
	Ai = 2.3956d3
	ki = 0.8724d0
	ni = 5.9734d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! dbar quark
	case(4)
	
	Ai = 2.2814d3
	ki = 0.8807d0
	ni = 5.9859d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! gluon
	case(5)
	
	Ai = 6.2d4
	ki = 0.98d0
	ni = 6.22d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! s and sbar quarks		
	case(6, 7)
	
	Ai = 933.358d0
	ki = 1.0536d0
	ni = 6.666d0
	
	fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
	
	case default
	
	print *, "Need to check the parton IDH in fik"
	stop
	
	end select
!--------------------------------	
	endif


!========================
! for PbPb at 5.02 TeV system
! Logarithmic interpolation of the parameters between 200 and 5.5 TeV in PRC 67, 034903(2003).
!========================

    if(ISys.eq.3)then

!    print *, 'system: PbPb502'
!----------------------------
    select case(IDH)
! u quark
    case(1)

    Ai = 2.0232d4
    ki = 0.5882d0
    ni = 5.3072d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! d quark
    case(2)

    Ai = 2.2790d4
    ki = 0.5774d0
    ni = 5.2902d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! ubar quark
    case(3)

    Ai = 4.2041d3
    ki = 0.7535d0
    ni = 5.5227d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! dbar quark
    case(4)

    Ai = 3.9673d3
    ki = 0.7625d0
    ni = 5.5339d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! gluon
    case(5)

    Ai = 1.1216d5
    ki = 0.8d0
    ni = 5.683d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! s and sbar quarks
    case(6, 7)

    Ai = 1.5397d3
    ki = 0.9329d0
    ni = 5.6333d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni

    case default

    print *, "Need to check the parton IDH in fik"
    stop

    end select
!--------------------------------
    endif

!========================
! for XeXe at 5.44 TeV system
! Logarithmic interpolation of the parameters between 200 and 5.5 TeV in PRC 67, 034903(2003).
!========================

    if(ISys.eq.4)then

!    print *, 'system: XeXe'
!----------------------------
    select case(IDH)
! u quark
    case(1)

    Ai = 2.1858d4
    ki = 0.5665d0
    ni = 5.2481d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! d quark
    case(2)

    Ai = 2.4663d4
    ki = 0.5552d0
    ni = 5.2311d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! ubar quark
    case(3)

    Ai = 4.5340d3
    ki = 0.7282d0
    ni = 5.4473d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! dbar quark
    case(4)

    Ai = 4.2734d3
    ki = 0.7377d0
    ni = 5.4583d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! gluon
    case(5)

    Ai = 1.2156d5
    ki = 0.775d0
    ni = 5.61d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! s and sbar quarks
    case(6, 7)

    Ai = 1.6468d3
    ki = 0.9096d0
    ni = 5.65582d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni

    case default

    print *, "Need to check the parton IDH in fik"
    stop

    end select
!--------------------------------
    endif


!============================================
! for dAu system 0-20% centrality
! R.C. Hwa and C.B. Yang, PRL 93, 082302(2004)
!=============================================
    
    if(ISys.eq.5)then
!-----------------------------
!    print *, 'system: dAu'

    select case(IDH)
! u quark
    case(1)

    Ai = 12.371d0
    ki = 1.440d0
    ni = 7.673d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! d quark
    case(2)

    Ai = 12.888d0
    ki = 1.439d0
    ni = 7.662d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! ubar quark
    case(3)

    Ai = 2.638d0
    ki = 1.768d0
    ni = 8.574d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! dbar quark
    case(4)

    Ai = 2.613d0
    ki = 1.766d0
    ni = 8.586d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
    ! gluon
    case(5)

    Ai = 63.116d0
    ki = 1.718d0
    ni = 8.592d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! s, sbar quark
    case(6, 7)

    Ai = 1.144d0
    ki = 1.935d0
    ni = 8.721d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni

    case default

    print *, "Need to check the parton IDH in fik"
    stop

    end select
!----------------------------
    endif

!========================
! for pPb system
! fit the data given by Hanzhong Zhang
!========================

    if(ISys.eq.6)then

!    print *, 'system: pPb'
!----------------------------
    select case(IDH)
! u quark
    case(1)

    Ai = 1.40393d0
    ki = 3.10961d0
    ni = 4.55291d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! d quark
    case(2)

    Ai = 1.3858d0
    ki = 3.19304d0
    ni = 4.61317d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! ubar quark
    case(3)

    Ai = 1.20505d0
    ki = 3.55904d0
    ni = 4.89068d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! dbar quark
    case(4)

    Ai = 1.23995d0
    ki = 3.51454d0
    ni = 4.86533d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! gluon
    case(5)

    Ai = 21.3404d0
    ki = 5.04642d0
    ni = 5.29d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! s and sbar quarks
    case(6, 7)

    Ai = 1.02512d0
    ki = 3.75027d0
    ni = 4.92425d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni

    case default

    print *, "Need to check the parton IDH in fik"
    stop

    end select
    endif
    !========================
! for AuAu39 system
! fit the data 2021/5/31
!========================

    if(ISys.eq.7)then

!    print *, 'system: AuAu39'
!----------------------------
    select case(IDH)
! u quark
    case(1)
        
    Ai = 1.89115d2
    ki = 1.90071d0
    ni = 8.88206d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! d quark
    case(2)

    Ai = 1.92446d2
    ki = 1.91823d0
    ni = 8.86506d0

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! ubar quark
    case(3)

    Ai = 4.36722d1
    ki = 2.28108d0
    ni = 1.00795d1

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! dbar quark
    case(4)

    Ai = 4.43760d1
    ki = 2.26444d0
    ni = 1.01040d1

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! gluon
    case(5)

    Ai = 8.67373d2
    ki = 2.26153d0
    ni = 1.00947d1

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni
! s and sbar quarks
    case(6, 7)

    Ai = 2.64301d1
    ki = 2.34232d0
    ni = 1.01681d1

    fik = 2.5d0*Ai/(1.0d0+k/ki)**ni

    case default

    print *, "Need to check the parton IDH in fik"
    stop

    end select
!--------------------------------
    endif
	
	return			
	
	end	


