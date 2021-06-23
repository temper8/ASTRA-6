C======================================================================|
	subroutine	NCTJ2
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include 'for/status.inc'
        integer	 j,jrw,jiw,jwork(21),lrho
        double	 precision  yt,ydt,yy,ydot,ywork(36),NUE,NUEE,COULG
        double	 precision  YGE,YGI,YCO(4,4),YNE,YNI,YTE,YTI,YH
	data	jiw/21/	jrw/36/  
	common	 /A_NCTJ2/ lrho,yge,ygi,yco
C----------------------------------------------------------------------|
 10	continue
	j  = 1
	include 'fml/nue'
        ydt = 5./NUE
	do	j=1,NA
	   lrho = j
	   yt = 0.
	   yy = ER(j)+1.d-3
	   call     markloc("in ESolver"//char(0))
 11	   continue
	   call	ESolver(yt,ydt,1,yy,jrw,ywork,jiw,jwork)
	   call	ncflux (1, yt, yy, ydot)
C Note! The limit for ydot should be consistent (at least by an 
C       order of magnitude larger than) with yatol in ESolver
	   if (abs(ydot) .gt. 1.d-1)	goto	11
	   ER(j) = yy
C----------------------------------------------------------------------|
C The rest of this subroutine is used exclusively for drawing. Exchange
C     with the Astra main can be made via C-arrays (CAR*) or WORK array.
C Beware that WORK arrays can be used also by other subroutines.
C     It is user's reponsibility to avoid a possible overlapping.
C The following can be omitted when not needed any more. 
C---
C Particle fluxes:
C Electrons (total)
	   CAR7(j) =-YGE*G11(j)
C Ions (total)
	   CAR8(j) =-YGI*ZMJ*G11(j)
C Relative discrepancy (G_e-G_i)/G_e
	   include 'fml/nue'
	   CAR6(j) = (YGE-YGI*ZMJ)/YGE
	   YH = HRO
	   if (j .eq. NA)	YH = HROA
	   YNE = (NE(j)+NE(j+1))*0.5
	   YNI = (NI(j)+NI(j+1))*0.5
	   YTE = (TE(j)+TE(j+1))*0.5
	   YTI = (TI(j)+TI(j+1))*0.5
C Partial fluxes:
C Diffusive (electrons)
	   CAR20(j) =-DN(J)*(NE(J+1)-NE(J))/YH*G11(j)
C Convective (electrons)
	   CAR21(j) = G11(j)*(CN(J)-
     1		(HN(j)*(TE(J+1)-TE(J))/YTE-
     2		 XN(j)*(TI(J+1)-TI(J))/YTI)/YH)*YNE
C Diffusive (ions)
	   CAR22(j) =-YCO(3,3)*(NI(J+1)-NI(J))/YH*G11(j)
C Convective (ions)
	   CAR23(j) = ((YCO(3,1)*(NE(J+1)-NE(J))/YNE+
     2			YCO(3,2)*(TE(J+1)-TE(J))/YTE+
     3			YCO(3,4)*(TI(J+1)-TI(J))/YTI)/YH-
     4	   1.d-3*ER(j)*(YCO(3,1)/YTE-YCO(3,3)/YTI*ZMJ))*YNI*G11(j)
C	   YGE = (DN(j)*(NE(J+1)-NE(J))+
C     1		 (HN(j)*(TE(J+1)-TE(J))/YTE+
C     2		  XN(j)*(TI(J+1)-TI(J))/YTI)*YNE)/YH-
C     3		  CN(J)*YNE
C	   YGI = (YCO(3,3)*(NI(J+1)-NI(J))+
C     1	         (YCO(3,1)*(NE(J+1)-NE(J))/YNE+
C     2		  YCO(3,2)*(TE(J+1)-TE(J))/YTE+
C     3		  YCO(3,4)*(TI(J+1)-TI(J))/YTI)*YNI)/YH-
C     4	         (YCO(3,1)/YTE1-YCO(3,3)*ZMJ/YTI1)*YER*YNI
	   WORK(j,1)  = YCO(1,1)	! DN
	   WORK(j,2)  = YCO(1,2)	! HN
	   WORK(j,3)  = YCO(1,3)
	   WORK(j,4)  = YCO(1,4)	! XN
	   WORK(j,5)  = YCO(2,1)	! DE
	   WORK(j,6)  = YCO(2,2)	! HE
	   WORK(j,7)  = YCO(2,3)
	   WORK(j,8)  = YCO(2,4)	! XE
	   WORK(j,9)  = YCO(3,1)
	   WORK(j,10) = YCO(3,2)
	   WORK(j,11) = YCO(3,3)
	   WORK(j,12) = YCO(3,4)
	   WORK(j,13) = YCO(4,1)	! DI
	   WORK(j,14) = YCO(4,2)	! HI
	   WORK(j,15) = YCO(4,3)
	   WORK(j,16) = YCO(4,4)	! XI
	enddo
	CAR6(NA1)  = CAR6(NA)
	CAR7(NA1)  = CAR7(NA)
	CAR8(NA1)  = CAR8(NA)
	CAR20(NA1) = CAR20(NA)
	CAR21(NA1) = CAR21(NA)
	CAR22(NA1) = CAR22(NA)
	CAR23(NA1) = CAR23(NA)
C----------------------------------------------------------------------|
C Edge values are not used by the code 
C but for drawing it is better to define them
	ER(NA1) = ER(NA)
	DE(NA1) = DE(NA)
	XE(NA1) = XE(NA)
	CE(NA1) = CE(NA)
	HE(NA1) = HE(NA)
	DI(NA1) = DI(NA)
	HI(NA1) = HI(NA)
	CI(NA1) = CI(NA)
	XI(NA1) = XI(NA)
	HN(NA1) = HN(NA)
	XN(NA1) = XN(NA)
	DN(NA1) = DN(NA)
	CN(NA1) = CN(NA)
	return
	end
C======================================================================|
	subroutine ncflux (neq, yt, yy, ydot)
C----------------------------------------------------------------------|
C     dy/dt = f(t,y) ,  or, in component form,
C     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include 'for/status.inc'
	double precision yt, yy(*), ydot(*), NUE, NUEE, COULG
	integer	j,neq,lrho
	double precision YCO(4,4),YGE,YGI,YH,YNE1,YNI1,YTE1,YTI1
	double precision YBT,YNE,YNI,YTE,YTI,YER,YEH,YET,YMU,YRO,YRT,YZI
	common	/A_NCTJ2/ lrho,yge,ygi,yco
C----------------------------------------------------------------------|
	call	markloc("in ncflux"//char(0))
	YBT = BTOR
	YZI = ZMJ
	YEH = 5.d-2
	YET = 5.d-2
	YRT = RTOR
	j = lrho
	YRO = j*HRO
	YER = YY(1)
	YNE = (NE(j)+NE(j+1))*0.5
	YNI = (NI(j)+NI(j+1))*0.5
	YTE = (TE(j)+TE(j+1))*0.5
	YTI = (TI(j)+TI(j+1))*0.5
	YMU = (MU(j)+MU(j+1))*0.5
	YRO = j*HRO
	YNE1 = YNE*1.d19
	YNI1 = YNI*1.d19
	YTE1 = YTE*1.d3
	YTI1 = YTI*1.d3
	call	markloc("in w7coeff"//char(0))
	call	W7COEFF
     +	       (YCO,YBT,YNE1,YNI1,YTE1,YTI1,YER,YEH,YET,YMU,YRO,YRT,YZI)
	call	markloc("in ncflux"//char(0))
	YH = HRO
	if (j .eq. NA)	YH = HROA
	DE(j) = YCO(2,1)
	HE(j) = YCO(2,2)
	XE(j) = YCO(2,4)
	CE(j) = YCO(2,1)*YER/YTE1-
     +		YCO(2,3)*(ZMJ*YER/YTI1+(NI(J+1)-NI(J))/YNI/YH)
	DI(j) = YCO(4,1)
	HI(j) = YCO(4,2)
	XI(j) = YCO(4,4)
	CI(j) = YCO(4,1)*YER/YTE1-
     +		YCO(4,3)*(ZMJ*YER/YTI1+(NI(J+1)-NI(J))/YNI/YH)
	DN(j) = YCO(1,1)
	HN(j) = YCO(1,2)
	XN(j) = YCO(1,4)
	CN(j) = YCO(1,1)*YER/YTE1-
     +		YCO(1,3)*(ZMJ*YER/YTI1+(NI(J+1)-NI(J))/YNI/YH)
C Test version:
C	DE(j) = 0.
C	XE(j) = 0.
C	CE(j) = 0.
C	HE(j) = 1.
C	DI(j) = 0.
C	HI(j) = 0.
C	CI(j) = 0.
C	XI(j) = 1.
C	HN(j) = 0.
C	XN(j) = 0.
	YGE = (DN(j)*(NE(J+1)-NE(J))+
     1	      (HN(j)*(TE(J+1)-TE(J))/YTE+
     2	       XN(j)*(TI(J+1)-TI(J))/YTI)*YNE)/YH-
     3	       CN(J)*YNE
	YGI = (YCO(3,3)*(NI(J+1)-NI(J))+
     1	      (YCO(3,1)*(NE(J+1)-NE(J))/YNE+
     2	       YCO(3,2)*(TE(J+1)-TE(J))/YTE+
     3	       YCO(3,4)*(TI(J+1)-TI(J))/YTI)*YNI)/YH-
     4	      (YCO(3,1)/YTE1-YCO(3,3)*ZMJ/YTI1)*YER*YNI
	include 'fml/nue'
	ydot(1) = (YGE-YGI*ZMJ)*NUE*NE(j)*NE(j)
	return
	end
C======================================================================|
	subroutine	ESolver(yt,ydt,neq,yy,jrw,ywork,jiw,jwork)
C----------------------------------------------------------------------|
	implicit none
	external ncflux,jacdummy
	integer  neq,jrw,jiw,jtol,jstat,jopt,jt,jtask,jwork(jiw)
	double   precision yy(*),time,yt,ytout,tau,ydt
	double   precision yrtol,yatol,ywork(jrw)
C----------------------------------------------------------------------|
	call	markloc("in lsode"//char(0))
	jtask = 1
	jtol  = 1
	jstat = 1
	jopt  = 0
	yrtol = 0.
	yatol = 1.d-6
	ytout = yt+ydt
	jt    = 10					! lsode
	call	lsode(ncflux,neq,yy,yt,ytout,jtol,yrtol,yatol,jtask,
     +			jstat,jopt,ywork,jrw,jwork,jiw,jacdummy,jt)
C	jt    = 2					! dlsoda
C	call	dlsoda(ncflux,neq,yy,yt,ytout,jtol,yrtol,yatol,jtask,
C     1            	jstat,jopt,ywork,jrw,jwork,jiw,jacdummy,jt)
C jstat = 2  if DLSODA was successful, negative otherwise.
C        -1 means excess work done on this call (perhaps wrong JT).
C        -2 means excess accuracy requested (tolerances too small).
C        -3 means illegal input detected (see printed message).
C        -4 means repeated error test failures (check all inputs).
C        -5 means repeated convergence failures (perhaps bad Jacobian
C           supplied or wrong choice of JT or tolerances).
C        -6 means error weight became zero during problem. (Solution
C           component i vanished, and ATOL or ATOL(i) = 0.)
C        -7 means work space insufficient to finish (see messages).
	if (jstat .eq. 2)	return
	write(*,*)"Time =",yt,"   LSODA Exit code =",jstat
	end
C======================================================================|
	subroutine jacdummy (NEQ, YT, YY, ML, MU, PD, NROWPD)
C----------------------------------------------------------------------|
	double precision YT, YY(*), PD(NROWPD,*)
	end
C======================================================================|
c programm to estimate approximated monoenergetic neoclassical coefficients.
      subroutine w7coeff(co,b,ene,eni,te,ti,efi,eh,et,aiota,r,r0,z)
c
c output:4x4 coeff matrix (similar to eq 60 of the manual)
c input:
c     b  : magnetic field (T)
c     ene: electron density (m-3)
c     eni: ion density (m-3)
c     te : electron temperature (eV)
c     ti : ion temperature (eV)
c     efi: Radial electric field (V/m)
c     eh : helical ripple
c     et : toroidal ripple
c     aiota: rot transform /(2*pi)=1/q
c     r  : minor radius (m)
c     r0 : major radius (m)
c     z  : atomic number
c
      implicit real*8 (a-h,o-z)
      double precision co(4,4)
c
      xe=1.
      xi=1.
      e=1.602e-19
      me=9.11e-31
      mi_me=1836.
      e_me=0.1758e12
      pi=acos (-1.)
      wce=2.*pi*28.*b*1.e9        ! e-cyclotron freq
      wci=wce/mi_me               ! i-cyclotron freq
      wE=abs(efi/(r*b))           ! drift frequency
      vde=xe*Te*aiota/(r0*b)
      vthe=sqrt(2.*e_me*te)
      vdi=xi*Ti*aiota/(r0*b*z)
      vthi=sqrt(2.*e_me*ti/mi_me)
      ele=24.-log(sqrt(ene*1.e-6)/te)
      anue=2.9e-12*z*ene*ele/te**1.5
      eli=ele
      anui=4.8e-14*eni*eli/ti**1.5
      fbl=sqrt(et+2.*eh)-sqrt(2.*et)

c Axisymmetric coefficients (automatically ambipolar: does not depend on efi)
c electrons
      dpse=1.4*vde*r0*anue/(wce*aiota**2)
      dpe=0.4*vde*vthe/(wce*aiota)
      dbe=vthe*r0*anue/(et**1.5*wce*aiota**2)
      dbpe=dbe*dpe/(dbe+dpe)
c ions
      dpsi=1.4*vdi*r0*anui/(wci*aiota**2)
      dpi=0.4*vdi*vthi/(wci*aiota)
      dbi=vthi*r0*anui/(et**1.5*wci*aiota**2)
      dbpi=dbi*dpi/(dbi+dpi)
c Axisymmetric coefficients (automatically ambipolar)
      dae=(dpse**1.5+dbpe**1.5)**0.6666666666
      dai=(dpsi**1.5+dbpi**1.5)**0.6666666666

c non-axisymmetric coeff. (dpends on electric field)
c electrons      
      d1nue=1.257*eh**1.5*vde**2/(pi*anue)
      dnu12e=0.628*vde**2*anue**0.5/(pi*wE**1.5)
      dnue=0.5*vde**2*anue/(wE**2*fbl)
c ions
      d1nui=1.257*eh**1.5*vdi**2/(pi*anui)
      dnu12i=0.628*vdi**2*anui**0.5/(pi*wE**1.5)
      dnui=0.5*vdi**2*anui/(wE**2*fbl)
c non-axisymmetric coeff.
      dhe=1./(1./d1nue+1./dnu12e+1./dnue)
      dhi=1./(1./d1nui+1./dnu12i+1./dnui)
c total monoenergetic coeff.
      dde=dhe+dae
      ddi=dhi+dai

c approximated coeffs:
      d1e=dde
      d2e=1.5*dde
      d3e=3.75*dde
      d1i=ddi
      d2i=1.5*ddi
      d3i=3.75*ddi
      
c neoclassical fluxes
c      gge=-de1*ene* (grne/ene+q*efi/te)-(d2e-1.5*d1e)*ene*grte/te
c      qqe=-de2*te*ene*(grne/ene+q*efi/te)-(d3e-1.5*d2e)*ene*grte
c change to ASTRA code system (eq. 60 of the manual)
      co(1,1)=d1e             ! Dn_e
      co(1,2)=(d2e - 1.5*d1e) ! De_e
      co(1,3)=0.
      co(1,4)=0.
      co(2,1)=d2e
      co(2,2)=(d3e - 1.5*d2e)   ! Chi_e
      co(2,3)=0.
      co(2,4)=0.
      co(3,1)=0.
      co(3,2)=0.
      co(3,3)=d1i             ! Dn_i
      co(3,4)=(d2i - 1.5*d1i) ! De_i
      co(4,1)=0.
      co(4,2)=0.
      co(4,3)=d2i
      co(4,4)=(d3i - 1.5*d2i) ! Chi_i
      return
      end
C======================================================================|
