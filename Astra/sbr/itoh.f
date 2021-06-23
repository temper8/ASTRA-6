C======================================================================|
	subroutine ITOH(YNIMP,YCHI,YCOR)
C----------------------------------------------------------------------|
C		Transport coefficients for Itoh model.
C		More commentaries in the original core.
C	  Interface to ASTRA (G.V.Pereverzev 09-Feb-1999)
C----------------------------------------------------------------------|
C Parameters:
C   YNIMP 	 number of impurity species (total number: NSM=YNIMP+2)
C   YCHI(1:NA1)  thermal diffusivities (2/3?) (currently only one)
C   YCOR(1:NA1)  a correction factor due to rotational shear
C 	Example call:
C   ITOH(CIMP1,CAR20,CAR10):;	XI=CAR20*CAR10+...;
C----------------------------------------------------------------------|
C WARNING!    
C	This implementation assumes
C-(1)
C     Z [PZ(NSM)], nuclear charge of each species (including electrons ?) 
C	is measured in proton charge units
C-(2)
C	RW(NRM,NFM) (Beam energy density of each species) is understood
C       as different beam species: H, D, T, ..., not E, E/2, E/3. 
C-(3)
C     Correspondence to existing ASTRA variables
C      Itoh:	       Astra:
C	NRM		NRD
C	NRMAX		NA1
C	RR		RTOR
C	RA		ABC
C	RMP		AMETR (?)
C	RAVG		RHO   (?)
C	BB		BTOR
C	ZEFF		ZEF
C-(4)
C     Current version provides 
C	the same AKDW (anomalous thermal diffusivity) for all species. 
C	Neoclassic coefficients are not implemented.
C	Therefore, only one output array is enabled.
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	 JSMAX,JSM,JFMAX,JFM,JFAIL,J,J1
	parameter  (JSMAX=10, JFMAX=3)
	double precision YNIMP,YCHI(*),YCOR(*),YEB,WE1,TAUAp
	double precision BETI,ALMHD,CS,G1
	real	PA(JSMAX),PZ(JSMAX),PTS(JSMAX),PNSS(JSMAX),RN(NRD,JSMAX)
	real	RT(NRD,JSMAX),RW(NRD,JFMAX),ZEFF(NRD),QP(NRD),BP(NRD)
	real	AKNC(NRD,JSMAX),AKDW(NRD,JSMAX),AK(NRD,JSMAX)
	real	RR,RA,BB,RAVG(NRD),RMP(NRD)
C----------------------------------------------------------------------|
	equivalence	(QP(1),WORK(1,1)),(BP(1),WORK(1,2)),
     >			(RN(1,1),WORK(1,3)),(RT(1,1),WORK(1,JSMAX+3)),
     >	    (RW(1,1),WORK(1,4*JSMAX+3)),(AKNC(1,1),WORK(1,2*JSMAX+3)),
     >	    (AK(1,1),WORK(1,5*JSMAX+3)),(AKDW(1,1),WORK(1,3*JSMAX+3))
C (NRD)x(6*JSMAX+JFMAX+2) elements of WORK used
C----------------------------------------------------------------------|
C Number of beam species:
	if(QNBI.lt.1.E-3)	then
	    JFM = 0
	else
	    JFM = 1
	endif

C Number of species:
	JSM = YNIMP+2.5
	if (JSM .gt. JSMAX)	then
	   write(*,*)">>> Itoh model >>> Too many impurity species"
	   return
	endif
C.    species order: 1) e, 2) dominant ion, 3) another ion, etc.
	PA(1) = 5.4463e-4
	PA(2) = AMJ
	PA(3) = AIM1
	PA(4) = AIM2
	PA(5) = AIM3
	PZ(1) = -1.
	PZ(2) = ZMJ
	PZ(3) = ZIM1(1)
	PZ(4) = ZIM2(1)
	PZ(5) = ZIM3(1)
	PTS(1) = TE(NA1)
	PTS(2) = TI(NA1)
	PTS(3) = TI(NA1)
C	PTS(4) = TI(NA1)
C	PTS(5) = TI(NA1)
	PNSS(1) = 0.1*NE(NA1)
	PNSS(2) = 0.1*NI(NA1)
	PNSS(3) = 0.1*NIZ1(NA1)
	YEB = 1.
C	YEB = EBEAM*(DBM1+0.5*DBM2+0.333*DBM3)/(DBM1+DBM2+DBM3)
	RR = RTOR
	RA = ABC
	BB = BTOR
	do	j=1,NA1
	    RN(j,1) = 0.1*NE(j)
	    RN(j,2) = 0.1*NI(j)
	    RN(j,3) = 0.1*NIZ1(j)
	    RT(j,1) = TE(j)
	    RT(j,2) = TI(j)
	    RT(j,3) = TI(j)
C  RW(NRM,NFM)  Beam energy density of each species (z.c.) (keV * 1.E20 /m3).
C  ** RW : 2/3 * (energy density of fast ions) measured by 10^20 m^-3 keV
	    RW(j,1) = .0667*NIBM(j)*YEB
	    QP(j) = 1./MU(j)
	    BP(j) = BTOR*j*HRO*MU(j)/RTOR
	    ZEFF(j) = ZEF(j)
	    RAVG(j) = RHO(j)
	    RMP(j) = AMETR(j)
	enddo
C----------------------------------------------------------------------|
	call	ITOH_ORIG(NRD,JSM,JFM,NA,RR,RA,RMP,RAVG,BB,
     >		PA,PZ,PTS,PNSS,RN,RT,RW,ZEFF,QP,BP,AKNC,AKDW,AK,JFAIL)
C----------------------------------------------------------------------|
C    Output: (2*YNIMP+4) arrays (currently only one!)
	do  10	j=1,NA
	   YCHI(j) = AKDW(j,1)
	   if (j .eq. NA)	goto	10
C The current version does not distinguish between species 
C		      and returns zeros for neoclassics.
C	    do	j1=1,JSM
C		YCHI(j+(j1-1)*NRD) = AKDW(j,j1)
C	    enddo
C	    do	j1=1,JSM
C		YCHI(j+(JSM+j1-1)*NRD) = AKNC(j,j1)
C	    enddo

C added 28-05-99 to include rotational shear according to S.-I.Itoh et al.,
C		  		Phys.Rev.Let., Vol.72, Nu.8, 1200 (1994)
C   				PPPC, Vol.38, 1743
!----------------------------------------------------------------------!
	   include 'fml/cs'
	   include 'fml/beti'
	   include 'fml/almhd'
	   TAUAp = RTOR/CS*sqrt(BETI*TE(j)/TI(j))
	   TAUAp = TAUAp*ROC/(j*HRO*MU(j))		! ROC <-> ABC ?
	   WE1 = (ER(j+1)-ER(j))/HRO
	   WE1 = TAUAp*WE1/SHEAR(j)/RHO(j)*ROC/BTOR	! ROC <-> ABC ?
	   YCOR(j) = 1./(1.+max(0.d0,G1(SHEAR(j),ALMHD))*WE1**2)
c	   YCOR(j) = 1./(1.+20.*WE1**2)
 10	continue
	YCHI(NA1) = YCHI(NA)
	end
C======================================================================|
	double precision FUNCTION G1(S,ALFA)
	double precision S,ALFA,SA,SA2
C----------------------------------------------------------------------|
	G1 = 0.
	SA = S-ALFA
	SA2 = SA*SA
	G1 = (0.375*(1-SA2)+0.5*(1+6*SA2)/ALFA)**2
	G1 = G1+0.375+0.5*(0.25*(1-5*SA2)-2*SA2/ALFA)**2
	G1 = 8*ALFA*(1-2*SA)*G1/(2+6*SA2/(1-2*SA))
	return
	end
C======================================================================|
C   A somewhat different physics based model is the Current Diffusive
C Ballooning Mode (CDBM) model [2.8.26].  This is based on a one point
C renormalization of pressure driven 'resistive MHD' turbulence but with the
C important difference that a self-consistent turbulent electron viscosity
C due to electron inertia replaces collisional resistivity in Ohm's Law and
C sustains the turbulent transport.  In this theory the turbulence has a
C sub-critical nature, which is supported by direct numerical simulations
C [2.8.27] and the transport is not particularly dependent on the linear
C instability criterion.  The model incorporates effects of a large Shafranov
C shift in the equilibrium and reflects favorable aspects of ideal MHD
C ballooning stability: reduced transport for low (or negative) and high
C magnetic shear and high pressure gradients; transport reductions due to
C sheared radial electric fields can also be included [2.8.28].  The theory
C involves one undetermined numerical coefficient which is chosen once and
C for all to optimize the fit to a dataset.  The model has captured
C satisfactorily the essential features of the Ohmic, L-mode, the internal
C transport barrier for the high bp mode of JT-60U [2.8.29] and current
C profile control by LHCD [2.8.30].
C 
C [2.8.26]	ITOH K, ITOH S-I, FUKUYAMA A, YAGI M, AZUMI M,
C Self-sustained turbulence and L-mode confinement in toroidal plasma I,
C Plasma Phys Control Fus 36 (1994) 279
C 
C [2.8.27]	YAGI M, ITOH S-I, ITOH K, FUKUYAMA A, AZUMI M,
C Self-sustained plasma turublence due to current diffusion, Physics of
C Plasmas 2 (1995), pp 4140-4148
C 
C [2.8.28]	ITOH S-I, ITOH K, FUKUYAMA A, YAGI M, Theory of anomalous
C transport in H-mode plasmas, Phys Rev Lett 72 (1994) pp 1200-1203
C 
C [2.8.29]	FUKUYAMA A, ITOH K, ITOH S-I, YAGI M, AZUMI M, Theory of
C improved confinement in high-bp tokamaks, Plasma Phys Contr Fus 36 (1994),
C pp1385-1390
C 
C [2.8.30]	FUKUYAMA A, ITOH K, ITOH S-I, YAGI M, AZUMI M, Transport
C simulation on L-mode and improved confinement associated with current
C profile modification, Plasma Phys Contr Fus 37 (1995) pp611-631
C======================================================================|
       SUBROUTINE ITOH_ORIG(NRM,NSM,NFM,NRMAX,RR,RA,RMP,RAVG,BB,
     1 PA,PZ,PTS,PNSS,RN,RT,RW,ZEFF,QP,BP,
     2 AKNC,AKDW,AK,KFAIL)
C  Calculate transport coefficients for Itoh model.  Modified from
C  subroutine supplied to ITER profile database by A. Fukuyama.
C
       IMPLICIT REAL*4 (A-F,H,O-Z)
C
C.    Begin index of input arguments:
CI NRM          Dimension of radial index (first index) of RN, RT, etc.
CI NSM          Number of species in RN, RT, etc.
CI NFM          Number of species in RW.
CI NRMAX        Index of last radial zone.
CR RR           Ro, major radius (m)
CR RA           a, minor radius (m)
CR RMP(NRM)     Midplane minor radius of zone center (m).
CR RAVG(NRM)    'Average' minor radius (m).
CR BB           Bo, toroidal magnetic field (T)
C.    species order: 1) e, 2) dominant ion, 3) another ion, etc.
CR PA(NSM)      A, atomic mass of each species (H=1, D=2, etc.).
CR PZ(NSM)      Z, nuclear charge of each species.
CR PTS(NSM)     Edge temperature of each species.
CR PNSS(NSM)    Edge density of each species (1.E20 /m3).
CR RN(NRM,NSM)  Density of each species at each zone center (1.E20 /m3).
CR RT(NRM,NSM)  Temperature of each species at each zone center (keV).
CR RW(NRM,NFM)  Beam energy density of each species (z.c.) (keV * 1.E20 /m3).
CR ZEFF(NRM)    Zeff at each zone center.
CR QP(NRM)      q, safety factor at outer boundary of each zone.
CR BP(NRM)      Poloidal magnetic field at outer boundary of each zone (T).
C.    Output arguments:
CR AKNC(NRM,NSM)  Neoclassical thermal diffusivity (m2/sec).
CR AKDW(NRM,NSM)  Turbulent thermal diffusivity (m2/sec).
CR AK(NRM,NSM)    Total thermal diffusivity (AKNC + AKDW) (m2/sec).
C.    Internal constants:
CR CNC          Coefficient of turbulent thermal diffusivity.
CR CK0          Coefficient of neoclassical thermal diffusivity.
CI KFAIL        Error flag: 0 ok; 1 NSM to big.
C.    End of variable index.
C.    Begin code generated by VNDX.
       INTEGER NRM,NSM,NFM,NRMAX
       REAL RR,RA,RMP(NRM),RAVG(NRM),BB,PA(NSM),PZ(NSM),PTS(NSM),
     1 PNSS(NSM),RN(NRM,NSM),RT(NRM,NSM),RW(NRM,NFM),BP(NRM),ZEFF(NRM),
     2 QP(NRM),AK(NRM,NSM),AKNC(NRM,NSM),AKDW(NRM,NSM),CNC,CK0
C.    End code generated by VNDX.
C
       PARAMETER (NSMAX=10) ! Maximum # species in local arrays:
       DIMENSION ZMASS(NSMAX),ZTEMP(NSMAX),ZVTH(NSMAX),ZTAU(NSMAX),
     1 ZRHO2(NSMAX),ZRNU(NSMAX),ZRK2(NSMAX)
C  ****** Local constants ******
       LOGICAL LOC_NEO
       DATA PI,AME,AMM,AEE,VC,AMYU0,AEPS0,RKEV/
     1 3.141592,9.10953E-31,1.67290E-27,1.60219E-19,2.99792E8,
     2 1.256637E-6,8.85419E-12,1.60219E-16/
       DATA CNC,CK0/ 1., 12./ LOC_NEO/.FALSE./
       DATA RK22,RA22,RB22,RC22/2.55,0.45,0.43,0.43/
       DATA RK2 ,RA2 ,RB2 ,RC2 /0.66,1.03,0.31,0.74/
C 1-April-96 Skip over the neoclassical calculation (no Bpol available)
C
C  ** RN : the density profile measured by 10^20 m^-3
C  ** RT : the temperature profile measured by keV
C  ** RW : 2/3 * (energy density of fast ions) measured by 10^20 m^-3 keV
C
C  ** PNSS : surface density [10^20 m-3]
C  ** PTS  : surface temerature [keV]
C
C  ***** DPP  is the derivative of total pressure (dnT/dr) *****
C  ***** ANE  is the electron density (ne) *****
C  ***** DQ   is the derivative of safety factor (dq/dr) *****
C
C       KFAIL=0
       IF(NSM .GT. NSMAX) THEN
C        Insufficient room in internal arrays so send a message:
         write(*,*)' Need longer arrays in ITOH'
         KFAIL=1
         DO JS=1,NSM
           DO JR=1,NRMAX
             AK(JR,JS)=0.
             AKDW(JR,JS)=0.
             AKNC(JR,JS)=0.
           ENDDO
         ENDDO
         RETURN
       ENDIF
C
       ZMASS(1)=AME
       DO JS=2,NSM
         ZMASS(JS)=PA(JS)*AMM
       ENDDO
C
       DO 100 NR=1,NRMAX
         IF(NR.EQ.NRMAX) THEN
           ANE=PNSS(1)
           DO JS=1,NSM
             ZTEMP(JS)=PTS(JS)
           ENDDO
C
           RPP=0.
           RPM=0.
           IF(NFM.GT.0) THEN
             DO JF=1,NFM
               RPM = RPM + RW(NR-1,JF)+ RW(NR,JF)
             ENDDO
           ENDIF
           DO JS=1,NSM
             RPP= RPP + PNSS(JS)*PTS(JS)
             RPM= RPM + RN(NR-1,JS)*RT(NR-1,JS)+RN(NR,JS)*RT(NR,JS)
           ENDDO
           RPM = 0.5*RPM
C
           ZEFFL=ZEFF(NR)
           DPP = (RPP-RPM)/(RMP(NR)-RMP(NR-1))
C
         ELSE
           ANE = 0.5*(RN(NR+1,1)+RN(NR,1))
           DO JS=1,NSM
             ZTEMP(JS)= 0.5*(RT(NR+1,JS)+RT(NR,JS))
           ENDDO
C
           RPP=0.
           RPM=0.
           IF(NFM.GT.0) THEN
             DO JF=1,NFM
               RPP = RPP + RW(NR+1,JF)
               RPM = RPM + RW(NR,JF)
             ENDDO
           ENDIF
           DO JS=1,NSM
             RPP= RPP + RN(NR+1,JS)*RT(NR+1,JS)
             RPM= RPM + RN(NR,JS)*RT(NR,JS)
           ENDDO
C
           ZEFFL  = 0.5*(ZEFF(NR+1)+ZEFF(NR))
           DPP = (RPP-RPM)/(RMP(NR+1)-RMP(NR))
         ENDIF
C
         IF(NR.EQ.1) THEN
           DR_Q=RAVG(NR+1)-RAVG(NR)
           DQ = (QP(NR+1)-QP(NR))/(1.5*DR_Q)
         ELSEIF(NR.EQ.NRMAX) THEN
           DR_Q=RAVG(NR)-RAVG(NR-1)
           DQ = (QP(NR)-QP(NR-1))/DR_Q
         ELSE
           DR_Q=RAVG(NR+1)-RAVG(NR-1)
           DQ = (QP(NR+1)-QP(NR-1))/(DR_Q)
         ENDIF
         QL = QP(NR)
         RL = RAVG(NR)
C
         EPS  = RL/RR
C
C        High-n ballooning mode section
C
         S=RL*DQ/QL
C  Use parentheses to avoid underflows in single precision.
         WPE2=(ANE*1.E20*AEE)*(AEE/AME)/AEPS0
         DELTA2=VC**2/WPE2
         AMD=ZMASS(2)
         VA=SQRT(BB**2/(AMYU0*(ANE*1.E20*AMD)))
         DBDR=DPP*(1.E20*RKEV)*RA/(BB**2/(2*AMYU0))
         ALFA=-QL*QL*DBDR*RR/RA
         RKCV=-EPS*(1.-1./(QL*QL))
         FS=TRCOFS(S,ALFA,RKCV)
C
         AKDWL=CK0*FS*SQRT(ABS(ALFA))**3*(DELTA2*VA)/(QL*RR)
         DO JS=1,NSM
           AKDW(NR,JS)=AKDWL
         ENDDO
C
         IF(.NOT.LOC_NEO) THEN
           CNC=0.
           GO TO 90
         ENDIF
C        ***** Neoclassical transport (Hinton, Hazeltine) *****
         EPSS=SQRT(EPS)**3
         DO JS=1,NSM
           ZRHO2(JS)=2.*(ZMASS(JS)/AEE)*ZTEMP(JS)*(RKEV/AEE)/
     1     (PZ(JS)*BP(NR))**2
C	write(*,*)NR,JS,ZMASS(JS),ZTEMP(JS),(RKEV/AEE),PZ(JS)*BP(NR)
         ENDDO
C
         COEF_NEW = 12.*PI*SQRT(PI)*SQRT(RKEV)*(SQRT(AMM)/AEE)*
     1   (AEPS0*SQRT(RKEV)/AEE)**2/(ANE*(1.E20*AEE)*ZEFFL*15.)
         ZTAU(1) = COEF_NEW*SQRT(PA(1))*(ZTEMP(1))**1.5/SQRT(2.)
         DO JS=2,NSM
           ZTAU(JS) = COEF_NEW*SQRT(PA(JS))*ZTEMP(JS)**1.5/PZ(JS)**2
         ENDDO
C        Moved this down from above to where it is needed:
         DO JS=1,NSM
           ZVTH(JS) = SQRT(ZTEMP(JS)*(RKEV/ZMASS(JS)))
         ENDDO
C
         DO JS=1,NSM
           ZRNU(JS) = ABS(QP(NR))*RR/(ZTAU(JS)*ZVTH(JS)*EPSS)
         ENDDO
         RNUE=ZRNU(1)
         ZRK2(1)=RK22*(1./(1.+RA22*SQRT(RNUE)+RB22*RNUE)
     1   +(EPSS*RC22)**2/RB22*RNUE/(1.+RC22*RNUE*EPSS))
         DO JS=2,NSM
           RNUI=ZRNU(JS)
           ZRK2(JS) =RK2 *(1./(1.+RA2 *SQRT(RNUI)+RB2 *RNUI)
     1     +(EPSS*RC2 )**2/RB2 *RNUI/(1.+RC2 *RNUI*EPSS))
         ENDDO
C
         DO JS=1,NSM
           AKNC(NR,JS) = SQRT(EPS)*ZRHO2(JS)*ZRK2(JS)/ZTAU(JS)
         ENDDO
C
   90  CONTINUE
         DO JS=1,NSM
           AK(NR,JS) = AKDW(NR,JS)+CNC*AKNC(NR,JS)
         ENDDO
  100  CONTINUE
C
       RETURN
       END
C/ MODULE TRCOFS
C
C.              %%%%%%%%%%%%%%%%%%%%%%%%%
C.              %                       %
C.              %       TRCOFS          %
C.              %                       %
C.              %%%%%%%%%%%%%%%%%%%%%%%%%
C
       FUNCTION TRCOFS(S,ALFA,RKCV)
C
       IMPLICIT REAL*4 (A-F,H,O-Z)
C
C added 15-02-00 to exclude division by zero shear (G.V.Pereverzev)
       S2 = S*S
       if (S2 .lt. 1.E-3)	S2 = 1.E-3

       IF(ALFA.GE.0.) THEN
         SA=S-ALFA
         IF(SA.GE.0.) THEN
           FS1=(1.+9.*SQRT(2.)*SA**2.5)
     1     /(SQRT(2.)*(1.-2.*SA+3.*SA*SA+2.*SA*SA*SA))
         ELSE
           FS1=1./SQRT(2.*(1.-2.*SA)
     1     *(1.-2.*SA+3.*SA*SA))
         ENDIF
         IF(RKCV.GT.0.) THEN
           FS2=SQRT(RKCV)**3/(S2)
         ELSE
           FS2=0.
         ENDIF
       ELSE
         SA=ALFA-S
         IF(SA.GE.0.) THEN
           FS1=(1.+9.*SQRT(2.)*SA**2.5)
     1     /(SQRT(2.)*(1.-2.*SA+3.*SA*SA+2.*SA*SA*SA))
         ELSE
           FS1=1./SQRT(2.*(1.-2.*SA)
     1     *(1.-2.*SA+3.*SA*SA))
         ENDIF
         IF(RKCV.LT.0.) THEN
           FS2=SQRT(-RKCV)**3/(S2)
         ELSE
           FS2=0.
         ENDIF
       ENDIF
       TRCOFS=MAX(FS1,FS2)
       RETURN
       END
