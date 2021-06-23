C======================================================================|
	subroutine TOSCOPE(JSBP)
C----------------------------------------------------------------------|
C Call of this subroutine is created by the Astra compiler
C The C-function to_neut_ fills shared memory segment Nr. (JSBP+1)
C Parameter JSBP is the ordinal number of the (&:) request in a model
C     to the C-function to_neut_
C No variable transfer to Astra
C----------------------------------------------------------------------|
	implicit none
	integer JSBP,j
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision mu0
	data   mu0/1.25663706d0/
C----------------------------------------------------------------------|
	do	j=1,NA1
	   work1(j,1) = RHO(j)			! [m]
	   work1(j,2) = TE(j)			! [keV]
	   work1(j,3) = TI(j)			! [keV]
	   work1(j,4) = NE(j)*1.d19		! [1/m^3]
	   work1(j,5) = NI(j)*1.d19		! [1/m^3]
	   work1(j,6) = CC(j)			! [MS/m]
	   work1(j,7) = CD(j)			! [MA/m^2]
	   work1(j,8) = 1.6d-3*mu0*		! [MPa?]
     .	 	(NE(J)*TE(J)+NI(J)*TI(J)+0.5*NB2EQL*(PBLON(J)+PBPER(J)))
	enddo
C	write(*,*)" Before SCoPE"
C	write(*,'(1P,2(6E13.6/))')TIME,RHO(1),RHO(NA1),
C     &        TE(1),CD(NA1),work1d(5)
	call	to_scope(JSBP,NA1,NRD,TIME,work1)
        end
C======================================================================|
	subroutine OTSCOPE(JSBP,JSBR)
C----------------------------------------------------------------------|
C Call of this subroutine is created by the Astra compiler
C Parameter JSBP is the ordinal number of the (&:) request in a model
C Parameter JSBR is the ordinal number of the (:) request in a model
C No variable exchange with Astra except for CPUse
C----------------------------------------------------------------------|
	implicit none
	integer JSBP,JSBR,j,j5,jot,jn,jcall
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
	double precision SC_TIME
	save	jcall
	data	jcall/0/
C Warning! Be carefull with using "work1" here.
C          Other parallel processes can spoil it.
C          WORK1 has to be replaced by actual arguments.
C Presently SCoPE returns
C    work1d(1) - time
C    jn        - size of the SCoPE radial grid
C    work1d(2:jn+1) - rho grid
C    work1d(2+NB1:2+NB1+jn) - T_e 
C    work1d(2+2*NB1:2+2*NB1+jn) - T_i
C    work1d(2+3*NB1:2+3*NB1+jn) - n_e 
C    work1d(2+4*NB1:2+4*NB1+jn) - n_i 
C    work1d(2+5*NB1:2+5*NB1+jn) - CC 
C    work1d(2+6*NB1:2+6*NB1+jn) - CD 
C    work1d(2+7*NB1:2+7*NB1+jn) - pressure 
C    work1d(1:1+8*NB1)=work1d(1:jot) - total range used in this group
C    jot = 1+8*NB1
C    work1d(jot+1) = RTOR 
C    work1d(jot+2) = BTOR 
C    work1d(jot+3) = IPL
C    work1d(jot+4) = R_magn.ax. 
C    work1d(jot+5) = Phi_edge
C    work1d(jot+6:) - output arrays
C Remember:   WORK1D(1) coincides with WORK(1:1)
C----------------------------------------------------------------------|
	call	ot_scope(JSBP,CPTSBR(JSBR),work,jot,jn)
	SC_TIME = work1d(1)
C	write(*,*)(work1d(j+1),j=1,jn)
C	write(*,*)RHO(1),RHO(2),RHO(NA),RHO(NA1),ROC
	do j=1,jn
	   work1(j,1) = work1d(j+1)/work1d(jn+1)
	enddo
	do j=1,NA1
	   work1(j,2) = RHO(j)/ROC
	enddo
C	write(*,*)"ASTRA",jn
C	write(*,'(1P,6(E12.4))')(work1d(j+1),j=1,jn)
C	write(*,'(1P,6(E12.4))')(work1(j,1),j=1,jn)

C The rest is not used
C	j5 = 2+NB1					! T_e
C	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR23,work1(1,2))
C	j5 = j5+NB1					! T_i
C	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR24,work1(1,2))
CC      write(*,*)"SCoPE_Ti -> ASTRA", jn,NA1,TIME
CC      write(*,'(4(2F8.3,2X))')(work1d(j5+j),j=0,jn-1)
CC      write(*,'(4(2F8.3,2X))')(work1d(j5+j),j=0,jn-1)
CC      write(*,'(4(2F8.3,2X))')(CAR24(j),j=1,NA1,2)
C	j5 = j5+NB1					! n_e
C	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR25,work1(1,2))
C	j5 = j5+NB1					! n_i
C	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR26,work1(1,2))
C	j5 = j5+NB1					! T_i
C	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR27,work1(1,2))
C	j5 = j5+NB1					! T_i
C	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR28,work1(1,2))
C	j5 = j5+NB1					! pressure
C	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR29,work1(1,2))

	j5 = jot+6					! Psi
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR1,work1(1,2))
	j5 = j5+NB1					! j
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR2,work1(1,2))
C	write(*,'(1P,6(E12.4))')(work1d(j5+j),j=0,jn-1)
	j5 = j5+NB1					! q
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR3,work1(1,2))
C	write(*,'(1P,6(E12.4))')(work1d(j5+j),j=0,jn-1)
C	write(*,'(1P,6(E12.4))')(1./work1d(j5+j),j=0,jn-1)
	j5 = j5+NB1					! a
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR4,work1(1,2))
	j5 = j5+NB1					! Sh shift
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR5,work1(1,2))
	j5 = j5+NB1					! elongation
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR6,work1(1,2))
	j5 = j5+NB1					! triangularity
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR7,work1(1,2))
	j5 = j5+NB1					! up-down shift
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR8,work1(1,2))
	j5 = j5+NB1					! V'
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR9,work1(1,2))
	j5 = j5+NB1					! G11
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR10,work1(1,2))
	j5 = j5+NB1					! G22
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR11,work1(1,2))
C      write(*,*)jn
C      write(*,*)(work1d(j5+j),j=0,jn-1)
C      write(*,*)(CAR11(j),j=1,NA1)
	j5 = j5+NB1					! G33
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR12,work1(1,2))
	j5 = j5+NB1					! I
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR13,work1(1,2))
	j5 = j5+NB1					! S_lat
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR14,work1(1,2))
	j5 = j5+NB1					! <B/B_0>
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR15,work1(1,2))
	j5 = j5+NB1					! <(B/B_0)^2>
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR16,work1(1,2))
	j5 = j5+NB1					! <(B_0/B)^2>
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR17,work1(1,2))
	j5 = j5+NB1					! B_max
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR18,work1(1,2))
	j5 = j5+NB1					! B_min
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR19,work1(1,2))
	j5 = j5+NB1					! trapped fraction
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR20,work1(1,2))
	j5 = j5+NB1					! Z_eff
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR21,work1(1,2))
	j5 = j5+NB1					! a
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR22,work1(1,2))
	j5 = j5+NB1					! J
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR23,work1(1,2))
	j5 = j5+NB1					! particles
C	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR24,work1(1,2))
	j5 = j5+NB1					! sigma
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR25,work1(1,2))
	j5 = j5+NB1					! j_ext
	call	SMOOTH(1.d-3,jn,work1d(j5),work1,NA1,CAR26,work1(1,2))
        CV1 = work1d(jot+1)			! R_0
        CV2 = work1d(jot+2)			! B_0
        CV3 = work1d(jot+3)*1.d-6		! Ipl
        CV4 = work1d(jot+4)			! Rma
        CV5 = work1d(jot+5)			! Phi_edge
C        RTOR = work1d(1)			! R_0
C        BTOR = work1d(2)			! B_0
C        IPL = work1d(3)*1.d-6		! Ipl
C        CV5 = work1d(5)			! Phi_edge
	do j=NA1,1,-1
	   car1(j) = car1(1)-car1(j)		! psi
	   car2(j) = 1.d-6*car2(j)		! j     [MA/m^2]
	   car5(j) = SHIFT-car5(j)		! Shift
	   car25(j) = 1.d-6*car25(j)		! sigma [MS/m]
	   car26(j) = 1.d-6*car26(j)		! j_ext [MA/m^2]
	enddo
        CV6 = work1d(23*NB1+jot+6)		! Total number of species
        CV7 = work1d(23*NB1+jot+7)/9.1094d-31	! 1st species mass [Au]
        CV8 = work1d(23*NB1+jot+8)/1.6022d-19	! 1st species mass [Au]
        CV9 = work1d(23*NB1+jot+9)/9.1094d-31	! 2nd species mass [Au]
        CV10= work1d(23*NB1+jot+10)/1.6022d-19	! 2nd species mass [Au]
	CV11= SC_TIME

C	write(*,*)" After  SCoPE:  t = ", SC_TIME,jn
C	write(*,'(1P,2(6E13.6/))')TIME,RHO(1),RHO(NA1),
C     &        TE(1),CD(NA1),work1d(5)
C	write(*,'(2I5,5F10.3)')NB1,NA1,work(1,1),work(2,1),
C     &		1.d-6*work(3,1),work(4,1),work(5,1)
C	write(*,*)
        end
C======================================================================|
