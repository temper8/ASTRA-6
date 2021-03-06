C======================================================================|
	subroutine TOGLF23P(JS,JE,JSBP)
C----------------------------------------------------------------------|
C Call of this subroutine is created by the Astra compiler
C The C-function to_neut_ fills shared memory segment Nr. (JSBP+1)
C Parameter JSBP is the ordinal number of the (&:) request in a model.
C           JSBP is generated by the Astra compiler.
C----------------------------------------------------------------------|
	implicit none
	integer JS,JE,JSBP,j
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
C	double precision ALMHD,CS,ROTSH,YAL
C	allocatable YAL(:,:)
C----------------------------------------------------------------------|
C	allocate (YAL(NA1,3))
C	do j=1,NA1
C	   include 'fml/almhd'
C	   YAL(j,1) = ALMHD
C	   include 'fml/cs'
C	   include 'fml/rotsh'
C	   YAL(j,2) = CS
C	   YAL(j,3) = ROTSH
C	enddo
	call	to_glf23p(JS,JE,NA1,ER,NN,NIBM,SHEAR,JSBP)
C	deallocate (YAL)
	end
C======================================================================|
	subroutine OTGLF23P(JS,JE,JSBP,JSBR)
C----------------------------------------------------------------------|
C Call of this subroutine is created by the Astra compiler
C Example call from a model:
C   glf23p(1,25)&<:; 	! the call is placed     in detvar.tmp
C   glf23p(26,NA1)&<:;	! the response is taken  in init.tmp @ eqns.tmp
C	less economical:
C   glf23p(1,25)&:; 	! both the call and the response are put in
C   glf23p(26,NA1)&:;	! init.tmp @ eqns.tmp i.e. after all
C
C Input: JS is the left grid point
C        JE is the right grid point
C Parameter JSBP is the ordinal number of the (&:) request in a model
C Parameter JSBR is the ordinal number of the (:) request in a model
C           JSBR is needed to store the CPU times used by all 
C       	 processes and subroutines in a homogeneous way.
C           JSBP and JSBR are generated by the Astra compiler.
C----------------------------------------------------------------------|
	implicit none
	integer JS,JE,JSBP,JSBR
	include	'for/parameter.inc'
	include 'for/status.inc'		! 
	include 'for/outcmn.inc'		! CPTSBR is used
C	double precision YAL
C	allocatable YAL(:,:)
C----------------------------------------------------------------------|
C	allocate (YAL(JE-JS,15))		! temporary work array
C	deallocate (YAL)
C----------------------------------------------------------------------|
	call	ot_glf23p(JS,JE,JSBP,CPTSBR(JSBR),
C glf output:
C here "21" is an arbitrary offset, i.e. 21 -> 21+14 NRD arrays are used
     &		work(1,21) )	! Employs work(1:NRD,21) -> work(1:NRD,35)
	end
C======================================================================|
