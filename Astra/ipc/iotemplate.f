C======================================================================|
	subroutine TOTEMPLATE(JSBP)
C----------------------------------------------------------------------|
C Call of this subroutine is created by the Astra compiler
C The C-function to_neut_ fills shared memory segment Nr. (JSBP+1)
C Parameter JSBP is the ordinal number of the (&:) request in a model
C     to the C-function to_neut_
C No variable transfer to Astra
C----------------------------------------------------------------------|
	implicit none
	integer JSBP
        end
C======================================================================|
	subroutine OTTEMPLATE(JSBP,JSBR)
C----------------------------------------------------------------------|
C Call of this subroutine is created by the Astra compiler
C Parameter JSBP is the ordinal number of the (&:) request in a model
C Parameter JSBR is the ordinal number of the (:) request in a model
C No variable exchange with Astra except for CPUse
C----------------------------------------------------------------------|
	implicit none
	integer JSBP,JSBR
	include	'for/parameter.inc'
	include 'for/outcmn.inc'
	call	ot_template(JSBP,CPTSBR(JSBR))
        end
C======================================================================|
