	subroutine	myout
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	j,IFKEY
C This command produces an action corresponding to pressing key 'G'
	j = IFKEY(ichar('G'))
	end
