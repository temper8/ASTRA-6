	subroutine	mykey
	implicit none
	integer	j,IFKEY
C A call of this subroutine produces the action caused by pressing 
C	key 'G' (output to PS file) in interactive regime. 
C (Enter "Astra -h" from UNIX prompt to find other key bindings)
	write(*,*)"MYKEY called"
	j = IFKEY(ichar('G'))
	j = IFKEY(ichar('N'))
	j = IFKEY(ichar('G'))
	end
