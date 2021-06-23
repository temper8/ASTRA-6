C======================================================================|
	subroutine	MYDRAW
C----------------------------------------------------------------------|
C The subroutine allows to display user's drawing
C----------------------------------------------------------------------|
	implicit 	none
	character	STRI*132
	integer	lonlen
C	call 	erasrw		! Uncomment this line to erase frame
	call	colovm(2)
	STRI = "This mode is reserved for user's drawing"
	call	textvm(150,150,STRI,lonlen(STRI))
	STRI = 'See example in the file "AWD/sbr/mydraw.f"'
	call	textvm(150,170,STRI,lonlen(STRI))
	call	rectvm(0,100,100,100,100)
	call	redraw(0)
C	call	sleep(1)
C	call	colovm(1)
C	call	rectvm(1,150,150,100,100)
c	STRI = 'Done'
c	call	textvm(150,130,STRI,lonlen(STRI))
C	call	redraw(0)
c	call	makeAUGcnf
c	call	drawconf
c	call	IFKEY(32)
	end
C=======================================================================
