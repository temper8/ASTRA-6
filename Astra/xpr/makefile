#XLBR=/usr/X11R6/lib/libX11.a
#LLIB= -L$(HOME)/work/f77 -lmyX11 -lm $(XLBR)
LLIB=-lm
OBJECTS=xpr/ipc.o
#AFL=$(myFC) -I. -w -o		#
#AFC=$(myFC) -I. -c		# Use Astra 
#ACC=cc -c			# environment
#ACL=cc -o			#
ALIB=.lbr/$(ObjectCode)/

#---------------
xpr/neutex:	xpr/neutex.f xpr/neut.o $(OBJECTS)
		$(AFL) xpr/neutex xpr/neutex.f xpr/neut.o $(OBJECTS) 
# For unclear reason xpr/neutex.f cannot be linked with $(ALIB)for.a 
# in place of xpr/ipc.o as it is done for xpr/glf23p.f
#
# At the same time xpr/ipc.o can be created making use of ipc/ipc_control.c
# and then used in the linker.
# The functions   a_stop(freeshm), swatch, write_aipc  and  SP_stamp 
# are used from ipc/ipc_control.c
#
xpr/neut.o:	xpr/neutex.c
		$(ACC) xpr/neutex.c -o xpr/neut.o
#---------------
$(OBJECTS):	ipc/ipc_control.c
		$(ACC) ipc/ipc_control.c -o $(OBJECTS)
#---------------
xpr/template:	xpr/template.f xpr/temp.o $(OBJECTS)
		$(AFL) xpr/template xpr/template.f xpr/temp.o $(OBJECTS)
#		rm xpr/temp.o

xpr/temp.o:	xpr/template.c
		$(ACC) xpr/template.c -o xpr/temp.o
#---------------
xpr/scope:	xpr/callscope.f xpr/scope.o $(OBJECTS)
		$(AFL) xpr/scope xpr/callscope.f $(OBJECTS) xpr/scope.o

xpr/scope.o:	xpr/scope.c
		$(ACC) xpr/scope.c -o xpr/scope.o
#---------------
# Use this to compile C-main
#scope:		xpr/scope.c xpr/callscope.o
#		$(ACL) xpr/scope xpr/callscope.o $(LLIB) -lifcore xpr/scope.c \
#		-L/opt/intel/fc/9.1.039/lib ipc/ipc_control.c
##		-L/opt/intel/Compiler/11.0/074/lib/intel64 ipc/ipc_control.c
#
#xpr/callscope.o: xpr/callscope.f
#		 $(AFC) scope.f -o xpr/callscope.o
#---------------
example: 	xpr/example.c
		$(ACL) xpr/example $(LLIB) -lm -lifcore \
		-L/opt/intel/fc/9.1.039/lib ipc/ipc_control.c xpr/example.c

