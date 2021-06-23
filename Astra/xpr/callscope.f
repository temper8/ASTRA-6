C======================================================================|
C Making SCoPE for ASTRA (Intel Fortran compiler)
C
C	cd a621
C	ifort -c xpr/callscope.f -o xpr/callscope.o
C	cc -c xpr/scope.c -o xpr/scope.o
C	cc -c for/ipc_control.c -o xpr/ipc.o
C	cd ../SRC_SCoPE/
C	ifort -c scope2astra.for
C	cd ../SCoPE/
C	rm -f 1* 2* [A-Z]*
C	ifort -o scope ../SRC_SCoPE/*.o ../a621/xpr/*.o
C
C The last line can be replaced with
C
C	ifort -o scope ../SRC_SCoPE/*.o \
C	   ../a621/xpr/scope.o ../a621/xpr/ipc.o ../a621/xpr/callscope.o
C----------------------------------------------------------------------|
	program   main
	implicit  none
	integer	  j,iargc,getpid,mampid,mamkey,eignr
	character*132 STRING,eigpath,mampath
	if (iargc() .ne. 4)	then
	   write(*,*)"Error"
	   call a_stop
	endif
	call	getarg(0,STRING)
	eigpath = STRING(1:len_trim(STRING))//char(0)
	call	getarg(1,STRING)
	mampath = STRING(1:len_trim(STRING))//char(0)
	call	getarg(2,STRING)
	read(STRING,*)mampid
	call	getarg(3,STRING)
	read(STRING,*)mamkey
	call	getarg(4,STRING)
	read(STRING,*)eignr
	write(*,*)
	call sbp2shm(eigpath, mampath, mampid, mamkey, eignr)
	end
C======================================================================|
	subroutine	callscope(
C Input:
     &			nmax,		! Maximum array dimension
     &			n1,		! Current array dimension
     &			time,		! Current time
     &			ro,		! rho_tor (see Eq.(24))
     &			te,		! T_e
     &			ti,		! T_i
     &			ne,		! n_e
     &			ni,		! n_i
     &			cc,		! Current conductivity
     &			cd,		! noninductive Current Density
     &			pr,		! plasma PRessure
C Output: length of the array can be changed by Astra 
     &			out)		! SCoPE output (currently 21048 bytes)
C----------------------------------------------------------------------|
	implicit none
        integer	 nmax,n1,j
        double precision time,out(*)
        double precision ro(*),te(*),ti(*),ne(*),ni(*),cc(*),cd(*),pr(*)
C----------------------------------------------------------------------|
C      call scope(
CC Input:
C     &                  nmax,           ! Maximum array dimension
C     &                  n1,             ! Current array dimension
C     &                  time,           ! Current time
C     &                  ro,             ! rho_tor (see Eq.(24))
C     &                  te,             ! T_e
C     &                  ti,             ! T_i
C     &                  ne,             ! n_e
C     &                  ni,             ! n_i
C     &                  cc,             ! Current conductivity
C     &                  cd,             ! noninductive Current Density
C     &                  pr,             ! plasma PRessure
CC Output: length of the array can be changed by Astra
C     &                  out)            ! SCoPE output (currently 21048 bytes)
C        out(3) = 			! Ipl
C        out(4) = 			! Rma
C        out(5) = 			! Phi_edge
C        do j=1,n1
C           out(5+j) = 			! Psi(1:n1)
C           out(nmax+5+j) = 		! current density
C           out(2*nmax+5+j) = 		! safety factor
C           out(3*nmax+5+j) = 		! r_min
C           out(4*nmax+5+j) = 		! Shafranov shift
C           out(5*nmax+5+j) = 		! elongation
C           out(6*nmax+5+j) = 		! triangularity
C           out(7*nmax+5+j) = 		! up-down shift
C           out(8*nmax+5+j) = 		! V'
C           out(9*nmax+5+j) = 		! G11
C           out(10*nmax+5+j) = 		! G22
C           out(11*nmax+5+j) = 		! G33
C           out(12*nmax+5+j) = 		! I (B_tor=I/r)
C           out(13*nmax+5+j) = 		! S_lat = V'<\nabla(rho)>
C           out(14*nmax+5+j) = 		! <B/B_0>
C           out(15*nmax+5+j) = 		! <(B/B_0))^2>
C           out(16*nmax+5+j) = 		! <(B_0/B)^2>
C           out(17*nmax+5+j) = 		! B_max
C           out(18*nmax+5+j) = 		! B_min
C           out(19*nmax+5+j) = 		! trapped fraction
C           out(20*nmax+5+j) = 		! Reserved
C           out(21*nmax+5+j) = 		! Reserved
C           out(22*nmax+5+j) = 		! Reserved
C           out(23*nmax+5+j) = 		! Reserved
C           out(24*nmax+5+j) = 		! Reserved
C           out(25*nmax+5+j) = 		! Reserved
C        enddo
        end
C======================================================================|
