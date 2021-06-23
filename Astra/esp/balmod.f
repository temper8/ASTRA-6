	subroutine balmod(yeqpf,ypsopt)
C...    BALLOONING AND MERCIER CRITERION
C...    PERIODICITY ASSUMED
C...    2*NPER PERIODS IN POLOIDAL DIRECTION SYMMETRIC AROUND J0
C...    NB > 2*NPER*M
C       double precision version for ASTRA 6.x.x ==18.01.09
c====================================================
c-3        IQAE       !q 0 axis 1 boundary switch
c2.75      QAE        !q axis or boundary
c5         NPER       !number of periods for balloon 
c-0.25                !theta0 in arcl
c0         INCR       !increment compute switch: 1 to compute
c100       IOPT       !local opt. iterations, 0 no opt.
c.01       EPSOPT     !local opt. PP accuracy
c--------------------------------------------------------
c	yeqpf = eqpf out
c	ypsopt = eqpf opt
c=======================================================
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer j
	real*8 ypsin(nrd),yeqpf(*),ypsopt(*),yfpc,yfpa

		yfpc=fp(1)-mu(1)*btor*gp*(rho(2)-rho(1))*(rho(2)+rho(1))
		yfpa=fp(na1)-yfpc

	do j=1,na1
	 	ypsin(j)=(fp(j)-yfpc)/yfpa
	enddo

	call	balm1(yeqpf,ypsopt,ypsin,na1)

c	write(*,*) 'KOHEU,'
	return
	end	
