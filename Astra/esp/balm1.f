	subroutine balm1(yeqpf,ypsopt,ypsin,na1)
c(yeqpf,ypsopt,yxtr)
C...    BALLOONING AND MERCIER CRITERION
C...    PERIODICITY ASSUMED
C...    2*NPER PERIODS IN POLOIDAL DIRECTION SYMMETRIC AROUND J0
C...    NB > 2*NPER*M
C
c====================================================
c-3        IQAE       !q 0 axis 1 boundary switch
c2.75      QAE        !q axis or boundary
c5         NPER       !number of periods for balloon 
c-0.25                !theta0 in arcl
c0         INCR       !increment compute switch: 1 to compute
c100       IOPT       !local opt. iterations, 0 no opt.
c.01       EPSOPT     !local opt. PP accuracy
c--------------------------------------------------------
c	OUT:
c		WORK(1,J) = fp norm
c		WORK(2,J) = eqpf out = pressure gradient
c		WORK(3,J) = eqpf opt = Mercier limit for pres. grad.
c		iplas = number of points in spider equilibrium
c=======================================================
!        include 'double.inc'
!	include 'for/parameter.inc'
!	include 'for/const.inc'
!	real*4	WORK
!	common /WORKAR/	WORK(NRD,2*NRD+7)
	   implicit real*8(a-h,o-z)
         implicit integer(j-n)
 	parameter (NPER=5, INCR=0, IOPT=100, KPG=1)
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
	   parameter(nbz=nper*2*ntp)
!	   include 'implic.inc'

           include 'compol.inc'
	DIMENSION yeqpf(*),ypsopt(*),ypsin(*),ypsout(nrp),yout(nrp)
        DIMENSION AP(1:NRP),YFR(1:NRP),VK(1:NTP),UK(1:NTP),
     ,	VM(1:NTP),UM(1:NTP)
        DIMENSION DMERC(1:NRP),RMERC(1:NRP),RV(1:NRP),SV(1:NRP)
        DIMENSION DINCR(1:NRP),PPOPT(1:NRP)
c=============================================
        DIMENSION W(1:NRP,1:3),WT(0:NTP,1:7)
        DIMENSION WJ(1:NRP,1:NTP),WQ(1:NRP,1:NTP)
        DIMENSION WB(1:NBZ,1:5)
	EPSOPT =1.d-2
!iplas=N
!nt1=M
		J0=1
	j1=6*iplas+2
	j2=j1+(nt1+1)*14
	j3=j2+2*iplas*nt1
	j4=j3+2*iplas*nt1
!nbz=nper*2*nt1
	jtot=j4+10*nbz
	if(jtot.gt.nrd*(2*nrd+7)) 
     ,	write(*,*) 'balmod: jtot =',jtot,' out of range'	
	do j=1,iplas
		AP(j)=dsqrt(1.d0-psia(j))
	do jt=1,nt1
	ro(j,jt)=ro(j,jt)/ro(iplas,jt)
	enddo
	enddo
		y=0.d0
	do j=1,nt1

	if(y.lt.r(iplas,j)) then
		y=r(iplas,j)
		J0=j
	endif
	UM(j)=r(1,j)
	VM(j)=z(1,j)
	UK(j)=r(iplas,j)
	VK(j)=z(iplas,j)

	enddo
C...    F OUTER
	F0=fvac/(0.4d0*pi)
	PSM=(psim-psibon)/(0.4d0*pi)
          YFR(iplas)=F0**2
          DO  I=iplas-1,1,-1
            PFC=dfdpsi(I)+dfdpsi(I+1)
            HPS=-PSM*(AP(I+1)**2-AP(I)**2)
            YFR(I)=YFR(I+1)-PFC*HPS
 	enddo	
           do I=1,iplas
	      IF(YFR(I).LE.0.d0) WRITE(6,*) '!!! F<0'
	      Y=DABS(YFR(I))
            YFR(I)=DSQRT(Y)
       enddo
c        if(nequil.gt.42.) 
	call BALM0(NPER,J0,INCR,IOPT,EPSOPT,
     &                  KPG,psim/(0.4d0*pi),psibon/(0.4d0*pi),2.d0*pi,
     &                  iplas,AP,dfdpsi,dpdpsi,YFR,
     &                  nt1,UM,VM,UK,VK,
     &                  nrp,ntp,RO,
     &                  DMERC,RMERC,RV,SV,
     &                  DINCR,PPOPT,
     &                  W,WT,WJ,WQ,
     &                  NBZ,WB)
!     &                  work(1,1),work(j1,1),work(j2,1),work(j3,1),
!     &                  NBZ,work(j4,1))
!	open(33,file='bal-o.dat')
!        WRITE(33,*) '     PP11'
!        WRITE(33,*) (dpdpsi(I),I=1,iplas)
!        WRITE(33,*) '     OPT PRESSURE'
!        WRITE(33,*) (ppopt(I),I=2,iplas)
!	close(33)
	do j=1,iplas
		ypsout(j)=(1.d0-psia(j))/(1.d0-psia(iplas))
		yout(j)=dpdpsi(j)*rtor
	enddo
	call	TRANSF(     iplas,yout,ypsout,NA1,yeqpf,ypsin)
	do j=1,iplas
		yout(j)=ppopt(j)*rtor
	enddo

	call	TRANSF(     iplas,yout,ypsout,NA1,ypsopt,ypsin)

!		jplas=iplas
c	write(*,*) 'KOHEU,'
	return
	end	
