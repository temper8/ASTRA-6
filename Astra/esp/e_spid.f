	subroutine e_spid
c(itrel,relerr)
c========================================= version 1-FEB-09 for Astra 6.2.xx
c	Out: 	itrel- number of iterrations to coverge
c		relerr - relerr for convergency
	 include 'for/parameter.inc'
	parameter(NBNDX=125, NBNTX=200)
	 include 'for/const.inc'
	 include 'for/status.inc'
	 include 'for/outcmn.inc'
	 integer  ngapp,npfc0,nrp,ntp,nr1p,nt1p,nr2p,nt2p,NJLIM,NpLIM
	include 'dimpl1.inc'
!         parameter(nrp=256,ntp=256,nr1p=nrp-1,nt1p=ntp-1,
!     *             nr2p=nrp-2,nt2p=ntp-2)
c_aai
          include 'parevo.inc'
            !PARAMETER ( NJLIM=1550 )
            !PARAMETER ( NPFC0 = 40 )
            !PARAMETER ( NpLIM = 400 )
          include 'dimcp.inc'
!       parameter(ngapp=10)
c_aai
	 integer ni_p,nj_p,nbnd_p,nbnd2_p
          parameter(ni_p=nrp,nj_p=ntp)
          parameter(nbnd_p=125,nbnd2_p=nbnd_p*2)
         character*40 prename
         character*40 eqdfn

c==================================================================
        integer ksym,knstep,keyctr,k_grid,k_auto,key_ini,k_con,kpr,
     ,  	key_dmf,nstep,k_fixfree,k_kinx,icount,
     ,  	n_gap,n_pfc,NBDX,NBTX
	integer jnstep,jni_p0,jnj_p0,j,ji,jj,jneql,jnteta,jna1,jdna1
c==================================================================
	real*8 	yrout(nrp,ntp),yzout(nrp,ntp),yrho(nrd),
     *		ya(nrd),ytria(nrd),yelon(nrd),yshif(nrd),yshiv(nrd),
     *		XTR(nrd),yxeq(NRD),yXTR(NRD),eqff1(NRD),yrr,
     *          yrmax,yrmin,yzmax,yzmin,yli3,yreler,platok,dpsdt,
     *          yrzmin,yrzmax,R_ref,Z_ref,tok_ref,alfa,ydrhodt,
     *		yrtor,ybtor,yroc1,yroc,yrocnew,yipl,yyipl
	real*8 	
     *  yg11(nrd),yg22(nrd),yg33(nrd),yvr(nrd),yvrs(nrd),
     * 	yslat(nrd),ygradro(nrd),ymu(nrd),yipol(nrd),
     *  ybmaxt(nrd),ybmint(nrd),ybdb02(nrd),yb0db2(nrd),ybdb0(nrd),
     *	yeqpf(nrd),yeqff(nrd),yvolum(nrd),yfpo(nrd),
     *	ycc(nrd),yTe(nrd),ycubs(nrd),ycd(nrd),ycu(nrd),yfp(nrd)
	real*8		d_t,t_ime,roc1,rocnew1,rocnew,NEQUIL1
	real*8		ypres(nrd),W_Dj(nrd)
!	real*8		pres(nrd),yp	!a620 only
                


	real*8		yneold,yerr,rzbnd(nbnd2_p)
!rzbnd(2*NBNDX)
	real*8		R_ax,Z_ax,tokpl,vz
       

        real*8 currpf_ref(npfc0),Dgap_ref(ngapp)
        real*8 cur_pf(npfc0),gap(ngapp)
       save R_ref,Z_ref,tok_ref
       save Dgap_ref,currpf_ref

	
	save 	jnstep,   jni_p0,   jnj_p0,   icount,   yneold,   yerr
	save 	nstep
	data	jnstep/0/ jni_p0/0/ jnj_p0/0/ icount/0/ yneold/0.d0/
	data	yerr/1.d-5/ nstep/0/


!=======================================BPEMEHHO read vvvvvvvvvvvvv      
	Open(1,file='exp/equ/spidatin.dat')
	 Read(1,*) kpr   	! k=1   !print in spider
	 Read(1,*) k_grid	! k_grid= 0   rect. grid ! k_grid= 1   adap. grid
	 Read(1,*) k_auto	!= 0   ! k_auto= 1->   full initialization 
!	 Read(1,*) keyctr	!= 0 
	 Read(1,*) k_fixfree	!= 0 fixed boundary, = 1 free boundary 
	 Read(1,*) k_kinx       != 0 no output for kinx, 1 output for kinx
	 Read(1,*) key_dmf	!= -2 dpi_ext/dt+ diff.mag.field
				! !=-3-> diff.mag.field, =0-> without
	 Read(1,*) prename	!='exp/equ/ITER/'
	 Read(1,*) ksym		!= 13 length of prenamne
	 Read(1,*) k_con		!0= no controller
	 Read(1,*) dpsdt		!dpsi_ext/df (if key_dmf = -2)
	 Read(1,*) key_ini	! =1 astra profiles, =0 start from EQDSK and SPIDER profiles	
	 Read(1,*) eqdfn		!eqdsk file name
        close(1)
c=======================================================================
C========================Input parameters
		yrtor	=rtor
		ybtor	=btor
		yroc	=roc
		yrocnew	=roc
		PLATOK	=IPL
		jna1	=na1
	  	d_t 	=TAU
	   	t_ime 	=TIME 
	do j=1,nrd
		ymu(j)		=mu(j)
		yeqpf(j)	=eqpf(j)
		yeqff(j)	=eqff(j)
		yrho(j)		=rho(j)
		ycc(j)		=cc(j)
		yTe(j)		=Te(j)
		ycubs(j)	=cubs(j)
		ycd(j)		=cd(j)
		ycu(j)		=cu(j)
		yfp(j)		=fp(j)
		yfpo(j)		=fp(j)
	if(key_ini.ne.0) 
     *		pres(j)		= 1.6d-3*(PFAST(j)+
     *	 	NE(J)*TE(J)+NI(J)*TI(J)+0.5*NB2EQL*(PBLON(J)+PBPER(J)))
	! read ASTRA pressure	
!		car2x(j)	=pres(j)
		ypres(j)	=pres(j)
	if(key_dmf.eq.-2) dpsdt=UEXT	
		yyipl=yfpo(na)+	(yfpo(na1)-yfpo(na))*hro/hroa !old value in a regular point
	roc1=roc
	enddo

c=========================================================================
 
	if(nstep.eq.0) then
	   write(*,*)'equlibrium input from ', prename
	if(key_dmf.eq.-2) then
	
	write(*,*) 'external voltage UEXT(t) for current diffusion'
	endif
	endif
         call  kpr_calc(kpr)
         call  put_name(prename,ksym)
         call put_key_con(k_con)

!=======================================BPEMEHHO read ^^^^^^^^^^^^ 
	 if(yneold.ne.0.d0) then
	if(nequil.ne.yneold.and.key_dmf.ne.0) then
	   nequil = yneold
	   write(*,*) 
     . 'exp/equ/spiderin.dat: ',
     . 'for current diffusion (key_dmf=-2,-3) NEQUIL can not be changed' 
	    endif
	    endif
!aai
	if(nequil.ne.yneold) then 
!aai		icount = 3
C...Initial iterrations start every time when NEQUIL > = 43 was changed
	write(*,*) 'Attention! 	Metric for SPIDER solver:'
	write(*,*) 
     .'boudary points Rij,Zij must be in contr clock wise direction'
     
	write(*,*) 'new mesh point numbers for SPIDER equilibrium solver'
		jnj_p0	=NEQUIL
		jni_p0	=(100*(NEQUIL-jnj_p0))
	if(jni_p0.le.43) then
		jni_p0 = 44
	write(*,*) 'for spider equilibrium is changed:'
	write(*,*) 'OLD NEQUIL = ', NEQUIL,'NEW NEQUIL =',jnj_p0 + 0.44d0	
		NEQUIL = jnj_p0 + 0.44d0

	endif
		yneold	=NEQUIL

	write(*,*) 'N poloidal = Int(NEQUIL) 		= ', jnj_p0
	write(*,*) 'N radial = 100*(NEQUIL - Int(NEQUIL)) 	= ', jni_p0

	else
		jnstep	=1
	endif
		jnteta=jnj_p0
		jneql=jni_p0

	if(icount.gt.0) then
		jnstep	= 0
	 write(*,*) 'E_spid: iteration N ',icount
	 icount = icount - 1
	endif

cw	write(*,*) NBND,NEQUIL,jnteta,jneql

c=======================================
	if(NBND.lt.4.or.NBND.gt.NBNDX ) then
	write(*,*) 'e_spid: wrong number of boundary points NBND = ', NBND
	write(*,*) 'NBND must be: 4 < NBND < ',NBNDX
 	write(*,*) 'contunue with 3M Equilibrium solver' 
		NEQUIL = 41.d0
		return
	else

		call BNDRY1(rzbnd)


	endif


c===== write 'dat/spidat.dat' for SPIDER stand alone analysis========
	if(kpr.eq.2) then
	include 'wrtspidat.inc'
	endif
c=====================================================================
c===================================== remember data for old grid 

	       NA1O =NA1
               HROA =ROC-RHO(NA1-1)
	       YIPL = (FP(NA1)-FP(NA))/HROA
	do	J = 1,NA1O-1
	   	yXTR(J) = RHO(J)/ROC	
	enddo
		yXTR(NA1O)=1.d0
c======================================================================
	
          if(jnstep.eq.0) then
	      nstep =jnstep
          else
	      nstep =nstep+1
          endif


c======================================================================
c============================================= start of spider

        call aspid_flag(1)
c========================================= interpolation to Spider grid
        call astra2spider(jneql,jnteta,nbnd,rzbnd,key_dmf,
     *                    jna1,yeqpf,yeqff,yfp,PLATOK,
     *                    yrtor,ybtor,yrho,yroc,nstep,yreler,ymu,
     *                    ycc,yTe,ycubs,ycd,key_ini,eqdfn )


!!!!!controllable variables:{gaps(1:6),Jpl,pfc(1:11)}

c=============================================== gets data for controller
	if(k_fixfree*k_con.ne.0) 
     *	  call put_conp(R_ref,Z_ref,tok_ref,Dgap_ref,currpf_ref)
c========================================= equil + current diffusion
      	  call spider(nstep,t_ime,d_t,key_dmf,k_grid,k_auto,k_fixfree,
     *            dpsdt,key_ini)  
c============ takes data for controller and contr. values to ASTRA for output
	if(k_fixfree*k_con.ne.0) then
	  call get_conp(R_ref,Z_ref,tok_ref,Dgap_ref,currpf_ref)
	  call get_cpar(R_ax,Z_ax,tokpl,gap,cur_pf,vz,n_gap,n_pfc)
	endif
c===========================================================================
	call get_roc(yROCNEW,btor) !takes new toroidal flux from SPIDER
		rocnew=yrocnew
		ydrhodt=(yrocnew-yroc)/TAU
cp	write(*,*) 'roc, rocnew = ',roc, rocnew
		ROC = yrocnew

c== changes nimber of radial grid points NA1 for new ROC for ASTRA if necessary
	call	NEWGRD	
!
	if(HROA.le.0.) write(*,*) 'na1,na1o,hroa',na1,na1o,hroa

		jna1=na1
cp	write(*,*) 'NA1O,NA1,jna1 =', NA1O,NA1,jna1
!====================================================== new normalised grid
	do	J = 1,NA1
	   	XTR(J)  = RHO(J)/ROC	
	enddo
		XTR(NA1) =1.d0
	do      j=1,nrd
	        yrho(j)  =rho(j)
	enddo
c=========================== interpolation of SPIDER output to ASTRA grid 
        call spider2astra(yrout,yzout,yrtor,ybtor,yrho,yrocnew,jna1,
     *                    yg11,yg22,yg33,yvr,yvrs,yslat,ygradro,yroc,
     *                 ymu,yipol,ybmaxt,ybmint,ybdb02,yb0db2,ybdb0,yxeq,
     *			  yreler,yli3,ni_p,nj_p,platok,ycu,yfp,ypres,W_Dj)
cp	write(*,*) 
cp     *  'end spider2astra NA1,NA1O,roc,rocnew',NA1,NA1O,roc,rocnew

	NA=NA1-1
C Volume (on the shifted grid) is calculated using VR:
		VOLUM(1) = HRO*yVR(1)
	do	J=2,NA
	   	VOLUM(J) = VOLUM(J-1)+HRO*yVR(J)
	enddo
		VOLUM(NA1) = VOLUM(NA)+(HROA-0.5*HRO)*yVR(NA)
		VOLUME = VOLUM(NA1)

c======================================== write li3 from SPIDER to dat/lin3
	open(33,file='dat/lin3')
		write(33,*) yli3
	close(33)

c================================ output for ASTRA ================================

		ipl	=PLATOK

c		volume	=yvolume
	do j=1,nrd
		g11(j)		=yg11(j)
		g22(j)		=yg22(j)
		g33(j)		=yg33(j)
		vr(j)		=yvr(j)
		vrs(j)		=yvrs(j)
		slat(j)		=yslat(j)
		gradro(j)	=ygradro(j)
		ipol(j)		=yipol(j)
		bmaxt(j)	=ybmaxt(j)
		bmint(j)	=ybmint(j)
		bdb02(j)	=ybdb02(j)
		b0db2(j)	=yb0db2(j)
		bdb0(j)		=ybdb0(j)
c		volum(j)	=yvolum(j)
	if(key_dmf.ne.0.or.key_ini.eq.0) then
		cu(j)		=ycu(j)
		fp(j)		=yfp(j)
		fpo(j)		=yfpo(j)
		mu(j)		=ymu(j)
	if(key_ini.eq.0) pres(j)=ypres(j) !read Spider pressure from EQDSK file
	endif
!		car1x(j)	=ypres(j)
	enddo
c=================================== dpsi/drho = 0.4*pi*Ipl/G2
	YIPL=0.4*GP*IPL*RTOR/(IPOL(NA1)*g22(NA))

c1
	if(key_dmf.ne.0) then
c1
		fp(na1)=yipl*hroa+fp(na)
c	write(*,*) 'fp(na1),yipl,hroa,fp(na)'
c	write(*,*) fp(na1),yipl,hroa,fp(na)
c	write(*,*) 'fp(na),yyipl,hro,fp(na-1)'
c	write(*,*) fp(na),yyipl,hro,fp(na-1)
c	write(*,*) 'Dpsi ASTRA',fp(na1)-fp(1)
c	write(*,*) 'Dpsi SPIDER',yfp(na1)-yfp(1)
c1
	endif
c 997	continue
c====================================================================
	if(NA1O.ne.NA1) then
c====================== extention of arrays to the new grid
	  if(NA1O.gt.NA1) then
	     jdna1=NA1O-NA1
	  else
	     jdna1=NA1-NA1O
	  endif
	  if(jdna1.gt.1) then
		d_t=d_t/jdna1
	write(*,*) 'Warning: too fast volume changes,',
     >		' time step is reduced'
	write(*,*) 'OLD grid: ',NA1O,' NEW grid: ',NA1
	write(*,*) 'OLD step,s: ',TAU,' NEW step: ',d_t
	  endif

	if(jnstep.eq.0.or.jdna1.gt.1) then
c=  remapping of input parameters to consistent grid at initial itterations
	do	J = 1,NA1
	   	XTR(J) = RHO(J)/ROC	
	enddo
		XTR(NA1)=1.d0
	
		ALFA=1.d-6 

	call	SMOOTH(ALFA,NA1O,TI,yXTR,NA1,TI,XTR)
	call	SMOOTH(ALFA,NA1O,TE,yXTR,NA1,TE,XTR)
	call	SMOOTH(ALFA,NA1O,NI,yXTR,NA1,NI,XTR)
	call	SMOOTH(ALFA,NA1O,NIZ1,yXTR,NA1,NIZ1,XTR)
	call	SMOOTH(ALFA,NA1O,NIZ2,yXTR,NA1,NIZ2,XTR)
	call	SMOOTH(ALFA,NA1O,NIZ3,yXTR,NA1,NIZ3,XTR)
	call	SMOOTH(ALFA,NA1O,NHE3,yXTR,NA1,NHE3,XTR)
	call	SMOOTH(ALFA,NA1O,NALF,yXTR,NA1,NALF,XTR)
	call	SMOOTH(ALFA,NA1O,NHYDR,yXTR,NA1,NHYDR,XTR)
	call	SMOOTH(ALFA,NA1O,NDEUT,yXTR,NA1,NDEUT,XTR)
	call	SMOOTH(ALFA,NA1O,NTRIT,yXTR,NA1,NTRIT,XTR)
	call	SMOOTH(ALFA,NA1O,NE,yXTR,NA1,NE,XTR)
	call	SMOOTH(ALFA,NA1O,F0,yXTR,NA1,F0,XTR)
	call	SMOOTH(ALFA,NA1O,F1,yXTR,NA1,F1,XTR)
	call	SMOOTH(ALFA,NA1O,F2,yXTR,NA1,F2,XTR)
	call	SMOOTH(ALFA,NA1O,F3,yXTR,NA1,F3,XTR)
	call	SMOOTH(ALFA,NA1O,F4,yXTR,NA1,F4,XTR)
	call	SMOOTH(ALFA,NA1O,F5,yXTR,NA1,F5,XTR)
	call	SMOOTH(ALFA,NA1O,F6,yXTR,NA1,F6,XTR)
	call	SMOOTH(ALFA,NA1O,F7,yXTR,NA1,F7,XTR)
	call	SMOOTH(ALFA,NA1O,F8,yXTR,NA1,F8,XTR)
	call	SMOOTH(ALFA,NA1O,F9,yXTR,NA1,F9,XTR)
	call	SMOOTH(ALFA,NA1O,pfast,yXTR,NA1,pfast,XTR)
	call	SMOOTH(ALFA,NA1O,pblon,yXTR,NA1,pblon,XTR)
	call	SMOOTH(ALFA,NA1O,pbper,yXTR,NA1,pbper,XTR)
	if(key_dmf.eq.0) then
	call	SMOOTH(ALFA,NA1O,fp,yXTR,NA1,fp,XTR)
	call	SMOOTH(ALFA,NA1O,fpo,yXTR,NA1,fpo,XTR)
	call	SMOOTH(ALFA,NA1O,cu,yXTR,NA1,cu,XTR)
	call	SMOOTH(ALFA,NA1O,mu,yXTR,NA1,mu,XTR)
	   endif
	call    EDCEL(NA1O,YIPL)
	else
	if(key_dmf.ne.0)   then
		fpo(na1)=fpo(na)
	if(na1.gt.na1o)	fpo(na)=yyipl	!map to a regular point
	endif
	call	EDCELL(NA1O,YIPL)
	call    EDCEL(NA1O,YIPL)
	endif
	endif

ccc 
 777	continue

c===========================================================================
	 if(key_dmf.ne.0) then

	 do j=1,NA
	         UPL(J)=(FP(j)-FPO(j))/TAU
		 ULON(J)=IPOL(J)*G33(J)*UPL(J)
	      enddo
	         UPL(NA1)=UPL(NA) + gp2*btor*roc*mu(na1)*ydrhodt
		 ULON(NA1)=IPOL(NA1)*G33(NA1)*UPL(NA1)
	endif
C================== calculation of 3M from SPIDER output geometry

	do ji=jneql,2,-1
		yrmax=-99999.d0
		yrmin=99999.d0
		yzmax=-99999.d0
		yzmin=99999.d0
	do jj=1,jnteta
		if(yzmin.ge.yzout(ji,jj)) then 
			yzmin=yzout(ji,jj)
			yrzmin=yrout(ji,jj)
		endif
		if(yzmax.le.yzout(ji,jj)) then 
			yzmax=yzout(ji,jj)
			yrzmax=yrout(ji,jj)
		endif
		if(yrmin.ge.yrout(ji,jj)) then 
			yrmin=yrout(ji,jj)

		endif
		if(yrmax.le.yrout(ji,jj)) then 
			yrmax=yrout(ji,jj)
		endif
c		yrmin=min(yrout(ji,jj),yrmin)
c		yrmax=max(yrout(ji,jj),yrmax)
	enddo
		yrr=.5d0*(yrmax+yrmin)

	if(ji.eq.jneql) SHIFT=yrr-RTOR
		ya(ji)=.5*(yrmax-yrmin)

		yshif(ji)=yrr-RTOR
		yshiv(ji)=.5d0*(yzmin+yzmax)
		yelon(ji)=(yzmax-yzmin)/(yrmax-yrmin)
		ytria(ji)=(yrr-0.5d0*(yrzmin+yrzmax))/ya(ji)

	enddo
		ya(1)=0.d0
		yelon(1)=yelon(2)
		ytria(1)=0.d0
		yshif(1)=yrout(1,1)-RTOR
		yshiv(1)=yzout(1,1)
C=========== output of R,Z coodinates for mag. surf. to file dat/spidat-RZ.dat
	if(kpr.eq.1) then
	   open(1,file='dat/spidat-RZ.dat')
	   write(1,*) 'jR,jteta = ', jNEQL,jnteta
	   write(1,*) 'axis: Z,R = ', yzout(1,1),yrout(1,1)
	   write(1,*) 'jR, jTet, Z(jR,jTeta), R(jR,jTet)'
	   do ji=jneql,2,-1
	      do jj=1,jnteta
		 write(1,*) ji,jj, yzout(ji,jj),yrout(ji,jj) 
	      enddo
	   enddo
	   close(1)
	endif
c============================================================= write  out end

	do	J = 1,NA1
	   	XTR(J) = RHO(J)/ROC
	enddo
		XTR(NA1)=1.d0
		yxeq(jneql)=1.d0
		yxeq(1)=0.d0
		alfa=1.d-6	
C BPEMEHHO 
cp	write(*,*) 'jneql,NA1,XTR,yxeq =',jneql,NA1,XTR(NA1),yxeq(jneql)
        call	TRANSF(     jneql,ya,yxeq,NA1,AMETR,XTR)
	call	SMOOTH(ALFA,jneql,ySHIF,yxeq,NA1,SHIF,XTR)
        call	SMOOTH(ALFA,jneql,ySHIV,yxeq,NA1,SHIV,XTR)
        call	SMOOTH(ALFA,jneql,yELON,yxeq,NA1,ELON,XTR)
        call	SMOOTH(ALFA,jneql,yTRIA,yxeq,NA1,TRIA,XTR)
c 

	do       j=NA1,NB1
	   shif(j) = shif(na1)
	   ametr(j)= ametr(na1)+(AB-ametr(na1))*(j-na1)/(nb1-na1)
	   elon(j)=elon(na1)+(elonm-elon(na1))*(j-na1)/(nb1-na1)
	   tria(j)=tria(na1)+(trich-tria(na1))*(j-na1)/(nb1-na1)
	enddo

		updwn=yzout(1,1)
		SHIF(NA1)=SHIFT
		tria(1)=0.
		droda(1) =rho(1)/AMETR(1)
	do	j=2,NA1
		droda(j) =(rho(j)-rho(j-1))/(AMETR(j)-AMETR(j-1))
	enddo
		roc	=rocnew
  	        ABC=AMETR(NA1)
		elong=ELON(NA1)
		TRIAN=TRIA(NA1)           

		jnstep=1+jnstep
c======================================== save data for KINX analysys
		if(k_kinx.ne.0) call wr_spik

cp	write(*,*) 'jneql,NA1,XTR,yxeq =',jneql,NA1,XTR(NA1),yxeq(jneql)
cp	write(*,*) 'end e_spid NA1=',NA1
c==================================================== new timestep
		TAU=d_t
	return
	end
C=========================================== remapping for cold neutrals
	subroutine	EDCEL(JNA1O,YIPL)
C NN defined in this subroutine will be assigned to FPO in OLDNEW
C
c	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	JNA1O,j,jdel
	real*8	YIPL
	if (NA1 .gt. JNA1O)	then	! NA1 increases: NA1O >= NA
		jdel=NA1-JNA1O
	do j=1+jdel,NA1
		NN(j)=NN(j-jdel)
	enddo
	do j=1,jdel
		NN(j)=NN(1)
	enddo
	endif
	if (NA1 .lt. JNA1O)	then		! NA1 decreases


		jdel=JNA1O-NA1
	do j=1,NA1
		NN(j)=NN(j+jdel)
	enddo
	endif
	end
C======================================================================|
	

