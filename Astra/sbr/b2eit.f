c************************************************************************	
        subroutine B2EIT(YGPUFF,YGDTIS,YSCORE,YPFUS,YPSOL,
     >                            YPI,YSDT,YCAR,
     >                     YQPK,YNES,YGDTS,YNNWM,YNHES,YGHES,
     >                     YTES,YTIS,YNCS,YGCT,YPIMP, YENWM,YEHE)
c================================================Version 18-MAR-2004
c         Calculation of Core/SOL/DIV compatible boundary conditions for 
c ITER on the basis of B2/EIRENE CAlculations for inductive reference 
c scenario 
c=========================================================== References:	
c     [1]Kukushkin A.S. et al, NF 43 (2003) 716-723
c     [2] A.S.Kukushkin,H.D.Pacher,D.Coster et al, 
c         "Effect of Carbon Redeposition on the 
c     Divertor Performance", EPS 2003, Snct.Petersburg
c============================================ programmed by A.Polevoi
c============================================          and A.Zolotukhin
c======================================================================
c	INPUT:
c               YGPUFF 10^19 atoms/s    gas puffing flux
c               YGDTIS 10^19 ions/sec   radial ion net flux 
c                                       across the separatrix
c		YSTOT  10^19 atoms/s	total particle throughput
c		YSCORE 10^19 atoms/s	core fuelling (< STOT)
c		YPFUS, MW		total fusion power
c		YPSOL, MW		total power to SOL
c		YPI, MW			total power to SOL with ions (< YPSOL)
c		YSDT, m3/s 		pumping speed
c 		YCAR 0 = metallic wall; ne.0 full carbon wall
c	OUTPUT:
c		
c		YQPK,'	= qpk, MW/m2 out	divertor load
c		YNES,'	= ne_s, 10^19/m3	ne(s)
c               YGDTS,' = Gamma_DT_sep, 10^19/s'DT neut gas puff to the core
c		YNNWM,   = n_neut, 10^19/m^-3	Neutral density at the separatrix
c		YNHES,'	= nHe_s, 10^19/m3	nHe(s)
c		YGHES,'	= GHe_n_s, 10^19/s	He neut. recycling to the core
c		YTES,'	= Te_s, keV 		Te(s)
c		YTIS,'	= Ti_s, keV 		Ti(s)
c		YNCS,'	= nC_s, 10^19/m3	nC(s)
c		YGCT,'	= GC_t,  10^19/s	C flux at target
c		YPIMP,'	= Pimp, MW		Pimp_rad, (radiation in DIV?)
c		YENWM,'	= En_s, keV		EDT,neut(s)
c		YEHE,'	= EHe_s, keV		EHe,neut(s)
C====== YGDT,YGCOR Pa m3/s (molecular species) total and core prtcl throughput: 
C=====	STOT 10^19 atoms/s (1 Pa m3/s = 54 (27) 10^19 at/s for DT(molecules) 
C=====                                              (He, C  atoms))
c=======================================================================
c   Some formulas and definition of input/output arrays in ASTRA
c
c   YSTOT = Gamma_puff + Gamma_{DT,ion} - Gamma_{DT,neutral},
c
c   where   Gamma_puff <= 10^4 [10^19 part/sec] ,
c
c                            /      dn
c           Gamma_{DT,ion} = \ (-D --- + V n) dS = QNB - (Z_1*QF1B+Z_2*QF2B+...) 
c                            /      dr
c
c           Gamma_{DT,neutral} = G_{DT,n,in} - G_{DT,n,out} = gamma*n_DT               
c
c           /
c         = \dV n_e*(S_ion - S_rec) ,
c           /
c          where S_ion = [<sigmaV>_ie + <sigmaV>_ii*n_i/n_e]*NNWM*NN,
c                S_rec = <sigmaV>_rec*n_i
c
c                                      /
c                 Gamma_{DT,neutral} + \dV n_e*S_rec
c         NNWM = ----------------------/---------------------
c                 /
c                 \dV n_e*[<sigmaV>_ie + <sigmaV>_ii*n_i/n_e]*NN
c                 /
c          Here the density and the energy of incoming neutrals 
c          are determined as "warm neutals", NNWM and ENWM, respectively)
c       ------------------------------------------------------------
c
c   YSCORE = Gamma_{NB,pellet} - 2*Gamma_fusion ,
c          
c   YSCOR = Gamma_{DT,neutral} + Gamma_{NB,pellet} - 2*Gamma_fusion ,
c          
c   where
c            Gamma_{NB,pellet} = ... + vint(SNEBM,ABC)
c
c            Gamma_fusion = vint(CARxx, ABC),
c                           CARxx = NDEUT*NTRIT*SVDT
c       ------------------------------------------------------------
c   YPFUS = 5.* QDTB
c   YPSOL = QIB + QEB
c   YPI   = QIB
c   YSDT =< 50 ~ 20
c   YCAR = 0 if metallic wall, ne.0 if full carbon wall
c=======================================================================
        implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
c	Input
        double precision YGPUFF, YGDTIS, YSCORE 
	double precision YPFUS, YPSOL, YPI, YSDT, YCAR
c	Output
        double precision YQPK,YNES,YGDTS,YNNWM
	double precision YNHES,YGHES,YTES,YTIS,YNCS
        double precision YGCT, YPIMP, YENWM, YEHE
c	Intrinsic variables
        double precision YSCOR, YSTOT,YGDT, YGCOR, YETAC, YFF
        double precision YPE, YFHE, YXEI
	double precision YGDTSold, YDELTA, Yeps
	double precision Y, SVIE, SVII, SVREC
	double precision YSNI(NRD),YSNR(NRD)
	double precision YVINTREC,YVINTION
	double precision vint
	double precision SNNI,SNNR,SNNEU
	integer j
	
	parameter (Yeps=.01)
		
	save YGDTSold
	data YGDTSold /0./
	
		if (NN(NA).le.0.) return	

		
 100		YSTOT   = YGPUFF + YGDTIS - YGDTSold
                YSCOR   = YSCORE  + YGDTSold
		
                YGDT    = YSTOT/54.        !10^19 atoms/s -> 1 Pa*m3/s
                YGCOR   = YSCOR/54.        !10^19 atoms/s -> 1 Pa*m3/s
		
        if(YGCOR.lt.0.) then
                YGCOR =0.
        write(*,*) 'from B2ITER:    
     >              negative core fueling? GCORE= ',YGCOR
        endif
	
        if(YGCOR.gt.YGDT)       YGDT    = YGCOR
	
                YETAC   = YGCOR/YGDT
                YFF     = 1.+0.18*YETAC
        
	if(YPSOL.le.0.)         YPSOL   = 0.1
        
	if(YPI.gt.YPSOL) then
           write(*,*) 'from b2iter: warning Pi = ',
     >              YPI,'> Psol ', 'Suggest Psol=Pi'
                YPSOL = YPI
        endif
        
	if(YPI.le.0.)   YPI=.0001*YPSOL
                YPE     = YPSOL - YPI
                YFHE    = 0.21*YPFUS/YPSOL
                YXEI    = YPE/YPI
        
	call B2EM
     >                 (YFHE,YFF,YSDT,YPSOL,YGDT,YXEI,
     >                  YQPK,YNES,YGDTS,YNHES,YGHES,
     >                  YTES,YTIS,YNCS,YGCT,YPIMP,YCAR)
     
     
                YENWM    =.001*YTIS         ! En_S keV
                YEHE    =.001*min(YTIS,60.) !EHe_s keV
                YGDT    = YGDT*54.          !1Pam3/s -> 10^19atoms/s 
                                            !           (molecules)
                YGHES   = YGHES*27.         !1Pam3/s -> 10^19atoms/s
                YGDTS   = YGDTS*54.         !1Pam3/s -> 10^19atoms/s 
                                            !           (molecules)	
		
                YTES    = .001*YTES         !eV -> keV
                YTIS    = .001*YTIS         !eV -> keV
                YGCT    = YGCT*27.          !1Pam3/s -> 10^19atoms/s
		
		YDELTA = ABS((YGDTS-YGDTSold)/max(1.e-5,YGDTS))
		if (YDELTA.gt.Yeps) then
		    YGDTSold = YGDTS
		    go to 100
		endif
	
	YGDTSold = YGDTS
	YVINTREC=0.
	YVINTION=0.
	do j=1,NA1
	        INCLUDE	'fml/svie'
	        INCLUDE	'fml/svii'
                YSNI(j) = (SVIE+SVII*NI(J)/NE(J))*NN(J)*NE(j)
	        INCLUDE	'fml/svrec'
	        YSNR(j)	=-SVREC*NI(J)*NE(j)
	enddo
        YVINTREC = vint(YSNR,ROC)
	YVINTION = vint(YSNI,ROC)

	YNNWM=(YGDTS-YVINTREC)/YVINTION 
	
cp
         return
c        open(1,file='B2EOUT.txt')
         write(*,*) YFHE,'	= fHe (= 0.21*PFUS/PSOL) 	input'
         write(*,*) YFF,'	= ff (= 1+.18 GCORE/GDT,tot) 	input'
         write(*,*) YSDT,'	= SDT m3/s 			input' 
         write(*,*) YPSOL,'	= PSOL, MW 			input'
         write(*,*) YGDT,'	= GDT, 10^19/s			input'
c======================================================================
         write(*,*) YQPK,'	= qpk, MW/m2 out		out'
         write(*,*) YNES,'	= ne_s, 10^19/m3		out'
         write(*,*) YGDTS,'	= GDT_n_s, 10^19/s		out'
         write(*,*) YNNWM,'	= NNWM_s, 10^19/m^-3            out'
         write(*,*) YNHES,'	= nHe_s, 10^19/m3		out'
         write(*,*) YGHES,'	= GHe_n_s, 10^19/s		out'
         write(*,*) YTES,'	= Te_s, keV 			out'
         write(*,*) YTIS,'	= Ti_s, keV 			out'
         write(*,*) YNCS,'	= nC_s, 10^19/m3		out'
         write(*,*) YGCT,'	= GC_t,  10^19/s		out'
         write(*,*) YPIMP,'	= Pimp, MW			out'
         write(*,*) YENWM,'	= En_s, keV			out'
         write(*,*) YEHE,'	= EHe_s, keV			out'
c        close(1)
        return
        end
c======================================================================
        subroutine B2EM
     >          (YFHE,YFF,YSDT,YPSOL,YGDT,YXEI,
     > YQPK,YNES,YGDTS,YNHES,YGHES,YTES,YTIS,YNCS,YGCT,YPIMP,
     >  YCAR)
c-------------------------------------------------------------------------------
c   Kukushkin e.a., NF 43	(2003) 716-723 and EPS (2003)
c	if: Ycar=0. then: metallic wall	else: full carbon wall 
        implicit none
c       Input	
	double precision YFHE,YFF,YSDT,YPSOL,YGDT,YXEI,YCAR
c	Output
	double precision YQPK,YNES,YGDTS,YNHES,YGHES
	double precision YTES,YTIS,YNCS,YGCT,YPIMP
c       Intrinsic variables
        double precision PSOL,SDT,GDT,YGDTC,YMUCR,PDT
	double precision YQPK0,YNS0,yGDT0,YTSE0,YTSI0
	double precision YGC0,YGC1,YGC2,YGC3
	double precision YNC0,YNC1,YNC2,YNC3
	double precision YP0,YP1,YP2,YP3
	
        PSOL    = YPSOL/100.
        SDT     = YSDT/20.
        GDT     = YGDT/124.
        YGDTC   = 124.* YFF**2*PSOL**.87*SDT    !GDT from Table 1
        YMUCR   = 1.1   *1.2  ! 20% higher from Kukushkin 
                              !                 et al EPS 2003?
        if(YCAR.eq.0.) then
c= metallic wall
                YQPK0   =7.55
                YNS0    =3.89
                YGDT0   =16.4
                YTSE0   =162.
                YTSI0   =270.

                YGC0    =342.
                YGC1    =.66
                YGC2    =-.11
                YGC3    =.18
                YNC0    =2.3e-2
                YNC1    =1.2
                YNC2    =.6
                YNC3    =-.5
                YP0     =39.
                YP1     =-.075
                YP2     =.15
                YP3     =1.42
        else
c= carbon wall
                YQPK0   =5.36
                YNS0    =3.
                YGDT0   =12.3
                YTSE0   =178.
                YTSI0   =324.

                YGC0    =181.
                YGC1    =.44
                YGC2    =.33
                YGC3    =.37
                YNC0    =7.6e-2
                YNC1    =.78
                YNC2    =1.44
                YNC3    =-.14
                YP0     =59.
                YP1     =.36
                YP2     =-.72
                YP3     =1.04
        endif

        if(YGDT.ge.(YMUCR*YGDTC)) then
	write(*,*)'>>> 
     >  B2EIT: parameters are beyond the critical point'
        write(*,*) 'Gas puff decreasing is suggested!'
	goto 20
c Saturation
           YQPK  = YQPK0*YFF**(-1.7)*PSOL**1.26
           YNES  = YNS0*YFF**1.25*PSOL*.55*YXEI**.05
           YGDTS = YGDT0*YFF**(-2.5)*SDT*.3
           YNHES = 3.06e-2*YFHE/YFF**5/SDT*PSOL**.7*YXEI**(-.1)
           YGHES = 0.512*YFHE/YFF**5.42/SDT*PSOL**.52
           YTES  = YTSE0*YFF**(-.4)*SDT**(-.02)*PSOL**.32*YXEI**.049
           YTIS  = YTSI0*YFF**(-.9)*SDT**(-.04)*PSOL**.36*YXEI**(-.115)
           YNCS  = YNC0*YFF**3.*PSOL**.54*YXEI**(-.13)
           YGCT  = YGC0*YFF**1.21*PSOL**0.75
           YPIMP = YP0*PSOL**1.35
           YGDT  = YGDTC
        else
 20        PDT   = YGDT/YSDT/6.2
           YQPK  = YQPK0*PSOL**2.*PDT**(-.85)
           YNES  = YNS0*YFF**.53*PSOL**.24*YXEI**.05*PDT**.36
           YGDTS = YGDT0*YFF**(-3.)*SDT**.3*PSOL**(-.22)*PDT**.25
           YNHES = 3.06e-2*YFHE/(YFF*SDT)*PSOL**2.44*YXEI**(-.1)
     >                    /PDT**2.
           YGHES = 0.512*YFHE/(YFF*SDT)*PSOL**2.44/PDT**2.21
        YTES = YTSE0*YFF**(-.06)*SDT**(-.02)*PSOL**.47*YXEI**.05
     >                    *PDT**(-.17)
        YTIS = YTSI0*YFF**(-.32)*SDT**(-.04)*PSOL**.61*YXEI**(-.116)
     >                    *PDT**(-.29)
            YNCS = YNC0*PDT**YNC1*YFF**YNC2*PSOL**YNC3*YXEI**(-.13)
            YGCT = YGC0*PDT**YGC1*YFF**YGC2*PSOL**YGC3
            YPIMP = YP0*PDT**YP1*YFF**YP2*PSOL**YP3
        endif
        return
        end

