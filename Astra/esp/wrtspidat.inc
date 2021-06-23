	do j=1,na1
	eqff1(j)=eqff(j)
	enddo
	EQFF1(NA1)=EQFF(NA)
!+(EQFF(NA)-EQFF(NA-1))/HRO*HROA
	open(1,file='dat/spidat.dat')
	write(1,*) 'jNEQL,jnteta,jNBND'
	write(1,*) jNEQL,jnteta,jNBND
	write(1,*) '(rzbnd(j),j=1,2*jNBND)'
	write(1,*) (rzbnd(j),j=1,2*jNBND)
	write(1,*) 'na1'
	write(1,*) na1
	write(1,*) '(eqpf(j),j=1,na1)'	 
	write(1,*) (eqpf(j),j=1,na1)	 
	write(1,*) '(eqff(j),j=1,na1)'	 
	write(1,*) (eqff1(j),j=1,na1)	 
	write(1,*) '(fp(j),j=1,na1)'	 
	write(1,*) (fp(j),j=1,na1)	 
	write(1,*) '(rho(j),j=1,na1)'	 
	write(1,*) (rho(j),j=1,na1)	 
	write(1,*) 'ipl,rtor,btor,roc,jnstep,'
	!write(1,*) ipl,yrtor,ybtor,roc,jnstep
	write(1,*) ipl,rtor,btor,roc,jnstep
	write(1,*) '(fp(j),j=1,na1)'	 
	write(1,*) (fp(j),j=1,na1)	 
	write(1,*) '(g11(j),j=1,na1)'	 
	write(1,*) (g11(j),j=1,na1)	 
	write(1,*) '(g22(j),j=1,na1)'	 
	write(1,*) (g22(j),j=1,na1)	 
	write(1,*) '(g33(j),j=1,na1)'	 
	write(1,*) (g33(j),j=1,na1)	 
	write(1,*) '(vr(j),j=1,na1)'	 
	write(1,*) (vr(j),j=1,na1)	 
	write(1,*) '(vrs(j),j=1,na1)'	 
	write(1,*) (vrs(j),j=1,na1)	 
	write(1,*) '(slat(j),j=1,na1)'	 
	write(1,*) (slat(j),j=1,na1)	 
	write(1,*) '(gradro(j),j=1,na1)'	 
	write(1,*) (gradro(j),j=1,na1)	 
	write(1,*) '(mu(j),j=1,na1)'	 
	write(1,*) (mu(j),j=1,na1)	 
	write(1,*) '(ipol(j),j=1,na1)'	 
	write(1,*) (ipol(j),j=1,na1)	 
!aai 030108
	write(1,*) '(cc(j),j=1,na1)'	 
	write(1,*) (cc(j),j=1,na1)
	write(1,*) '(Te(j),j=1,na1)'	 
	write(1,*) (Te(j),j=1,na1)
	write(1,*) '(cubs(j),j=1,na1)'	 
	write(1,*) (cubs(j),j=1,na1)
	write(1,*) '(cd(j),j=1,na1)'	 
	write(1,*) (cd(j),j=1,na1)
!aai 030108
	close(1)