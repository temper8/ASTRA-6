	subroutine	BALST
	double precision aBal(256)
	integer kBal(256),nBal
	integer getbalst
	external getbalst
	j=getbalst(aBal,kBal,nBal)
C	write(*,'(a,i3)')'nBal=',nBal
	if(j .eq. 0) then
	   do j=1,nBal
	      if(kBal(j) .ne. 0) write(*,'(a,1p1e10.3,i2)')'Unstable a'
     >			,aBal(j),kBal(j)
	   enddo
	endif
        end
