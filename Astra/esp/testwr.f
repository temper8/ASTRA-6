	SUBROUTINE testwr
      open(1,file='dat/test.wr')	
        write(1,*) 'write check'      
      close(1)
      	RETURN
	END
