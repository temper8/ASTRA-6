       subroutine MOD(period,duty_off,percpower_off,switch)
c Square modulation; duty cycle and "off" power in %
       implicit none
       include 'for/parameter.inc'
       include 'for/const.inc'
       include 'for/status.inc'
       double precision duty,duty_off,phase,switch,percpower_off,period
       if (duty_off .eq. 0) then
         switch=1
         return
         endif
       duty=1-duty_off
       switch=percpower_off
c phase = (t - n*period)/period  .le. 1
       phase = time/period -INT(time/period)
       if (phase .le. duty)  switch=1. 
       RETURN
       END
       
