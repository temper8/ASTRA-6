/**************************************************************************/
/* This block is used in 
 */
int a_stop__ ()
{   a_stop_(); }
int a_stop ()
{   a_stop_(); }
int a_stop_()
{  char stri[64];
//   a_postmortem_();
#ifdef IPCACTIVE
   freeshm();
#endif
   system("rm tmp/*.ppm >& /dev/null");
   /* system("if ( -e ${HOME}/bin/cln ) ${HOME}/bin/cln"); */
   exit(0);
}
/* The function returns the total CPU time [sec] spent by the calling process.
      It adds the CPU time spent between two successive calls to the argument.
   The  function "times" returns the number of clock ticks that have elapsed
      since the moment the system was booted. 
   The  "tms_utime"  field contains the CPU time spent executing instructions
      of the calling process.
   The  "tms_stime"  field contains the CPU time spent in the system while 
      executing tasks on behalf of the calling process.                     */
double swatch_(secs)		double	*secs;
{  double	runsec;
   clock_t	  cpu_time, run_time;
   static clock_t time0=0, prev_time;
   static double  secs_per_tick;
   struct   tms	  buf;
   if (time0 == -1) { return -1.; }		/* Overflow range of clock_t*/
   if (time0 ==  0) { 				/* Set time0 at start */
     secs_per_tick = 1./sysconf(_SC_CLK_TCK);
     time0 = times(&buf);
     prev_time = buf.tms_utime+buf.tms_stime;
     return 0.;
   }
   run_time = times(&buf)-time0;		/* Set time difference */
   cpu_time = buf.tms_utime+buf.tms_stime;
   *secs += (cpu_time-prev_time)*secs_per_tick;
   /* printf("Process: %f  %f\n", cpu_time*secs_per_tick, */
   /* (cpu_time-prev_time)*secs_per_tick); */
   prev_time = cpu_time;	
   runsec = run_time*secs_per_tick;   /* printf("Process: %f\n", runsec); */
   return runsec;
}
double swatch (secs)		double	*secs;
{  double	runsec;
   runsec = swatch_(secs);
   return runsec;
}
/****************************************************************************/
