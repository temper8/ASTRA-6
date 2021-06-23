/* Transfer two time strings "hh:mm:ss" into time difference  */
#include	<stdio.h>
#include	<math.h>
#include	<stdlib.h>
#include	<string.h>
int   main(argc, argv, envir)
int   argc;
char *argv[], *envir[];
{char	buf[1];	int	j;	long	n, m;
 if (argc >= 4)
 {	
     printf("Start error: too many arguments\n");
     exit(0);
 }
 if (argc <= 1)
 {
     printf("Start error: no arguments\n");
     exit(0);
 }
 ++argv;	
 for (j=n=0; j<=2; j++)
 {
     strncpy(buf,*argv+3*j,2);
     n+=atoi(buf); n*=60;
 }
 ++argv;
 for (j=m=0; j<=2; j++)
 {
     strncpy(buf,*argv+3*j,2); 
     m+=atoi(buf); 
     m*=60;
 }
 m=m/60;
 n=n/60;
 j=(m-n)/3600;
 n=(m-n)%3600;
 m=n/60;
 n=n%60;
 printf("%2.2d:%2.2d:%2.2d\n",j,m,n);
}
