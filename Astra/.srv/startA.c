#include 	<stdio.h>
#include 	<stdlib.h>
#include 	<string.h>
int  main(argc, argv, envir)
int  argc;
char *argv[], *envir[];
{
char	stri1[80];
FILE	*fp1, *fp2;

if(argc >= 4)
  {	printf("Start error: too many arguments\n");
	exit(0);
   }
if(argc <= 1)
  {	printf("Start error: no arguments\n");
	exit(0);
   }

fp1 = fopen("tmp/astra.log","r");
fp2 = fopen("tmp/Astart.tmp","w");

++argv;
while(fgets(stri1,80,fp1))
  {
	if(strncmp(*argv,stri1,3)) fputs(stri1,fp2);
	else
	{
	if(!strncmp(*argv,"WHO",3))
		{fputs("WHOME:",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"ARO",3))
		{fputs("AROOT:",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"AWD",3))
		{fputs("AWD:  ",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"VAR",3))
		{fputs("VAR:  ",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"MOD",3))
		{fputs("MOD:  ",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"MACH",4))
		{fputs("MACH: ",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"SHOT",4))
		{fputs("SHOT: ",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"XFI",3))
		{fputs("XFILE:",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"RUNN",4))
		{fputs("RUNN: ",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"TIM",3))
		{fputs("TIME: ",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"TEN",3))
		{fputs("TEND: ",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"RSN",3) || !strncmp(*argv,"PRA",3))
		{fputs("RSNAM:",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"PRO",3))
		{fputs("PROFT:",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"RTY",3))
		{fputs("RTYPE:",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"DAT",3))
		{fputs("DAT:  ",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"RES",3))
		{fputs("RES:  ",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"TSK",3))
		{fputs("TSK:  ",fp2);
		goto	next;
		}
	if(!strncmp(*argv,"TOT",3))
		{fputs("TOT:  ",fp2);
		goto	next;
		}
	printf("Unrecognized 2nd argument \n");
 	fputs(stri1,fp2);
	goto out;
	next:
	if(argc == 3)	fputs(*++argv,fp2);
		fputs("\n",fp2);
	out: ;
	}
   }
fclose(fp1);
fclose(fp2);
system("mv tmp/Astart.tmp tmp/astra.log");
}
