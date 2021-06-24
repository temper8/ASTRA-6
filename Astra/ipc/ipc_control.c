#include	<unistd.h>
#include	<time.h>
#include	<string.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<unistd.h>
#include	<errno.h>
#include	<sys/times.h>
#include	<sys/types.h>
#include	<sys/ipc.h>
#include	<sys/sem.h>
#include	<sys/shm.h>
double	swatch (double*);		double	swatch_(double*);
//char  *AWD, *PLATF, *MOD, *DATA;
char   AWD[40], MOD[20], PLATF[20], DATA[20];
//extern char *A_ipc_file;          /* The name is defined in ipc_lib.c */
char A_ipc_file[64];
char ASTRA_task[132];
const  char *A_log_file = "./tmp/astra.log";
key_t  my_key;
/*#if defined(_SEM_SEMUN_UNDEFINED) */
union semun {
  int val;                    /* value for SETVAL */
  struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
  unsigned short int *array;  /* array for GETALL, SETALL */
  struct seminfo *__buf;      /* buffer for IPC_INFO */
};
/*#endif */
pid_t A_PID = 0;
int A_NB1 = 0;
int A_SemID = 0;
int A_Nsems = 0;       /* the number of semaphores */
int A_ShmNum = -1;
#define NC1 A_NB1
#define	A_ShmShift 2
#define	A_Nsemx	   20
char A_ChNa[A_ShmShift+A_Nsemx][64]; /* Child process name (not used) */
int  A_ChID[A_ShmShift+A_Nsemx] = {0,0,0};        /* Child process ID */
int  A_ShmL[A_ShmShift+A_Nsemx];	/* Child Shmem segment length */
int  A_ShmID[A_ShmShift+A_Nsemx] = {0,0};
void *A_ShmAdr[A_ShmShift+A_Nsemx];
            /* {sem_num,sem_op,sem_flag}; */
// struct sembuf buf0 = {0, 0, SEM_UNDO};
struct sembuf buf0 = {0, 0, ~SEM_UNDO&~IPC_NOWAIT};
ushort semarr[A_Nsemx], semar0[A_Nsemx];  /* Is used in "checkeach()" */
/**********************************************************************/
/* Check existence of executable files listed in subs
   Reads tmp/astra.log and fills external variables AWD, PLATF, MOD, DATA
 */
/* *Nsub - total number of files_names/strings in subs,
   *Lstr - length of an element of the character ARRAY "subs",
           maximum length of the subprocess_name,
   *subs - character ARRAY, described in a calling Fortran routine as
           character*(*Lstr) ARRAY(max_length)
	   each element includes a name of external process to be called
 */
int checkexec (int* Nsub, int *Lstr, char *subs)
{   checkexec_(Nsub, Lstr, subs); }
int checkexec_(int* Nsub, int *Lstr, char *subs)
{
  char stri[164], name[32], path[64];
  int  j, i;
  if (A_Nsems != 0)	return(A_Nsems);	/* Do check only once */
  /* Check existence of subprocess executable files */
  /* Read tmp/astra.log and store run info */
  FILE *A_LOG;
  A_LOG = fopen(A_log_file,"r");
  if (!A_LOG) 
  { printf("Cannot open Astra log file: \"%s\"\n",A_log_file);
    exit(0);
  }
  while (fgets(stri,7,A_LOG) != NULL)
  { if ( strncmp("AWD:  ",stri,6) == 0 ) fscanf(A_LOG,"%s",AWD);
    if ( strncmp("PROFT:",stri,6) == 0 ) fscanf(A_LOG,"%s",PLATF);
    if ( strncmp("MOD:  ",stri,6) == 0 ) fscanf(A_LOG,"%s",MOD);
    if ( strncmp("VAR:  ",stri,6) == 0 ) fscanf(A_LOG,"%s",DATA);
  }
  fclose(A_LOG);
  strcpy(stri,PLATF);
  i = strchr(stri,'/')-stri;
  PLATF[i+1] = '\0';

  for(j=0; j < *Nsub; j++)
  {  if ( strlen(&subs[*Lstr*j]) == 0 )   goto	Error1;
     strcpy(path,&subs[*Lstr*j]);
//   printf("Input: \"%s\"\n",path);
     if ( strchr(path,'~') != NULL )
     {  if ( getenv("HOME") == NULL ) goto Error2;
        strcpy(stri,getenv("HOME"));
        strcat(stri,&path[1]);
        strcpy(path,stri);
        i = strrchr(path,'/')-&path[0];
        strcpy(name,&path[i+1]);
        if ( strlen(path) > 62 ) goto Error3;
        path[i+1] = '\0';
     }
     else if ( strrchr(path,'/') != NULL ) 
     {  i = strrchr(path,'/')-&path[0];
        path[i+1] = '\0';
        strcpy(name,&subs[*Lstr*j+i+1]);
     }
     else
     {  strcpy(name,path);
        path[0] = '\0';
     }
     i = strlen(path);
     /*  printf("path  = [%s]  %d\n", &path[0], i);
	 printf("PLATF = [%s]  %d\n", &PLATF[0], strlen(PLATF));
	 printf("name  = [%s]  %d\n", &name[0], strlen(name));       */
     if ( i == 0 )
     {  strcpy(stri,"test -x .tsk/");
        strcat(stri,PLATF);
        strcat(stri,name);
     }
     else
     {
        strcpy(stri,"test -x ");
        strcat(stri,path);
        strcat(stri,name);
//   printf("Check: [%s]\n%d\n",stri,system(stri));

/*            Slow option:
        strcpy(stri,"cd ");		strcat(stri,path);
        strcat(stri," ; test -x ");	strcat(stri,name);
        printf("Check: [%s]\n");
*/
     }
/*     printf("The executable file \"%s\"\n",stri);	*/

     if ( system(stri) == 0 )
     {
	++A_Nsems;
     }
     else
     {  if ( i == 0 )
        printf("The executable file \"%s\" (#%d) does not exist\n",
                  &stri[8],j+1);
        else
        printf("The executable file \"%s%s\" (#%d) does not exist\n",
                  path,name,j+1);
        exit(j);
     }  /* printf("%d child processes found\n",A_Nsems); */
  }
  return(++A_Nsems);
  Error1:
  printf(" >>> Xroutine call string (#%d) \"%s\" error\n",j+1,&subs[*Lstr*j]);
  exit(0);
  Error2:
  printf(" >>> Xroutine call string \"%s\" error:\n",&subs[*Lstr*j]);
  printf(" >>> Symbol \"~\" is not allowed.\n");
  exit(0);
  Error3:
  printf(" >>> Xroutine call string \"%s\" error:\n",&subs[*Lstr*j]);
  printf(" >>> Absolute path is too long.\n");
  exit(0);
}
/**********************************************************************/
/* "call initipc(NB1)" is placed in init.inc
   Get PID and key for the Astra main process
   Create and initialize a set of A_Nsems semaphores
   Assign NB1 (= *Ngrid) to A_NB1 (alias NC1)
   Allocate two shared memory segments for Astra datasets 
 */
int initipc (int* Ngrid)
{   initipc_(Ngrid); }
int initipc_(int* Ngrid)
{
  int   l, is=0, ds, j, *k;
  FILE  *A_PDF;
  char	hostname[32];
  size_t namlen;
  time_t hold_time;
  static union semun Mysemun;

  if (A_Nsems == 0)	return(0); /* Remove this line if ESC is enabled */
  if (A_NB1 != 0)	return(0);	/* Initialize only once */

/* Collecting data */
  A_PID = getpid();		/* printf("Apid = %d\n",A_PID); */
   /* Define the absolute path name of Astra executable ASTRA_task */
  strcpy(ASTRA_task,AWD);     /* printf("ASTRA_task = %s\n",ASTRA_task); */
  strcat(ASTRA_task,".tsk/"); /* printf("ASTRA_task = %s\n",ASTRA_task); */
  strcat(ASTRA_task,PLATF);   /* printf("ASTRA_task = %s\n",ASTRA_task); */
  strcat(ASTRA_task,MOD);     /* printf("ASTRA_task = %s\n",ASTRA_task); */
  strcat(ASTRA_task,".exe");  /* printf("ASTRA_task = %s\n",ASTRA_task); */
  my_key = ftok( ASTRA_task, (int)A_PID);    /* Get System V IPC key */
  /* printf("MyKey = %d  %d\n",my_key,(int)my_key); */
  if (A_Nsems != 0)
  {			/* Create a set of A_Nsems semaphores */
      A_SemID = semget(my_key, A_Nsems, 0660|IPC_CREAT);
      /* Initialize all semaphores in the set A_SemID as {0,0,...} */
      Mysemun.val = 0;	for(j=0; j < A_Nsems; j++)
      {   semctl(A_SemID, j, SETVAL, Mysemun);
      //  WhatSem();	printf("\n");
      }   /* printf("\nA set of %d semaphores is created,  SemID = %d\n",
	                                            A_Nsems,A_SemID); */
  }		/* Write file tmp/astra.ipc */
  /* Form absolute path for A_ipc_file */
  getcwd(A_ipc_file,(size_t)64);
  strcat(A_ipc_file,"/tmp/astra.ipc");
  A_PDF = fopen(A_ipc_file,"w");
  if (!A_PDF) 
  { printf("Cannot open Astra IPC file: \"%s\"\n",A_ipc_file);
    exit(0);
  }
  fprintf(A_PDF," Astra task:  \"%s\"\n",ASTRA_task);
  fprintf(A_PDF," Astra files:  \"%s\",  \"%s\"\n",DATA,MOD);
  gethostname(hostname,(size_t)32);     
  // printf("hostname = %s\n",hostname);
  hold_time=time(NULL);
  fprintf(A_PDF," Astra@%s started on:  %s",&hostname,ctime(&hold_time));
//  getdomainname(hostname,(size_t)32);  printf("domainname = %s\n",hostname);
  /* The next line will be used if ESC is enabled */
  if (A_Nsems == 0)	{   fclose(A_PDF);	return(0);  }
  fprintf(A_PDF," Astra(main):  PID = %d,  SemID = %d\n",
	  (int)A_PID,A_SemID);
  A_NB1 = *Ngrid;
#include "A_vars.h"
  l = sizeof(struct A_vars);
  AllocateShmem(l);
  fprintf(A_PDF," Astra_inout: ShmID(A_vars):%12d%12d\n", A_ShmID[0], l);
#include "A_arrs.h"
  l = sizeof(struct A_arrs);
  AllocateShmem(l);
  fprintf(A_PDF," Astra_inout: ShmID(A_arrs):%12d%12d\n", A_ShmID[1], l);
  fprintf(A_PDF,"        PID      ShmID       ShmSize     Process\n");
  fclose(A_PDF);

  AVARS = (struct A_vars *)A_ShmAdr[0];
  AVARS->c_size[0] = sizeof(struct A_vars);
  AVARS->c_size[1] = sizeof(int);
  AVARS->c_size[2] = sizeof(double);
  AVARS->c_size[3] = sizeof(key_t);
  AVARS->c_size[4] = sizeof(pid_t);
  AVARS->CheckWord = 314159265;

  AARRS = (struct A_arrs *)A_ShmAdr[1];
  AARRS->Size = sizeof(struct A_arrs);
  AARRS->CheckWord = 314159265;

  /*
  printf("l = %d, %d, %d, %d, %d, %d, %d\n",
	 sizeof(struct A_vars), 
         AVARS->c_size[0],
         AVARS->c_size[1],
         AVARS->c_size[2],
         AVARS->c_size[3],
         AVARS->c_size[4],
         AVARS->CheckWord
        );
  printf("  sizeof(key_t) = %d,   ",sizeof(key_t));
  printf("  sizeof(pid_t) = %d\n",  sizeof(pid_t));
  printf("  sizeof(int)   = %d,   ",sizeof(int));
  printf("  sizeof(double) = %d\n", sizeof(double));
  printf("  sizeof(*AVARS) = %d,  %d\n",sizeof(*AVARS),(AVARS->c_size[0]));
  printf("\n  AARRS:   sizeof = %d = %d = %d\n",
       sizeof(*AARRS), sizeof(struct A_arrs), AARRS->Size);
  printf("CheckWord = %d\n", AARRS->CheckWord);
  */
}
/*---------------------------------------------------------------------*/
/* Get ShMemIDs for Astra datasets (const.inc) and (status.inc) */
int AllocateShmem (int l)
{
  if (++A_ShmNum > A_ShmShift+A_Nsemx)
  {  printf(">>> ERROR >>> Too many shared memory segments requested\n");
     a_stop();
  }
  /* Allocate a shared memory segment starting at A_ShmAdr[A_ShmNum] */
  A_ShmID[A_ShmNum] = shmget((key_t)(my_key+A_ShmNum),l,0660|IPC_CREAT);
  /* Attach shared memory to the process */
  A_ShmAdr[A_ShmNum] = shmat(A_ShmID[A_ShmNum], NULL, 0);
/*  printf("SegNum = %d, Shmem ID = %d,  Key_t = %d,  Adr = %d\n", 
    A_ShmNum, A_ShmID[A_ShmNum], my_key+A_ShmNum,A_ShmAdr[A_ShmNum]); */
  return(0);
}
/*---------------------------------------------------------------------*/
/* Set (lock) the primary semaphore to -(Number_of_processes)
   Launch subordinate processes
   Those in turn do the following
     (1) get all PIDs, Keys, ShMems
     (2) create own ShMem segment
     (3) append file tmp/astra.ipc
     (4) increment primary semaphore
     (5) wait until the primary opens track
*/
int inikids (int* Nsub, int *Lstr, char *subs)
{   inikids_(Nsub, Lstr, subs); }
int inikids_(int* Nsub, int *Lstr, char *subs)
{
  if (A_ChID[A_ShmShift] != 0)	return(0);	/* Initialize only once */
//   printf("inikids called Nsub = %d\n",*Nsub);
  char stri[164], name[32], path[64];
  int  j, i;
// printf("inikids called Nsub = %d,   A_Nsems = %d\n",*Nsub, A_Nsems);
  if (A_Nsems <= *Nsub)
  {  printf(" >>> ERROR >>> Incomplete semaphore set\n");
     return(1);
  }
// printf("inikids called Nsub = %d,   A_Nsems = %d\n",*Nsub, A_Nsems);
  for(j=0; j < *Nsub; j++) 
  {  if ( strlen(&subs[*Lstr*j]) == 0 )   goto	Error1;
     strcpy(path,&subs[*Lstr*j]);
//   printf("Input: \"%s\"\n",path);
     if ( strchr(path,'~') != NULL )
     {  strcpy(stri,getenv("HOME"));
        strcat(stri,&path[1]);
        strcpy(path,stri);
        i = strrchr(path,'/')-&path[0];
        strcpy(name,&path[i+1]);
        path[i+1] = '\0';
     }
     else if ( strrchr(path,'/') != NULL ) 
     {  i = strrchr(path,'/')-&path[0];
        path[i+1] = '\0';
        strcpy(name,&subs[*Lstr*j+i+1]);
     }
     else
     {  strcpy(name,path);
        path[0] = '\0';
     }
     i = strlen(path);
/*   printf("\npath = [%s]  %d\n", &path[0], i);
     printf("name = [%s]  %d\n", &name[0], strlen(name));	*/
     if ( i == 0 )
     {	sprintf(stri,"%s%s%s %s %d %d %d &",
	  "./.tsk/", PLATF, name, ASTRA_task, A_PID, (int)my_key, j+1);
	/* printf("\"%s\" launched.\n",stri);	     WhatSem(); */
        i = system(stri);
     }
     else
     {
/* The following OS command ("cd xxx") is time consuming.
   Below it is replaced with C command chdir();
	sprintf(stri,"cd %s >/dev/null; ./%s %s %d %d %d &",
		path, name, ASTRA_task, A_PID, (int)my_key, j+1);
   Note! chdir() does not recognize ~ as home directory.
*/
/* Get and print current working directory:
      getcwd(stri,164);   printf("%s\n", stri);		*/
      chdir(path); 
      sprintf(stri,"./%s %s %d %d %d &",
		name, ASTRA_task, A_PID, (int)my_key, j+1);
      i = system(stri);
      chdir(AWD);
//    getcwd(stri,164);  printf("%s\n", stri);  printf("\"%s\"\n", AWD);
     }
// printf("Process No.%d is launched as:  \"%s\"\n",j+1, stri);      

/* Wait until proc No.(j+1) opens the PRIMARY semaphore (sem_num=0)    */
/*                          incrementing its initial value (0) by 1    */
struct sembuf bufj = {0,-1,~IPC_NOWAIT}; /* {sem_num, sem_op, sem_flag}*/
semop(A_SemID, &bufj, 1);
/* printf("\nProcess No.%d: is started\n",j+1);		WhatSem();     */
     if ( i == -1 ) return(j+1);
  }
  /* All secondaries all launched and the file A_ipc_file is completed */
  /* Now the data from A_ipc_file have to be retrieved by the main     */
  if ( read_aipc(Nsub, Lstr, subs) ) goto Error2;
  return(0);

 Error1:
  printf("Error in input SBP string [%s]\n",&subs[*Lstr*j]);
  return(j);
 Error2:
  printf("Error in input data interpretation (function read_aipc)\n");
  exit(j);
}
/*---------------------------------------------------------------------*/
void whatsem()
{ WhatSem(); }
void whatsem_()		/* Callable from FORTRAN as */
{ WhatSem(); }		/*      call whatsem()      */
int WhatSem()
{   int j;
  if (A_Nsems == 0)	return(0);
  ushort semarray[A_Nsems];
  union semun Mysemun;
  Mysemun.array = &semarray[0]; 
  semctl(A_SemID, 0, GETALL, Mysemun);
//  printf(" Semaphore set of %d\n",Mysemun.buf.sem_nsems);
  printf(" Semaphore set = {");
  for (j=0; j < A_Nsems-1; j++) printf("%d,",semarray[j]);
  printf("%d}\n",semarray[A_Nsems-1]);
  return (0);
}
/*--------------------- Check if ipc is activated --------------------*/
int ifipc ()
  { ifipc_(); }
int ifipc_()
{ if ( A_NB1 == 0 )  return(0);   return(1); 
}
/*--------------------- Unlock subprocess ----------------------------*/
int letsbp (int* n)
  { letsbp_(n); }
int letsbp_(int* n)
{
  auto struct sembuf bufN = {*n, 1, IPC_NOWAIT};
  --buf0.sem_op;   /* Each call decrements Sem0 value by 1 */
  /* Increment semval # sem_num=*n, Open subprocess */
  semop(A_SemID, &bufN, 1);
  /* printf("Process %d let,  buf0.sem_op = %d, Total %d\n",
        *n,buf0.sem_op,A_Nsems);			WhatSem();     */
}
/*---------------------------------------------------------------------*/
/* Compare Sem0 value with buf0.sem_op ( == -A_Nsems) and
	wait until all subprocesses increment Sem0 by 1 
        so that Sem0 reaches value A_Nsems                       ------*/
int wait4all_()
{   wait4all(); } //{   checkeach(); } 
int wait4all()
{ int j;
  if (A_ShmNum < 0) return(0);	/* Do check only after initialization   */
//  static struct timespec timeout = {0,1000000};   /* timeout = 1 msec */
  static struct timespec timeout = {0,100000000};   /* timeout = .1 sec */

 MinorLoop:
  {/* Here the primary process can do limited actions e.g. analyze keys */
      (void) AstraEvent();
  }
  /* Go on if [Sem0value+buf0.sem_op==0], goto Minorloop after timeout  */
  /* Here: buf0 = {(ushort_t) buf0.sem_num = 0 = Const,
                -1 <= (short) buf0.sem_op <= Number_of_active_subprocesses, 
		      (short) buf0.sem_flg = SEM_UNDO = Const)} 
     Note! Each call increments value of semadj by 1.
           Overflow occurs when the number of successive calls exceeds 32767
  */
  //  if (semop(A_SemID, &buf0, 1))
  if (j=semop(A_SemID, &buf0, 1))
    {                   /*  WhatSem();  printf("j = %d  %d\n",j,errno); */
    switch( errno )
    { case EAGAIN:			/* Timed out */
	/* printf("Timed out\n"); */
	  goto  MinorLoop;
      case EIDRM:
	  printf("The semaphore set was removed\n");
	  goto  Exit;
      case EINTR:
	  printf("While blocked in this call, the process caught a signal\n");
	  goto  Exit;
      case ERANGE:	printf("\n");		WhatSem();
	  printf("ERANGE error:  sem_op = %d,   SEMVMX = ?\n\n",buf0.sem_op);
	  goto  Exit;
      case EINVAL:	printf("EINVAL:\t");	goto Exit;
      case EFAULT:	printf("EFAULT:\t");	goto Exit;
      case E2BIG:	printf("E2BIG:\t");	goto Exit;
      case EACCES:	printf("EACCES:\t");	goto Exit;
      case EFBIG:	printf("EFBIG:\t");	goto Exit;
      case ENOMEM:	printf("ENOMEM:\t");	goto Exit;
      default:
	  printf(">>> PRIMARY >>> Unrecognised semtimedop error %d\n",errno);
	  goto Exit;
     }
  }                       	     /* printf("PRIMARY unlocked\n"); */
  buf0.sem_op = 0;
  return(0);
 Exit:
  printf(">>> PRIMARY >>> semtimedop error = %d\n",errno);
  a_stop();
}
/*---------------------------------------------------------------------*/
/*    Read IPC file after launching all processes and                  */
/* (1) fill arrays ofChild_process_IDs, Child_shmem_lengths/IDs,       */
/* (2) attach child process shmem segments to the main process memory. */
int read_aipc (int* Nsub, int *Lstr, char *subs)
  { read_aipc_(Nsub, Lstr, subs); }
int read_aipc_(int* Nsub, int *Lstr, char *subs)
{
  FILE *A_PDF;
  char stri[132], name[32];
  int  j, i, k, ID, ShmID, kS; 

  A_PDF = fopen(A_ipc_file,"r");
  if (!A_PDF) 
  { printf("Cannot open Astra IPC file: \"%s\"\n",A_ipc_file);
    exit(0);
  }
  for (j=0; j++ < 5+A_ShmShift;) fgets(stri, 132, A_PDF); /* Skip lines */
  i = A_ShmShift;
  while (EOF != fscanf(A_PDF,"%12d%12d%12d%s", &ID, &ShmID, &kS, stri) )
  {
    for (j = i-A_ShmShift; j < A_Nsems; j++)
    {
       if ( strrchr(&subs[*Lstr*j],'/') != NULL )
       {  k = strrchr(&subs[*Lstr*j],'/')-&subs[*Lstr*j];
          strcpy(name,&subs[*Lstr*j+k+1]);
       }
       else
       {  strcpy(name,&subs[*Lstr*j]);
       }
       if ( j == i-A_ShmShift ) break;  /* SBP name matches the record */
       else goto Out1;             /* Inconsistency in file A_ipc_file */
    }
    if (++A_ShmNum <= A_ShmShift+A_Nsemx)
    {
       A_ChID[A_ShmNum] = ID;
       A_ShmL[A_ShmNum] = kS;   /* A_ShmNum == ProcNumber+A_ShmShift-1 */ 
       A_ShmID[A_ShmNum] = ShmID;    /* Here (A_ShmNum == Shmem_index) */ 
       A_ShmAdr[A_ShmNum] = shmat(ShmID, NULL, 0);
       strcpy(A_ChNa[A_ShmNum],stri);
       /* fprintf(stdout,"%d %12d%12d%12d  %s\n\n",
          i-A_ShmShift, A_ShmL[i], A_ShmID[i], A_ShmAdr[i], A_ChNa[i]); */
    }
    else
    {
       printf(">>> ERROR >>> Too many shmem segments\n");
       a_stop();
    }
    i++;
  }		/* End while */
  fclose(A_PDF);
  if ( A_ShmNum+1-A_ShmShift == *Nsub )
      return(0);

 Out1:
  printf(" >>> File \"%s\" error >>> Missing SBP(%d) name: [%s]\n",
	 i+1-A_ShmShift,name);
  goto Out;
 Out2:
  printf(" >>> File \"%s\" error >>> Incomplete record.\n", A_ipc_file);
  goto Out;
 Out:
  printf("\n  File \"%s\" contents:\n", A_ipc_file);
  system("cat /home/grp/MAstra/tmp/astra.ipc");
  printf("\n  File processing:\n");
  printf("A_Nsems = %d,  A_ShmShift = %d,  A_ShmNum = %d,  Nsub = %d\n",
	  A_Nsems, A_ShmShift, A_ShmNum, *Nsub);
  for (j=A_ShmShift-1; j++ <= *Nsub; )
  {     fprintf(stdout,"%d %12d%12d%12d%12d  %s\n",
	  j-A_ShmShift+1, A_ChID[j], A_ShmL[j], A_ShmID[j],
          A_ShmAdr[j], A_ChNa[j]);
  }
  return(1);
}
/*---------------------------------------------------------------------*/
/****************** Get ID of ShMem for const.inc **********************/
/* Used in ipc/setipcdata.f90 only */
int getshmadr (const int* n, int* fpoint)
{   getshmadr_(fpoint); }
int getshmadr_(const int* n, int* fpoint)
{
  //#include "A_vars.h"
  if (A_Nsems == 0)	return(0);
  if (A_ShmNum < 0)	return(0);

  printf("\n Before:  &fpoint = %d,    fpoint = %d\n",fpoint, *fpoint);
  //  AVARS = (struct A_vars *)A_ShmAdr[*n];
  //  *fpoint = AVARS;
  //  *fpoint = (double *)A_ShmAdr[*n];
  //	printf("\n  C:     &cvars = %d,    cvars = %g = %g\n",
  //	 AVARS, AVARS->ab, *AVARS);

  switch( *n )
    { case 0:
	*fpoint = (double *)A_ShmAdr[*n];
	break;
      case 1:
	*fpoint = (int *)A_ShmAdr[*n];
	/* *fpoint = (struct A_arrs *)A_ShmAdr[*n]; */
	break;
      default:
        printf(">>> GetShmAdr >>> Illegal call:  n = %d\n",*n);
	return(1);
    }
  printf(" After:   &fpoint = %d,    fpoint = %d\n",fpoint, *fpoint);
  return(0);
/*
  AARRS = (struct A_arrs *)A_ShmAdr[1];
  *fpoint = AARRS;
*/
}
/****************** Get ID of ShMem for const.inc **********************/
/* Used in ipc/setipcdata.f90 only */
int getshmadrp (const int* n, int* fpoint)
{   getshmadrp_(fpoint); }
int getshmadrp_(const int* n, int* fpoint)
{
  //#include "A_vars.h"
  if (A_Nsems == 0)	return(0);
  if (A_ShmNum < 0)	return(0);

  printf("\n Before:  &fpoint = %d,    fpoint = %d\n",fpoint, *fpoint);
  //  AVARS = (struct A_vars *)A_ShmAdr[*n];
  //  *fpoint = AVARS;
  //  *fpoint = (double *)A_ShmAdr[*n];
  //	printf("\n  C:     &cvars = %d,    cvars = %g = %g\n",
  //	 AVARS, AVARS->ab, *AVARS);

  switch( *n )
    { case 0:
	fpoint = (double *)A_ShmAdr[*n];
	break;
      default:
        printf(">>> GetShmAdr >>> Illegal call:  n = %d\n",*n);
	return(1);
    }
  printf(" After:   &fpoint = %d,    fpoint = %d\n",fpoint, *fpoint);
  return(0);
/*
  AARRS = (struct A_arrs *)A_ShmAdr[1];
  *fpoint = AARRS;
*/
}
/****************** Get ID of ShMem for const.inc **********************/
int setvars (double* AB, double* T, double* H, int* N)
{   setvars_(AB, T, H, N); }
int setvars_(double* AB, double* T, double* H, int* N)
{ /* First active only after "initipc", i.e. after "init.inc" */
  int *I, j;
#include "A_vars.h"
  if (A_Nsems == 0) {
     printf(" >>> setvars >>> Illegal call: Semaphores are not created\n");
     return(0);     }
  if (A_ShmNum < 0)	return(0);

  AVARS = (struct A_vars *)A_ShmAdr[0];
  /*
  printf("AVARS->ab  = %g   %d\n",AVARS->ab, &(AVARS->ab));
  printf("AVARS->abc = %g   %d\n",AVARS->abc,&(AVARS->abc));
  printf("AVARS->tim = %g   %d\n",AVARS->time,&(AVARS->time));
  printf("AVARS->tau = %g   %d\n",AVARS->tau,&(AVARS->tau));
  printf("AVARS->na1 = %d   %d\n",AVARS->na1,&(AVARS->na1));
  printf("AVARS->nb1 = %d   %d\n",AVARS->nab,&(AVARS->nab));
  printf("AVARS->nb1 = %d   %d\n",AVARS->nb1,&(AVARS->nb1));
  printf("AVARS->ChW = %d   %d\n",AVARS->CheckWord,&(AVARS->CheckWord));
  */
  AVARS->ab    = *AB;
  AVARS->abc   = *(AB+1); 
  AVARS->aim1  = *(AB+2); 
  AVARS->aim2  = *(AB+3); 
  AVARS->aim3  = *(AB+4); 
  AVARS->amj   = *(AB+5); 
  AVARS->btor  = *(AB+7); 
  AVARS->elong = *(AB+8); 
  AVARS->encl  = *(AB+10);
  AVARS->enwm  = *(AB+11);
  AVARS->ipl   = *(AB+18);
  AVARS->nncl  = *(AB+20);
  AVARS->nnwm  = *(AB+21);
  AVARS->rtor  = *(AB+27);
  AVARS->shift = *(AB+28);
  AVARS->trian = *(AB+29);
  AVARS->updwn = *(AB+32);
  AVARS->zmj   = *(AB+36);
  AVARS->time = *T;
  AVARS->hro  = *H;
  AVARS->hroa = *(H+1);
  AVARS->roc  = *(H+4);
  AVARS->tau  = *(H+7);
  I = (int*)(H+38);
  AVARS->na1 = *(I+1);
  AVARS->nab = *(I+2);
  AVARS->nb1 = *(I+3);
  AVARS->nrd = *N;
  return(0);
}
/****************** Get ID of ShMem for status.inc  ***********************/
int setarrs (double* TE, double* AM, double* NN, int* NRD)
{   setarrs_(TE, AM, NN, NRD); }
int setarrs_(double* TE, double* AM, double* NN, int* NRD)
{ 
  int	j, i, i1, j1;
  double *D;
#include "A_arrs.h"
  if (A_Nsems == 0)	return(0);
  if (A_ShmNum < 0)	return(0);

  AARRS = (struct A_arrs *)A_ShmAdr[1];
  i1 = 0;
  for (j=0; j < NC1; j++)  AARRS->te[j]  = *(TE+j+i1);  i1 =  2*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->ti[j]  = *(TE+j+i1);  i1 =  4*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->ne[j]  = *(TE+j+i1);  i1 =  6*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->fp[j]  = *(TE+j+i1);  i1 =  8*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->mu[j]  = *(TE+j+i1);  i1 = 10*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->cu[j]  = *(TE+j+i1);  i1 = 20*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->zef[j] = *(TE+j+i1);  i1 = 44*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->upl[j] = *(TE+j+i1);  i1 = 47*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->ni[j]  = *(TE+j+i1);  i1 = 49*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->zmain[j]= *(TE+j+i1); i1 = 50*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->amain[j]= *(TE+j+i1); i1 = 52*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->vr[j]   = *(TE+j+i1); i1 = 53*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->shif[j] = *(TE+j+i1); i1 = 55*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->elon[j] = *(TE+j+i1); i1 = 56*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->tria[j] = *(TE+j+i1);
  i1 = 0;
  for (j=0; j < NC1; j++)  AARRS->ametr[j] = *(AM+j+i1);  i1 = *NRD;
  for (j=0; j < NC1; j++)  AARRS->rho[j] = *(AM+j+i1);
  i1 = 0;
  for (j=0; j < NC1; j++)  AARRS->nn[j]   = *(NN+j+i1); i1 =     *NRD;
  for (j=0; j < NC1; j++)  AARRS->tn[j]   = *(NN+j+i1); i1 =  2*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->niz1[j] = *(NN+j+i1); i1 =  3*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->niz2[j] = *(NN+j+i1); i1 =  4*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->niz3[j] = *(NN+j+i1); i1 =  5*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->nalf[j] = *(NN+j+i1); i1 =  8*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->nhydr[j]= *(NN+j+i1); i1 =  9*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->ndeut[j]= *(NN+j+i1); i1 = 10*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->ntrit[j]= *(NN+j+i1); i1 = 11*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->nhe3[j] = *(NN+j+i1); i1 = 12*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->zim1[j] = *(NN+j+i1); i1 = 13*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->zim2[j] = *(NN+j+i1); i1 = 14*(*NRD);
  for (j=0; j < NC1; j++)  AARRS->zim3[j] = *(NN+j+i1);
  return(0);
		/* Alternatively the same can be done as */
  D = (double*)A_ShmAdr[1];
  /* Common block A_VECTORS: TE, TI, NE, FP, MU, CU */
  j1 = i1 = 0;
  for (i=0; i <= 5; i++) 
    {  for (j=0; j <= NC1; j++)  *(D+j+j1) = *(TE+j+i1); 
       i1 += (*NRD);    i1 += (*NRD);  j1 += NC1;
    }
  /* Common block A_VECTORS: ZEF 20, UPL 44, NI 47, ZMAIN 49, AMAIN 50 */
  i1 = 20*(*NRD); for (j=0; j < NC1; j++) *(D+j+j1)=*(TE+j+i1); j1+=NC1;
  i1 = 44*(*NRD); for (j=0; j < NC1; j++) *(D+j+j1)=*(TE+j+i1); j1+=NC1;
  i1 = 47*(*NRD); for (j=0; j < NC1; j++) *(D+j+j1)=*(TE+j+i1); j1+=NC1;
  i1 = 49*(*NRD); for (j=0; j < NC1; j++) *(D+j+j1)=*(TE+j+i1); j1+=NC1;
  i1 = 50*(*NRD); for (j=0; j < NC1; j++) *(D+j+j1)=*(TE+j+i1); j1+=NC1;
  /* Common block A_VECTORS: VR 52, SHIF 53, ELON 55, TRIA 56 */
  i1 = 52*(*NRD); for (j=0; j < NC1; j++) *(D+j+j1)=*(TE+j+i1); j1+=NC1;
  i1 = 53*(*NRD); for (j=0; j < NC1; j++) *(D+j+j1)=*(TE+j+i1); j1+=NC1;
  i1 = 55*(*NRD); for (j=0; j < NC1; j++) *(D+j+j1)=*(TE+j+i1); j1+=NC1;
  i1 = 56*(*NRD); for (j=0; j < NC1; j++) *(D+j+j1)=*(TE+j+i1); j1+=NC1;
  i1 = 0;
/*
  printf("Zef(0) = %g,\tni(0)= %g,\t  V\'(0) = %g\n",
  *(TE+0*(*NRD)),*(TE+2*(*NRD)),*(TE+4*(*NRD))); */
  printf("\nU(0) = %g,\tni(0)= %g,\t  A(0) = %g\n",
	 *(D+7*NC1),*(D+8*NC1),*(D+10*NC1));
  /* Common block A_EQUIL: AMETR, RHO */
  for (i=0; i <= 1; i++) 
    {  for (j=0; j < NC1; j++) *(D+j+j1) = *(AM+j+i1);
       i1 += *NRD;  j1 += NC1;
    }  i1 = 0;
  /* Common block A_IONS: NN, TN, NIZ1, NIZ2, NIZ3, NALF */
  for (i=0; i <= 5; i++) 
    {  for (j=0; j < NC1; j++) *(D+j+j1) = *(NN+j+i1);
       i1 += *NRD;  j1 += NC1;
    }  i1 = 8*(*NRD);
  for (i=8; i <= 14; i++) 
    {  for (j=0; j < NC1; j++) *(D+j+j1) = *(NN+j+i1);
       i1 += *NRD;  j1 += NC1;
    }
  return(0);
}
/*---------------------------------------------------------------------*/
int freeshm ()
{   freeshm_(); }
int freeshm_()
{ struct shmid_ds Myshmid_ds;
  int j;
  if (A_Nsems != 0)	semctl(A_SemID, 0, IPC_RMID);
  if (A_ShmNum < 0)	return(0);
  for (j=0; j <= A_ShmNum; j++)
  {
      /* Detach and remove all shared memory segments */
      if (shmdt(A_ShmAdr[j]) < 0)
      {  switch( errno )
         { case EINVAL:
	     printf(">>>  Invalid ShmAdr[j] value\n",j);
	 }
      }
      else if (shmctl(A_ShmID[j], IPC_RMID, &Myshmid_ds) < 0)
      {  switch( errno )
         { case EPERM:
	     printf(">>> User cannot remove shared memory segment #%d\n",j);
	     break;
           default:
	     printf(">>> Process # %d >>>  unknown shmctl error\n",j);
	 }
      }
   }
  /* printf(">>> IPC cleaning done\n"); */
  return(0);
}
/**************************************************************************/
/* Call as:  i = write_aipc_(My.OrdNr, My.ShMid, lS);                     */
#include "A_vars.h"
int write_aipc (const struct A_proc_info My, char* AWD, int* lS)
{ FILE *A_PDF;
 char A_IPC[96]; 
  strcpy(A_IPC,AWD);
  strcat(A_IPC,"tmp/astra.ipc");
  A_PDF = fopen(A_IPC,"a");
  if (!A_PDF) 
  { printf("Cannot open existing Astra IPC file: \"%s\"\n",A_ipc_file);
    exit(0);
  }
  fprintf(A_PDF,"%12d%12d%12d   %s\n", getpid(), My.ShMid, *lS, My.Path);
  fclose(A_PDF);
/*fprintf(stdout,"%12d%12d%12d  %s\n", getpid(), My.ShMid, *lS, My.Path);
  printf("\n\n  After Proc No. %d\n",My.OrdNr);
  system("cat /home/grp/MAstra/tmp/astra.ipc");                           */
  return(0);
}
/**************************************************************************/
/*    The function can be (optionally) called from each subprocess        */
/* to submit its communication report                                     */
/* Call as:  SP_stamp(Mama, My, SemID);                                   */
void SP_stamp(const struct A_proc_info Mama, 
              const struct A_proc_info My,
 int SemID)
{
/*printf("I am (#%d) \"%s\",\tMy_PID %d,  SemID %d\n  A_PID %d,  A_Key %d\n",
         My.OrdNr, My.Path, My.Pid, SemID, Mama.Pid, Mama.Key);
printf("I am (#%d) \"%s\",\tMy_PID %d,  SemID %d\n",
	 My.OrdNr, My.Path, My.Pid, SemID);
*/
}
/**************************************************************************/
/* Compare Sem0 value with buf0.sem_op and
	wait until all subprocesses increment Sem0 by 1
   This is the same as "wait4all" but additionally it reports
   all changes in a status of the semaphore set                    -------*/
int checkeach ()                /* Normally, this function is not used    */
{   checkeach_(); } 
int checkeach_()
{ 
  if (A_ShmNum < 0) return(0);	  /* Do check only after initialization  */
//                                  {sec,nanosec}
//  static struct timespec timeout = {0,1000000};   /* timeout = 1 msec  */
//  static struct timespec timeout = {0,1000};      /* timeout = 1 mksec */
  static struct timespec timeout = {0,100000000};   /* timeout = 1 mksec */

  int j;
  static int i=0;
  /*  union semun { 
     int val; 
     struct semid_ds *buf; 
     ushort *array; 
  } arg; 
  */
  union semun arg; 
  //  ushort semarr[A_Nsems], semar0[A_Nsems]; 
  arg.array = &semarr[0]; 

 MinorLoop:
  {/* Here the primary process can do limited actions e.g. analyze keys */
      /* (void) AstraEvent(); */
      if (semctl(A_SemID, 0, GETALL, arg) == 0) 
      {               /* Report any change in the semaphore set status */
	if (i == 0)
	    { 	 					  /* 1st check */
	       for (j=0; j < A_Nsems; j++) semar0[j] = *(arg.array+j);
	       i = 1;
	       goto MinorLoop;
	    }
	 else
	    {  for (j=0; j < A_Nsems; j++)		  /* 2nd check */
		  if (semar0[j] != *(arg.array+j) ) i = 2;
	    }
	 if (i == 2) printf("Semval = {%d,%d,%d,%d},    buf0.sem_op = %d\n",
			      *arg.array, *(arg.array+1), *(arg.array+2),
			      *(arg.array+3), buf0.sem_op);
	 i = 1;
         for (j=0; j < A_Nsems; j++) semar0[j] = *(arg.array+j);
      }
      else
      {  printf("semctl error >>> errno %d  EFAULT %d EIDRM %d\n",
		errno,EFAULT,EIDRM);
      }
  }
  /* Go on if [Sem0value+buf0.sem_op==0], goto Minorloop after timeout */
  if (semop(A_SemID, &buf0, 1))
  {  switch( errno )
    { case EAGAIN:			/* Timed out */
	/* printf("Timed out\n"); */
	  goto  MinorLoop;
      case EIDRM:		     /* The semaphore set was removed */
	  printf("The semaphore set was removed\n");
	  goto  Exit;
      case EINTR:/*While blocked in this call, the process caught a signal*/
	  printf("While blocked in this call, the process caught a signal\n");
	  goto  Exit;
      case EINVAL:		     /*  */
	  printf("EINVAL\n");
	  goto  Exit;
      case EFAULT:		     /*  */
	  printf("EFAULT\n");
	  goto  Exit;
      case E2BIG:		     /*  */
	  printf("E2BIG\n");
	  goto  Exit;
      case EACCES:		     /*  */
	  printf("EACCES\n");
	  goto  Exit;
      case EFBIG:		     /*  */
	  printf("EFBIG\n");
	  goto  Exit;
      case ENOMEM:		     /*  */
	  printf("ENOMEM\n");
	  goto  Exit;
      case ERANGE:		     /*  */
	  printf("\n");
	  WhatSem();
	  printf("ERANGE error:  sem_op = %d,   SEMVMX = ?\n",
		 buf0.sem_op);
	  printf("\n");
	  goto  Exit;
      default:
	  printf(">>> PRIMARY >>> Unrecognised semtimedop error %d\n",errno);
	  goto Exit;
     }
  }                       	     /* printf("PRIMARY unlocked\n"); */
  //  printf("Check passed  "); WhatSem();
  buf0.sem_op = 0;
  return(0);
 Exit:
  a_stop();
}                  /* End checkeach */
/*--------------------------------------------------------------------*/
