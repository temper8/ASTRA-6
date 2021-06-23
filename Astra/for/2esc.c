#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>
#include <string.h>
void	escvout (double*,double*,double*);
void	escvout_(double*,double*,double*);
void	escmout (double*,double*,double*,double*,double*);
void	escmout_(double*,double*,double*,double*,double*);
void	escaout (double*,double*,double*,double*,double*,double*,
		 double*,double*,double*,double*,double*);
void	escaout_(double*,double*,double*,double*,double*,double*,
		 double*,double*,double*,double*,double*);
void	get3ma (double*,double*,double*,double*,double*,double*);
void	get3ma_(double*,double*,double*,double*,double*,double*);
void	drawconfig(int*,int*,double*,int*);	
void	drawconfig_(int*,int*,double*,int*);
void	SemCheck();
const	double	greekp=3.1415926536;
const	double	gp24=39.4784176;
extern	int	Xmode;		/* Default: NoX - Batch mode */
extern	int ESNp;
extern	double ESgt[],*ESsr,*ESsz,*ESaR0,*ESaZ0;
static	int isw[128];
static	double sr[128],sra[128],srq[128] ,sz[128],sza[128],szq[128]
	,aB[128],aBa[128],aBq[128],gH[128],gHa[128],gHq[128];
static	double F[128],Fa[128],gFa[128],gFaa[128],gYa[128],gYaa[128]
	,T[128],Ta[128],P[128],Pa[128];
static	int Ifesilinked=1;
/*#if defined(_SEM_SEMUN_UNDEFINED) */
union semun {
  int val;                    /* value for SETVAL */
  struct semid_ds *buf;       /* buffer for IPC_STAT, IPC_SET */
  unsigned short int *array;  /* array for GETALL, SETALL */
  struct seminfo *__buf;      /* buffer for IPC_INFO */
};
/*#endif */

#define Sem0 0
#define Sem1 1

int ESNa,ESNa1,ESNp,ESNp1,ESFp,ESFp1;
double *ESaR0,*ESaZ0;
double *ESsr,*ESsz,*ESaB;
double *ESsb,*ESsb2a,*rcT,*rcT2a,*rsT,*rsT2a;

double ESgt[129];

static int Pid;
static int SemID;
static int ShmID;
static union semun Mysemun;
static struct sembuf bufP	={Sem0,-1,SEM_UNDO};
static struct sembuf bufV	={Sem1, 1,IPC_NOWAIT};
static char *lShm;
static int kD,kI;

static int *ShMemI;
static double *ShMemT;
static double *ShMemB;
static double *ShMemP;
static double *ShMemE;
static double *ShMBal;
static double *aBal;
static int *kBal;
static int nBal,NBal=0x100;

/* Addition for spline reconstruction */

static double cr3,cr6,cr4,cr2;
static double cr3ha,cha,crha,c2rha,cHa,crHa,chha,ch2ha,ch2h1,ch4h0,ch4h1;
static int iSplA=0,iSplA1;
static double SplAx,SplAxx,SplAX,SplAXX;
static double ESsa[0x200];
/*---------------------------------------------------------------------*/
int ESFirstSetSplA()
{
  int i;
  double h;

  cr3	=1./3.;
  cr6	=1./6.;
  cr4	=0.25;
  cr2	=0.5;
  h	=1./ESNa;
  for(i=0; i < ESNa1; i++) ESsa[i]=h*i;

  h	=ESsa[1]-ESsa[0];
  cha	=cr6*h;
  cr3ha	=cr3*h;
  crha	=1./h;
  c2rha	=2.*crha;
  chha	=cr2*h;
  ch2ha	=chha*cha;
  cHa	=cr6*h*h;
  crHa	=6.*crha*crha; 
  ch2h1	=cr3*h*h;
  ch4h1	=-cr3*h*h*h*h/30.;
  ch4h0	=7.*cr4*ch4h1;
  ch4h1	*=2.;
  return(0);
}

int ESSetSplA(double A)
{
  int iA;
  int i,j,kT1;
  
  if(A > 1.) A=1.;
  if(A < 0.) A=0.;
  while(iSplA < ESNa && ESsa[iSplA+1] < A) iSplA++;
  while(iSplA > 0 && ESsa[iSplA] > A) iSplA--;
  iSplA1	=iSplA+1;
  SplAx		=(A-ESsa[iSplA])*crha;
  SplAxx	=SplAx*SplAx;
  SplAX		=1.-SplAx;
  SplAXX	=SplAX*SplAX;
  return(0);
}

int splRA(double*f, double*df, double*g, double*d2g)
{
  *f	=SplAX*g[iSplA]+SplAx*g[iSplA1]
    +(SplAX*(SplAXX-1.)*d2g[iSplA]+SplAx*(SplAxx-1.)*d2g[iSplA1])*cHa;
  if(df != NULL){
    *df	=(g[iSplA1]-g[iSplA]
	  +((3.*SplAxx-1.)*d2g[iSplA1]-(3.*SplAXX-1.)*d2g[iSplA])*cHa)*crha;
  }
  return(0);
}
/* End of additions */

int LaunchEsc(int ifX)
{
  int i;
  char elem;
  char *lc,*ls,ln[64];
  double t,dt;
  const char *ObjCode = "ObjectCode";

  kD	=sizeof(double);
  kI	=sizeof(int);

  /* Initialize */
  Pid	=getpid();

  /* Init semaphores */
  /* Two semaphores (Empty and Full) are created in one set */
  SemID	=semget((key_t)Pid, 2, 0666|IPC_CREAT);

  /* Init Empty to number of elements in shared memory */
  Mysemun.val	=1;
  semctl(SemID, Sem0, SETVAL, Mysemun);		/* Set Sem0 to 1 */

  /* Init Full to zero, no elements are produced yet */
  Mysemun.val	=0;
  semctl(SemID, Sem1, SETVAL, Mysemun);		/* Set Sem1 to 0 */
  if (getenv("XESC") == NULL) /* Check if $XESC is set, then use X graphics */
/*if ( ifX == 0 ) */        /* The same can be controlled via parameter ifX */
  {
/*sprintf(ln,"cd ESC;Obj/esc2.%s %d > /dev/null &",getenv(ObjCode),Pid); */
  printf(" ESC is started in no X mode\n");
  sprintf(ln,"cd ESC;Obj/esc2.%s %d &",getenv(ObjCode),Pid);
  }
  else
  {
/*sprintf(ln,"cd ESC;Obj/esc2.X%s %d > /dev/null &",getenv(ObjCode),Pid); */
  printf(" ESC is started in X mode\n");
  sprintf(ln,"cd ESC;Obj/esc2.X%s %d &",getenv(ObjCode),Pid);
  }
/*  printf(" Calling string: \"%s\"\n",ln); */
/*  (void) SemCheck(); */
  semop(SemID, &bufP, 1);   /*Set Sem1=0, Lock ASTRA*/
  i = system(ln);
  semop(SemID, &bufV, 1);   /*Set Sem0=1, Open ESC*/
  /*  printf("Open ESC passed\n"); */
  (void) SemCheck();
  semop(SemID, &bufP, 1);   /*Wait until Sem0++ by ESC: Lock ASTRA*/
  /*  printf("Lock ASTRA passed\n"); */
  
  /* Init Shared memory */
  ShmID	=shmget((key_t)(Pid+1),0, 0666 | IPC_CREAT);
  /* attach shared memory to process */
  lShm=shmat(ShmID,NULL,0);

  /*		Put this group in esc2 
  FILE *A_PDF;
  A_PDF = fopen("tmp/astra.ipc","a");
  if (A_PDF == NULL){}
  fprintf(A_PDF,"%-16s%12d%12d%12d\n", getpid(), ShmID, SemID);
  fclose(A_PDF);"ESC: PID, ShmID, SemID"
  */

  ShMemI	=(int*)lShm;
  ShMemT	=(double*)(ShMemI+32);
  ShMemB	=ShMemT+16;
  ShMemP	=ShMemB+480;
  ShMemE	=ShMemP+0x600;

  ESNa	=ShMemI[1];
  ESNp	=ShMemI[2];
  ESFp	=ShMemI[3];
  ESaR0	=ShMemT+1;
  ESaZ0	=ESaR0+1;

  ESNa1	=ESNa+1;
  ESNp1	=ESNp+1;
  ESFp1	=ESFp+1;

  ShMBal	=ShMemE+ESNp1+(11+18*ESNp1)*ESNa1;
  aBal		=ShMBal;
  kBal		=(int*)(aBal+NBal);

  ESISetMem(ShMemE,ESNa1,ESNp1);
  ESsr	=ShMemE+ESNp1+11*ESNa1;
  i	=ESNp1*ESNa1;
  ESsz	=ESsr+4*i;
  ESaB	=ESsz+4*i;

  ESsb	=(double*)(kBal+NBal);
  ESsb2a=ESsb+ESNa1;
  rcT	=ESsb2a+ESNa1;
  i	=ESNa1*ESFp1;
  rcT2a	=ESsb2a+i;
  rsT	=rcT2a+i;
  rsT2a	=rsT+i;

  dt	=8.*atan(1.)/ESNp;
  for(i=0; i < ESNp1; i++){
    t	=dt*i;
    ESgt[i]	=t;
  }
  ESFirstSetSplA();

  semop(SemID, &bufV, 1);   /*Open ESC*/   		   
  /*  printf("LaunchESC exit: Sem0=%d, Sem1=%d\n",
      semctl(SemID, Sem0, GETVAL),semctl(SemID, Sem0, GETVAL)); */
  return(0);
}

void esc0_(ifX) int *ifX;
{
  LaunchEsc(*ifX);
  return;
}
void esc0(ifX) int *ifX;
{ esc0_(ifX); }

void esc1_()
{
  int k;
  struct shmid_ds Myshmid_ds;

  esifree_();
  (void) SemCheck();
  semop(SemID, &bufP, 1);   
  k	=4;
  memcpy(lShm,(void*)&k,kI);
  semop(SemID, &bufV, 1);   		   
  (void) SemCheck();
  semop(SemID, &bufP, 1);   
  /* Remove semaphores */
  semctl(SemID,0, IPC_RMID, Mysemun);
  /* Remove shared memory */
  shmctl(ShmID, IPC_RMID, &Myshmid_ds);
  return;
}

void esc_(
	  int *Ffail 		/* Flag of failure:
				   0 - normal operation;
				   1 - some problems;
				   2 - problems near the boundary;
				   4 - problems near the axis;
				   8 - problems in the middle;
				   ....;
				   */
	  ,int *Fjob 		/* Flag for job assignment:
				   0 - nothing special was requested;
				   1 - save the geometry;
				   2 - take the geometry;
				   4 - write ESI data files;
				   ...;
				   */ 
	  ,int *Mpol		/* Working number of Fourier harmonics
				   in $\Psi$ */
	  ,double *sTime	/* time */
	  )
{
  int k;
  double *ld;
  (void) SemCheck();
  semop(SemID, &bufP, 1);
  ShMemT[0]	=*sTime;
  ShMemI[8]	=*Fjob;
  ShMemI[10]	=0;
  *sTime	=ShMemT[15];

  semop(SemID, &bufV, 1);
  (void) SemCheck();
  semop(SemID, &bufP, 1);
  *Mpol	=ShMemI[7];
  *Ffail=ShMemI[9];
  ESIInit();
  semop(SemID, &bufV, 1);
  return;
}
void esc(
	  int *Ffail
	  ,int *Fjob
	  ,int *Mpol
	  ,double *sTime
	  )
{
  esc_(Ffail,Fjob,Mpol,sTime);
}

int getbalst_(
	      double *a /* array of ESC radii (0< a <= 1) of tested surfaces */
	      ,int *ind /* array of 0 - stable, 1 - ballooning unstable,
			      2 - Mercier unstable */
	      , int *n /* number of points in return arrays */
	      )
{
  *n	=ShMemI[10];
  if(*n == 0) return(1);
  memcpy((void*)a,(void*)aBal,ShMemI[10]*sizeof(double));
  memcpy((void*)ind,(void*)kBal,ShMemI[10]*sizeof(int));
  return(0);
}

void escin_(double	*pProf	/* p[] - pressure related profile */
	    ,double	*jProf	/* j[] - current related profile */
	    ,double	*aProf	/* a[] - normalized sqrt(Phi) square root 
				   of toroidal flux.
				   If a=NULL (dropped in FORTRAN) - uniform 
				   grid from 0 till 1,
				   otherwise it is considered as an array of 
				   grid points */
	    ,int	*nProf	/* number <= 129 of profile points including
				   magnetic axis and plasma edge */
	    ,int	*nrhH	/* First digit (kAst)
				   n - 0 - normalized radial coordinate (0-1)
				       1 - absolute radial coordinate;
				   r - 0 - b (vertical semiaxis)
				           as the radial coordinate,
				       1 - V (volume),
				       2 - gF (toroidal flux),
				       3 - gY (poloidal flux, DO NOT use),
				   h specifies the meaning
				   of pProf:
				   0 - j_p;
				   1 - P =dp/dpsi;
				   2 -p [MPa].
				   Second digit H specifies the meaning
				   of jProf:
				   0 - j_s [MA/m^2];
				   1- j_|| [MA/m^2]=(j dot B)/(R_0 B grad phi),
				                   R_0 = magnetic axis;
				   2- j_||R_0 [MA/m] = (j dot B)/(B grad phi);
				   3 - T=FF';
				   6 - q;
				   7 - 1/q;
				   8 - \gY;
				   9 - \sigma, conductivity [MS}
				   For Example:
				   26 - p[] and q[] profiles are supplied;
				   21 - p[] and j||[] profiles are supplied;
				   0 -  jp[] and js[] profiles are supplied.

				   Possible combinations are limited to:
				    0, 1, 2, 6, 7
				   10,11,12,16,17
				   21,22,26,27
				 */
	    ,double	*Rpv	/* R[m]- plasma boundary */
	    ,double	*Zpv	/* Z[m]- plasma boundary */
	    ,int	*Npv	/* number <= 257 of the plasma-vacuum 
				   points.
				   If *Npv >12, the first and last points
				   coincide */
	    ,double	*RBtor	/* RBtor [m Tesla] outside the 
				   plasma */
	    ,double	*Rext	/* Reference major radius [m] */
	    )
{
  int i,k;
  void *lc;
  (void) SemCheck();
  semop(SemID, &bufP, 1);   /*Lock ASTRA*/
  ShMemI[4]	=*nrhH;
  ShMemI[5]	=*Npv;
  ShMemI[6]	=*nProf;
  ShMemT[3]	=*RBtor;
  ShMemT[4]	=*Rext;
  k	=(*Npv)*kD;
  memcpy((void*)ShMemB,(void*)Rpv,k);
  memcpy((char*)ShMemB+k,(void*)Zpv,k);

  k	=(*nProf)*kD;
  if(aProf != NULL) memcpy((void*)ShMemP,(void*)aProf,k);
  memcpy((char*)ShMemP+k,(void*)pProf,k);
  memcpy((char*)ShMemP+2*k,(void*)jProf,k+2*sizeof(double));
  semop(SemID, &bufV, 1);   /*Open ESC for return confirmation*/
  (void) SemCheck();
  semop(SemID, &bufP, 1);   /*Waits until ESC does Sem0++ (Locks ASTRA)*/
  semop(SemID, &bufV, 1);   /*Open ESC*/   		   
  /*  printf("ESCIN exit: Sem0=%d, Sem1=%d\n",
      semctl(SemID, Sem0, GETVAL),semctl(SemID, Sem0, GETVAL)); */
  return;
}
void escin(double	*pProf
	    ,double	*jProf
	    ,double	*aProf
	    ,int	*nProf
	    ,int	*nrhH 
	    ,double	*Rpv  
	    ,double	*Zpv  
	    ,int	*Npv  
	    ,double	*RBtor
	    ,double	*Rext 
	    )
{    escin_(pProf,jProf,aProf,nProf,nrhH,Rpv,Zpv,Npv,RBtor,Rext);    }

/* See also "man sigwaitinfo"  "man sigtimedwait" */
void SemCheck()
{
  const int i=1;
  if (Xmode) goto Loop;
  return;
 Loop:
  if ( semctl(SemID, Sem0, GETVAL, Mysemun) > 0 )	return;
  usleep(i);
  (void) AstraEvent();
  goto Loop;
}

#include	<errno.h>
int SemTimedCheck()		/*!!!!!!!!!!!! To be checked */
{  static struct timespec timeout = {0,1000};   /* timeout = 1 mksec */
MinorLoop:
  (void) AstraEvent();
  if (semtimedop(SemID, &bufP, 1, &timeout))
  {  switch( errno )
    { case EAGAIN:			/* printf("Timed out\n"); */
	  goto  MinorLoop;
      case EIDRM:		     /* The semaphore set was removed */
	  goto  Exit;
      default:
	  printf(">>> PRIMARY >>> Unrecognised semtimedop error \n");
	  goto Exit;
     }
  }                       	     /* printf("PRIMARY unlocked\n"); */
Exit:
  /* If needed put here 
  bufP.sem_op = 0;
     check if semtimedop has already done it */
  return(0);
}
/* Astra definitions: */
/* static struct sembuf bufP	= {Sem0,-1, SEM_UNDO};	*/ /* Lock */
/* static struct sembuf bufV	= {Sem1, 1, IPC_NOWAIT};*/ /* Open */
/* Esc definitions: */
/* static struct sembuf bufP	= {Sem1,-1, SEM_UNDO};	*/ /* Lock */
/* static struct sembuf bufV	= {Sem0, 1, IPC_NOWAIT};*/ /* Open */
/*---------------------------------------------------------------------*/
/************************** Return average quantities ***************/
void escvout ( double	*Rmaxis
	      ,double	*Zmaxis
	      ,double	*gFt
	      )
{    escvout_(Rmaxis,Zmaxis,gFt); 
}
void escvout_( double	*Rmaxis
	      ,double	*Zmaxis
	      ,double	*gFt
	      )
{
  double gFtor,gY,aq;
  if(Ifesilinked){
    esilink2c_(F,Fa,gFa,gFaa,gYa,gYaa,T,Ta,P,Pa
	       ,sr,sra,srq,sz,sza,szq,aB,aBa,aBq,gH,gHa,gHq,isw);
    Ifesilinked	=0;
  }
  aq = 1.;
  *Rmaxis = ESaR0[0];
  *Zmaxis = ESaZ0[0];
  esigetmagfluxes_(&gFtor, &gY, &aq);
  *gFt = gFtor;
  return;
}
/************************** Return average quantities ***************/
void escmout ( double	*a
	      ,double	*Va
	      ,double	*G33
	      ,double	*Fpol
	      ,double	*sa
	      )
{    escmout_( a, Va, G33, Fpol, sa );
}
void escmout_( double	*a
	      ,double	*Va
	      ,double	*G33
	      ,double	*Fpol
	      ,double	*sa
	      )
{
  int j;
  double d,r,ra,rt,za,zt,rh;
  double aq[128];

  if(Ifesilinked){
    esilink2c_(F,Fa,gFa,gFaa,gYa,gYaa,T,Ta,P,Pa
	       ,sr,sra,srq,sz,sza,szq,aB,aBa,aBq,gH,gHa,gHq,isw);
    Ifesilinked	=0;
  }
  rh	=1./ESNp;
  for(j=0; j < ESNp; j++){
    aq[j]	=a[0];
  }
  esiget2dfunctions_(aq, ESgt, &ESNp);
  Fpol[0] = F[0];
  Va[0]	  = 0.;
  G33[0]  = 0.;
  sa[0]  = 0.;
  if(a[0] != 0.){
    for(j=0; j < ESNp; j++){
      r	=sr[j];
      ra	=sra[j];
      rt	=srq[j];
      za	=sza[j];
      zt	=szq[j];
      d	=ra*zt-rt*za;
      Va[0]	+=r*d;
      G33[0]	+=d/r;
    }
    G33[0]	/=Va[0];
    Va[0]	*=gp24*rh;
    sa[0]	=-a[0]*(gYaa[0]/gYa[0]-gFaa[0]/gFa[0]);
  }
  return;
}
/*************** Map ESC configuration to 3M geometry *******************/
void escaout (double	*a
	     ,double	*Vas
	     ,double	*Slat
	     ,double	*g_11
	     ,double	*g_22
	     ,double	*gM
	     ,double	*Bmin
	     ,double	*Bmax
	     ,double	*Bavr
	     ,double	*B2av
	     ,double	*rB2a
	     )
{   escaout_ (a,Vas,Slat,g_11,g_22,gM,Bmin,Bmax,Bavr,B2av,rB2a);
}
void escaout_(double	*a
	     ,double	*Vas
	     ,double	*Slat
	     ,double	*g_11
	     ,double	*g_22
	     ,double	*gM
	     ,double	*Bmin
	     ,double	*Bmax
	     ,double	*Bavr
	     ,double	*B2av
	     ,double	*rB2a
	     )
{
  int j;
  double B2,sqg,d,g22,r,ra,rt,za,zt,rh;
  double aq[128];

  if(Ifesilinked){
    esilink2c_(F,Fa,gFa,gFaa,gYa,gYaa,T,Ta,P,Pa
	       ,sr,sra,srq,sz,sza,szq,aB,aBa,aBq,gH,gHa,gHq,isw);
    Ifesilinked	=0;
  }
  rh	=1./ESNp;
  for(j=0; j < ESNp; j++){
    aq[j] = a[0];
  }
  esiget2dfunctions_(aq, ESgt, &ESNp);
  gM[0]	=-gYa[0]/gFa[0];
  Vas[0] = 0.;
  Slat[0] = 0.;
  g_11[0] = 0.;
  g_22[0] = 0.;
  Bmin[0] = 1.e5;
  Bmax[0] = 0.;
  Bavr[0] = 0.;
  B2av[0] = 0.;
  rB2a[0] = 0.;
  if(a[0] == 0.){
    printf(">>> ERROR >>> Average on the magnetic axis is requested\n");
    return;
  }
  for(j=0; j < ESNp; j++){
    r	= sr[j];
    ra	= sra[j];
    rt	= srq[j];
    za	= sza[j];
    zt	= szq[j];
    d	= ra*zt-rt*za;
    g22	= rt*rt+zt*zt;
    sqg = r*d;
    Vas[0] += sqg;				/* VRS*ROC */
    Slat[0] += sqrt(g22)*r;			/* SLAT */
    g_11[0] += g22*r/d;				/* G11/ROC */
    g_22[0] += g22/(sqg*F[0]);			/* G22/(ROC*BTOR*RTOR*RTOR) */
    if (Bmin[0] > aB[j]) Bmin[0] = aB[j];
    if (Bmax[0] < aB[j]) Bmax[0] = aB[j];
    B2 = aB[j]*aB[j];
    Bavr[0] += aB[j]*sqg;
    B2av[0] += B2*sqg;
    rB2a[0] += sqg/B2;
  }
  Slat[0] *= gp24*rh;
  g_11[0] *= gp24*rh;
  g_22[0] *= rh;
  Bavr[0] /= Vas[0];
  B2av[0] /= Vas[0];
  rB2a[0] /= Vas[0];
  Vas[0] *= gp24*rh;
  /*    printf("%f,   %f\n",a[0],F[0]); */
  return;
}
/*************** Map ESC configuration to 3M geometry *******************/
void get3ma (gr,a,sh,el,tr,ud)
     double	*gr;
     double	*a, *sh, *el, *tr, *ud;
{
     get3ma_(gr,a,sh,el,tr,ud);
}
void get3ma_(gr,a,sh,el,tr,ud)
     double	*gr;
     double	*a,*sh,*el,*tr,*ud;
{
  int j,n;
  double aq[4],gq[4];

  gq[0]	=0;
  gq[1]	=0.25*ESgt[ESNp];
  gq[2]	=2.*gq[1];
  gq[3]	=3.*gq[1];
  
  for(j=0; j < 4; j++){
    aq[j]	=gr[0];
  }
  j	=4;
  esiget2dfunctions_(aq,gq,&j);
  if(gr[0] != 0.){
    a[0]	=0.5*(sr[0]-sr[2]);
    sh[0]	=0.5*(sr[0]+sr[2]);
    el[0]	=0.5*(sz[1]-sz[3])/a[0];
    tr[0]	=0.5*(sr[0]+sr[2]-sr[1]-sr[3])/a[0];
    ud[0]	=0.5*(sz[1]+sz[3]);
  }
  else{
    a[0]	=0.;
    sh[0]	=0.;
    el[0]	=(sza[1]-sza[3])/(sra[0]-sra[2]);
    tr[0]	=0.;
    ud[0]	=0.;
  }
  return;
}
/************** Plot ESC configuration in 8th output  mode **************
   j0     r-axis position,
   ysc    horizontal/vertical scale for the picture,
   is     array for saving/erasing the previous picture,
   *ESsr  1-st element of 21x65 array r(\rho_N,\tau),  \rho_N = \rho/\rho_bnd
   *ESsz  1-st element of 21x65 array z(\rho_N,\tau),
************************************************************************/
/**************** Draw equilibrium in ESC mode *************************/
/* id (0) - Pixmap id */
/* j0  - mid-plane position */
/* ysc - scale factor */
/* is  - old configuration to be erased */
/* Note! When OUTDSP in 8th mode is called from review then this function is missing
         therefore a dummy subroutine is created when .exe/makereview is running  */
void drawconfig (id,j0,ysc,is)	int *id,*j0,*is;	double	*ysc;
{	drawconfig_(id,j0,ysc,is);	}
void drawconfig_(id,j0,ysc,is)
     int *id,*j0,*is;
     double	*ysc;
{
  extern int ESNa,ESNa1,ESNp1;
  int i, j,ji, k, x10, y10, x1, x2, y1, y2, EraseColor=31;
  colovm(&EraseColor);			/* Erase old */
  drawvm(id, is, is+1, is, is+1);
  for (k=2, i=2; i <ESNa1; i+=2){
    for (j=0; j < ESNp1;  j++){
      drawvm(id, is+k, is+k+1, is+k+2, is+k+3);
      k += 4;
    }
  }
  i = 14;	colovm(&i);
  k = 0;	x1 = ESsr[0]*(*ysc);	y1 = *j0-ESsz[0]*(*ysc);
  is[k++] = x1;	is[k++] = y1;		drawvm(id, &x1, &y1, &x1, &y1);
  for (i=2; i < ESNa1; i+=2){
    if(i == ESNa ) { x1=2;	colovm(&x1);	}
    x1 = x10 = ESsr[ESNp1*i]*(*ysc);
    y1 = y10 = *j0-ESsz[ESNp1*i]*(*ysc);
    for (j=0; j <ESNp1;  j++){
      ji=ESNp1*i+j;
      x2 = ESsr[ji]*ysc[0];
      y2 = *j0-ESsz[ji]*ysc[0];
      if (j == ESNp1-1) {
	x2 = x10;
	y2 = y10;
      }
      *(is+k++) = x1;
      *(is+k++) = y1;
      *(is+k++) = x2;
      *(is+k++) = y2;
      /*if(i==10) printf("%d  %d  %d  %f  %f\n",j,ji,x2,ESsr[ji],ESsz[ji]); */
      drawvm(id, &x1, &y1, &x2, &y2);
      x1 = x2;
      y1 = y2;
      if(k > (ESNa/2+1)*ESNp1*4-4)	return;
    }
  }
  return;
}
/**********************************************************************/
