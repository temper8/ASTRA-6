struct A_proc_info
{
  char	 Path[96];      /* Path of a process */
  key_t	 Key;           /* Key of a process */
  pid_t	 Pid;           /* pid of a process */
  int	 OrdNr;         /* Ordinal number of a process */
  int	 ShMid;         /* ID of ShMem for a process */
  double CPUse;         /* CPU time uzage of a process */
  int	 Size;          /* Size of A_proc_info */
  int    CheckWord;     /* Control number */
};
struct A_vars
{   /* Common block A_VARIABLES */
  double ab;   		/* AB   */
  double abc;  		/* ABC  */
  double aim1;   	/* AIM1 */
  double aim2;   	/* AIM2 */
  double aim3;   	/* AIM3 */
  double amj;    	/* AMJ,   AWALL*/
  double btor;   	/* BTOR */
  double elong;  	/* ELONG, ELONM*/
  double encl;   	/* ENCL */
  double enwm;   	/* ENWM, FECR,FFW,FICR,FLH,GN2E,GN2I */
  double ipl;    	/* IPL,  LEXT*/
  double nncl;   	/* NNCL */
  double nnwm;   	/* NNWM, QECR,QFW,QICR,QLH,QNBI*/
  double rtor;   	/* RTOR */
  double shift;  	/* SHIFT*/
  double trian;  	/* TRIAN, TRICH, UEXT*/
  double updwn;  	/* UPDWN, WNE,WTE,WTI*/
  double zmj;	  	/* ZMJ  */
  double MoreDbls[24];  /* free */
	/* Common block A_GRIDS */
  double time;  	/* TIME */
  double hro;   	/* HRO  */
  double hroa;  	/* HROA */
  double roc;   	/* ROC  */
  double tau;   	/* TAU  */	   
  int    na1; 	  	/* NA1 */ 
  int    nab;	  	/* NAB */
  int    nb1;	  	/* NB1 */
  int    nrd;	  	/* NRD */
  int	 c_size[12];  /* sizeof: int, double, key_t, pid_t, etc */
  int    MoreInts[63];	/* free */
  int    CheckWord;     /* Control number */
} *AVARS;
