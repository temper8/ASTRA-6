#include 	<stdio.h>
#include	<stdlib.h>
#include	<ctype.h>
#include	<string.h>
#include <sys/types.h>
#include <dirent.h>
#define		MxLen	6	/* 	Maximal_(fml/fnc)_name_length */

/*
   In order to test FileList:
   1. Move the next line outside this comment
#define DEB
   2. Type
cc -o ss Src/Sort.c <Ret>
   3. Type the command
ss <name of Directory>  <Ret>
 */

int FileList(char *FList[],int NList, char*DirNm)
{
  int nFList=0;
  int i,k,k1,k2;
  char *p, *lc1,*lc2;
  DIR *dirp;
  struct dirent *direntp;

  nFList	=0;
  if((dirp	=opendir(DirNm)) != NULL){
    while((direntp=readdir(dirp)) != NULL){
      p	=direntp->d_name;
      i = *(p+strlen(p)-1);		/* printf("<%c>\n",i); */
      if(*p != '.' && i != '~' && i != '%' && i != '#' && nFList < NList){
	/* Sorting */
	i	=0;
	while(i < nFList){
	  lc1	=FList[i];
	  lc2	=p;
	  while(*lc1 == *lc2 && *lc1 != '\0'){
	    lc1++;
	    lc2++;
	  }
	  k1	=*lc1 != '.' ? *lc1 : 1;
	  k2	=*lc2 != '.' ? *lc2 : 1;
	  if(k2 < k1){
	    break;
	  }
	  i++;
	}
	k	=nFList;
	while(k > i){
	  strcpy(FList[k],FList[k-1]);
	  k--;
	}
	strcpy(FList[i],p);
	nFList++;
	/*        printf("<%s>  is processed\n",p); */
      }			/* End sorting */
      /*    else printf("<%s>  is skipped\n",p); */
    }
    closedir(dirp);
  }
#ifdef DEB
  for(i=0; i < nFList; i++){
    printf("%3d <%s>\n",i,FList[i]);
  }
#endif
  return(nFList);
}

#ifdef DEB
int main(int argc, char ** argv)
{
  if(argc < 2){
    printf("Usage:\nss <DirName>\n");
    return(0);
  }
  FileList(argv[1]);
  return(0);
}
#endif

void getfml (ierr,nfml,fmlnames,lnam)
   int *ierr, *nfml, *lnam;   char fmlnames[];
{  void	getfml_(int*,int*,char*,int*);	
	getfml_(ierr,nfml,fmlnames,lnam);	}
void    getfml_(ierr,nfml,fmlnames,lnam)
   int *ierr, *nfml, *lnam;   char fmlnames[];	
{ int	i, j, k, l, n, offset=4;	/* Cmp. with "fmlname" declaration */
  char	stri80[80], fmlname[MxLen+5];	/* MxLen + offset(fml/) + 1('\n') */
  FILE	*fp1, *fp2;
	n = 0;		if ( *lnam != MxLen )	{
	    printf(">>> ERROR >>> Inconsistency in array FMLNAM length\n");
			      *ierr=-1; return; }
	strcpy (stri80,"/bin/ls -1 fml/* > tmp/fml.list");
	j = system(stri80);	if ( j != 0 )	{
	    printf(">>> ERROR >>> Directory \"fml\" not found or empty\n");
			      *ierr=-1; return; }
	fp1 = fopen("tmp/fml.list","r");	if ( fp1 == NULL ) {
	    printf(">>> ERROR >>> Cannot write to tmp/\n");
						 *ierr=-1; return; }
	fp2 = fopen("tmp/declar.fml","w");		if ( fp2 == NULL ) {
	    printf(">>> ERROR >>> Cannot write to the working directory\n");
						 *ierr=-1; return; }
	fputs("      double precision ",fp2);	l = 23;
	while (fgets(fmlname,80,fp1))
	{ i = strlen(fmlname);
 	  if ( isalnum(fmlname[i-2]) )
	     {  i -=(offset+1);		k = l+i+1;
		strncpy(fmlnames+n*MxLen,fmlname+offset,i);
		if ( i < MxLen )
			strncpy(fmlnames+n*MxLen+i,"      ",MxLen-i);
		if ( l > 23 && k <= 71 )  fputc(',',fp2);
		if ( k > 71 ) { fputs("\n      double precision ",fp2); l = 23;
				k = l+i+1; }
		for( j=0; j < i; j++ ) fputc(toupper(fmlname[j+offset]),fp2);
		l = k; n = n+1;		if ( n > *nfml )  {
		    printf(">>> ERROR >>> Total FML number > %d\n",*nfml);
						*ierr=-2; return; }
	     }
	}
	*ierr = 0;	*nfml = n;	fputc('\n',fp2);
	fclose(fp1);	fclose(fp2);	remove("tmp/fml.list");
}
/*	printf("Total: %d\n",n);
	strncat(stri80,stri50,j);
	printf("Returned code: %d\n",j);
	strcat (stri80," tmp/model.txt");
	printf("\"%s\"\n",stri80);		 j = system(stri80);
*/
void getfnc (ierr,nfnc,fncnames,lnam)
     int *ierr, *nfnc, *lnam;   char fncnames[];
{
  void	getfnc_(int*,int*,char*,int*);	
  getfnc_(ierr,nfnc,fncnames,lnam);
}

void getfnc0_(ierr,nfnc,fncnames,lnam)
     int *ierr, *nfnc, *lnam;   char fncnames[];	
{
  int	i, j, k, l, n, offset=4;	/* Cmp. with "fncname" declaration */
  char	stri80[80], fncname[MxLen+7];	/* MxLen + offset(fnc/) +
						 + ext(.f) + 1('\n') */
  FILE	*fp1, *fp2;
	n = 0;		if ( *lnam != MxLen )	{
	    printf(">>> ERROR >>> Inconsistency in array FNCNAM length\n");
			      *ierr=-1; return; }
	strcpy (stri80,"/bin/ls -1 fnc/* > tmp/fnc.list");
	j = system(stri80);	if ( j != 0 )	{
	    printf(">>> ERROR >>> Directory \"fnc\" not found or empty\n");
			      *ierr=-1; return; }
	fp1 = fopen("tmp/fnc.list","r");	if ( fp1 == NULL ) {
	    printf(">>> ERROR >>> Cannot write to tmp/\n");
						 *ierr=-1; return; }
	fp2 = fopen("tmp/declar.fnc","w");		if ( fp2 == NULL ) {
	    printf(">>> ERROR >>> Cannot write to the working directory\n");
						 *ierr=-1; return; }
  fputs("      double precision VINT,IINT,GRAD,GRADS,FRMAX,FRMIN,RFMIN,RFMAX\n",
      fp2);
  fputs("      double precision RFVAL,AFVAL,RFVEX,AFVEX,RFVIN,AFVIN,RFA\n",
       fp2);
  fputs("      double precision RFAN,XFA,XFAN,AFR,AFX,RECR,ATR,ATX\n",fp2);
  fputs("      double precision TIMINT,TIMDER,TIMAVG,GAUSS\n",fp2);
  fputs("      double precision RADIAL,RADINT,ASTEP,RSTEP,XSTEP,STEP,CUT\n",
       fp2);
  fputs("      double precision FTBOX,FXBOX,FABOX\n",fp2);
  fputs("      double precision FIXVAL,FTAV,FTMIN,FTMAX,FRAMP,FJUMP\n",fp2);
  fputs("      external IINT\n",fp2);
	fputs("      double precision ",fp2);	l = 23;
	while (fgets(fncname,80,fp1))
	{ i = strlen(fncname);
	if ( isalnum(fncname[i-2]) ){
	       i -=(offset+3);		k = l+i+2;
/*	        printf("%d \"%s\"\n",i,fncname);               */
		strncpy(fncnames+n*MxLen,fncname+offset,i);
		if ( i < MxLen )
			strncpy(fncnames+n*MxLen+i,"      ",MxLen-i);
		if ( l > 23 && k <= 71 )	fputc(',',fp2);
		if ( k > 71 ) { fputs("\n      double precision ",fp2); l = 23;
				k = l+i+2; }
		for( j=0; j < i; j++ ) fputc(toupper(fncname[j+offset]),fp2);
		fputc('R',fp2);	l = k; n = n+1;	if ( n > *nfnc )  {
		    printf(">>> ERROR >>> Total FNC number > %d\n",*nfnc);
						*ierr=-2; return; }
	     }
	}
	*ierr = 0;	*nfnc = n;	fputc('\n',fp2);
	fclose(fp1);	fclose(fp2);	remove("tmp/fnc.list");
}

void getfnc_(int *ierr,int *nfnc, char *fncnames,int *lnam)
{

  char *FList[512];
  char FBuff[512][16];
  char *lc1,*lc2;
  int i,k,L,n,nFList,NList=512;
  FILE *fp2;


  if(*lnam != MxLen){
    printf(">>> ERROR >>> Inconsistency in array FNCNAM length\n");
    *ierr	=-1;
    return;
  }
  for(i=0; i < NList; i++){
    FList[i]	=FBuff[i];
  }
  nFList	=FileList(FList,NList,"fnc");
  if(*nfnc < nFList){
    printf(">>> ERROR >>> Total number of functions > %d\n",*nfnc);
    *ierr=-2;
    return;
  }

  fp2	=fopen("tmp/declar.fnc","w");
  if(fp2 == NULL){
    printf(">>> ERROR >>> Cannot write to the working directory\n");
    *ierr	=-1;
    return;
  }
  fputs("      double precision VINT,IINT,GRAD,GRADS,FRMAX,FRMIN,RFMIN,RFMAX\n",
      fp2);
  fputs("      double precision RFVAL,AFVAL,RFVEX,AFVEX,RFVIN,AFVIN,RFA\n",
       fp2);
  fputs("      double precision RFAN,XFA,XFAN,AFR,AFX,RECR,ATR,ATX\n",fp2);
  fputs("      double precision TIMINT,TIMDER,TIMAVG,GAUSS\n",fp2);
  fputs("      double precision RADIAL,RADINT,ASTEP,RSTEP,XSTEP,STEP,CUT\n",
       fp2);
  fputs("      double precision FTBOX,FXBOX,FABOX\n",fp2);
  fputs("      double precision FIXVAL,FTAV,FTMIN,FTMAX,FRAMP,FJUMP\n",fp2);
  fputs("      external IINT\n",fp2);
  n	=0;
  L	=0;
  lc2	=fncnames;
  for(i=0; i < nFList; i++){
    lc1	=FList[i];    k	=0;
    /*    printf("Nm[%2d]=<%s>\n",i,FList[i]);	*/
    while(*lc1 != '.' && *lc1 != '\0' && k < MxLen){
      *lc2++	=toupper(*lc1);
      lc1++;
      k++;
    }
    if(*lc1 != '.'){
      printf("%s function name is longer than %d\n",FList[i],MxLen);
      exit(0);
    }
    if(L+k+1+n < 72){
      if(n){
	putc(',',fp2);
	L++;
      }
    } 
    else{
      putc('\n',fp2);
      n	=0;
    }
    if(n == 0){
      fputs("      double precision ",fp2);
      L	=23;
    }
    lc1	=FList[i];
    while(*lc1 != '.' && *lc1 != '\0'){
      putc(toupper(*lc1),fp2);
      lc1++;
      L++;
    }
    putc('R',fp2);
    L++;
    n	=1;
    while(k < MxLen){
      *lc2++	=' ';
      k++;
    }
  }
  if(n){
    putc('\n',fp2);
  }
  fclose(fp2);
  *ierr =0;
  *nfnc =nFList;
  return;
}
