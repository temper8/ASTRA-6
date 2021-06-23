/* Here are C-functions:
	viewmn		If_empty	GetKey		GetKeySym
	asrumn		pollevent	asktab
	waitevent	nextevent	askcol		askgrf
	Cursor_in_Box	Menu_table 	Put_button	asklis	
	askuna		GetName		setcol		xaxis
	GetValue	FindBoxNum	GetEsc		Root_window_event
	Del_Menu	Del_button	getatr				*/
#include 	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
/*#include 	<ctype.h> */
#include 	<unistd.h>
#include	<X11/Xlib.h>
#include	<X11/Xutil.h>
#include	<X11/keysym.h>
#include	<X11/keysymdef.h>
/*#include	<signal.h>*/
/*#include	<math.h>*/
/*#include 	<X11/cursorfont.h>*/
void	asrumn(int*);			void	asrumn_(int*);
void	xaxis(int*);			void	xaxis_(int*);
void	viewmn(int*);			void	viewmn_(int*);
void	mvcursor(int*,int*,int*);	void	mvcursor_(int*,int*,int*);
void	ProcessRootWindowEvent (XEvent*);
int	setcol(int*);			int	setcol_(int*);
int	pollevent(int*);		int	pollevent_(int*);
int	waitevent(int*,int*,int*);	int	waitevent_(int*,int*,int*);
int	Cursor_in_Box();
int	Root_window_event(int*);
/*int	GetKeySym (XKeyEvent); */
extern	Display	*theDisplay;
extern	Window	theRootWindow;
extern	Cursor	theRootCursor, theMenuCursor;
extern	GC	theGCA, hghGC, hgh_menuGC, hintGC;
extern	int	theScreen, theDepth, iconState, ColorNum;
extern	int	XWX, XWY, XWW, XWH, XCorrection, YCorrection, theHeight;
int	XRW, YRW;
extern	int	Xmode;
extern	int	AstraColorNum[];
static	int	iact = -1, ibcursor = -1;
		      /*#define R7 65496
			#define R9 65498
			#define R13 65502
			#define R15 65504*/
#define ROOT_EV_MASK  (  ButtonPressMask | KeyPressMask | EnterWindowMask |\
  StructureNotifyMask |  FocusChangeMask | ExposureMask | LeaveWindowMask )
#define POLL_EV_MASK (ButtonPressMask | KeyPressMask | ExposureMask |\
  StructureNotifyMask | FocusChangeMask | EnterWindowMask | LeaveWindowMask)
#define GET_EV_MASK (POLL_EV_MASK | PointerMotionMask | ButtonMotionMask |\
  ButtonReleaseMask)
#define WRITE_ XDrawImageString(theDisplay,theWindow,
#define MVPOINTER_ XWarpPointer(theDisplay,None,theWindow,0,0,0,0,
#define RETURNPOINTER_ XWarpPointer(theDisplay, None, theRootWindow, 0,0,0,0, Kevent.xcur,Kevent.ycur)
#define WRITER_ XDrawImageString(theDisplay,theRootWindow,
#define DLINER_ XDrawLine(theDisplay,theRootWindow,
#define LINGCA_ XDrawLine (theDisplay,theWindow,theGCA,
#define CLREC_ XClearArea(theDisplay,theRootWindow,
#define	F1sw	8			/*font 1 symbol width		*/
#define	F1sh	13			/*font 1 symbol height		*/
#define	BTyof	0			/*button text vertical offset	*/
#define	BTxof	0			/*button text horizontal offset	*/
#define	BVers	4			/*button vertical separation	*/
#define NMRMB   36			/*main Review menu table	*/
int	nmamb=NMRMB;			/*main Astra  menu table       	*/
char	*MAMT[NMRMB];			/*   Dimension is defined as	*/
char	MAMK[NMRMB];			/*	max(nmamb,NMRMB)	*/
#define Button struct BUTTON
Button 	{int mexl, meyd, mexr, meyu, melen; char keysym, *mename;};
Button	MAMB[NMRMB];
#define NDW     12			/*number of dialog windows	*/
char	*DWTitle[NDW] =	
  {	"Variable control",	/*No. 1  "V"*/
	"Constant control",	/*No. 2  "C"*/
	"Times & Grids",	/*No. 3  "D"*/
	"Window control",	/*No. 4  "W"*/
	"Scale control",	/*No. 5  "S"*/
	"Y-shift",		/*No. 6  "Y"*/
	"Time interval",	/*No. 7  "M" in the mode 7*/
	"Mark times:  < 0 - skip,  0 - dim,  > 0 - color #", /*"M" modes 4,5*/
	"Equilibrium control",	/*No. 9  call from metric */
	"1D_Ufile",		/*No. 10 Not used */
	"2D_Ufile",		/*No. 11 Not used */
	"NBI const for beam No",/*No. 12 call from nbiext */
   };
int	DWnamlen[NDW] =	{
	    6,	/*No. 1  PRNAME*/	6,	/*No. 2  CFNAME*/
	    6,	/*No. 3  DTNAME*/	4,	/*No. 4  NAME[TR]*/
	    4,	/*No. 5  NAME[TR]*/	4,	/*No. 6  NAME[TR]*/
	    6,	/*No. 7  NAM7*/		6,	/*No. 8  NAMEP*/
	    6,	/*No. 9  DTNAME*/	6,	/*No. 10 1D_Ufi*/
	    6,	/*No. 11 2D_Ufi*/	6,	/*No. 12 NBI control*/
			};
struct	{int xcur, ycur;}	Kevent;
/**********************************************************************/
void xaxis (modex) int *modex;		/* Change of X-axis */
{    xaxis_(modex);	}
void xaxis_(modex) int *modex;		/* Change of X-axis */
{  
  switch (*modex)
    {	case 0:	MAMT[0]="16*f(a)  "; MAMT[1]=" 8*f(a)    ";
		break;
	case 1:	MAMT[0]="16*f(a)  "; MAMT[1]=" 8*f(a)    ";
		break;
	case 2:	MAMT[0]="16*f(rho)"; MAMT[1]=" 8*f(rho)  ";
		break;
	case 3:	MAMT[0]="16*f(psi)"; MAMT[1]=" 8*f(psi)  ";
		break;
    }
  return;
}
/********** Draw lower part of the root window ************************/
/* The function defines menu table in the Astra run mode */
void asrumn (modex)	int *modex;
{    asrumn_(modex);	}
void asrumn_(modex)	int *modex;
{   int	Xx, Xy, i, ixx=8, iyy;
    int	BHors=5;		/*button horizontal separation	*/
    char *BUTEXT[NMRMB] = 
/*	{" 16*f(a)   "," 8*f(a)   "," Refresh  ",
	 " 2*f(a,t)  "," 8*f(psi) ","User graph",
	 " 2*f(R,t)  "," 8*f(t)   ","  Next    ",
	 "Phase space","Eqilibrium"," Backward ",
*/	{"16*f(a)  "," 8*f(a)    "," Refresh  ",
	 "2*f(a,t) "," 8*f(t)    ","User graph",
	 "2*f(R,t) ","Phase space","  Next    ",
	 "8*f(psi) ","Eqilibrium "," Backward ",
			"  Scales  ","Variables "," Type data ","Port_PS",
			" Windows  ","Constants ","Save tuning","Land_PS",
			"  Select  ","  Grids   ","Write data ","U-files",
			"  Style   ","Type model","What X-axis","Y-shift",
			"Run ","Step","Quit","Help"};
	char	BUTKEY[NMRMB]={	'1','2','R',
				'4','6','9',
				'5','7','N',
				'3','8','B',
				'S','V','T','G',
				'W','C','I','Q',
				'M','D','F','U',
				'.','L','X','Y',
				'\015','\040','\057','H'};
	nmamb = 32;
	/* Draw separating lines between plots and menu */
	Change_Color (theGCA, AstraColorNum[2], AstraColorNum[3]);
	DLINER_ theGCA, 0, XWH-128, XWW-1, XWH-128);
	DLINER_ theGCA, 0, XWH-110, XWW-1, XWH-110);
	DLINER_ theGCA, 0, XWH-109, XWW-1, XWH-109);
	for ( i=0; i<nmamb; i++ )
	    { MAMT[i] = BUTEXT[i];	MAMK[i] = BUTKEY[i];	}
	xaxis (modex);	iyy=XWH-74;
	Menu_table (ixx, iyy, 12, 3, MAMK, MAMT, MAMB, BHors);
	Menu_table (ixx+272, iyy, 16, 4, MAMK+12, MAMT+12, MAMB+12, BHors);
	Menu_table (ixx+610, iyy, 4, 1, MAMK+28, MAMT+28, MAMB+28, BHors);
	if (ibcursor >= 0) Put_button(MAMB[ibcursor],hgh_menuGC);
/* Menu titles */
	Change_Color (hghGC, AstraColorNum[2], AstraColorNum[3]);
	Xx=ixx+40;	Xy=iyy-5;
	XDrawImageString(theDisplay,theRootWindow,hghGC,Xx,Xy,
	"Graphic mode           Presentation       Control",49); Xx=Xx+475;
	XDrawImageString(theDisplay,theRootWindow,hghGC,Xx,Xy,
	"In/Out     Status",17);
}
/**********************************************************************/
void viewmn (modex)	int *modex;
{	viewmn_ (modex);	}
void viewmn_(modex)	int *modex;
{	int	Xx, Xy, i;
   int	ixx=8, iyy, BHors=8;	 /*button horizontal separation	*/
   char	*BUTEXT[NMRMB]={"16*f(a)  "," 8*f(a)    ","Refresh","Scales ",
			"2*f(a,t) "," 8*f(t)    ","Y-shift","Windows",
			"2*f(R,t) ","Phase space"," Next  ","Select ",
			"8*f(psi) ","Eqilibrium "," Back  ","Style  ",
			"PS-portr. ","Type data",
			"PS-landsc."," U-files ",
			"Write data","Variables",
			"Model info","Constants",
			" >"," End ",
			"< ","Start",
			">>","Load>",
			"<<","Load<",
			"Grids ","X-axis"," Help "," Quit "};
   char	BUTKEY[NMRMB]={	'1','2','R','S',
			'4','6','Y','W',
			'5','7','N','M',
			'3','8','B','.',
			'G','T','Q','U',
			'F','V','L','C',
			'\000','\000','\000','\000',
			'\000','\000','\000','\000',
			'D','X','H','\033'};
   /* Draw separating lines between plots and menu */
	Change_Color (theGCA, AstraColorNum[2], AstraColorNum[3]);
	DLINER_ theGCA, 0, XWH-128, XWW-1, XWH-128);
	DLINER_ theGCA, 0, XWH-110, XWW-1, XWH-110);
	DLINER_ theGCA, 0, XWH-109, XWW-1, XWH-109);
	for ( i=0; i<nmamb; i++ )
	    { MAMT[i] = BUTEXT[i];	MAMK[i] = BUTKEY[i];	}
	xaxis (modex);	iyy=XWH-74;
	Menu_table (ixx, iyy, 16, 4, MAMK, MAMT, MAMB, BHors);
	Menu_table (ixx+311, iyy, 8, 2, MAMK+16, MAMT+16, MAMB+16, BHors);
	Menu_table (ixx+490, iyy, 8, 2, MAMK+24, MAMT+24, MAMB+24, BHors);
	Menu_table (ixx+575, iyy, 4, 1, MAMK+32, MAMT+32, MAMB+32, BHors);
	Change_Color (hghGC, AstraColorNum[2], AstraColorNum[3]);
	Xx=ixx+35;	Xy=iyy-5;
	XDrawImageString(theDisplay,theRootWindow,hghGC,Xx,Xy,
	"Graphic mode       Presentation         Data output",51); Xx=Xx+454;
	XDrawImageString(theDisplay,theRootWindow,hghGC,Xx,Xy,
	"Time scan   Info",16);
}
/**********************************************************************/
Menu_table (xm, ym, nbutt, nclmn, mek, met, butt, BHors)
		 int nclmn, xm, ym, nbutt, BHors;
		 char mek[], *met[]; Button butt[];
{	int		i, leng, lenm, ic = 1;
	for (ic=0; ic<nclmn; ic++)			{ 
	    if(ic) { lenm=0;	for (i=ic-1; i<nbutt; i=i+nclmn) 	
			{ leng=strlen(met[i]); if(leng>lenm) lenm=leng; } }
	    for (i=ic; i<nbutt; i=i+nclmn)	{ 
		leng=strlen(met[i]);	
		if(ic) butt[i].mexl=butt[i-1].mexl+lenm*F1sw+2*BTxof+3+BHors;
		else	butt[i].mexl=xm;
		butt[i].mexr=butt[i].mexl+leng*F1sw+2*BTxof+3;
		butt[i].meyu=ym+((i-ic)/nclmn)*(F1sh+2*BTyof+1+BVers);
		butt[i].meyd=butt[i].meyu+F1sh+2*BTyof+1;
		butt[i].melen=leng;
		butt[i].mename=met[i];
		butt[i].keysym=mek[i];	}	}
	for (i=0; i<nbutt; i++)	Put_button(butt[i],theGCA);
}
/*************** Button drawing ************************************/
Put_button (but, aGC) Button but; GC aGC;
{	int	yd, xd, yu, xu;
	yd = but.meyd;	xd = but.mexl;
	yu = but.meyu;	xu = but.mexr;
	Change_Color(aGC,1,0);
	CLREC_ xd-2, yu-2, xu-xd+4, yd-yu+4, False);
	WRITER_ aGC,xd+BTxof+2,yd-BTyof-3,but.mename,but.melen); 
	DLINER_ aGC,xd+2,yd,xu-2,yd);	DLINER_ aGC,xu-2,yd,xu,yd-2);
	DLINER_ aGC,xu,yd-2,xu,yu+2);	DLINER_ aGC,xu,yu+2,xu-2,yu);
	DLINER_ aGC,xu-2,yu,xd+2,yu);	DLINER_ aGC,xd+2,yu,xd,yu+2);
	DLINER_ aGC,xd,yu+2,xd,yd-2);	DLINER_ aGC,xd,yd-2,xd+2,yd);
	/* No round caps:
	DLINER_ aGC,xd,yd,xu,yd);	DLINER_ aGC,xu,yd,xu,yu);
	DLINER_ aGC,xu,yu,xd,yu);	DLINER_ aGC,xd,yu,xd,yd);*/
}
/**********************************************************************/
void AstraEvent_ ()
{
  XEvent theEvent;
    if (!Xmode) return;
    XSelectInput (theDisplay, theRootWindow, ROOT_EV_MASK);
    if( XEventsQueued(theDisplay, QueuedAfterReading) != 0 )
    { 
	XNextEvent (theDisplay, &theEvent);
	if( theEvent.xany.window == theRootWindow )
	    ProcessRootWindowEvent (&theEvent);
	XFlush(theDisplay);
    }
}
/**********************************************************************/
void AstraEvent ()
{ if (!Xmode) return;	(void)AstraEvent_ (); 
}
/**********************************************************************/
void astraevent_ ()
{ (void)AstraEvent_ (); }
/**********************************************************************/
int pollevent (key) int *key;
{	pollevent_ (key);	}
int pollevent_(key) int *key;
{ 	*key = 0;	return Root_window_event(key);
}
/***************** Polling for events *********************************/
int Root_window_event(key)	int *key;
{	XEvent		theEvent;
	XKeyEvent	theKeyEvent;
	XComposeStatus	theComposeStatus;
	KeySym		theKeySym;
	int		i, count, mode = QueuedAfterReading;
	int		theKeyBufferMaxLen = 4;
	char		ch, theKeyBuffer[5];

   XSelectInput (theDisplay, theRootWindow, ROOT_EV_MASK);
   /*	printf("\nXSelect done,   "); */
   count = XEventsQueued(theDisplay, mode);
   /*	printf("EventsQueued:  %d\n",count); */
   if ( count == 0 )
     { XFlush(theDisplay); *key = 0;	return 0; }
XNextEvent (theDisplay, &theEvent);
   /*	printf("Event.type:  %d\n",theEvent.type); */
   switch( theEvent.type )
     {	case ButtonPress:				/* Linux 4 */
	     Kevent.xcur = theEvent.xbutton.x;	
	     Kevent.ycur = theEvent.xbutton.y;
	     i=Cursor_in_Box ();
	     if (i>=0 && i<= nmamb)	*key=(int)MAMK[i];
	     if (*key=='\033') { *key='c'; return 1; /* <Ctrl>+C */}
	     if (theDepth==1 &&  *key=='A') *key=0;
	     return 0;
	case KeyPress:					/* Linux 2 */
	     Kevent.xcur = theEvent.xbutton.x;	
	     Kevent.ycur = theEvent.xbutton.y;
	     theKeyEvent = theEvent.xkey;
	  /* count = XLookupString (&theEvent.xkey, theKeyBuffer,
		theKeyBufferMaxLen, &theKeySym, &theComposeStatus);*/
	     count = XLookupString (&theKeyEvent, theKeyBuffer,
		theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
	     if ( count > theKeyBufferMaxLen ) 
	        { printf("Keyboard error:  Bufferlength =%d\n",count);
			break; }
	     if ( theKeyEvent.state & Mod1Mask)
	       { if ( theKeySym <= 127 )	{ ch = theKeySym;
	       /*printf("<Alt>+%s is pressed,   ASCII code =%d\n",
		 theKeyBuffer,theKeySym);	*/
		 *key = ch;	return 2;	}
	       }
	     if ( theKeyEvent.state & ControlMask)
	       { if ( isascii(theKeySym) )	{ ch = theKeySym;
	       /* printf("<Ctrl>+%c is pressed,   ASCII code =%d\n",
		 ch,theKeySym);			*/
		 *key = ch;	return 1;	}
	       }
	     theKeyBuffer[count]='\0';		*key = theKeyBuffer[0];
	     if (theDepth==1 && (*key=='A' || *key=='a') ) *key=0;
	       /* printf("%d  %d  %d\n",theKeyEvent,count,theKeySym); */
	       /* printf("%d  %d  %d  %d  %d\n",Mod1Mask,Mod2Mask
				,Mod3Mask,Mod4Mask,&theEvent); */
	     return 0;
	case Expose:
	     if(theEvent.xexpose.count == 0)
	     *key=(int)'R';	return 0;
	case UnmapNotify:
	  /*	     printf("UnmapNotify  %d\n",theEvent.type); */
	     *key = 0;		return 65005;
	case FocusIn:
	  /*	     printf("FocusIn     %d\n",theEvent.type); */
	     *key = 0;		return 65006;
	case EnterNotify:
	  /*	     printf("EnterNotify %d\n",theEvent.type); */
	     *key = 0;		return 65006;
	case MapNotify:
	  /*	     printf("MapNotify   %d\n",theEvent.type); */
	     *key = 0;		return 0;
	case FocusOut:
	  /*         printf("theRootWindow is out of focus\n\n"); */
	     *key = 0;		return 0;
	case ConfigureNotify:				/* Linux 22 */
	  /* printf("Configure Notify Event\n"); */ /* Window movement */
	     *key = 0;		return 0;
       	default:
	  /*	     printf("Event.type %d\n",theEvent.type); */
	     *key = 0;		return 0;
     }
/*    if( XCheckWindowEvent( theDisplay, theRootWindow,
	ButtonPressMask, &theEvent) == True )
	   { Kevent.xcur = theEvent.xbutton.x;	
	     Kevent.ycur = theEvent.xbutton.y;
	     return 0;}
    if( XCheckWindowEvent( theDisplay, theRootWindow,
	KeyPressMask, &theEvent) == True )
	   { XLookupString (&theEvent.xkey, theKeyBuffer, theKeyBufferMaxLen,
		 			&theKeySym, &theComposeStatus);
	     Kevent.kbdk = theKeyBuffer[0];
	     theKeyBuffer[1]='\0';
	     return 0; }
    if( XCheckWindowEvent( theDisplay, theRootWindow,
	ExposureMask, &theEvent) == True )	
	   { if(theEvent.xexpose.count == 0)
	     return 0; }
*/
	*key = 0;			return 0;
}
/************************* presently not used **************************/
GetKeySym (theKeyEvent) XKeyEvent theKeyEvent;
{	XComposeStatus	theComposeStatus;
	KeySym		theKeySym;
	int		theKeyBufferMaxLen = 4, count, key;
	char		ch, theKeyBuffer[5];
	count = XLookupString (&theKeyEvent, theKeyBuffer,
		theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
	if ( count > theKeyBufferMaxLen ) 
	   { printf("Keyboard error:  Bufferlength =%d\n",count);
	     return 0; }
	if ( theKeyEvent.state & Mod1Mask)
	   { if ( theKeySym <= 127 )
		  printf("<Alt>+%s is pressed,   ASCII code =%d\n",
			theKeyBuffer,theKeySym);  return 0;
		}
	if ( theKeyEvent.state & ControlMask)
	        { if ( isascii(theKeySym) )	{ ch = theKeySym;
		     printf("<Ctrl>+%c is pressed,   ASCII code =%d\n",
			ch,theKeySym);		  return 0;	}
		}
	theKeyBuffer[count]='\0';
	key = theKeyBuffer[0];
	/*     printf("key = %d,    theKeySym = %d\n",key,theKeySym); */
	if (theDepth==1 && (key == 'A' || key == 'a') )	return 0;
	return key;
}
/*  switch(theKeySym)
      {	case	XK_Escape:			return (-1);
	case	XK_BackSpace:
	case	XK_Delete:			return (1);
	case	XK_Return:
	case	XK_KP_Enter:			return (2);
	case	XK_Tab:				return (3);
	case	XK_R7:				return (11);
	case	XK_Left:			return (21);
	case	XK_Up:				return (12);
	case	XK_Right:			return (23);
	case	XK_Down:			return (32);
	case	XK_R9:				return (13);
	case	XK_R13:				return (31);
	case	XK_R15:				return (33);
	case	XK_Shift_L:
	case	XK_Shift_R:  			return (34);
	case	XK_Control_L:
	case	XK_Control_R:			return (35);
	case	XK_Meta_L:
	case	XK_Meta_R:   		  	return (36);
	case	XK_Alt_L:
	case	XK_Alt_R:		   	return (37);
	default:
	  if ( !isprint(theKeySym) )		return (99);
	   printf("[%1.1s]%d\n",theKeyBuffer,theKeySym); 
		 str[*pos]=theKeySym;		return (0);
      } */
/**********************************************************************/
void GetRWgeometry (XRW,YRW)	int *XRW, *YRW;
{
   Window	        theRW, theCW;
   int		Wx, Wy, iX, iY;
   unsigned int	Ww, Wh, Wb, Wd;
   /*   Status		theS; */
   /*   XWindowAttributes theAttributes; */
   /* printf("%d %d (x,y)=(%d,%d)\n",*x,*y,XWX,XWY); */
   XGetGeometry(theDisplay,theRootWindow,
			       &theRW, &Wx, &Wy,  &Ww, &Wh,  &Wb, &Wd);
   XTranslateCoordinates(theDisplay, theRootWindow,
			theRW, 0, 0, &iX, &iY, &theCW);
   *XRW = iX-Wx;	*XRW -= XCorrection;
   *YRW = iY-Wy;	*YRW -= YCorrection;
   /* printf("Translate: (x,y)=(%d,%d) -> (x,y)=(%d,%d) get\n",Wx,Wy,iX,iY); */
   /* printf(" -> (x,y)=(%d,%d)\n",Wx+XRW,Wy+YRW); */
   /* printf("{%d,%d}\n",theRW,theCW); */
   /* printf("New position (x,y)=(%d,%d)\n",Wx+XRW,Wy+YRW); */
   /* printf("(x,y)=(%d,%d);   dest_(x,y)=(%d,%d);\n",Wx,Wy,iX,iY); */
   /* theS = XGetWindowAttributes(theDisplay,theRootWindow, */
   /*			 &theAttributes); */
   /* printf("theEvent.type %d\n",theEvent.type); */
   /* printf(" status %d\n",theS); */
   /* printf(" (x,y)=(%d,%d)\n",theAttributes.x,theAttributes.y); */
   /* printf("(width,heigth)=(%d,%d);   ",Ww,Wh); */
   /* printf("(border,depth)=(%d,%d)\n",Wb,Wd); */
   return;
}
/**********************************************************************/
waitevent (theKey, xCursor, yCursor)	int *theKey, *xCursor, *yCursor;
{	int i;
	i = waitevent_ (theKey, xCursor, yCursor);
	return i;
}
/* Waiting events from Astra in WAIT mode & and from dialog windows */
waitevent_ (theKey, xCursor, yCursor)	int *theKey, *xCursor, *yCursor;
{	int i;
  Wait:	i=nextevent(theKey, xCursor, yCursor, nmamb, MAMB, MAMK, MAMT);
	if (i==0 && *theKey==0 ) goto Wait;
	/*printf("next_event: %d   %d\n",i,*theKey); */
	return i;
}
/*********************** Waiting events ********************************/
nextevent(theKey, xCursor, yCursor, nmbt, mbb, mbk, mbt) 
		int *theKey, *xCursor, *yCursor; int nmbt;
		Button mbb[]; char mbk[]; char *mbt[];
	/* ibcursor    -  current box # or -1 */
	/* iact     - highlighted box # or -1 */
{	XEvent		theEvent;
	XKeyEvent	theKeyEvent;
	XComposeStatus	theComposeStatus;
	KeySym		theKeySym;
	int		length, ibox, theKeyBufferMaxLen = 4;
	char		theKeyBuffer[5];
	int		longKey;
XSelectInput (theDisplay, theRootWindow, GET_EV_MASK);
XNextEvent (theDisplay, &theEvent);
/* printf("next_event: %d   %d\n",theEvent.type, theEvent.xany.window); */

switch( theEvent.type )
  { case MapNotify:
	*theKey = 0;
	return 65004;

    case UnmapNotify:
	*theKey = 0;
	return 65003;

    case KeyPress:
 	Kevent.xcur = theEvent.xbutton.x;   *xCursor=Kevent.xcur;
	Kevent.ycur = theEvent.xbutton.y;   *yCursor=Kevent.ycur;
	theKeyEvent = theEvent.xkey;

	/* *theKey = GetKeySym(theKeyEvent);
	   printf("%d,%d,%d\n",*theKey,theKeyEvent,&theKeyEvent);
	   return; */
	length = XLookupString (&theKeyEvent, theKeyBuffer,
		theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
		/*  printf("%d,  %d\n",length,theKeySym); */
	if ( length > theKeyBufferMaxLen || length < 0 ) 
	  { printf("String translation error:  Bufferlength =%d\n",length);
	    break; }
	if ( length == 0 )
	  { if( theKeySym >= 65000 )
	    /* A special key (as <Ctrl>, <Alt>, <F1>, Home, etc) is pressed */
	    { *theKey = theKeySym;
	      switch(theKeySym)
	       {  case	XK_Left:	return	65361;
	          case	XK_Up:		return	65362;
	          case	XK_Right:	return	65363;
	          case	XK_Down:	return	65364;
	          default:		return	0;
	       }
	     } 
	  else
	    { *theKey = theKeySym;	return 65000;     }
	  }
	if ( theKeyEvent.state & Mod1Mask)
	  { if ( theKeySym <= 127 )	{ 
	    /*	 printf("<Alt>+%s is pressed,   ASCII code =%d\n",
		 theKeyBuffer,theKeySym);	*/
	    *theKey = theKeySym;	return 2;	}
	  }
	if ( theKeyEvent.state & ControlMask)
	  { if ( isascii(theKeySym) )	{ 
	    /*	 printf("<Ctrl>+%c is pressed,   ASCII code =%d\n",
		 toascii(theKeySym),theKeySym);		*/
	    *theKey =  theKeySym;	return 1;	}
	  }
	longKey	=theKeyBuffer[0];
	theKeyBuffer[1]='\0';
	if( (length > 0 ) )
          { if( (theKeySym >= ' ' ) &&			/* ASCII 32  */
	        (theKeySym <= '~' ) )			/* ASCII 126 */
	      {	*theKey	= longKey;
		if(theDepth==1 && (theKeySym=='A' || theKeySym=='a'))
		   *theKey=0;
		/*printf("Key   %d \n",*theKey);*/
	      }
	    else
	      {	/*printf("ASCII key was hit: [%s]\n",theKeyBuffer);
		  printf("longKey   %d \n",longKey);*/
		*theKey = longKey;
		if( longKey == 8 )  { return 65361; /* BS equiv <- */}
	      }
	    return 65000;
	  }
    case ButtonPress:
      /*printf("ButtonPressMask:  %d,  %d\n", ButtonPressMask, ButtonPress); */
 	Kevent.xcur = theEvent.xbutton.x;   *xCursor=Kevent.xcur;
	Kevent.ycur = theEvent.xbutton.y;   *yCursor=Kevent.ycur;
	ibox = Cursor_in_Box ();
	/* printf("ibox:  %d, %d, %s\n", ibox, (int)mbk[ibox], mbt[ibox]); */
	/* printf("Button  pressed  x = %d,  y = %d,  key = %d\n"
			,*xCursor,*yCursor,*theKey); */
	if (ibox == -1) return 65001;
	if (ibox >= 0) *theKey=(int)mbk[ibox];
	if (theDepth==1 && *theKey=='A') *theKey=0;    /*No color table*/
	if (*theKey == 015) Put_button(mbb[ibox],theGCA); /*Undo highlight*/
	if (*theKey != 0)   ibox = -1;
	if (ibox == 24)	return 65363;	/*   >    */
	if (ibox == 25)	return 65502;	/*  End   */
	if (ibox == 26)	return 65361;	/*   <    */
	if (ibox == 27)	return 65496;	/*  Home  */
	if (ibox == 28)	return 65364;	/*   >>   */
	if (ibox == 29)	return 65504;	/*  PgDn  */
	if (ibox == 30)	return 65362;	/*   <<   */
	if (ibox == 31)	return 65498;	/*  PgUp  */
	return 65001;

    case ButtonRelease:		
      /* printf("ButtonReleaseMask: %d, %d\n",ButtonReleaseMask,ButtonRelease); */
 	Kevent.xcur = theEvent.xbutton.x;   *xCursor=Kevent.xcur;
	Kevent.ycur = theEvent.xbutton.y;   *yCursor=Kevent.ycur;
	/* printf("Button released  x = %d,  y = %d  key = %d\n"
			,*xCursor,*yCursor,*theKey); */
	return 65002;

    case MotionNotify:
      /*     	printf("Case Motion\n"); */
 	Kevent.xcur = theEvent.xbutton.x;   *xCursor=Kevent.xcur;
	Kevent.ycur = theEvent.xbutton.y;   *yCursor=Kevent.ycur;
	/*printf("Cursor position   x=%d,  y=%d\n",*xCursor,*yCursor);*/
	ibcursor = Cursor_in_Box ();
	if (iact == ibcursor )  return 65000;	/* No changes */
	if (iact >= 0) 				/* Hgh_box -> std */
	   { Put_button(mbb[iact],theGCA);  iact = -1;}
	if (ibcursor >= 0)			/* Std_box -> hgh */
	   { Put_button(mbb[ibcursor],hgh_menuGC); iact = ibcursor; }
	*theKey = 0;
	return 65000;

    case Expose:
	if(theEvent.xexpose.count == 0)	longKey = 'R';	/* ASCII 82 */
	else				longKey = 0;
	*theKey	=longKey;
	return 65000;

    case ConfigureNotify:				/* Linux 22 */
      /* printf("Configure Notify Event\n");*//* Occurs when moving a window */ 
	*theKey = 0;	     return 0;

/*    case EnterWindow:
	printf(" EnterWindowMask #:  %d  \n", theEvent.type);
	return 65000;
    case LeaveWindow:
	printf(" LeaveWindowMask #:  %d  \n", theEvent.type);
	return 65000;
    case StructureNotify:
	printf(" StructureNotifyMask #:  %d  \n", theEvent.type);
	return 65000;
    case DestroyNotify:
	printf(" DestroyNotifyMask #:  %d  \n", theEvent.type);
	return 65000;
    default:
	longKey=theEvent.type;
	if(longKey == 18)	*theKey	=318;
	if(longKey == 19)	*theKey	=319;
	printf(" default  %d  \n", theEvent.type);
	if(longKey == 18 || longKey == 19 ) return 65000;
*/
  }
return 0;
}
/****** Returns:  menu box number, -1 if no match; 
        Makes use of global structures 
	Button	MAMB[NMRMB];	Event	Kevent;		******/
int Cursor_in_Box ()
{	int	i;
	for(i=0; i<nmamb; i++)
	{   if (Kevent.xcur>=MAMB[i].mexl && Kevent.xcur<=MAMB[i].mexr &&
		Kevent.ycur>=MAMB[i].meyu && Kevent.ycur<=MAMB[i].meyd  )
		return i;
	}
	return -1;
}
/**********************************************************************/
void mvcursor (key,ix,iy) int *key, *ix, *iy;
{	mvcursor_ (key,ix,iy);	}
/* Move cursor by one pixel */
void mvcursor_ (key,ix,iy) int *key, *ix, *iy;
{	if (*key == 361) *ix = *ix-1; /* <- */
	if (*key == 362) *iy = *iy-1; /* Up */
	if (*key == 363) *ix = *ix+1; /* -> */
	if (*key == 364) *iy = *iy+1; /* Dn */
	XWarpPointer(theDisplay, None, theRootWindow, 0,0,0,0, *ix,*iy);
}
/*******    ASKLIS is called by IFKEY (keys "S,Y,D,C,V,M[modes 4,5,7])
************************************************************************/
/*	Example call from FORTRAN
	if (NEQL .eq. -1)	then		! Option (2)
		 NEQL = int(NEQUIL)
		 CHAR6 = DTNAME(17)
		 YACC = DELOUT(17)
		 DELOUT(17) = ACC
		 DTNAME(17) = "Toler"
		 call ASKLIS(4,DELOUT(17),DTNAME(17),9)
		 DTNAME(17) = Char6
		 ACC = DELOUT(17)
		 DELOUT(17) = YACC
C		 write(*,*)"Time =",TIME
		 goto	3
	endif
*/
int  asklis (nofbox, array, theNames, id)
     int    *nofbox, *id;
     char   theNames[]; double *array;
{  int	i;
	i = asklis_(nofbox, array, theNames, id);
	return i;
}
int  asklis_ (nofbox, array, theNames, id)
    int  *nofbox, *id;
    char theNames[];
    double *array;
{	Window		theWindow;
	XEvent		theEvent;
	double		valn;
	int	UpLeftx =2, UpLefty =10,/* window corner location	*/
		Width =452, Height =480;/* window size			*/
	int	ix, iy, wbox, hbox, 	/* parameter box corner and size*/
		wsym = 8, hsym = 13,	/* symbol width and height	*/
		lvalue = 6, lname = 6,	/* length of value and name	*/
		xshif = 5, yshif = 3,	/* table corner 		*/
		nclmn = 4, nparam, 	/* # of columns and parameters	*/
		spos  = 0, vpos, 	/* edit & edit start positions	*/
		ixa0=0, iya0=0, ixa,iya,/* arrow (textcursor) position	*/
		namlen, ibox, iret, oldparam, ixold, iyold;
	int	lline, xButton, yButton, i, ii=-1, ind, ind1,
		ihelp, selalb=0,  contr, icol, dcol, irow, drow;
	float	param;
	char	value[10],  ovalue[10],  stri[10],  vsym = '=',title[70],
		stri50[50],stri40[40], stri256[256], theName[10], legend[70];
	i = *id-1;	namlen = DWnamlen[i];	strcpy(title,*(DWTitle+i));
/*	printf("theID = %d,  theTitle = %s[%d],   \t theNameLength = %d\n",
			*id, title, strlen(title),namlen);
*/	lline = lname + lvalue + 1; 
	if ((lvalue == 0) || (lname == 0)) lline = lname+lvalue;
	hbox = hsym+3;		wbox = wsym*(lline+1);	vpos= wsym*(lname+1);
 	selalb = 1;	nparam   = *nofbox;	contr=1;
	if( nparam==0 ) return (-1);
	if( nparam<0 ) { nparam = -nparam; selalb = 0; }
	if( nparam>1000 ) { nparam=nparam-1000;	contr=0; }
	if( *id == 3 )  contr=1;
	Width =2*xshif+nclmn*wbox-wsym; 
	Height =2*yshif+((nparam-1)/nclmn+2)*hbox;
	ibox = 1;	if (contr == 0 ) ibox = 0;
       	if( *id == 8 )  { ibox = 0;	selalb = -1; }

	GetRWgeometry (&XRW,&YRW);
	/*	UpLeftx = 2;     i = XWX-Width-10; */
	UpLeftx = 2;     i = XRW-Width-8;		UpLefty = YRW;
	if (i > UpLeftx)  UpLeftx = i;		i = YRW+Height+30;
	if (i > theHeight) UpLefty = theHeight-Height-30;
	theWindow 
	  = Open_Window (UpLeftx, UpLefty, Width, Height, 0, title, 0,
	     RootWindow(theDisplay,theScreen), theMenuCursor);
	XSelectInput (theDisplay, theWindow, POLL_EV_MASK);	ihelp = 0;
Create_table:
	Change_Color (theGCA, 1, 0);	/* white background, black foreground */
	for (i=1; i <= nparam; i++) 	{
	ix =wbox*((i-1)%nclmn)+xshif;	iy =hbox*((i-1)/nclmn)+yshif;
	if ( lname ) {	ind = (i-1)*namlen;
		WRITE_ theGCA, ix, hsym+iy, theNames+ind, namlen); }
	if ( lname*lvalue ) { 
		WRITE_ theGCA, ix+vpos-wsym, hsym+iy, &vsym, 1); }
	if ( lvalue ) { valn=*(array+i-1); numstrA(valn, value, lvalue); 
		WRITE_ theGCA, ix+vpos, hsym+iy, value, lvalue); }
	XFlush(theDisplay);		}
	oldparam = nparam;	ixold = ix;	iyold = iy; 
	for (i=0; i<lvalue; i++) 	ovalue[i]=value[i];
	if (ibox == 1 && selalb != -1 && ii < 0) MVPOINTER_ xshif+40,yshif+13);
	Change_Color (theGCA, 50, 0);
	for (i=1; i < nclmn; i++) { ind = xshif+i*wbox-wsym/2;
	LINGCA_ ind,0,ind,Height-hbox); }
	i = Height-hbox; LINGCA_ 0,i,Width,i);
	i--; 		 LINGCA_ 0,i,Width,i);
	WRITE_ hghGC, xshif, Height-3, "OK",2);
	if ( ibox == 0 && contr == 0 ) MVPOINTER_ xshif+5,Height-5);
       	if ( *id == 8 )  MVPOINTER_ xshif+Width-50,Height-5);
	i = nclmn*wbox/wsym-5;	ind = 50;	if(i<50)	ind=i;
	strcpy(legend,"/ Esc - done;   ");
	if (*id <= 6) {	strcat(legend,"Button / Tab - select;");
			strcat(legend,"  Return/Tab - enter;");
			strcat(legend,"  ? - quick help");		}
	else	      {	strcat(legend,"     Button / Tab - select;");
			strcat(legend,"     Return/Tab - enter;");	}
	Change_Color (hintGC, 42, 0);	ind = strlen(legend);
	if (contr == 1) WRITE_ hintGC, xshif+2*wsym, Height-3, legend, ind);
	Change_Color (hintGC, 1, 0);
	if (contr == 0) WRITE_ theGCA, xshif+4*wsym, Height-3, 
	"    Information table.   No changes permitted.     ", ind);
	if ( *id == 8 )	{ 
	      if (selalb == 0)	WRITE_ hintGC,
			xshif+8+(i-6)*wsym,Height-3," (Select all ) ",15);
	      if (selalb == -1)	WRITE_ hintGC, 
			xshif+8+(i-6)*wsym,Height-3,"(Unselect 0)",12);
			}
	Change_Color (theGCA, 1, 0);
	if (ibox) goto Newparam;
Table_control:
	XNextEvent (theDisplay, &theEvent);
	if ( theEvent.xany.window == theRootWindow )
           { /*  printf("theRootWindowEvent.type =%d\n",theEvent.type); */
	     ProcessRootWindowEvent (&theEvent);
	     goto Table_control;
	   }
	if(theEvent.type == Expose)	{ ii = 0;  goto	Create_table;}
	if(theEvent.type == ButtonPress)
	   {	xButton	= theEvent.xbutton.x;
		yButton = theEvent.xbutton.y;
		ind  = FindBoxNum(xButton-xshif, yButton-Height+hbox+1,
						 	4*wsym, hbox, 1, 1);
		if(contr)
		{
		if (ibox) { oldparam = ibox;	ixold = ix;	iyold = iy; 
			if (spos>0)  { 	sscanf(stri,"%6g",&param);
					*(array+ibox-1)=param; valn=param; 
					numstrA(valn, value, lvalue); spos=-1;
				     }
			for (i=0; i<lvalue; i++) 	ovalue[i]=value[i];
			  }
		ibox = FindBoxNum(xButton-xshif+wsym/2, yButton-yshif,
				  wbox, hbox, nclmn, nparam);
		ind1 = FindBoxNum(xButton-xshif-nclmn*wbox+10*wsym,
				yButton-Height+hbox+1, 10*wsym, hbox, 1, 1);
		if(ind1!=0 && selalb!=1)
		   { for (i=0; i<nparam; i++) 
			{ if( *(array+i)<=0 ) *(array+i)=selalb; }
			{ if(selalb)  selalb=0; else selalb=-1; }
			  goto Create_table;
		   }
		}
		if( ind ) { iret = -1;	ibox = oldparam;  goto Escend; }
		if (ibox == 0)	goto Table_control;
Newparam:
		ix =wbox*((ibox-1)%nclmn)+xshif;
		iy =hbox*((ibox-1)/nclmn)+yshif;
		if ( lname ) {ind = (oldparam-1)*namlen; ind1 = (ibox-1)*namlen;
		    WRITE_ theGCA, ixold, hsym+iyold, theNames+ind, namlen);
		    Change_Color(hghGC,AstraColorNum[8],AstraColorNum[9]);
 		    WRITE_ hghGC, ix, hsym+iy, theNames+ind1, namlen);
		    Change_Color(hghGC,AstraColorNum[2],AstraColorNum[3]); }
		if ( lname*lvalue ) { 
		    WRITE_ theGCA, ixold+vpos-wsym, hsym+iyold, &vsym, 1); 
		    Change_Color(hghGC,AstraColorNum[8],AstraColorNum[9]);
		    WRITE_ hghGC, ix+vpos-wsym, hsym+iy, &vsym, 1);
		    Change_Color(hghGC,AstraColorNum[2],AstraColorNum[3]); }
		if ( lvalue ) {valn=*(array+ibox-1); numstrA(valn,value,lvalue);
		    WRITE_ theGCA, ixold+vpos, hsym+iyold, ovalue, lvalue); 
		    Change_Color(hghGC,AstraColorNum[8],AstraColorNum[9]);
		    WRITE_ hghGC, ix+vpos, hsym+iy, value, lvalue); 
		    Change_Color(hghGC,AstraColorNum[2],AstraColorNum[3]);
		    spos = -1; for (i=0; i<lvalue; i++) stri[i]=' '; 
		    ixa = ix+vpos;	iya = iy+hsym+2;
		    MoveArrow(theWindow,ixa0,iya0,ixa,iya);
		    ixa0 = ixa;		iya0 = iya; }
		XFlush(theDisplay);
		goto Table_control;
	   }
	if ((ibox == 0)||(lvalue == 0)) {
			if(GetEsc(theEvent.xkey)) goto EndDialog;
				else	goto Table_control; }
	iret = GetValue(theEvent.xkey, stri, lvalue, "0123456789.-+eE", &spos);
Escend:	if ( contr==0 ){/*printf("iret %d\n",iret);*/	goto EndDialog;}
	if ( spos >= 0 )
	   {	ixa = ix+vpos+spos*wsym;	iya = iy+hsym+2;
		MoveArrow(theWindow,ixa0,iya0,ixa,iya);
	        ixa0 = ixa;			iya0 = iya; }
	if ( iret && (spos>0) )
	   { 	sscanf(stri,"%6g",&param);
		/*printf("input =%s,    double =%g \n",stri,param);*/
		*(array+ibox-1) = param; valn=param;
		numstrA(valn, value, lvalue); spos=-1;	}
	if ( iret == -1 )  	goto EndDialog;
	if ( iret == -2 )
	   { i = 0;	ii = namlen;
Fstblanc:    strncpy(theName,theNames+(ibox-1)*namlen+i,ii);
	     if ( theName[0]  == ' ' )	{ii--;	i++;  goto Fstblanc;}
	     theName[ii] = '\0';
	     while ( theName[--ii] == ' ' ) theName[ii] = '\0';	ii++;
/*	     printf("theName: [%s] ii = %d\n",theName,ii);	
	     for (i=0; i<ii; i++) theName[i] = toupper(theName[i]);
	     printf("theName: [%s] ii = %d\n",theName,ii);	*/
/*	     printf("[%s]\n",stri50); printf(" %s\n",stri256);  */

	     if (*id == 1)		/*Variables <- for/const.in*/
		{ strcpy(stri50,"grep -i \"C ");   strncat(stri50,theName,ii);
	          strcat(stri50," \" for/const.inc");   system(stri50);  }
	     if (*id == 2)		/*Constant table*/
		{ strcpy (stri256,"grep -i -w ");  strncat(stri256,theName,ii);
		  getnames_(stri40,stri50);
		  /* printf("[%s],  %d\n",stri40,strlen(stri40)); */
		  /* printf("[%s],  %d\n",stri50,strlen(stri50)); */
		  /* printf("[%s],  %d\n",stri256,strlen(stri256)); */
   		  strcat (stri256," equ/txt/");
		  strncat(stri256,stri40,strlen(stri40));	  printf(
		    "\nControl constant \"%s\" -> using in the model \"%s\"\n",
		     theName, stri40);		 ii = system(stri256);
		  if (ii == 256) printf("The constant is not used.\n");   }
	     if (*id == 3)		/*Time, grid control*/
		{  strcpy(stri50,"grep -i \"C ");  strncat(stri50,theName,ii);
	           strcat(stri50," \" for/const.inc");   system(stri50);  }
	     if (*id == 4)
		{  strcpy (stri256,"grep -w ");    strncat(stri256,theName,ii);
		   strcat (stri256," tmp/model.txt");
		   if (!ihelp) printf(" No. Scale  Name  Output expression\n");
		   system(stri256);	ihelp = 1;	}
	     if (*id == 5)		/*Scale control*/
		{  strcpy (stri256,"grep -w ");    strncat(stri256,theName,ii);
		   strcat (stri256," tmp/model.txt");
		   if (!ihelp) printf(" No. Scale  Name  Output expression\n");
		   system(stri256);	ihelp = 1;	}
	     if (*id == 6)		/*Vertical shift*/
		{  strcpy (stri256,"grep -w ");    strncat(stri256,theName,ii);
		   strcat (stri256," tmp/model.txt");
		   if (!ihelp) printf(" No. Scale  Name  Output expression\n");
		   system(stri256);	ihelp = 1;	}
	     if (*id == 7)
		{  printf("Curve presentation\n");	}
	     if (*id == 8)
		{  printf("Mark time slices\n");	}
		/* Control print:	printf("[%s]\n",stri50);  */
	   }
	if ( spos < 0 )
	   {	spos = -1; for (i=0; i<lvalue; i++) stri[i]=' '; 
		Change_Color(hghGC,AstraColorNum[8],AstraColorNum[9]);
		WRITE_ hghGC, ix+vpos, hsym+iy, value, lvalue);
		Change_Color(hghGC,AstraColorNum[2],AstraColorNum[3]); }
	else
	   {    Change_Color(hghGC,AstraColorNum[8],AstraColorNum[9]);
		WRITE_ hghGC, ix+vpos, hsym+iy, stri, lvalue);
		Change_Color(hghGC,AstraColorNum[2],AstraColorNum[3]); }
	if ( iret == 2 )
	   {	oldparam = ibox;	ixold = ix;	iyold = iy;
		for (i=0; i<lvalue; i++) 	ovalue[i]=value[i]; 
		ibox++; if(ibox > nparam) ibox=1; goto Newparam; }
	if ( iret > 2 )
	   {	oldparam = ibox;	ixold = ix;	iyold = iy;
		for (i=0; i<lvalue; i++) 	ovalue[i]=value[i];
 		dcol=iret%10-2; 	icol =(ibox-1)%nclmn+dcol;
		if (icol<0) icol=0;	if (icol>=nclmn) icol--;
		drow=iret/10-2;		irow =(ibox-1)/nclmn+drow;
		if (irow<0) irow=0;	if (irow>(nparam-1)/nclmn) irow--;
		ibox =irow*nclmn+icol+1;
		if(ibox == oldparam && icol == nclmn-1 && drow != -1) ibox++;
		if(ibox == oldparam && icol == 0       && drow != 1)  ibox--;
		if (ibox < 1) ibox=1;	if (ibox > nparam) ibox=nparam;
		goto Newparam; }
				goto Table_control;
EndDialog:
	RETURNPOINTER_;
	XSetInputFocus(theDisplay,theRootWindow,
	         RevertToPointerRoot,theEvent.xkey.time);
	XFlush(theDisplay);
	XDestroyWindow(theDisplay, theWindow);
	ibcursor = -1;
	XFlush(theDisplay);
	return (0);
}
/************************************************************************/
int  caution ()
{    int i;	i = caution_();	return i;	}
int  caution_ ()
{
    Window		theWindow;
    int	UpLeftx = 2, UpLefty = 10,/* window corner location	*/
	Width, Height;		/* window size			*/
    int	wsym = 8, hsym = 13,	/* width and height, font 8x13bold */
	margin = 6, i, charLength; 	/* Margin */
    char	stri256[256], title[70] = "Caution!";

    strcpy (stri256,"Writing movie");
    charLength = strlen(&stri256[0]);
    Width = wsym*charLength+2*margin;
    Height =hsym+2*margin;

    GetRWgeometry (&XRW,&YRW);
    UpLeftx = 2;     i = XRW-Width-8;		UpLefty = YRW;
    if (i > UpLeftx)  UpLeftx = i;		i = YRW+Height+30;
    if (i > theHeight) UpLefty = theHeight-Height-30;
/*  printf("stri256 = %s[%d],\t\n",stri256, strlen(&stri256[0]));			
    printf("theTitle = %s[%d],\t\n", title, strlen(title));
    printf("Width = %d,\tHeight = %d,\t%d,\t%d,\t\n", Width, Height);
*/
    printf("XRW = %d,\tYRW = %d,\tUpLeftX = %d,\tUpLeftY = %d,\t \n"
	   ,XRW,YRW, UpLeftx, UpLefty);

    theWindow = Open_Window (UpLeftx, UpLefty, Width, Height, 1, title, 0,
		       RootWindow(theDisplay,theScreen), theMenuCursor);
    i = 30; colorb(&i);
//  Change_Color(hghGC,AstraColorNum[14],AstraColorNum[20]);/* Magenta on Yellow */
    XDrawImageString(theDisplay, theWindow, hghGC, margin, margin+10,
		     stri256, charLength);
    XFlush(theDisplay);
    sleep(5);
EndDialog:
    XDestroyWindow(theDisplay, theWindow);
    XFlush(theDisplay);
    return (0);
}
/********** ASKTAB is called from ASXWIN and ASTWIN (file surv.f)
            that in turn are invoked by key "M" from IFKEY	*********/
int  asktab (title, template, array, len, nrows, ngroup, morow)
     int     *len, *nrows, *ngroup, *morow;
     char    title[], template[], array[];
/*
Input:	title	- Title of the table 
	template  string defining a structure of the table and its 1st line
		  1st line does not appear if all non-'|' symbols are spaces
	array	- data (numbers or strings) for input and output
	len	- length of "array" element according to description 
			in the calling routine (80 in the example below)
	nrows	- number of rows in a table to be created
	ngroup	- if > 0 distance (in rows) between blue separating lines
	morow	- if > 0 separates bottom of the table with a fat blue line

Returned value:
	-1  - Error (window was not created)
	 0  - Normal exit (table updated)
	>0  - Temporary exit (the returned value goes to the calling routine)
Call from FORTRAN
    Example 1: 
	integer		nrow, len
	character*80	title, template, array(10)
	len = 80
	nrow = 5	! must be <= 10
	title = "Example"//char(0)
	template = " field1 |field2|    | f4 |       "//char(0)
	array(1) = " CumBol &  10.1& 2  &    &-1.E+5 "
	array(2) = "                                 "	! Empty string
	array(3) = "                                 "	! Empty string
	array(4) = "                                 "	! Empty string
	array(5) = "                                 "	! Empty string
C symbols in the positions "|" will not appear in the table
C 4th parameter must be 80
C 5th parameter must be <= 10
	call	asktab(title,template,array,80,3,0,0)

More examples in the file "for/surv.f" subroutine "ASXWIN" or "ASTWIN"

Add features:
   If (Width=2*xshif+wsym*(last_pos) < 37*wsym)
       or (Height=2*yshif+(nlines+3)*hbox > theHeight) double Width
*/
{  int	i;
	i = asktab_(title, template, array, len, nrows, ngroup, morow);
	return i;
}
int asktab_ (title, template, array, len, nrows, ngroup, morow)
				int *len, *nrows, *ngroup, *morow;
				char title[], template[], array[];
{	XEvent	theEvent;
	static	Window	theWindow;
	static	int	xshif = 5, yshif = 3,	/* table corner 	*/
		ind, ibold, hbox,	/* current/old box No. & height	*/
		irow0=0, 	wbox, 	/* number of chars in box 	*/
		Width =452, Height =480,/* defaut window sizes		*/
		inold, ixold, iyold, icold, ixa0, iya0, arrdim,
		nclmn, nwid[20], nsta[20], mode=0, ii=0,
		icol, isym, irow, ixa, iya, /* arrow (textcursor) position */
		spos=0,	/* abs. and rel. position of symbol in table*/
		lvalue,		/* 1 if box was changed, 0 otherwise	*/
		wsym = 8, hsym = 13;	/* symbol width and height	*/
	int	UpLeftx, UpLefty =10,	/* window corner location	*/
		nend, ibox, iret, ix, ix1, iy, i, j,
		mxfields = 20;		/* max / actual # of columns	*/
		/* int	count, mode = QueuedAfterReading; */
	static	char	stri[132];
	if ( mode == 1 )  goto	Table_control;		mode = 0;
	j = strlen(title);	if ( j > 132 )			{
		printf("%s\n\"%s\"\n%s\n","ASKTAB >>> Title",title,
		     "           is too long");		return -1;	}
	if  ( *len > 132 )				{
		printf("%s \"%s\" %s\n","ASKTAB >>> Table",title,
		     ". Input string is too long");	return -1;	}
	nend = strlen(template);
	if ( nend > 132 )	{
		printf("%s \"%s\" %s\n","ASKTAB >>> Table",title,
		     ". Requested width is too large");	return -1;	}
	for ( j=0; j < nend ; j++)
	      if ( template[j] != '|'  &&  template[j] != ' ' ) irow0 = 1;

		/****   Split "template" in blocks  ****/
	for ( j=nclmn=nsta[0]=0; j < nend ; j++ )
	    { if ( template[j] == '|' ) 
		 { nwid[nclmn] = j-nsta[nclmn]; nsta[++nclmn] = j+1;}
	      if ( nclmn > mxfields ) {
		printf("%s %d\n","ASKTAB >>>  Too many input fields ",nclmn);
			return -1;	}
	    }
	nwid[nclmn] = nend-nsta[nclmn];		nclmn++;	arrdim = *len;

	ixa0  = ibold = 0;	UpLeftx = 2;   	hbox   = hsym+3;
	Width = 2*xshif+nend*wsym;   Height = 2*yshif+(*nrows+4+irow0)*hbox;
	if (Kevent.xcur+Kevent.ycur < 10) 
	  { Kevent.xcur = 330;
	    Kevent.ycur = 430;
	  }
       	/*  j = XWX-Width-10;	if (j > UpLeftx)  UpLeftx = j; */
	GetRWgeometry (&XRW,&YRW);
	UpLeftx = 2;	   j = XRW-Width-8;	UpLefty = YRW;
	if (j > UpLeftx)   UpLeftx = j;		j = YRW+Height+30;
	if (j > theHeight) UpLefty = theHeight-Height-30;
	theWindow = Open_Window(UpLeftx, UpLefty, Width, Height, 0, title, 0,
		    RootWindow(theDisplay,theScreen), theMenuCursor);
	MVPOINTER_ xshif+wsym*nwid[0]/2,yshif+irow0*hbox+hsym);
	XSelectInput (theDisplay, theWindow, POLL_EV_MASK);

Create_table:
	Change_Color (theGCA, 3, 0); /* white background, blue  foreground */
	if  ( irow0 ) 
	    {	/*********  Type template in blue  *********/
	    for ( j=0, ix=xshif; j < nclmn ; j++)
		{ WRITE_ theGCA, ix, hsym, template+nsta[j], nwid[j]);
		  if (j+1 == nclmn) break;	ix += wsym*(nwid[j]+1);
		}
		/*********  Draw blue horizontal lines  *********/
	    j = hbox;			LINGCA_ 0,j+4,Width,j+4);
	    }
	if  ( *ngroup > 0 )
	    for ( j = (*ngroup+irow0)*hbox; j+4+yshif < Height-2*hbox ;
		  j += (*ngroup)*hbox )	  LINGCA_ 0,j+4,Width,j+4);
	if  ( *morow > 0 )
	    {	j = Height-(*morow+2)*hbox-3;	LINGCA_ 0,j,Width,j);
		++j;				LINGCA_ 0,j,Width,j);
	    }
		/*********  Draw blue vertical lines  *********/
	for ( j = 0, ix=xshif-0.5*wsym; j < nclmn-1; j++)
	    { ix += wsym*(nwid[j]+1);	ix1 = ix;
	      if (nwid[j] == 0)		ix1 = ix1-2;
	      if (nwid[j+1] == 0)	ix1 = ix1+2;
	      LINGCA_ ix1,0,ix1,Height-2*hbox);
	    }
		/*********  Draw bottom line & comments  *********/
	WRITE_ theGCA, xshif, Height-3-hbox,
		"Button - select,  <Tab>, <Ret>, Arrows - move", 45);
	WRITE_ theGCA, xshif+4*wsym, Height-4, "/<ESC> - done", 13);
	Change_Color (theGCA, 50, 0);  /* white background, red foreground */
	j = Height-2*hbox; 	LINGCA_ 0,j,Width,j);
	j++; 		   	LINGCA_ 0,j,Width,j);
	Change_Color (hghGC, 50, 0);
	WRITE_ hghGC,  xshif, Height-4, " OK ",4);
	XDrawRectangle (theDisplay,theWindow,hghGC,5L,Height-17L,30L,16L);
	Change_Color (theGCA, 1, 0);	Change_Color (hghGC, 1, 0);
					/* white background, black foreground */
		/*********  Draw input data  *********/
	for (   i=0; i < *nrows; i++ )				/* row loop */
 	    {	iy =yshif+hsym+hbox*(i+irow0);
	    for ( j=isym=0, ix=xshif; j < nclmn ; j++)	     /* column loop */
	    	{ if (j > 0) { isym = nsta[j]; ix += wsym*(nwid[j-1]+1); }
	    	  WRITE_ theGCA, ix, iy, array+i*arrdim+isym, nwid[j]);
		}
	    }
	if ( ibold > 0 )				/* after Expose event */
	   {	isym = 0;	if ( icol > 1)  isym = nsta[icol-1];
		ix=xshif+wsym*isym;	iy=yshif+hsym+hbox*(irow+irow0-1);
		WRITE_ hghGC, ix, iy, stri, wbox);
		MoveArrow(theWindow, 0, iya0, ixa, iya);
	   }
	if ( ii == 0 )		/* Select Upper Left box on entry */
	   { ibox = ibold = irow = icol = icold =ii = 1;
		inold = ind = isym = lvalue = 0;	wbox = nwid[0];
		ix = ixold = ixa = ixa0 = xshif+wsym*isym;
		iy = iyold = iya = iya0 = yshif+hsym+hbox*irow0;
		WRITE_ hghGC, ix, iy, array, wbox);
		strncpy(stri,array,wbox);
		MoveArrow(theWindow, 0, iya0, ixa, iya);
	   }
	XFlush(theDisplay);

Table_control:
	XNextEvent (theDisplay, &theEvent);
	if ( theEvent.xany.window == theRootWindow )
           { ProcessRootWindowEvent (&theEvent);
	     goto Table_control;
	   }
	if (theEvent.type == Expose)	goto	Create_table;
	if (theEvent.type == ButtonPress)
	{   ix = theEvent.xbutton.x;	iy = theEvent.xbutton.y;
	    irow = iy-yshif-hbox*irow0-2;
	    if ( irow >= 0) 
	    {	irow = irow/hbox+1;		i=(ix-xshif)/wsym;
		if ( ibold > 0 )	strncpy(array+ind,stri,wbox);
		if ( irow > *nrows )		  /* No selection or Exit */
		   { if ( i > 3 || iy < Height-4-hsym ) goto Unselect;
			goto EndDialog;			/* OK was pressed */
		   }
		for ( j=0; j < nclmn; )			/* Selection is made */
		    { i -= nwid[j]+1;	icol=++j;   if ( i < 0 ) break;	}
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
				|| nwid[icol-1] == 0 )	goto	Unselect;
Fillbox:	isym = spos = 0;	if ( icol > 1)  isym = nsta[icol-1];
		ix=ixa=xshif+wsym*isym;	iy=iya=yshif+hsym+hbox*(irow+irow0-1);
		ind = (irow-1)*arrdim+isym;	wbox = nwid[icol-1];
		ibox = (irow-1)*nclmn+icol;
		if ( ibold != 0 )
		     WRITE_ theGCA, ixold, iyold, stri, nwid[icold-1]);
		strncpy(stri,array+ind,wbox);
		WRITE_ hghGC, ix, iy, stri, wbox);
		MoveArrow(theWindow, ixa0, iya0, ixa, iya);
		icold = icol;	inold = ind;	ixold = ix;	iyold = iy;
		ixa0 = ixa;	iya0 = iya;	ibold = ibox;	lvalue = 0;
	    }
	    else 
Unselect:   {	if ( ibold > 0 )			{
		WRITE_ theGCA, ixold, iyold, array+inold, nwid[icold-1]);
		MoveArrow(theWindow,ixa0,iya0,0,iy);	ibold=spos=0;	}
	    }
	    XFlush(theDisplay);
	    goto Table_control;
	}

	if ( ibold > 0 )
	{
	   iret = GetKey (theEvent.xkey, stri, &spos);
	   if ( iret == -1 && ibold > 0 )			/* <Esc> */
	      { strncpy(array+ind,stri,wbox);   goto EndDialog;
	      }
	   if ( iret == 0 )
	      {	if (spos == 0 ) { 		/* 1-st entry in the box */
		for (j=1; j < 132; j++) stri[j]=' '; stri[wbox] = '\0';  }
		if (++spos > wbox) spos = wbox;		lvalue = 1;
		WRITE_ hghGC, ix, iy, stri, wbox);	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 1 )				/* <Del> or <BS> */
	      {	if (--spos >= 0) { strncpy(stri+spos,stri+spos+1,wbox-spos);
				  stri[wbox-1] = ' ';	}
		if (spos < 0) { spos = 0;  strncpy(stri,array+ind,wbox); }
		WRITE_ hghGC, ix, iy, stri, wbox);	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( (iret == 2 || iret == 3) && ibold > 0 )	/* <Tab> or <Ret> */
	      {	Move_right:
		do { if ( ++icol > nclmn )	icol = 1; }
		while ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
				|| nwid[icol-1] == 0 );
		if ( icol == 1 )  goto	Move_down;
		if ( spos != 0 )  goto	Insert;		goto Fillbox;
	      }
	   if ( iret == 12 )				/* case	XK_Up:	  */
	      { Move_up:
		if ( --irow == 0 ) irow = *nrows;
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
			 && irow >= 0)		goto	Move_up;
		if ( spos != 0 )  goto	Insert;	  goto	Fillbox;
	      }
	   if ( iret == 21 )				/* case	XK_Left:  */
	      { Move_left:
		if ( --spos < 0 )
		{ spos = 0;
		     if (--icol == 0 ) { icol = nclmn;   --irow;}
		     if ( irow == 0 )	 icol = irow = 1;
		     if (If_empty(icol, irow, nclmn, nsta, nwid, array, arrdim)
				|| nwid[icol-1] == 0 )	goto	Move_left;
		     if ( lvalue != 0 )	goto  Insert;	goto	Fillbox;
		}	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 23 )				/* case	XK_Right: */
	      {	if (spos == 0 && lvalue == 0)	strncpy(stri,array+ind,wbox);
		if (++spos > wbox)	goto	Move_right;
		ixa = ix+spos*wsym;		lvalue = 1;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 32 )				/* case	XK_Down:  */
	      { Move_down:
		if ( ++irow > *nrows)	irow = 1;
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
			 && irow <= *nrows)	goto	Move_down;
		if (spos != 0)	goto Insert;	goto Fillbox;
 	      }
	}
	if ( GetEsc(theEvent.xkey) ) goto EndDialog;
	goto Table_control;
Insert:
	strncpy(array+ind,stri,wbox);   goto Fillbox;

EndDialog:
	RETURNPOINTER_;
	XFlush(theDisplay);
	XDestroyWindow(theDisplay, theWindow);
	XFlush(theDisplay); 
	ibcursor = -1;
	mode = 0;
	return 0;
}
/********** The function is called from ASKXGR (file for/surv.f) 
            that in turn is called from IFKEY (key "O")		*********/
int askgrf  (title, template, array, len, nrows, ngroup, morow, modex)
   int	*len, *nrows, *ngroup, *morow, *modex;
   char title[], template[], array[];
{  int	i;
   i = askgrf_ (title, template, array, len, nrows, ngroup, morow, modex);
	return i;
}
/* The same as "asktab", but the 1st column and the column next to "||"
   		are drawn in blue and closed for access
Input:	title	- Title of the table 
	template  string defining a structure of the table and its 1st line
		  1st line does not appear if all non-'|' symbols are spaces
	array	- data (numbers or strings) for input and output
	len	- length of "array" element according to description 
			in the calling routine (80 in the example below)
	nrows	- number of rows in a table to be created
	ngroup	- if > 0 distance (in rows) between blue separating lines
	morow	- if > 0 separates bottom of the table with a fat blue line
		  this part of the table represents switched off windows,
		  i.e. those with negative number.
   ! Note: this parameter is used for output
*/
int askgrf_ (title, template, array, len, nrows, ngroup, morow, modex)
			int *len, *nrows, *ngroup, *morow, *modex;
			char title[], template[], array[];
{	XEvent	theEvent;
	Window	theWindow;
		int	UpLeftx =2, UpLefty =10,/* window corner location */
		Width =452, Height =480,/* defaut window sizes		*/
		ibox, ibold, hbox,	/* current/old box No. & height	*/
		wsym = 8, hsym = 13,	/* symbol width and height	*/
		wbox, lbox,		/* number of chars in box 	*/
		xshif = 5, yshif = 3,	/* table corner 		*/
		mxfields = 20,	nclmn,	/* max / actual # of columns	*/
		nsta[20],nwid[20],	/* ! <= mxfields=20 are allowed	*/
		nlock[10]={1,0},	/* no more than 10 can be locked*/
		ixa0, iya0, ixa, iya,   /* arrow (textcursor) position	*/
		ind,  inold,	/* address of the selected word in "array" */
		isym, spos=0,	/* abs. and rel. position of symbol in table*/
		lvalue,	ixbox=1,/* 1 if box was changed, 0 otherwise	*/
		iret, ixold, iyold, icold, i, j, nend,
		icol, irow, irow0=0, arrdim, ix, ix1, iy, ihelp ;
	char	stri[132], test='0', str[80];
	char	xlab[6][14] = { {"  \"a\"   [m]  "},{" \"a_N\"  [d/l]"},
				{"\"rho_N\" [d/l]"},{" \"psi\"  [Vs] "},
				{"\"rho_V\" [m]  "},{"\"rho_p\" [d/l]"} };
	strncpy(xlab[0],"  \"a\"   [m]  ",13);
	strncpy(xlab[1]," \"a_N\"  [d/l]",13);
	strncpy(xlab[2],"\"rho_N\" [d/l]",13);
	strncpy(xlab[3]," \"psi\"  [Vs] ",13);
	strncpy(xlab[4],"\"rho_V\" [m]  ",13);
	strncpy(xlab[5],"\"rho_p\" [d/l]",13);
	j = strlen(title);	if ( j > 132 )			{
		printf("%s\n\"%s\"\n%s\n","ASKCOL >>> Title",title,
		     "           is too long");		return (-1);	}
	if  ( *len > 132 )				{
		printf("%s \"%s\" %s\n","ASKCOL >>> Table",title,
		     ". Input string is too long");	return (-1);	}
	nend = strlen(template);	if ( nend > 132 )	{
		printf("%s \"%s\" %s\n","ASKCOL >>> Table",title,
		     ". Requested width is too large");	return (-1);	}
	for ( j=0; j < nend ; j++)
	      if ( template[j] != '|'  &&  template[j] != ' ' ) irow0 = 1;

		/****   Split "template" in blocks  ****/
	for ( j=nclmn=nsta[0]=0; j < nend ; j++ )
	    { if ( template[j] == '|' ) 
		 { nwid[nclmn] = j-nsta[nclmn];
		   nsta[++nclmn] = j+1;	/*	nlock[nclmn] = 0;
		   if ( template[j+1] == '|' )  nlock[nclmn] = 1;*/
		 }
	      if ( nclmn > mxfields ) {
		printf("%s %d\n","ASKCOL >>>  Too many input fields ",nclmn);
			return (-1);	}
	    }
	nwid[nclmn] = nend-nsta[nclmn];		nclmn++;	arrdim = *len;

	ixa0  = ibold = ihelp = 0;	hbox = hsym+3;	UpLeftx = 2;
	Width = 2*xshif+nend*wsym;	lbox = 3*hbox;
	Height = 2*yshif+(*nrows+irow0)*hbox+lbox;
	if (Kevent.xcur+Kevent.ycur < 10) 
		{Kevent.xcur = 330; Kevent.ycur = 430;}
	/*  j = XWX-Width-20;	if (j > UpLeftx)  UpLeftx = j; */
	GetRWgeometry (&XRW,&YRW);
	UpLeftx = 2;	   j = XRW-Width-8;	UpLefty = YRW;
	if (j > UpLeftx)   UpLeftx = j;		j = YRW+Height+30;
	if (j > theHeight) UpLefty = theHeight-Height-30;
	theWindow = Open_Window(UpLeftx, UpLefty, Width, Height, 0, title, 0,
		    RootWindow(theDisplay,theScreen), theMenuCursor);
	MVPOINTER_ xshif+38*wsym,Height+0.5*hbox-lbox);
	XSelectInput (theDisplay, theWindow, POLL_EV_MASK);
	if ( *modex == 0 )	test = '0';
	if ( *modex == 1 )	test = '1';
	if ( *modex == 2 )	test = '2';
	if ( *modex == 3 )	test = '3';
	if ( *modex == 4 )	test = '4';
	if ( *modex == 5 )	test = '5';

Create_table:
	Change_Color (theGCA, 3, 0); /* white background, blue  foreground */
	if  ( irow0 ) 
	    {	/*********  Type template in blue  *********/
	    for ( j=0, ix=xshif; j < nclmn ; j++)
		{ WRITE_ theGCA, ix, hsym, template+nsta[j], nwid[j]);
		  if (j+1 == nclmn) break;	ix += wsym*(nwid[j]+1);
		}
		/*********  Draw blue horizontal lines  *********/
	    j = hbox;			LINGCA_ 0,j+4,Width,j+4);
	    }
	if  ( *ngroup > 0 )
	    for ( j = (*ngroup+irow0)*hbox; j+4+yshif < Height-lbox ;
		  j += (*ngroup)*hbox )	  LINGCA_ 0,j+4,Width,j+4);
	if  ( *morow > 0 )
	    {	j = Height-*morow*hbox-lbox-3;	LINGCA_ 0,j,Width,j);
		++j;				LINGCA_ 0,j,Width,j);
	    }
		/*********  Draw blue vertical lines  *********/
	for ( j = 0, ix=xshif-0.5*wsym; j < nclmn-1; j++)
	    { ix += wsym*(nwid[j]+1);	ix1 = ix;
	      if (nwid[j] == 0)		ix1 = ix1-2;
	      if (nwid[j+1] == 0)	ix1 = ix1+2;
	      LINGCA_ ix1,0,ix1,Height-lbox);
	    }
	    j = Height-lbox;	LINGCA_ 0,j,Width,j);
		++j;		LINGCA_ 0,j,Width,j);
	    if ( *modex >= 0 )	LINGCA_ ix1,j,ix1,j+hbox);
		/*********  Draw radial cooordinate box *********/
						/* Change of X-axis */
/*	if ( *modex >= 0 )
	   { strcpy(stri,"Your flux coordinate is ");
	     WRITE_ theGCA, xshif, Height-3+hbox-lbox, stri, 24);
	   }
	else
*/	   { strcpy(stri,"The current abscissa is ");
	     WRITE_ theGCA, xshif, Height-3+hbox-lbox, stri, 24);
	   }

		/*********  Draw bottom line & comments  *********/
	Change_Color (theGCA, 1, 0);	Change_Color (hghGC, 1, 0);
	WRITE_ theGCA, xshif, Height-3-hbox,
		"Button - select, <Tab>,<Ret>,Arrows - move", 42);
	WRITE_ theGCA, xshif+4*wsym, Height-4, "/<ESC> - done,", 14);
	WRITE_ hghGC,  xshif+28*wsym, Height-4, "?",1);
	WRITE_ theGCA, xshif+29*wsym, Height-4, " - help", 7);
	Change_Color (theGCA, 50, 0); 		/* red on white */
	j = Height-2*hbox; 	LINGCA_ 0,j,Width,j);
	j++; 		   	LINGCA_ 0,j,Width,j);
	Change_Color (hghGC, 50, 0);
	WRITE_ hghGC,  xshif, Height-4, " OK ",4);
	Change_Color (theGCA, 1, 0);	Change_Color (hghGC, 1, 0);
						/* black on white */
/*printf("&Root %d  &Win %d;",&theRootWindow,&theWindow);
  printf("   Root %d  Win %d\n\n",theRootWindow,theWindow);
  printf("Start: &theEvent %d  theEvent %d\n\n",&theEvent,theEvent);
*/
		/*********  Draw input data  *********/
	for (   i=0; i < *nrows; i++ )				/* row loop */
 	    {	iy =yshif+hsym+hbox*(i+irow0);	Change_Color(theGCA,3,0);
	    for ( j=isym=0, ix=xshif; j < nclmn ; j++)	     /* column loop */
	    	{ if (j > 0) { isym = nsta[j]; ix += wsym*(nwid[j-1]+1);
			if ( nwid[j-1] == 0)	Change_Color(theGCA,3,0);
			       else		Change_Color(theGCA,1,0);  }
	    	  WRITE_ theGCA, ix, iy, array+i*arrdim+isym, nwid[j]);
		}
	    }
	if ( ibold > 0 && ixbox == 0 )		/* after Expose event */
	   {	isym = 0;	if ( icol > 1)  isym = nsta[icol-1];
		ix=xshif+wsym*isym;	iy=yshif+hsym+hbox*(irow+irow0-1);
		strncpy(stri,array+ind,wbox);
		WRITE_ hghGC, ix, iy, stri, wbox);
		MoveArrow(theWindow, 0, iya0, ixa, iya);
	   }
	Change_Color (hghGC, 3, 0);	ix1 = xshif+24*wsym;
	if ( *modex >= 0 )
	   { WRITE_ hghGC, ix1, Height-3+hbox-lbox, xlab[*modex], 13);
	   }
	else
	   { WRITE_ hghGC, ix1, Height-3+hbox-lbox, "\"time\"", 6);
	   }
	Change_Color(hghGC,1,0);
	if ( ixbox > 0 )		/* Select Xbox on entry */
	   { irow = 1;	ibox = ibold = icol = icold = 0;    lvalue = 0;
	     inold = ind = isym = 0;		wbox = 1;
	     ix = ixold = ixa = ixa0 = xshif+41*wsym-3;
	     iy = iyold = iya = iya0 = Height-3+hbox-lbox;
	     goto	FillXbox;
	   }
	 else if ( *modex >= 0 )
	   { ix1 = xshif+41*wsym;	Change_Color(theGCA,1,0);
	     WRITE_ theGCA, ix1, Height-3+hbox-lbox, &test, 1);
	   }

Table_control:
	XNextEvent (theDisplay, &theEvent);
	if ( theEvent.xany.window == theRootWindow )
{ /* printf("theRootWindowEvent.type =%d\n",theEvent.type); */
	     ProcessRootWindowEvent (&theEvent);
	     goto Table_control;
	   }
	if (theEvent.type == Expose)		goto	Create_table;
	if (theEvent.type == ButtonPress)
	{						/* ButtonPressed */
	   ix = theEvent.xbutton.x;	iy = theEvent.xbutton.y;
	   irow = iy-yshif-hbox*irow0-2;	
	   irow = irow/hbox+1;		i=(ix-xshif)/wsym;

	   if ( irow < 0)	goto	Unselect;
	   if ( irow > *nrows+1 && i < 3 && iy > Height-4-hsym )
	      {	if ( spos != 0 )	goto	Escape;	/* OK was pressed */
					goto	EndDialog;
	      }
	   if ( irow == *nrows+1 && i >= 38 && i <= 44 && *modex >= 0 )
	      {	if ( ixbox > 0 )  goto	Table_control;	/* Staying in Xbox */
		ixbox = 1; if ( ibox == 0 ) goto FillXbox;
		goto	Unselect;			/* Entering Xbox */
	      }
NewSelection:						/* Obox selected */
	   for ( j=0; j < nclmn; )
	      { i -= nwid[j]+1;	icol=++j;   if ( i < 0 ) break;	}
	   if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
		   || icol == 1	|| nwid[icol-2] == 0 || nwid[icol-1] == 0 )
	      { ibox = 0;	goto	Unselect;   /* Empty Obox selected */
	      }
	   else
	      { i = (irow-1)*nclmn+icol;	/* Non-empty Obox selected */
		if ( i == ibox)	goto	Table_control;
		if ( i != ibox) { ibold = ibox;	ibox = i;  goto	Unselect; }
	      }
Unselect:
	   if ( ibold > 0 )
		{ strncpy(array+ind,stri,wbox);		ibold = spos = 0;
		  WRITE_ theGCA, ixold, iyold, array+inold, nwid[icold-1]);
		  MoveArrow(theWindow,ixa0,iya0,0,iy);	}
	   else if ( ixbox > 0 )
		{ ixbox = 0;				/* Leaving Xbox */
		  if ( *modex < 0 )	{ ixa0 = 0;    goto	Highlight; }
		  ix1 = xshif+41*wsym;	j = Height-3+hbox-lbox;
		  MoveArrow(theWindow, ixa0, iya0, 0, j);	ixa0 = 0;
		  Change_Color(theGCA,1,0);	WRITE_ theGCA,ix1,j,&test,1);
		  Change_Color (theGCA, 50, 0);	j = Height-2*hbox;
		  LINGCA_ 0,j,Width,j);	j++; 	LINGCA_ 0,j,Width,j);
		  Change_Color(theGCA,1,0);
		}
Highlight:
	   if ( ixbox > 0 )	goto	FillXbox;
	   if ( ibox == 0 )	goto	Table_control;

FillObox:  isym = spos = 0;	if ( icol > 1)  isym = nsta[icol-1];
	   ix=ixa=xshif+wsym*isym;	iy=iya=yshif+hsym+hbox*(irow+irow0-1);
	   ind = (irow-1)*arrdim+isym;	wbox = nwid[icol-1];
	   if ( ibold != 0 ) 
	        WRITE_ theGCA, ixold, iyold, stri, nwid[icold-1]);
	   if ( ixbox == 0 )
	      { strncpy(stri,array+ind,wbox);
		WRITE_ hghGC, ix, iy, stri, wbox);
		MoveArrow(theWindow, ixa0, iya0, ixa, iya);
	      }
	   icold = icol;	inold = ind;	ixold = ix;	iyold = iy;
	   ixa0 = ixa;	iya0 = iya;	ibold = ibox;	lvalue = 0;
	   XFlush(theDisplay);
	   goto Table_control;

FillXbox:  if ( *modex < 0 )  goto	Unselect;
	   ix1 = xshif+41*wsym;		iy = Height-3+hbox-lbox;
	   MoveArrow(theWindow, ixa0, iya0, ix1-2, iy);
	   ixa0 = ix1-2;		iya0 = iy;	ibox = 0;
	   sscanf(&test,"%d",modex);	/* Equivalent to atoi(&test) */
	   Change_Color (hghGC, 3, 0);	ix1 = xshif+24*wsym;
	   WRITE_ hghGC, ix1, Height-3+hbox-lbox, xlab[*modex], 13);
 	   ix1 = xshif+41*wsym;	Change_Color(hghGC,1,0);
	   WRITE_ hghGC, ix1, iy, &test, 1);
	   XFlush(theDisplay);
	   goto Table_control;
	}
/*printf("theEvent.type = %d,  KeyPress = %d\n",theEvent.type,KeyPress); */
   if (theEvent.type != KeyPress)	goto	Table_control;
	iret = GetKey (theEvent.xkey, stri, &spos);
	if ( ixbox > 0 )
	   {  					/*  Key pressed in Xbox */
	   if ( iret == 3  || iret == 12 || iret == 23 )
	      { test=test+1;	   if ( test > '5' )	test = '0';
		goto	FillXbox;
	      }
	   if ( iret == 21 || iret == 32 )
	      { test=test-1;	   if ( test < '0' )	test = '5';
		goto	FillXbox;
	      }
	   if ( iret == 0 )
	      {
	      if ( stri[0] == '?' )	goto	Printinfo;
	      if ( stri[0] < '0' )	{
		   stri[0] = '0';	goto	Printinfo;
					}
	      if ( stri[0] > '5')	{
		   stri[0] = '5';	goto	Printinfo;
					}
//	      if ( stri[0] < '0' || stri[0] > '5')	goto	FillXbox;
//              strncpy(&test,&stri[0],1);
//              test = stri[0];
              strcpy(&test,&stri[0]);
	      goto	FillXbox;
Printinfo:    printf(" The following flux labels are allowed:\n");
	      for (j=0; j<6; j++) printf("       %d   %.13s\n",j,xlab[j]);
	      }
	   }
	if ( ibold > 0 )
	   {  				/*  Key pressed in active Obox */
	   if ( iret == -1 )				/* <Esc> */
Escape:	      { strncpy(array+ind,stri,wbox);   goto EndDialog;
	      }
	   if ( iret == 7 )			   /* <Alt><Esc> */
	      { strncpy(array+ind,stri,wbox);
		*morow=-1;			goto EndDialog;
	      }
	   if ( iret == 0 )
	      {	if (spos == 0 ) { 		/* 1-st entry in the box */
		for (j=1; j < 132; j++) stri[j]=' '; stri[wbox] = '\0';  }
	        if ( stri[spos] == '?' )
		 { printf("%d  %.13s\n",spos,stri);
		printf("\"%.8s\"\n",array+ind);

		 /*		sscanf(&str,"%[^ ]%s",array+ind);*/
		printf("\"%s\"\n",str);
		 /*		sscanf(&str,"%s","grep -w \"");*/
		printf("\"%s\"\n",array+ind);
		for(j=0; j<=wbox; j++)
		if (array[ind+j] != ' ') strncat(str,array+ind+j,1);
		printf("\"%s\"\n",str);

		   strcpy (str,"grep -w \"");
		for(j=0; j<=wbox; j++)
		if ( array[ind+j] != ' ' ) strncat(str,array+ind+j,1);
		   strcat (str,"\" tmp/model.txt");
		   printf("%.50s \n",str);
		   if (!ihelp) printf(" No. Scale  Name  Output expression\n");
		   system(str);		ihelp = 1;
		   goto Table_control;
              }
		if (++spos > wbox) spos = wbox;		lvalue = 1;
		WRITE_ hghGC, ix, iy, stri, wbox);	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	   }
	   if ( iret == 1 )				/* <Del> or <BS> */
	      {	if (--spos >= 0) { strncpy(stri+spos,stri+spos+1,wbox-spos);
				  stri[wbox-1] = ' ';	}
		if (spos < 0) { spos = 0;  strncpy(stri,array+ind,wbox); }
		WRITE_ hghGC, ix, iy, stri, wbox);	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( (iret == 2 || iret == 3) && ibold > 0 )	/* <Tab> or <Ret> */
	      {	Move_right:
		do { if ( ++icol > nclmn )	icol = 2; }
		while ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
		 	|| nwid[icol-2] == 0 || nwid[icol-1] == 0 );
		if ( icol == 2 )  goto	Move_down;
		if ( spos != 0 )  goto	Insert;   goto FillObox;
	      }
	   if ( iret == 12 )				/* case	XK_Up:	  */
	      { Move_up:
		if ( --irow == 0 ) irow = *nrows;
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
			 && irow >= 0)	goto	Move_up;
		if ( spos != 0 )  goto	Insert;   goto FillObox;
	      }
	   if ( iret == 21 )				/* case	XK_Left:  */
	      { Move_left:
		if ( --spos < 0 )
		{ spos = 0;
		     if (--icol == 1 )	{ icol = nclmn;   --irow; }
		     if ( irow == 0 )	{ irow = 1;	icol = 2; }
		     if ( If_empty(icol, irow, nclmn, nsta, nwid, array, arrdim)
			|| nwid[icol-2]==0 || nwid[icol-1]==0) goto Move_left;
		     if ( lvalue != 0 )	goto  Insert;   goto FillObox;
		}	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 23 )				/* case	XK_Right: */
	      {	if (spos == 0 && lvalue == 0)	strncpy(stri,array+ind,wbox);
		if (++spos > wbox)	goto	Move_right;
		ixa = ix+spos*wsym;		lvalue = 1;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 32 )				/* case	XK_Down:  */
	      { Move_down:
		if ( ++irow > *nrows)	irow = 1;
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
			 && irow <= *nrows)	goto	Move_down;
		if (spos != 0)	goto Insert;	goto	FillObox;
 	      }
	   }
	if ( GetEsc(theEvent.xkey) )		goto	EndDialog;
	goto Table_control;
Insert:
	strncpy(array+ind,stri,wbox);
	goto	FillObox;
EndDialog:
	RETURNPOINTER_;
	XFlush(theDisplay);
	XDestroyWindow(theDisplay, theWindow);
	XFlush(theDisplay); 
	ibcursor = -1;
	return (0);
}
/***********************************************************************/
/* The function processes selected events in theRootWindow,
   e.g. redrawing in case of exposure. 
	changing mode
	   0 otherwise */
void ProcessRootWindowEvent (theEvent)   XEvent	*theEvent;
{  int	j, k;
   XKeyEvent	theKeyEvent;
   XComposeStatus	theComposeStatus;
   KeySym		theKeySym;
   char 		theKeyBuffer[5];
   switch( theEvent->type )
     {	case ButtonPress:
	   Kevent.xcur = theEvent->xbutton.x;
	   Kevent.ycur = theEvent->xbutton.y;
	   k = Cursor_in_Box ();		  /* Returns box number */
	   if (k >= 0 && k <= nmamb)	j = (int)MAMK[k];
			   /* printf("Button pressed,  No = %d\n",j); */
	   goto 	Call_ifkey;
	case KeyPress:
	   theKeyEvent = theEvent->xkey;
	   j = XLookupString (&theKeyEvent, theKeyBuffer, 4,
			      &theKeySym, &theComposeStatus);
	   if ( j == 0 )		/* <Ret>, <Alt>, etc. is pressed */
	   return; 
	   if ( j != 1 ) 
	      { printf("String translation error:  Bufferlength =%d\n", j);
	        return; }
	   j = theKeySym;
	      /* printf("Key pressed:  Key = \"%c\",  Key = %d\n",
					theKeySym, theEvent->xkey); */
	   if (  j > 96 && j < 123 )	j -= 32;
Call_ifkey:
	   if ( j == 47 ||				/* "/" Exit -> */
	        j == 13 || 				/* <Enter> */
	        j == 32 || 				/* <Space> */
	        j == 37 || 				/* "%" */
	       (j >= 48 && j <= 57  ) ||		/* "0" -> "9" */
	        j == 66 || 				/* "B" */
	        j == 67 || 				/* "C" */
	        j == 68 || 				/* "D" */
	        j == 72 ||  				/* "H" */
	        j == 78 ||  				/* "N" */
	        j == 79 ||  				/* "O" */
	        j == 82 ||  				/* "R" */
	        j == 83 ||  				/* "S" */
	        j == 87 ||  				/* "W" */
	        j == 88  				/* "X" */
	      )		ifkey_(&j);
	   j = 0; 
	   return;
	case Expose:		/* printf("Expose %d\n",theEvent->type); */
	   j = 82;	ifkey_(&j);	j = 0;
	   return;
/*	case MapNotify:		 printf("Mapping %d\n",theEvent->type);
	   return;
	case UnmapNotify:	 printf("UnmapNotify %d\n",theEvent->type);
	   return;
*/	default:	      /* printf("Unrecognized %d\n",theEvent->type); */
	   return;
      }
}
/***********************************************************************/
askcol  (title, template, array, len, nrows, ngroup, morow)
				int *len, *nrows, *ngroup, *morow;
				char title[], template[], array[];
/* The same as "asktab", but the 1st column 
   		is drawn in blue and closed for access
Input:	title	- Title of the table 
	template  string defining a structure of the table and its 1st line
		  1st line does not appear if all non-'|' symbols are spaces
	array	- data (numbers or strings) for input and output
	len	- length of "array" element according to description 
			in the calling routine (80 in the example below)
	nrows	- number of rows in a table to be created
	ngroup	- if > 0 distance (in rows) between blue separating lines
	morow	- if > 0 separates bottom of the table with a fat blue line

*/
{	askcol_(title, template, array, len, nrows, ngroup, morow);	}
askcol_ (title, template, array, len, nrows, ngroup, morow)
				int *len, *nrows, *ngroup, *morow;
				char title[], template[], array[];
{	Window	theWindow;
	XEvent	theEvent;
	int	UpLeftx =2, UpLefty =10,/* window corner location	*/
		Width =452, Height =480,/* defaut window sizes		*/
		ibox, ibold, hbox,	/* current/old box No. & height	*/
		wsym = 8, hsym = 13,	/* symbol width and height	*/
		wbox, 			/* number of chars in box 	*/
		xshif = 5, yshif = 3,	/* table corner 		*/
		mxfields = 20,	nclmn,	/* max / actual # of columns	*/
		nsta[20],nwid[20],nend,	/* ! <= mxfields=20 are allowed	*/
		ixa0, iya0, ixa, iya,   /* arrow (textcursor) position	*/
		ind,  inold,	/* address of the selected word in "array" */
		isym, spos=0,	/* abs. and rel. position of symbol in table*/
		lvalue,		/* 1 if box was changed, 0 otherwise	*/
		iret, ixold, iyold, icold, i, j,
		icol, irow, irow0=0, arrdim, ix, ix1, iy, ii=0;
	char	stri[132];
	j = strlen(title);	if ( j > 132 )			{
		printf("%s\n\"%s\"\n%s\n","ASKCOL >>> Title",title,
		     "           is too long");		return;	}
	if  ( *len > 132 )				{
		printf("%s \"%s\" %s\n","ASKCOL >>> Table",title,
		     ". Input string is too long");	return;	}
	nend = strlen(template);	if ( nend > 132 )	{
		printf("%s \"%s\" %s\n","ASKCOL >>> Table",title,
		     ". Requested width is too large");	return;	}
	for ( j=0; j < nend ; j++)
	      if ( template[j] != '|'  &&  template[j] != ' ' ) irow0 = 1;

		/****   Split "template" in blocks  ****/
	for ( j=nclmn=nsta[0]=0; j < nend ; j++ )
	    { if ( template[j] == '|' ) 
		 { nwid[nclmn] = j-nsta[nclmn]; nsta[++nclmn] = j+1;}
	      if ( nclmn > mxfields ) {
		printf("%s %d\n","ASKCOL >>>  Too many input fields ",nclmn);
			return;	}
	    }
	nwid[nclmn] = nend-nsta[nclmn];		nclmn++;	arrdim = *len;

	ixa0  = ibold = 0;		hbox   = hsym+3;	UpLeftx = 2;
	Width = 2*xshif+nend*wsym;	Height = 2*yshif+(*nrows+2+irow0)*hbox;
	if (Kevent.xcur+Kevent.ycur < 10) 
			{Kevent.xcur = 330; Kevent.ycur = 430;}
	/*  j = XWX-Width-20;	if (j > UpLeftx)  UpLeftx = j; */
	GetRWgeometry (&XRW,&YRW);
	UpLeftx = 2;	   j = XRW-Width-8;	UpLefty = YRW;
	if (j > UpLeftx)   UpLeftx = j;		j = YRW+Height+30;
	if (j > theHeight) UpLefty = theHeight-Height-30;
	theWindow = Open_Window(UpLeftx, UpLefty, Width, Height, 0, title, 0,
		    RootWindow(theDisplay,theScreen), theMenuCursor);
	MVPOINTER_ xshif+wsym*nwid[0]/2,yshif+irow0*hbox+hsym);
	XSelectInput (theDisplay, theWindow, POLL_EV_MASK);

Create_table:
	Change_Color (theGCA, 3, 0); /* white background, blue  foreground */
	if  ( irow0 ) 
	    {	/*********  Type template in blue  *********/
	    for ( j=0, ix=xshif; j < nclmn ; j++)
		{ WRITE_ theGCA, ix, hsym, template+nsta[j], nwid[j]);
		  if (j+1 == nclmn) break;	ix += wsym*(nwid[j]+1);
		}
		/*********  Draw blue horizontal lines  *********/
	    j = hbox;			LINGCA_ 0,j+4,Width,j+4);
	    }
	if  ( *ngroup > 0 )
	    for ( j = (*ngroup+irow0)*hbox; j+4+yshif < Height-2*hbox ;
		  j += (*ngroup)*hbox )	  LINGCA_ 0,j+4,Width,j+4);
	if  ( *morow > 0 )
	    {	j = Height-(*morow+2)*hbox-3;	LINGCA_ 0,j,Width,j);
		++j;				LINGCA_ 0,j,Width,j);
	    }
		/*********  Draw blue vertical lines  *********/
	for ( j = 0, ix=xshif-0.5*wsym; j < nclmn-1; j++)
	    { ix += wsym*(nwid[j]+1);	ix1 = ix;
	      if (nwid[j] == 0)		ix1 = ix1-2;
	      if (nwid[j+1] == 0)	ix1 = ix1+2;
	      LINGCA_ ix1,0,ix1,Height-2*hbox);
	    }
		/*********  Draw bottom line & comments  *********/
	WRITE_ theGCA, xshif, Height-3-hbox,
		"Button - select, <Tab>,<Ret>,Arrows - move", 42);
	WRITE_ theGCA, xshif+4*wsym, Height-4, "/<ESC> - done", 13);
	Change_Color (theGCA, 50, 0);  /* white background, red foreground */
	j = Height-2*hbox; 	LINGCA_ 0,j,Width,j);
	j++; 		   	LINGCA_ 0,j,Width,j);
	Change_Color (hghGC, 50, 0);
	WRITE_ hghGC,  xshif, Height-4, " OK ",4);
	Change_Color (theGCA, 1, 0);	Change_Color (hghGC, 1, 0);
				/* white background, black foreground */
		/*********  Draw input data  *********/
	for (   i=0; i < *nrows; i++ )				/* row loop */
 	    {	iy =yshif+hsym+hbox*(i+irow0);	Change_Color (theGCA, 3, 0);
	    for ( j=isym=0, ix=xshif; j < nclmn ; j++)	     /* column loop */
	    	{ if (j > 0) { isym = nsta[j]; ix += wsym*(nwid[j-1]+1); }
	    	  WRITE_ theGCA, ix, iy, array+i*arrdim+isym, nwid[j]);
		  Change_Color (theGCA, 1, 0);
		}
	    }
	if ( ibold > 0 )			/* after Expose event */
	   {	isym = 0;	if ( icol > 1)  isym = nsta[icol-1];
		ix=xshif+wsym*isym;	iy=yshif+hsym+hbox*(irow+irow0-1);
		WRITE_ hghGC, ix, iy, stri, wbox);
		MoveArrow(theWindow, 0, iya0, ixa, iya);
	   }
	if ( ii == 0 )		/* Select Upper Left box on entry */
	   { irow =ii = 1;	ibox = ibold = icol = icold = 2;    lvalue = 0;
		inold = ind = isym = nsta[1];		wbox = nwid[1];
		ix = ixold = ixa = ixa0 = xshif+wsym*isym;
		iy = iyold = iya = iya0 = yshif+hsym+hbox*irow0;
		WRITE_ hghGC, ix, iy, array+ind, wbox);
		strncpy(stri,array+ind,wbox);
		MoveArrow(theWindow, 0, iya0, ixa, iya);
	   }
	XFlush(theDisplay);

Table_control:
	XNextEvent (theDisplay, &theEvent);
	if ( theEvent.xany.window == theRootWindow )
           { ProcessRootWindowEvent (&theEvent);
	     goto Table_control;
	   }
	if (theEvent.type == Expose)	goto	Create_table;
	if (theEvent.type == ButtonPress)
	{   ix = theEvent.xbutton.x;	iy = theEvent.xbutton.y;
	    irow = iy-yshif-hbox*irow0-2;
	    if ( irow >= 0) 
	    {	irow = irow/hbox+1;		i=(ix-xshif)/wsym;
		if ( ibold > 0 )	strncpy(array+ind,stri,wbox);
		if ( irow > *nrows )		  /* No selection or Exit */
		   { if ( i > 3 || iy < Height-4-hsym ) goto Unselect;
			goto EndDialog;			/* OK was pressed */
		   }
		for ( j=0; j < nclmn; )			/* Selection is made */
		    { i -= nwid[j]+1;	icol=++j;   if ( i < 0 ) break;	}
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
		   || icol == 1	|| nwid[icol-1] == 0 )	goto	Unselect;
Fillbox:	isym = spos = 0;	if ( icol > 1)  isym = nsta[icol-1];
		ix=ixa=xshif+wsym*isym;	iy=iya=yshif+hsym+hbox*(irow+irow0-1);
		ind = (irow-1)*arrdim+isym;	wbox = nwid[icol-1];
		ibox = (irow-1)*nclmn+icol;
		if ( ibold != 0 )
		     WRITE_ theGCA, ixold, iyold, stri, nwid[icold-1]);
		strncpy(stri,array+ind,wbox);
		WRITE_ hghGC, ix, iy, stri, wbox);
		MoveArrow(theWindow, ixa0, iya0, ixa, iya);
		icold = icol;	inold = ind;	ixold = ix;	iyold = iy;
		ixa0 = ixa;	iya0 = iya;	ibold = ibox;	lvalue = 0;
	    }
	    else 
Unselect:   {	if ( ibold > 0 )			{
		WRITE_ theGCA, ixold, iyold, array+inold, nwid[icold-1]);
		MoveArrow(theWindow,ixa0,iya0,0,iy);	ibold=spos=0;	}
	    }
	    XFlush(theDisplay);
	    goto Table_control;
	}
	if ( ibold > 0 )
	{
	   iret = GetKey (theEvent.xkey, stri, &spos);
	   if ( iret == -1 && ibold > 0 )			/* <Esc> */
	      { strncpy(array+ind,stri,wbox);   goto EndDialog;
	      }
	   if ( iret == 0 )
	      {	if (spos == 0 ) { 		/* 1-st entry in the box */
		for (j=1; j < 132; j++) stri[j]=' '; stri[wbox] = '\0';  }
		if (++spos > wbox) spos = wbox;		lvalue = 1;
		WRITE_ hghGC, ix, iy, stri, wbox);	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 1 )				/* <Del> or <BS> */
	      {	if (--spos >= 0) { strncpy(stri+spos,stri+spos+1,wbox-spos);
				  stri[wbox-1] = ' ';	}
		if (spos < 0) { spos = 0;  strncpy(stri,array+ind,wbox); }
		WRITE_ hghGC, ix, iy, stri, wbox);	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( (iret == 2 || iret == 3) && ibold > 0 )	/* <Tab> or <Ret> */
	      {	Move_right:
		do { if ( ++icol > nclmn )	icol = 2; }
		while ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
				|| nwid[icol-1] == 0 );
		if ( icol == 2 )  goto	Move_down;
		if ( spos != 0 )  goto	Insert;		goto Fillbox;
	      }
	   if ( iret == 12 )				/* case	XK_Up:	  */
	      { Move_up:
		if ( --irow == 0 ) irow = *nrows;
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
			 && irow >= 0)		goto	Move_up;
		if ( spos != 0 )  goto	Insert;	  goto	Fillbox;
	      }
	   if ( iret == 21 )				/* case	XK_Left:  */
	      { Move_left:
		if ( --spos < 0 )
		{ spos = 0;
		     if (--icol == 1 )	{ icol = nclmn;   --irow; }
		     if ( irow == 0 )	{ irow = 1;	icol = 2; }
		     if (If_empty(icol, irow, nclmn, nsta, nwid, array, arrdim)
				|| nwid[icol-1] == 0 )	goto	Move_left;
		     if ( lvalue != 0 )	goto  Insert;	goto	Fillbox;
		}	ixa = ix+spos*wsym;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 23 )				/* case	XK_Right: */
	      {	if (spos == 0 && lvalue == 0)	strncpy(stri,array+ind,wbox);
		if (++spos > wbox)	goto	Move_right;
		ixa = ix+spos*wsym;		lvalue = 1;
		MoveArrow(theWindow,ixa0,iya0,ixa,iy);	ixa0 = ixa; iya0 = iy;
	      }
	   if ( iret == 32 )				/* case	XK_Down:  */
	      { Move_down:
		if ( ++irow > *nrows)	irow = 1;
		if ( If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
			 && irow <= *nrows)	goto	Move_down;
		if (spos != 0)	goto Insert;	goto Fillbox;
 	      }
	}
	if ( GetEsc(theEvent.xkey) ) goto EndDialog;
	goto Table_control;
Insert:
	strncpy(array+ind,stri,wbox);   goto Fillbox;

EndDialog:
	RETURNPOINTER_;
	XFlush(theDisplay);
	XDestroyWindow(theDisplay, theWindow);
	XFlush(theDisplay); 
	ibcursor = -1;
	return (ibold);
}
/***************  True if cursor is in empty box ***********************/
If_empty (icol, irow, nclmn, nsta, nwid, array, arrdim)
/*
    icol - current column
    irow - current row
    nclmn - total # of columns
    
*/
int	icol, irow, nclmn, arrdim, nsta[], nwid[];
char	array[];
{	int	j, imin, imax;
	imin = 1;	imax = nclmn;
R1:	for ( j=imin-1; j < nclmn ; j++)	{
	      if ( nwid[j] == 0 )
		 { if ( j >= icol )   {	imax = j;	goto	R2; }
					imin = j+2;	goto	R1;
		 }
	   					}
R2:	if ( strspn(array+(irow-1)*arrdim+nsta[imin-1]," ")
		> nsta[imax-1]+nwid[imax-1]-nsta[imin-1])	return (1);
	return (0);
}
/***********************************************************************/
GetValue (theEvent, str, lvalue, tstri, pos) XKeyEvent theEvent; 
			char str[], tstri[]; int *pos, lvalue;
{	int		i;
	XComposeStatus	theComposeStatus;
	KeySym		theKeySym;
	int		theBufferLength, theKeyBufferMaxLen = 4;
	char		theKeyBuffer[5];
	theBufferLength=XLookupString(&theEvent, theKeyBuffer,
		theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
	theKeyBuffer[theBufferLength]='\0';
	switch(theKeySym)
	{   case	XK_BackSpace:
	    case	XK_Delete:	if(*pos>0) {
	    	for (i=*pos; i<lvalue; i++) str[i-1]=str[i];
			str[lvalue-1] = ' '; } *pos = *pos - 1;
							return (0);
	    case	XK_Return:			return (1);
	    case	XK_KP_Enter:			return (1);
	    case	XK_Tab:				return (2);
	    case	XK_Escape:			return (-1);
	    case	XK_question:			return (-2);
	    case	XK_R7:				return (11);
	    case	XK_Left:			return (21);
	    case	XK_Up:				return (12);
	    case	XK_Right:			return (23);
	    case	XK_Down:			return (32);
	    case	XK_R9:				return (13);
	    case	XK_R13:				return (31);
	    case	XK_R15:				return (33);
	    default:	
		if(*pos<lvalue)
			{ for (i=0; tstri[i] != 0; i++) {  
			  if(theKeySym == tstri[i]) { if( *pos<0 ) *pos=0; 
			  str[*pos] = theKeyBuffer[0]; *pos=*pos+1; 
			}			    }	}
							return (0);
	}
}
/*******************************************************************/
GetKey (theKeyEvent, str, pos) XKeyEvent theKeyEvent; char str[]; 
					int *pos;
{	XComposeStatus	theComposeStatus;
	KeySym		theKeySym;
	int		theKeyBufferMaxLen = 4;
	char		theKeyBuffer[5];
    XLookupString(&theKeyEvent, theKeyBuffer,
		theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
	  if ( theKeyEvent.state & Mod1Mask)
	     { if ( theKeySym == XK_Escape )	{
		 /* printf("<Alt>+<Esc> is pressed\n"); */
	     return(7); }
	     }
/*	  if ( theKeyEvent.state & Mod1Mask)
	     { if ( theKeySym <= 127 )	{ ch = theKeySym;
	        printf("<Alt>+%s is pressed,   ASCII code =%d\n",
		theKeyBuffer,theKeySym);}
	     }
	  if ( theKeyEvent.state & ControlMask)
	     { if ( isascii(theKeySym) )	{ 
	        printf("<Ctrl>+%c is pressed,   ASCII code =%d\n",
		(char)theKeySym,theKeySym);          }
	     }
	  printf("%s %d  %d\n","GetKey >>>  theKeySym =",theKeySym,*pos);
*/
    switch(theKeySym)
      {	case	XK_Escape:			return (-1);
	case	XK_BackSpace:
	case	XK_Delete:			return (1);
	case	XK_Return:
	case	XK_KP_Enter:			return (2);
	case	XK_Tab:				return (3);
	case	XK_R7:				return (11);
	case	XK_Left:			return (21);
	case	XK_Up:				return (12);
	case	XK_Right:			return (23);
	case	XK_Down:			return (32);
	case	XK_R9:				return (13);
	case	XK_R13:				return (31);
	case	XK_R15:				return (33);
	case	XK_Shift_L:  /* Left shift */
	case	XK_Shift_R:  /* Right shift */	return (34);
	case	XK_Control_L:/* Left control */
	case	XK_Control_R:/* Right control*/	return (35);
	case	XK_Meta_L:   /* Left meta */
	case	XK_Meta_R:   /* Right meta */  	return (36);
	case	XK_Alt_L:    /* Left alt */
	case	XK_Alt_R:    /* Right alt */   	return (37);
	default:
	  if ( !isprint(theKeySym) )		return (99);
		/* printf("[%1.1s]%d\n",theKeyBuffer,theKeySym); */
		 str[*pos]=theKeySym;		return (0);
      }
}
/*********************************************************************/
GetEsc (theKeyEvent) XKeyEvent theKeyEvent;
{	XComposeStatus	theComposeStatus;
	KeySym		theKeySym;
	int		theKeyBufferMaxLen = 4;
	char		theKeyBuffer[5];
	XLookupString(&theKeyEvent, theKeyBuffer,
		theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
	switch(theKeySym) {
		case	XK_Escape:	return (1) ;
		default: 		return (0) ; }
}
/********************************************************************/
GetName (theKeyEvent, str, lvalue, pos) XKeyEvent theKeyEvent; char str[]; 
					int *pos, lvalue;
{   XComposeStatus	theComposeStatus;
    KeySym		theKeySym;
    int		theKeyBufferMaxLen = 4, i;
    char		theKeyBuffer[5], lsym;
    XLookupString(&theKeyEvent, theKeyBuffer,
		  theKeyBufferMaxLen, &theKeySym, &theComposeStatus);
    switch(theKeySym)
      {	case  XK_BackSpace: 
	case  XK_Delete:
	  /*    if( *pos > 0 ) {for (i=*pos-1; i<lvalue; i++) str[i]=' ';} */
	  *pos = *pos-1;	if( *pos >= 0 ) str[*pos]=' ';
      /* str[*pos]='\0'; */
		return (0);
	case  XK_Return:			return (1);
	case  XK_KP_Enter:			return (1);
	case  XK_Tab:				return (2);
	case  XK_Escape:			return (-1);
	case  XK_R7:				return (11);
	case  XK_Left:				return (21);
	case  XK_Up:				return (12);
	case  XK_Right:				return (23);
	case  XK_Down:				return (32);
	case  XK_R9:				return (13);
	case  XK_R13:				return (31);
	case  XK_R15:				return (33);
	default:
	  if ( theKeySym > 127 )	return(0);
	  if ( *pos < lvalue )
	    {  if ( *pos < 0 ) *pos=0;	/* Restore default name */
	    for ( i=*pos; i<lvalue; i++ ) str[i]=' ';
	    lsym=theKeySym;
	    if ( isalnum(lsym) || lsym=='/' || lsym=='.' || lsym=='_' )
	      { str[*pos]=lsym; *pos=*pos+1;  /*str[*pos]='\0';*/
	      }
	    }	 return (0);
      }
}
/*********************** U-file name setting ***************************/
int askuna (nofbox, una, theNames, unad)
    int *nofbox; char una[], theNames[], unad[];
{   int	i;
	i = askuna_(nofbox, una, theNames, unad);
	return i;
}
int askuna_ (nofbox, una, theNames, unad) int *nofbox; char una[], 
					theNames[], unad[];
{	Window		theWindow;
/*	Cursor		theTextCursor, theBoxCursor;*/
	XEvent		theEvent;
	int	UpLeftx =2, UpLefty =10,/* window corner location	*/
		Width =452, Height =480;/* window size			*/
	int	ix, iy, wbox, hbox, 	/* parameter box corner and size*/
		wsym = 8, hsym = 13,	/* symbol width and height	*/
		lvalue = 10, lname = 4,	/* length of value and name	*/
		xshif = 5, yshif = 3,	/* table corner and 		*/
		nclmn = 3, nparam, 	/* # of columns and parameters	*/
		ixa0=0, iya0=0, ixa,iya,/* arrow (textcursor) position	*/
		i, j, ibox, spos = 0, vpos, lsym=4,
		iret, oldparam, ixold, iyold, ii=-1,
		icol, dcol, irow, drow;
	int	lline, xButton, yButton, ind, ind1, ierr;
	char	Ufile_Name[40], vsym[5];
	/*	FILE	*cfilen; */
	strcpy(vsym," => ");		strcpy(Ufile_Name,"udb/");
	lline = 18; 	hbox = hsym+3;		wbox = wsym*(lline+1);
	vpos= wsym*(lname+4); 	nparam   = *nofbox;
	Width =2*xshif+nclmn*wbox-wsym; 
	Height =2*yshif+((nparam-1)/nclmn+3)*hbox;
	ibox = 0;                       UpLeftx = 2;
	/*  j = XWX-Width-10;	if (j > UpLeftx)  UpLeftx = j; */
	GetRWgeometry (&XRW,&YRW);
	UpLeftx = 2;	   j = XRW-Width-8;	UpLefty = YRW;
	if (j > UpLeftx)   UpLeftx = j;		j = YRW+Height+30;
	if (j > theHeight) UpLefty = theHeight-Height-30;
	theWindow = Open_Window(UpLeftx, UpLefty, Width, Height, 0,
		 "Save data in U-file format", 0,
		  RootWindow(theDisplay,theScreen), theMenuCursor);
	XSelectInput (theDisplay, theWindow, POLL_EV_MASK);
Create_table:
	Change_Color (theGCA, 1, 0);	/* white background, black foreground */
	for (i=1; i <= nparam; i++)
 	{
	ix =wbox*((i-1)%nclmn)+xshif;
	iy =hbox*((i-1)/nclmn)+yshif;
	ind = (i-1)*lname;
	WRITE_ theGCA, ix, hsym+iy, theNames+ind, lname);
	ind = (i-1)*lvalue;
	for (j=0; j<lvalue; j++) una[ind+j]=unad[j];
	WRITE_ theGCA, ix+lname*wsym, hsym+iy, vsym, lsym);
	WRITE_ theGCA,ix+vpos, hsym+iy, una+ind, lvalue);
	XFlush(theDisplay);
	}
	oldparam = nparam;	ixold = ix;	iyold = iy; 
	if ( ii < 0 ) MVPOINTER_ xshif+40,yshif+13);
	Change_Color (theGCA, 50, 0);
	for (i=1; i < nclmn; i++) { ind = xshif+i*wbox-wsym/2;
	LINGCA_ ind,0,ind,Height-2*hbox); }
	i = Height-2*hbox; 	LINGCA_ 0,i,Width,i);
	i--; 		 	LINGCA_ 0,i,Width,i);
	Change_Color (hghGC, 50, 0);
	WRITE_ hghGC, xshif, Height-4, " OK ",4);
	XDrawRectangle (theDisplay,theWindow,hghGC,5L,Height-17L,30L,16L);
	Change_Color (theGCA, 1, 0);
	i = nclmn*wbox/wsym-5;	WRITE_ theGCA, xshif+4*wsym, Height-4, 
	"/<ESC> - done;    Button, <TAB> or Arrow - select   ", i);
	WRITE_ theGCA, xshif, Height-4-hbox, 
	"           Select box and enter U-file name              ", i+4);
Table_control:
	XNextEvent (theDisplay, &theEvent);
	if ( theEvent.xany.window == theRootWindow )
           { ProcessRootWindowEvent (&theEvent);
	     goto Table_control;
	   }
	if (theEvent.type == Expose)	{ii = 0;  goto	Create_table;}
	if (theEvent.type == ButtonPress)
	   {	xButton	= theEvent.xbutton.x;
		yButton = theEvent.xbutton.y;
		if (ibox) { oldparam=ibox; ixold=ix; iyold=iy; }
		ibox = FindBoxNum(xButton-xshif+wsym/2, yButton-yshif,
				  wbox, hbox, nclmn, nparam);
		i = FindBoxNum(xButton-xshif, yButton-Height+hbox+1,
						 	4*wsym, hbox, 1, 1);
		if( i ) { iret = -1;	ibox = oldparam;  
				ixold = ix;	iyold = iy; goto Escend; }
		if (ibox == 0)	goto Table_control;
		ind = (oldparam-1)*lvalue;
Newparam:	i = ibox; 
		ix =wbox*((i-1)%nclmn)+xshif;	iy =hbox*((i-1)/nclmn)+yshif;
		ind = (oldparam-1)*lname; ind1 = (i-1)*lname;
		WRITE_ theGCA, ixold, hsym+iyold, theNames+ind, lname);
 		WRITE_ hghGC, ix, hsym+iy, theNames+ind1, lname);
		WRITE_ theGCA, ixold+lname*wsym, hsym+iyold, vsym, lsym);
		WRITE_ hghGC, ix+lname*wsym, hsym+iy, vsym, lsym);
		ind = (oldparam-1)*lvalue; ind1 = (i-1)*lvalue;
		WRITE_ theGCA, ixold+vpos, hsym+iyold, una+ind, lvalue);
		WRITE_ hghGC, ix+vpos, hsym+iy, una+ind1, lvalue); 
		spos = 0; 	/*MVPOINTER_ ix+vpos, iy+hsym);*/
		ixa = ix+vpos;	iya = iy+hsym+2;
		MoveArrow(theWindow,ixa0,iya0,ixa,iya);
		ixa0 = ixa;		iya0 = iya;
		XFlush(theDisplay);
		goto Table_control;
	   }
	if (ibox == 0)
	{   if(GetEsc(theEvent.xkey)) goto EndDialog;
	else
	    goto Table_control;
	}
	if (theEvent.type == KeyPress)
	{
	   ind = (ibox-1)*lvalue;
	   iret = GetName (theEvent.xkey, una+ind, lvalue, &spos);
/*	printf("File name: \"%.10s\"\n",una+ind); */
Escend:	   if ( spos >= 0 )
	      { ixa = ix+vpos+spos*wsym;	iya = iy+hsym+2;
		MoveArrow(theWindow,ixa0,iya0,ixa,iya);
		ixa0 = ixa;			iya0 = iya;
	      }
	   if ( spos < 0)
	      {	spos = -1;	strncpy(una+ind,unad,lvalue);
	      }
	   if ( strncmp(una+ind,unad,lvalue) ) 
	      {	stcopy(Ufile_Name+4,una+ind,spos);
		ierr = 0;	ierr = access(Ufile_Name,F_OK);
		if ( ierr )
		   { i = nclmn*wbox/wsym-5;
			WRITE_ theGCA, xshif+4*wsym, Height-4, 
		   "/<ESC> - done;    Button, <TAB> or Arrow - select   ", i);
		   }
		if ( !(ierr) )
		   { if ( !(strncmp(una+ind,"          ",lvalue)))
			  strncpy(una+ind,unad,lvalue);
		     else
			{ WRITE_ hghGC,xshif+19*wsym,Height-4,
					"              ",14);
			  WRITE_ hghGC,xshif+33*wsym,Height-4,
					" - file already exists",22);
			  WRITE_ hghGC,xshif+(33-4-spos)*wsym,Height-4,
					Ufile_Name,4+spos);
			  goto	HighlightName;
			}
		   }
	      }
	   i = nclmn*wbox/wsym-5;
	   WRITE_ theGCA,xshif+4*wsym,Height-4, 
		"/<ESC> - done;    Button, <TAB> or Arrow - select   ", i);
HighlightName:
/*	   if ( strncmp(una+ind,unad,lvalue) ) 
	      {	ierr = 0;	ierr = access(Ufile_Name,W_OK);
	      if ( ierr )  printf("Bad  name: \"%.14s\",%d\n",Ufile_Name,ierr);
	      if ( !(ierr) )printf("Good name: \"%.14s\"\n",Ufile_Name);
	      }
*/	   WRITE_ hghGC, ix+vpos, hsym+iy, una+ind, lvalue);
	   if ( iret == 0 )  	goto Table_control;
	   if ( iret == -1 )  	goto EndDialog;
	   if ( iret == 2 ) { oldparam = ibox;	ixold = ix;	iyold = iy;
		spos=0; ibox++; if(ibox > nparam) ibox=1; goto Newparam; }

	   if ( iret > 2 )
	   {    oldparam = ibox;	ixold = ix;	iyold = iy;
		spos=0; stcopy(Ufile_Name+4,una+ind,lvalue);
		dcol=iret%10-2; 	icol =(ibox-1)%nclmn+dcol;
		if (icol<0) icol=0;	if (icol>=nclmn) icol--;
		drow=iret/10-2;		irow =(ibox-1)/nclmn+drow;
		if (irow<0) irow=0;	if (irow>(nparam-1)/nclmn) irow--;
		ibox =irow*nclmn+icol+1;
		if(ibox == oldparam && icol == nclmn-1 && drow != -1) ibox++;
		if(ibox == oldparam && icol == 0       && drow != 1)  ibox--;
		if (ibox < 1) ibox=1;	if (ibox > nparam) ibox=nparam;
		goto Newparam;
	    }
	}
	goto Table_control;
EndDialog:
	RETURNPOINTER_;
	XFlush(theDisplay);
	XDestroyWindow(theDisplay, theWindow);
	ibcursor = -1;
	XFlush(theDisplay); 
	return (0);
}
/****************** Color table setting ********************************/
int setcol  (Atable)  int  *Atable;
{   int	i;
	i = setcol_ (Atable);
	return i;
}
int setcol_ (Atable)  int  *Atable;
{ Window	theWindow;
  XEvent	theEvent;
  int	UpLeftx=200,UpLefty=10,	/* Upper left window corner	*/
	Width =230, Height =220;/* window size			*/
  int	ibox1, wbox1, hbox1, 	/* # of parameter box and size	*/
	llin1 = 8, nclmn1 = 1, 	/* 1st str length & # of columns*/
	xshif1 =9, yshif1 =3,	/* 1st table corner and # of columns*/
	nparam1 = 14, 		/* # of parameters		*/
	ibox2, wbox2 = 27,	/* # of parameter box and width	*/
	nclmn2 = 5,		/* # of columns*/
	xshif2=80, yshif2 =5,	/* 2nd table corner and # of columns*/
	wsym = 8, hsym = 13,	/* font symbol width and height	*/
	i, ii=-1, ix, iy, xButton, yButton, ind;
  char	but[113];
    strcpy(but,"BackgrndCurve 1 Curve 2 Curve 3 Curve 4 Curve 5 Curve 0 ");
    strcpy(but+56,"Color 1 Message History TEXT frgTEXT bkgHIGH frgHIGH bkg");
    ibox1 = 0;	hbox1 = hsym+2;		wbox1 = wsym*(llin1+1);
    Height =2*yshif1+((nparam1-1)/nclmn1+2)*hbox1;
    /*  i = XWX-Width-10;     UpLeftx = 2;    if (i > UpLeftx)  UpLeftx = i; */
    GetRWgeometry (&XRW,&YRW);
    UpLeftx = 2;       i = XRW-Width-8;		UpLefty = YRW;
    if (i > UpLeftx)   UpLeftx = i;		i = YRW+Height+30;
    if (i > theHeight) UpLefty = theHeight-Height-30;
    theWindow = Open_Window(UpLeftx, UpLefty, Width, Height, 0, 
		"Set ASTRA colors",0,
		RootWindow(theDisplay,theScreen), theMenuCursor);
    XSelectInput (theDisplay, theWindow, POLL_EV_MASK);
Create_table:
    i = Height-hbox1; 	LINGCA_ 0,i,Width,i);
    i--;   		LINGCA_ 0,i,Width,i);
    WRITE_ hintGC, xshif1,    Height-2, " OK ",4);
    WRITE_ hintGC, xshif1+60, Height-2, " Click (1) type, (2) color box",30);
    if ( ii < 0 ) MVPOINTER_ xshif1+10,yshif1+8);
    for (i = 1; i < ColorNum; i++)
	{   Change_Color(theGCA, i, 0);	ix = i-1;
	    iy = yshif2+(ix-ix%nclmn2)*3;	ix = xshif2+30*(ix%nclmn2);
	    XFillRectangle(theDisplay,theWindow,theGCA,ix,iy,wbox2,hsym);
	}
Newcolors:
    Change_Color (theGCA, *Atable, *(Atable+1));
    Change_Color (hghGC, *(Atable+22), *(Atable+23));
    for (i=0; i < nparam1; i++)
	{ind = i*llin1;
	 ix =wbox1*(i%nclmn1)+xshif1;	iy =hbox1*(i/nclmn1)+yshif1;
	if (i==0  || i==9 || i==12 || i==13){
		Change_Color(theGCA,*(Atable+22),*(Atable+23));
/*printf("Box No. %d,   Color No.  %d   %d  %s\n",i,*(Atable+22),ind,but+ind); */
}
	 if (i>=1  && i<=8 ) {
		Change_Color(theGCA,*(Atable+2*i),*(Atable+2*i+1));
/*printf("Box No. %d,   Color No.  %d   %d  %s\n",i,*(Atable+2*i),ind,but+ind); */
}
	 if (i==10 || i==11) {
		Change_Color(theGCA,*(Atable+20), *(Atable+21));
/*printf("Box No. %d,   Color No.  %d   %d  %s\n",i,*(Atable+20),ind,but+ind); */
}
	 WRITE_ theGCA, ix, hsym+iy, but+ind, llin1);
	}
    XFlush(theDisplay);
Table_control:
	XNextEvent (theDisplay, &theEvent);
	if ( theEvent.xany.window == theRootWindow )
           { ProcessRootWindowEvent (&theEvent);
	     goto Table_control;
	   }
	if (theEvent.type == Expose)	{ii = 0;  goto	Create_table;}
	if (theEvent.type == ButtonPress)
	{  xButton	= theEvent.xbutton.x;
	   yButton = theEvent.xbutton.y;
	   i = FindBoxNum(xButton-xshif1+wsym/2,yButton-yshif1,
				  wbox1,hbox1,nclmn1,nparam1);
	   if(i) {ibox1 = i;
		  if (ibox1 >= 1 && ibox1<=10) i = *(Atable+2*ibox1-2);
		  if (ibox1 == 11)	i = *(Atable+20);
		  if (ibox1 == 12)	i = *(Atable+21);
		  if (ibox1 == 13)	i = *(Atable+22);
		  if (ibox1 == 14)	i = *(Atable+23); 
		  PutColorName(theWindow,80,199,Width-80,i);
		 }
	   ibox2 = FindBoxNum(xButton-xshif2+2,yButton-yshif2,
				  wbox2+3,hsym+2,nclmn2,ColorNum+5);
	   if (ibox2)
	      {	if (ibox2 >= ColorNum)   ibox2 = 0;
		PutColorName(theWindow,80,199,Width-80,ibox2);
		if(ibox1==1)
		   { *(Atable+1) = ibox2; 
		     for(i=0; i<10; i++) *(Atable+2*i+1) = ibox2; 
		   }
		if(ibox1>1 &&  ibox1<=10) *(Atable+2*ibox1-2) =ibox2;
		if(ibox1==11)	*(Atable+20) = ibox2;
		if(ibox1==12)	*(Atable+21) = ibox2;
		if(ibox1==13)	*(Atable+22) = ibox2;
		if(ibox1==14)	*(Atable+23) = ibox2; 
		goto Newcolors;
	      }
	   ind  = FindBoxNum(xButton-xshif1, yButton-Height+hbox1+1,
					 	4*wsym, hbox1, 1, 1);
	   if(ind)	goto Close_Win;
	   goto Table_control;
	}
    if (GetEsc(theEvent.xkey) != 1) goto Table_control;
Close_Win:
    for( i = 0; i < 24; i++ ) AstraColorNum[i]=*(Atable+i) ;
    RETURNPOINTER_;
    XFlush(theDisplay);
    XDestroyWindow(theDisplay,theWindow);
    ibcursor = -1;
    XFlush(theDisplay); 
    return (0);
}
/************** Returns box # or 0 ***********************************/
FindBoxNum (ix,iy,wBox,hBox,nclmn,nBox)
int	ix,iy;				/*current coordinates*/
int	nclmn,nBox;			/*# of columns and boxes*/
int	wBox,hBox;			/*box width & height*/
{
	int	icol, pnum;
	if ((ix <= 0) || (iy < 0))	return (0);
	icol	= ix/wBox;
	if (icol >= nclmn) 		return (0);
	pnum = (iy/hBox)*nclmn+icol+1;
	if (pnum > nBox) 		return (0);
	return (pnum);
}
/*********************** Sets Window Attributes **********************/
int getatr(j)	int j;
{   int i;   i = getatr_(j);	return (i);  }
int getatr_(j)	int j;
{  int i;
   XWindowAttributes	theAttributes;
   Status		theStatus;
	theStatus =
	  XGetWindowAttributes (theDisplay, theRootWindow, &theAttributes);
	printf("theStatus = %d\n",theStatus);
	i = theAttributes.map_state;
	printf("theAttributes.x = %d,   ",theAttributes.x);
	printf("theAttributes.map_state = %d,   ", i);
	if (i == IsUnmapped)	printf("IsUnmapped\n");
	if (i == IsUnviewable)	printf("IsUnviewable\n");
	if (i == IsViewable)	printf("IsViewable\n");
	return (i);
}
/*********************** Not Used ************************************/
void Del_Menu (NButt, B) int NButt; Button B[];
{	int i;
	for (i=0; i<NButt; i++) 	CLREC_ B[i].mexl, B[i].meyu, 
		B[i].mexr-B[i].mexl+1, B[i].meyd-B[i].meyu+1, False);
}
/*********************** Not Used ************************************/
void Del_button (B)  Button B;
{	CLREC_ B.mexl-2, B.meyu-2, B.mexr-B.mexl+4, B.meyd-B.meyu+4, False);
}
/**************************************************************************/
/*	printf("%d,<%.20s><%.20s>\n",spos,stri,array+ind);
	printf("<Tab>,<%.20s><%.20s>\n",stri,array+ind);
	printf("%d,<%.20s><%.20s>\n",spos,stri,array+ind);
	printf("theEvent.type = %d\n",theEvent.type);
	for (j=0; j < nclmn; j++) printf("%d, ",nwid[j]);	printf("  ");
	for (j=0; j < nclmn; j++) printf("%d, ",nsta[j]);	printf("\n");
for (j=0; j < nclmn; j++) printf("%d + %d,  ",nsta[j],nwid[j]);	printf("\n");
	for (j=0; j < nclmn; j++) printf("\n[%10.10s] ",stri+nwid[j]+1);
	printf("\n%d\n",nclmn);
	printf("Width = %d, Height = %d;  %d x %d\n",Width,Height,nlines,nclmn);
*/
