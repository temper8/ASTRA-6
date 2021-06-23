/*The file includes the C-functions:
	createGC	Change_Color 	Open_Screan	Close_Screan	
	Open_Window	MoveArrow	PutColorName	initDefaultColors	
	PSASetForeground PSADrawRectangle PSADrawLine	PSAMove	
	PSALine		PSADrawLString	PSADrawRString	PSADrawMString
and interfaces for FORTRAN calls
	initvm	endvm	drawvm	rectvm	erasrw	setlin	redraw	puto
	psopen	psclos	pcurso	rcurso	textvm	textbf  pscom	colovm
*/
#include 	<stdio.h>
#include	<string.h>
#include	<unistd.h>	/* Used for sleep(unsigned int) */
#include	<stdlib.h>	/* Used for system(const char*) */
/*#include	<signal.h>*/
#include	<X11/Xlib.h>
#include	<X11/Xutil.h>
#include	<X11/keysym.h>
#include	<X11/keysymdef.h>
#include 	<X11/cursorfont.h>
/*#include 	<X11/olcursor.h>*/
/*#include	"theInteger"*/
#define	INT	long
void	initvm (int*,int*,int*,int*,int*,char*,int*);
void	initvm_(int*,int*,int*,int*,int*,char*,int*);
void	wintitle (char*);		void	wintitle_(char*);
void	redraw(int*);			void	redraw_(int*);
void	erasrw();			void	erasrw_();
void	flushrw();			void	flushrw_();
void	pcurso();			void	pcurso_();
void	endvm();			void	endvm_();
void	pcurxy(int*,int*);		void	pcurxy_(int*,int*);
void	savepm(int*,int*);		void	savepm_(int*,int*);
void	textvm(int*,int*,char*,int*);  	void	textvm_(int*,int*,char*,int*);
void	textbf(int*,int*,char*,int*);  	void	textbf_(int*,int*,char*,int*);
void	textnb(int*,int*,char*,int*);  	void	textnb_(int*,int*,char*,int*);
void	puto(int*,int*,int*,int*);  	void	puto_(int*,int*,int*,int*);
void	createpixmap(int*);		void	createpixmap_(int*);
void	getcolor(int*);			void	getcolor_(int*);
void	colovm(int*);			void	colovm_(int*);
void	colorb(int*);			void	colorb_(int*);
void	testmarkers(int*,int*);		void	testmarkers_(int*,int*);
void	cnmark(int*,int*,int*);		void	cnmark_(int*,int*,int*);
void	cxmark(int*,int*,int*,int*);	void	cxmark_(int*,int*,int*,int*);
void	pscom(char*,int*);		void	pscom_(char*,int*);
void	psopen (char*,int*,int*);	void	psopen_(char*,int*,int*);
void	psclos ();			void	psclos_();
void	drawvm (int*,int*,int*,int*,int*);	
void	drawvm_(int*,int*,int*,int*,int*);
void	rectvm (int*,int*,int*,int*,int*);	
void	rectvm_(int*,int*,int*,int*,int*);
void	cleare (int*,int*,int*,int*,int*);	
void	cleare_(int*,int*,int*,int*,int*);
void	d1line (int*,int*,int*,int*,int*);
void	d1line_(int*,int*,int*,int*,int*);
void	d2line (int*,int*,int*,int*,int*);
void	d2line_(int*,int*,int*,int*,int*);
#define	BORDER_WIDTH	2
#define	NORMAL_WINDOW	0
#define	POP_UP_WINDOW	1
#define	DEFAULT_GEOMETRY	NULL
#define	DEFAULT_FONT1	"variable"
#define	STDfont	"8x13"
#define	HGHfont	"8x13bold"
#define F1sh 13
#define maxPixels 66
#define maxPixmaps 6
/* Pixmaps' allocation:
   Pixmaps[0]=theRootWindow		Pixmaps[3] ECR footprint
   Pixmaps[1] wall, etc.		Pixmaps[4] free
   Pixmaps[2] NBI footprint		
   Pixmaps[5] used for intermediate storage (normally empty)
*/
Display		*theDisplay;
Window		theRootWindow;
GC		theGCA, hghGC, hgh_menuGC, hintGC, iconGC;
Colormap	theColormap;
Drawable	Pixmaps[maxPixmaps] = {0,0,0,0,0,0}, marker[11];
Drawable	theIconPixmap;
Cursor		theRootCursor, thePauseCursor, theMenuCursor;
int		Xmode=0;		/* Default: NoX - Batch mode */
int		XWX,XWY,XWW,XWH,XCorrection=0,YCorrection=0;
int		theScreen, theDepth, theWidth, theHeight, iconState;
int		ColorNum, theCurrentColorNo;
int		AstraColorNum[64] = 
{
        0,0,	/* Color 0,   White (no change allowed) */
        1,0,	/* Color 1,   Black (no change allowed) */
        50,0,	/* Color 2,   Curve 1  (default Red) */
        3,0,	/* Color 3,   Curve 2  (default Blue) */
        62,0,	/* Color 4,   Curve 3  (default VioletRed) */
        37,0,	/* Color 5,   Curve 4  (default MediumSeaGreen) */
        5,0,	/* Color 6,   Curve 5  (default Brown) */
        15,0,	/* Color 7,   Curve 6  (default DarkTurquoise) */
        20,0,	/* Color 8,   Curve 7  (default Goldenrod) */
        29,0,	/* Color 9,   Curve 8  (default LimeGreen) */
        64,0,	/* Color 10,  Curve 9  (default Yellow) */
        30,0,	/* Color 11,  Curve 10 (default DarkGreen) */
        7,0,	/* Color 12,  Curve 11 (default Coral) */
        48,0,	/* Color 13,  Curve 12 (default Pink) */
        30,0,	/* Color 14,  Curve 13 (default Magenta) */
        1,0,	/* Color 15, */
        1,0,	/* Color 16, */
        1,0,	/* Color 17, */
        1,0,	/* Color 18, */
        1,0,	/* Color 19, */
        1,0,	/* Color 20, */
        1,0,	/* Color 21, */
        1,0,	/* Color 22, */
        1,0,	/* Color 23, */
        1,0,	/* Color 24, */
        1,0,	/* Color 25, */
        1,0,	/* Color 26, */
        1,0,	/* Color 27, */
        1,0,	/* Color 28, */
        1,0,	/* Color 29, */
        50,0,	/* Color 30  Message  (default Red) */
        26,0	/* Color 31  History  (default LightBlue) */
 };
unsigned long	theBlackPixel, theWhitePixel, theCurrentColor;
unsigned long	thePixels[maxPixels]; 
char  *theColorNames[maxPixels] = 
   {  "White",		"Black",	"Aquamarine",		/* 0-2 */
      "Blue",		"BlueViolet",	"Brown", 		/* 3-5 */
      "CadetBlue",	"Coral",	"CornflowerBlue",	/* 6-8 */
      "Cyan",		"DarkGreen",	"DarkOliveGreen",	/* 9-11*/
      "DarkOrchid",	"DarkSlateBlue","DarkSlateGrey",	/*12-14*/
      "DarkTurquoise",	"DimGrey",	"Firebrick",		/*15-17*/
      "ForestGreen",	"Gold",		"Goldenrod",		/*18-20*/
      "Grey",		"Green",	"GreenYellow",		/*21-23*/
      "IndianRed",	"Khaki",	"LightBlue",		/*24-26*/
      "LightGrey",	"LightSteelBlue","LimeGreen",		/*27-29*/
      "Magenta",	"Maroon",	"MediumAquamarine",	/*30-32*/
      "MediumBlue",	"MediumForestGreen","MediumGoldenrod",	/*33-35*/
      "MediumOrchid",	"MediumSeaGreen","MediumSlateBlue",	/*36-38*/
      "MediumSpringGreen","MediumTurquoise","MediumVioletRed",	/*39-41*/
      "MidnightBlue",	"Navy",		"Orange",		/*42-44*/
      "OrangeRed",	"Orchid",	"PaleGreen",		/*45-47*/
      "Pink",		"Plum",		"Red",			/*48-50*/
      "Salmon",		"SeaGreen",	"Sienna",		/*51-53*/
      "SkyBlue",	"SlateBlue",	"SpringGreen",		/*54-56*/
      "SteelBlue",	"Tan",		"Thistle",		/*57-59*/
      "Turquoise",	"Violet",	"VioletRed",		/*60-62*/
      "Wheat",		"Yellow",	"YellowGreen",		/*63-65*/
   };
double 	Reg_wA= 468., Reg_hA= 576., PS_xA= 72., PS_yA= 720., PSsc;
FILE    *PSAfile;
	/* FlagPSA = 0 (no hard copy), -1 (white), 1 (other colors) */
int 	FlagPSA=0, FigAcount=1; 
/************************************************************/
void createGC(theWindow, theNewGC, fontName) Window theWindow;
	GC *theNewGC; char fontName[];
{  XGCValues	theGCValues;
   XFontStruct	*fontStruct; 
   unsigned long	theValueMask;
   theValueMask	=0L;
   *theNewGC	= XCreateGC(theDisplay,theWindow,theValueMask,&theGCValues);
	if	(*theNewGC != 0)
	{  fontStruct = XLoadQueryFont (theDisplay, fontName);
	   if(fontStruct != 0) XSetFont(theDisplay,*theNewGC,fontStruct->fid);
	   XSetForeground(theDisplay,*theNewGC,theBlackPixel);
	   XSetBackground(theDisplay,*theNewGC,theWhitePixel);
	}
}
/************************************************************/
void getcolor (i)	int *i;	{      getcolor_(i);	   }
void getcolor_(i)	int *i;	{ *i = theCurrentColorNo; }
/************************************************************/
void Change_Color (anyGC, frgColor, bkgColor)
     int frgColor,bkgColor; GC anyGC;
{      	XSetForeground(theDisplay,anyGC,thePixels[frgColor]);
	theCurrentColor = thePixels[frgColor];
   	XSetBackground(theDisplay,anyGC,thePixels[bkgColor]);
}
/**************************************************************/
void Open_Screan (void)
{	theDisplay	=XOpenDisplay(NULL);
	if(theDisplay == NULL)	{printf(
	">>> ERROR: Cannot establish a connection to the X Server %s\n",
		 XDisplayName(NULL));	 return;}
	theScreen	= DefaultScreen(theDisplay);
	theDepth	= DefaultDepth(theDisplay,theScreen);
	iconState	= 0;
	theBlackPixel	= BlackPixel(theDisplay,theScreen);
	theWhitePixel	= WhitePixel(theDisplay,theScreen);
	theColormap	= XDefaultColormap(theDisplay,theScreen);
/*	theRootCursor	= XCreateFontCursor(theDisplay,47);*/
	theMenuCursor	= XCreateFontCursor(theDisplay,2);
	thePauseCursor	= XCreateFontCursor(theDisplay,130);/*97 or 34*/

/*The display information:*/
	theWidth	= DisplayWidth(theDisplay, theScreen);
	theHeight	= DisplayHeight(theDisplay, theScreen);
/*	printf("%s version %d of the X Window system, X%d R%d.\n",
		ServerVendor(theDisplay),
		VendorRelease(theDisplay),
		ProtocolVersion(theDisplay),
		ProtocolRevision(theDisplay) );			*/
/*	printf( "Display resolution %dx%d,  Display depth %d\n", 
		theWidth, theHeight, theDepth );		*/
}
/*************************************************************/
void Close_Screan (void)
{	/*XWarpPointer(theDisplay, None, theRootWindow, 0,0,0,0,0,420);*/
   int	i;
	for (i = 1; i < maxPixmaps-1; i++) 
	  if (Pixmaps[i]) XFreePixmap(theDisplay,Pixmaps[i]);
	XDestroySubwindows(theDisplay,theRootWindow);
	XDestroyWindow(theDisplay,theRootWindow);
	XFreeCursor(theDisplay, thePauseCursor);
	XFreeCursor(theDisplay, theMenuCursor);
/*	XFreeCursor(theDisplay, theRootCursor);  */
	XCloseDisplay(theDisplay);
/*	XFlush(theDisplay);			 */
}
/***************************************************************/
void createIconBitmap(theNewWindow)	Window theNewWindow;
{
//#define	afont "12x24bold"
//#define	afont "lucidasans-bolditalic-24"
#define	afont0 "9x18bold"
#define	afont "-cronyx-times-bold-r-normal--34-240-100-100-p-124-koi8-r"
#define	pfont "-cronyx-times-medium-r-normal--34-240-100-100-p-124-koi8-r"
  char A[]="A";
// Draw Icon Russian flag
  switch(1){
  case 0:			// Draw Icon Russian flag
    theIconPixmap = XCreatePixmap(theDisplay,theNewWindow,48,48,theDepth);  
    createGC (theNewWindow, &iconGC, afont0);
    XSetForeground(theDisplay,iconGC,theWhitePixel);
    XFillRectangle(theDisplay,theIconPixmap,iconGC,0,0,48,48);
    XSetForeground(theDisplay,iconGC,thePixels[50]);
    XDrawString(theDisplay,theIconPixmap,iconGC,1, 13,"ASTRA",5);
    XSetForeground(theDisplay,iconGC,thePixels[3]);
    XFillRectangle(theDisplay,theIconPixmap,iconGC,0,17,48,16);
    //    XSetForeground(theDisplay,iconGC,theWhitePixel);
    XDrawString(theDisplay,theIconPixmap,iconGC,1, 30,"ASTRA",5);
    XSetForeground(theDisplay,iconGC,thePixels[50]);
    XFillRectangle(theDisplay,theIconPixmap,iconGC,0,33,48,16);
    //    XSetForeground(theDisplay,iconGC,theWhitePixel);
    XDrawString(theDisplay,theIconPixmap,iconGC,1, 45,"ASTRA",5);
    break;
  case 1:			// Draw Icon "A-Icon"
    theIconPixmap = XCreatePixmap(theDisplay,theNewWindow,43,43,theDepth);  
    createGC (theNewWindow, &iconGC, afont);
    XSetLineAttributes(theDisplay,iconGC,3,LineSolid,CapNotLast,JoinMiter);
    XSetForeground(theDisplay,iconGC,theWhitePixel);
    XFillRectangle(theDisplay,theIconPixmap,iconGC,0,0,43,43);
    XSetForeground(theDisplay,iconGC,thePixels[50]);	//33 - MediumBlue
    XDrawArc(theDisplay,theIconPixmap,iconGC, 2, 6, 38, 30, 0, 64*360);
    XDrawString(theDisplay,theIconPixmap,iconGC,10, 31,"A",1);
    XWriteBitmapFile(theDisplay,"A_Icon",theIconPixmap,43,43,-1,-1);
    break;
  case 2:			 // Draw Icon "Phi-medium"
    theIconPixmap = XCreatePixmap(theDisplay,theNewWindow,21,20,theDepth);  
    createGC (theNewWindow, &iconGC, pfont);
    XSetLineAttributes(theDisplay,iconGC,3,LineSolid,CapNotLast,JoinMiter);
    XSetForeground(theDisplay,iconGC,theWhitePixel);
    XFillRectangle(theDisplay,theIconPixmap,iconGC,0,0,21,20);
    XSetForeground(theDisplay,iconGC,thePixels[50]);
    sprintf(A,"%c",0xe6); printf("A=\"%c\"\n",*A);
    XDrawString(theDisplay,theIconPixmap,iconGC,-2,21,A,1);
    XWriteBitmapFile(theDisplay,"PhiIcon",theIconPixmap,21,20,-1,-1);
    break;
  case 3:			// Draw Icon "Phi-bold"
    theIconPixmap = XCreatePixmap(theDisplay,theNewWindow,32,24,theDepth);  
    createGC (theNewWindow, &iconGC, afont);
    XSetLineAttributes(theDisplay,iconGC,3,LineSolid,CapNotLast,JoinMiter);
    XSetForeground(theDisplay,iconGC,theWhitePixel);
    XFillRectangle(theDisplay,theIconPixmap,iconGC,0,0,32,24);
    XSetForeground(theDisplay,iconGC,thePixels[50]);
    sprintf(A,"%c",0xe6);
    XDrawString(theDisplay,theIconPixmap,iconGC,3,23,A,1);
// Save theIcon:
    XWriteBitmapFile(theDisplay,"PhiIcon",theIconPixmap,32,24,-1,-1);
    break;
  }
}
/***************************************************************/
void makeIcon(theNewWindow)	Window theNewWindow;
{	/* Get pre-defined two-color Icon Pixmap */ 
#define PhiIcon_width 21
#define PhiIcon_height 20
static unsigned char PhiIcon_bits[] = {
   0xff, 0xf1, 0x1f, 0x1f, 0x00, 0x1f, 0xc7, 0x71, 0x1c, 0xe3, 0xf1, 0x18,
   0xf1, 0xf1, 0x11, 0xf1, 0xf1, 0x11, 0xf8, 0xf1, 0x03, 0xf8, 0xf1, 0x03,
   0xf8, 0xf1, 0x03, 0xf8, 0xf1, 0x03, 0xf8, 0xf1, 0x03, 0xf8, 0xf1, 0x03,
   0xf8, 0xf1, 0x03, 0xf8, 0xf1, 0x03, 0xf1, 0xf1, 0x11, 0xf1, 0xf1, 0x11,
   0xe3, 0xf1, 0x18, 0xc7, 0x71, 0x1c, 0x1f, 0x00, 0x1f, 0xff, 0xf1, 0x1f};
#define A_Icon_width 43
#define A_Icon_height 43
static unsigned char A_Icon_bits[] = {
   0xff, 0xff, 0xff, 0xff, 0xff, 0x07, 0xff, 0xff, 0xff, 0xff, 0xff, 0x07,
   0xff, 0xff, 0xff, 0xff, 0xff, 0x07, 0xff, 0xff, 0xff, 0xff, 0xff, 0x07,
   0xff, 0xff, 0xff, 0xff, 0xff, 0x07, 0xff, 0xff, 0x00, 0xf8, 0xff, 0x07,
   0xff, 0x1f, 0x00, 0xc0, 0xff, 0x07, 0xff, 0x07, 0x00, 0x00, 0xff, 0x07,
   0xff, 0x01, 0xfe, 0x03, 0xfc, 0x07, 0x7f, 0xc0, 0xcf, 0x1f, 0xf0, 0x07,
   0x3f, 0xf8, 0xcf, 0xff, 0xe0, 0x07, 0x1f, 0xfc, 0x8f, 0xff, 0xc1, 0x07,
   0x0f, 0xff, 0x87, 0xff, 0x87, 0x07, 0x8f, 0xff, 0x07, 0xff, 0x8f, 0x07,
   0x87, 0xff, 0x03, 0xff, 0x0f, 0x07, 0xc3, 0xff, 0x03, 0xff, 0x1f, 0x06,
   0xe3, 0xff, 0x03, 0xfe, 0x3f, 0x06, 0xe3, 0xff, 0x01, 0xfe, 0x3f, 0x06,
   0xf1, 0xff, 0x11, 0xfc, 0x7f, 0x04, 0xf1, 0xff, 0x18, 0xfc, 0x7f, 0x04,
   0xf1, 0xff, 0x38, 0xf8, 0x7f, 0x04, 0xf1, 0xff, 0x38, 0xf8, 0x7f, 0x04,
   0xf1, 0x7f, 0x3c, 0xf8, 0x7f, 0x04, 0xf1, 0x7f, 0x7c, 0xf0, 0x7f, 0x04,
   0xf1, 0x3f, 0x00, 0xf0, 0x7f, 0x04, 0xe3, 0x3f, 0x00, 0xe0, 0x3f, 0x06,
   0xe3, 0x3f, 0xfe, 0xe0, 0x3f, 0x06, 0xc3, 0x1f, 0xff, 0xc1, 0x1f, 0x06,
   0x87, 0x1f, 0xff, 0xc1, 0x0f, 0x07, 0x8f, 0x0f, 0xff, 0x80, 0x8f, 0x07,
   0x0f, 0x03, 0x3c, 0x00, 0x86, 0x07, 0x1f, 0xfc, 0xff, 0xff, 0xc1, 0x07,
   0x3f, 0xf8, 0xff, 0xff, 0xe0, 0x07, 0x7f, 0xc0, 0xff, 0x1f, 0xf0, 0x07,
   0xff, 0x01, 0xfe, 0x03, 0xfc, 0x07, 0xff, 0x07, 0x00, 0x00, 0xff, 0x07,
   0xff, 0x1f, 0x00, 0xc0, 0xff, 0x07, 0xff, 0xff, 0x00, 0xf8, 0xff, 0x07,
   0xff, 0xff, 0xff, 0xff, 0xff, 0x07, 0xff, 0xff, 0xff, 0xff, 0xff, 0x07,
   0xff, 0xff, 0xff, 0xff, 0xff, 0x07, 0xff, 0xff, 0xff, 0xff, 0xff, 0x07,
   0xff, 0xff, 0xff, 0xff, 0xff, 0x07};
 theIconPixmap = XCreatePixmapFromBitmapData (theDisplay, theNewWindow,
//	 PhiIcon_bits, PhiIcon_width, PhiIcon_height,	//Use PhiIcon
	 A_Icon_bits, A_Icon_width, A_Icon_height,	//Use A_Icon
//	 thePixels[3], theWhitePixel, theDepth);	// White on blue
//	 theWhitePixel, theBlackPixel, theDepth);	// Black on white
	 theWhitePixel, thePixels[50], theDepth);	// Red on white

//   createIconBitmap(theNewWindow);
 return;
}
/*************************************************************/
Window	Open_Window (x, y, width, height, flag, theTitle, 
		     iconicState, theParent, theCursor)
   int    x, y, width, height, flag, iconicState;
   char   theTitle[];
   Window theParent;
   Cursor theCursor;
{	XSetWindowAttributes 	theWindowAttributes;
	XSizeHints		theSizeHints;
	unsigned long 		theWindowMask;
	Window 			theNewWindow;
	XWMHints 		theWMHints;
	XClassHint		theClassHint;
	theWindowAttributes.border_pixel	= theBlackPixel;
	theWindowAttributes.background_pixel	= theWhitePixel;
	theWindowAttributes.cursor		= theCursor;
	if(flag == POP_UP_WINDOW) {
		theWindowAttributes.override_redirect	= True;
		theWindowAttributes.save_under		= True;
		theWindowMask 
		= CWBackPixel | CWBorderPixel | CWCursor | CWOverrideRedirect;}
	else {	theWindowAttributes.override_redirect	= False;
		theWindowMask 	= CWBackPixel | CWBorderPixel | CWCursor; }
	theNewWindow = XCreateWindow (theDisplay, theParent,
           	x, y, width, height, BORDER_WIDTH, theDepth, InputOutput,
		CopyFromParent, theWindowMask, &theWindowAttributes);
	theWMHints.input = True;
	if(iconicState == 0)	theWMHints.initial_state= NormalState;
		else		theWMHints.initial_state= IconicState;
	theWMHints.flags  = InputHint | StateHint;
	makeIcon(theNewWindow); theWMHints.icon_pixmap = theIconPixmap;
	theWMHints.flags  = theWMHints.flags | IconPixmapHint;
	XSetWMHints (theDisplay, theNewWindow, &theWMHints);
	XStoreName (theDisplay, theNewWindow, theTitle);
	theSizeHints.flags = USPosition | PSize | PMinSize | PMaxSize ;
	theSizeHints.x = x;
	theSizeHints.y = y;
	theSizeHints.width = width;
	theSizeHints.height = height;
	theSizeHints.min_width = width;		/* Resizing */
	theSizeHints.min_height = height;	/* forbidden */
	theSizeHints.max_width = width;
	theSizeHints.max_height = height;
	XSetNormalHints (theDisplay, theNewWindow, &theSizeHints);
	XMapWindow (theDisplay, theNewWindow);
	XFlush (theDisplay);
	if(flag == NORMAL_WINDOW) sleep(1);
	/* printf("newWindowID:  %d\n",theNewWindow); */
	return theNewWindow;
}
/*******************************************************************/
void wintitle_ (Title)	char Title[];
{   char	theTitle[120];
	strcpy(theTitle,Title);
	/* printf("The new title \"%s\",\n",theTitle); */
	XStoreName (theDisplay, theRootWindow, theTitle);
}
void wintitle (Title)
	char Title[];
{	wintitle_ (Title);	}
/*******************************************************************/
void MoveArrow (wind,fromx,fromy,tox,toy) 
	Window wind; int fromx, fromy, tox, toy;
{   int	Xx,Xy;
	if (fromx > 0)
	{ Xx = fromx;	Xy = fromy;
	  Change_Color (theGCA, 0, 0);
/*	  XCopyPlane (theDisplay, triarr, wind, theGCA, 
				0, 0, 7, 7, fromx, fromy, 0x01);*/
	  XDrawLine(theDisplay,wind,theGCA,Xx,Xy,Xx,Xy+5);
	  XDrawLine(theDisplay,wind,theGCA,Xx+1,Xy+1,Xx+1,Xy+5);
	  XDrawLine(theDisplay,wind,theGCA,Xx-1,Xy+1,Xx-1,Xy+5);
	  XDrawLine(theDisplay,wind,theGCA,Xx+2,Xy+3,Xx+2,Xy+5);
	  XDrawLine(theDisplay,wind,theGCA,Xx-2,Xy+3,Xx-2,Xy+5);
	  XDrawLine(theDisplay,wind,theGCA,Xx+3,Xy+5,Xx+3,Xy+5);
	  XDrawLine(theDisplay,wind,theGCA,Xx-3,Xy+5,Xx-3,Xy+5);
	}
	if (tox > 0)
	{ Xx = tox;	Xy = toy;
	  Change_Color (theGCA, 3, 0);	/* 3 - blue */
/*	  XCopyPlane (theDisplay, marker[7], wind, theGCA, 
				0, 0, 7, 7, tox, toy, 0x01);*/
	  XDrawLine(theDisplay,wind,theGCA,Xx,Xy,Xx,Xy+5);
	  XDrawLine(theDisplay,wind,theGCA,Xx+1,Xy+1,Xx+1,Xy+5);
	  XDrawLine(theDisplay,wind,theGCA,Xx-1,Xy+1,Xx-1,Xy+5);
	  XDrawLine(theDisplay,wind,theGCA,Xx+2,Xy+3,Xx+2,Xy+5);
	  XDrawLine(theDisplay,wind,theGCA,Xx-2,Xy+3,Xx-2,Xy+5);
	  XDrawLine(theDisplay,wind,theGCA,Xx+3,Xy+5,Xx+3,Xy+5);
	  XDrawLine(theDisplay,wind,theGCA,Xx-3,Xy+5,Xx-3,Xy+5);
	}
	Change_Color (theGCA, 1, 0);
}
/*=================================================================*/
/*	(  x >= 0,      y >= 0)  - upper left corner  of the rectangle
	(x+w <= XWW, y+h <= XWH) - lower right corner of the rectangle
	the function puts a mark 
			with a centre at the position (*x,*y),
			preliminary it saves an image under the mark,
			finally, it puts the image at the same place */
void puto (x,y,color,type)	int *x, *y, *color, *type;
{    puto_(x,y,color,type);	}
void puto_(x,y,color,type)	int *x, *y, *color, *type;
{  int	ix, iy, dst_x, dst_y, format;
   XImage *theImage;	unsigned int iw, ih;	unsigned long mask;
   iw = 7;	ih = 7;		ix = *x-iw/2;	iy = *y-ih/2;
   if(ix < 0)	ix = 0;		if(ix+iw > XWW) iw = XWW-ix;
   if(iy < 0)	iy = 0;		if(iy+ih > XWH) ih = XWH-iy;
/*   dst_x = ix;	dst_y = iy;
   mask = -1;	format = XYPixmap;
   theImage = XGetImage ( theDisplay, theRootWindow, ix, iy, iw, ih,
		mask, format);	*/
   ix = *x+7;		iy = *y+7;
   colovm(color);	format = *type-1;
   XCopyPlane (theDisplay, marker[format],  theRootWindow, theGCA, 
		0, 0, 7, 7, ix, iy, 0x01);
/*  ix = 0;	iy = 0;
   XPutImage(theDisplay, theRootWindow, theGCA, theImage,
		ix, iy, dst_x, dst_y, iw, ih);	*/
}
/*=================================================================*/
/*	(  x >= 0,      y >= 0)  - upper left corner  of the rectangle
	(x+w <= XWW, y+h <= XWH) - lower right corner of the rectangle
	the function puts a mark (circle) of the color 1 (red) 
			with a centre at the position (*x,*y),
			preliminary it saves an image under the mark,
			finally, it puts the image at the same place */
void testmarkers (x,y)	int *x, *y;
{	testmarkers_ (x,y);	}
void testmarkers_(x,y)	int *x, *y;
{  int	ix, iy, coord[2], dst_x, dst_y, format;
   XImage *theImage;	unsigned int iw, ih;	unsigned long mask;
	iw = 7	;	ih = 7;		ix = *x-iw/2;	iy = *y-ih/2;
	if(ix < 0)	ix = 0;		if(ix+iw > XWW) iw = XWW-ix;
	if(iy < 0)	iy = 0;		if(iy+ih > XWH) ih = XWH-iy;
/*	printf("x = %d; y = %d; w = %d; h = %d\n",ix,iy,iw,ih);	     	*/
/*	printf("x = %d -> %d; y = %d -> %d\n", *x, dst_x, *y, dst_y); 	*/
/*	dst_x = ix+20;		dst_y = iy+20;				*/
	dst_x = ix;		dst_y = iy;
	mask = -1;		format = XYPixmap;
	theImage = XGetImage (theDisplay, theRootWindow, ix, iy, iw, ih,
		   mask, format);		/* Store image in theImage*/
	coord[0]=*x-10;	coord[1] = *y-10;
	ix = 1;	iy = 8;	cxmark(&coord[0], &coord[1], &ix, &iy);
	ix = 0;	    iy = 0;
	XPutImage(theDisplay, theRootWindow, theGCA, theImage,
		  ix, iy, dst_x, dst_y, iw, ih);
}
/*=================================================================*/
void cxmark (x,y,color,type)	int *x, *y, *color, *type;		
{	cxmark_ (x,y,color,type);	}
void cxmark_(x,y,color,type)	int *x, *y, *color, *type;
{  int	i, ix, iy, ic;
/*	colovm(color);*/
	ic = *color;	ix = *x+7;	iy = *y+7;
	for (i = 0; i < 11; i++) { 	ic = 2;	colovm(&ic);
	XCopyPlane (theDisplay, marker[i],  theRootWindow, theGCA, 
				0, 0, 7, 7, ix, iy, 0x01);
				ic++; ix += 20; }
}
/*=================================================================*/
void cnmark (xy,color,type)	int *xy, *color, *type;		
{	cnmark_ (xy,color,type);	}
void cnmark_(xy,color,type)	int *xy, *color, *type;
/*type:	0 - diamond;	1 - circle;	2 - plus;	3 - aster;
	4 - cross;	5 - square;	6 - arrow;	7 - triarr; */
{  int	i, ix, iy, ic;
/*printf("     %d  %d,  color = %d, type = %d\n",*xy, *(xy+1),*color,*type); */
	ic = *color;	ix = *xy+7;	iy = *(xy+1)+7;
	/*	colovm(color);	colovm(&ic); */
	i = *type; XCopyPlane (theDisplay, marker[i], 
		    theRootWindow, theGCA, 0, 0, 7, 7, ix, iy, 0x01);
}
/*=================================================================*/
void PutColorName (wind,ix,iy,iw,num) Window wind; int ix, iy, iw, num;
{  int	len, i;
	len = strlen(theColorNames[num]);
	i = (iw-7*len)/2;
	XClearArea(theDisplay, wind, ix, iy, iw, 16, False);
	XDrawImageString(theDisplay, wind, hintGC, ix+i, iy+13,
		theColorNames[num],len);
}
/*=================================================================*/
void initDefaultColors () 
{   XColor theRGBColor, theHardwareColor;	int theStatus, i;
    if (theDepth > 1)
    {   ColorNum = maxPixels;
	for (i = 0; i < maxPixels; i++)
	{  theStatus = XLookupColor (theDisplay, theColormap,
		  theColorNames[i], &theRGBColor, &theHardwareColor);
	   if (theStatus != 0)
	   {	theStatus = XAllocColor (theDisplay,
	 			theColormap,&theHardwareColor);
		if (theStatus != 0) thePixels[i] = theHardwareColor.pixel; 
		else		    thePixels[i] = theBlackPixel; 
	   }
	}
    }
    else		/* Monochrome system */
    {	for (i = 0; i < maxPixels; i++)
	{   if (strcmp("White", theColorNames[i]) == 0)
			thePixels[i] = theWhitePixel; 
	     else	thePixels[i] = theBlackPixel;
	}
    }
}
/*********************  call initvm(XWX,XWY,XWW,XWH,COLTAB) ************/
void initvm_ (x,y,wid,hei,Atable,Title,titlen)
	int *x,*y,*wid,*hei,*Atable,*titlen; char Title[];
{    static char
       /* Use the unix command "bitmap" to create the data */
/*	aster1_bits[]  = {0x49, 0x2a, 0x1c, 0x7f, 0x1c, 0x2a, 0x49},
	cross1_bits[]  = {0x00, 0x22, 0x14, 0x08, 0x14, 0x22, 0x00}, */
	circle_bits[]  = {0x1c, 0x22, 0x41, 0x41, 0x41, 0x22, 0x1c},
	plus_bits[]    = {0x08, 0x08, 0x08, 0x7f, 0x08, 0x08, 0x08},
	diamond_bits[] = {0x08, 0x14, 0x22, 0x41, 0x22, 0x14, 0x08},
	square_bits[]  = {0x7f, 0x41, 0x41, 0x41, 0x41, 0x41, 0x7f},
	cross_bits[]   = {0x41, 0x22, 0x14, 0x08, 0x14, 0x22, 0x41},
	aster_bits[]   = {0x22, 0x14, 0x08, 0x7f, 0x08, 0x14, 0x22},
	arrow_bits[]   = {0x08, 0x1c, 0x1c, 0x2a, 0x2a, 0x08, 0x08},
	triarr_bits[]  = {0x08, 0x1c, 0x1c, 0x3e, 0x3e, 0x7f, 0x00},
	rtriar_bits[]  = {0x08, 0x18, 0x38, 0x78, 0x38, 0x18, 0x08},
	fsquare_bits[] = {0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f},
	fdiamon_bits[] = {0x08, 0x1c, 0x3e, 0x7f, 0x3e, 0x1c, 0x08},
	dies_bits[]    = {0x14, 0x14, 0x7f, 0x14, 0x7f, 0x14, 0x14},
	fcircle_bits[] = {0x1c, 0x3e, 0x7f, 0x7f, 0x7f, 0x3e, 0x1c};
    static char rct3_bits[] = {0x07, 0x05, 0x07};
    static char rct5_bits[] = {0x1f, 0x11, 0x11, 0x11, 0x1f};
    static char rct4_bits[] = {0x0f, 0x09, 0x09, 0x0f};
    Window	openWindow();	int i, ix, iy;
    char	theTitle[120];
	strncpy(theTitle,Title,*titlen);
	if ( strncmp(theTitle,"BGD",3) == 0 )   return; /* Background mode */
	Xmode = 1;
	/*printf("short int %d\n",sizeof(short));
	  printf("      int %d\n",sizeof(int));
	  printf("long  int %d\n",sizeof(long));
	  printf("    float %d\n",sizeof(float));
 	  printf("   double %d\n",sizeof(double));*/
        if (sizeof(int) != 4) 
	  printf(" >>> Warning >>> Incompatibility in color table\n");
	Open_Screan ();
	initDefaultColors();
	XWX = *x;      	XWY = *y;	XWW = *wid;	XWH=*hei;
	XWX = theWidth-XWW-2*(BORDER_WIDTH+5);	XWY = 5;
	theRootWindow	=Open_Window(XWX,XWY,XWW,XWH,0,theTitle, 
	   	iconState, RootWindow(theDisplay,theScreen),theRootCursor);
	GetRWgeometry (&ix,&iy);
/* Each bloody XWindow Manager tries to do it in its own unique manner */
	XCorrection = ix-XWX;		YCorrection = iy-XWY;
	createGC (theRootWindow, &hghGC, HGHfont);
	createGC (theRootWindow, &theGCA, STDfont);
	createGC (theRootWindow, &hgh_menuGC, "8x13bold");
	createGC (theRootWindow, &hintGC, "variable");
	XSetLineAttributes (theDisplay, hgh_menuGC, 2L, LineSolid,
						CapRound, JoinRound);
	XDrawRectangle (theDisplay,theRootWindow,theGCA,0L,0L,*wid-1L,*hei-1L);
	XWarpPointer(theDisplay, None, theRootWindow, 0,0,0,0,*wid/2,*hei-120);
/*        for( i = 0; i < 64; i++ ) AstraColorNum[i]=*(Atable+i); */

	Pixmaps[0] = theRootWindow;
	marker[0] = XCreateBitmapFromData (theDisplay, theRootWindow,
						fdiamon_bits, 7, 7);
	marker[1] = XCreateBitmapFromData (theDisplay, theRootWindow,
						circle_bits, 7, 7);
	marker[2] = XCreateBitmapFromData (theDisplay, theRootWindow,
						square_bits, 7, 7);
	marker[3] = XCreateBitmapFromData (theDisplay, theRootWindow,
						aster_bits, 7, 7);
	marker[4] = XCreateBitmapFromData (theDisplay, theRootWindow,
						cross_bits, 7, 7);
	marker[5] = XCreateBitmapFromData (theDisplay, theRootWindow,
						triarr_bits, 7, 7);
	marker[6] = XCreateBitmapFromData (theDisplay, theRootWindow,
						fcircle_bits, 7, 7);
	marker[7]= XCreateBitmapFromData (theDisplay, theRootWindow,
						fsquare_bits, 7, 7);
	marker[8] = XCreateBitmapFromData (theDisplay, theRootWindow,
						diamond_bits, 7, 7);
	marker[9] = XCreateBitmapFromData (theDisplay, theRootWindow,
						rct4_bits, 4, 4);
	marker[10] = XCreateBitmapFromData (theDisplay, theRootWindow,
						rct5_bits, 5, 5);
}
void initvm (x,y,wid,hei,Atable,Title,titlen)
	int *x,*y,*wid,*hei,*Atable,*titlen; char Title[];
{	initvm_ (x,y,wid,hei,Atable,Title,titlen);	}
/***************************************************************/
void createpixmap_(id)	int *id;
{   int i=0, j, k;
	if (Pixmaps[*id] == 0)	{ Pixmaps[*id] = 
	XCreatePixmap(theDisplay,theRootWindow,XWW,XWH,theDepth);
	j = XWW-1; k = XWH-1;  cleare(id,&i,&i,&j,&k); }
}
void createpixmap (id)	int *id;
{	createpixmap_(id);	}
/***************************************************************/
void endvm ()
{	Close_Screan();	}
void endvm_()
{	Close_Screan();	}
/**************************  Redraw pixmap  *******************/
void redraw (id)	int *id;
{	redraw_(id);	}
void redraw_(id)	int *id;
{   if (*id == 0)	goto Flushit;
    if (Pixmaps[*id] == 0)	return;
    XSetFunction(theDisplay,theGCA,GXxor);
    XCopyArea(theDisplay,Pixmaps[*id],theRootWindow,theGCA,
		0,0,XWW,XWH,0,0);
    XSetFunction(theDisplay,theGCA,GXcopy);
 Flushit:
    XFlush(theDisplay);
}
/**********************  Copy pixmap  **************************/
void savepm (source, destination)	int *source, *destination;
{	savepm_(source, destination);	}
void savepm_(source, destination)	int *source, *destination;
{
    XCopyArea(theDisplay,Pixmaps[*source],Pixmaps[*destination],theGCA,
		0,0,XWW,XWH,0,0);
    XFlush(theDisplay);
}
/***************** Clear the root window ******************/
void erasrw ()
{	XClearWindow(theDisplay, theRootWindow);	}
void erasrw_()
{	XClearWindow(theDisplay, theRootWindow);	}
/**********************************************************/
void cleare (id,ix,iy,iw,ih)	int *id, *ix, *iy, *iw, *ih;
{	cleare_ (id,ix,iy,iw,ih);	}
void cleare_ (id,ix,iy,iw,ih)	int *id, *ix, *iy, *iw, *ih;
{  if (Pixmaps[*id] == 0)	return;
   if (*id == 0) XSetForeground(theDisplay,theGCA, theWhitePixel);
   if (*id)      XSetForeground(theDisplay,theGCA, theBlackPixel);
   XFillRectangle(theDisplay,Pixmaps[*id],theGCA,*ix,*iy,*iw,*ih);
		 XSetForeground(theDisplay,theGCA, theCurrentColor);
}
/**********************************************************/
void rcurso ()
{	XUndefineCursor(theDisplay, theRootWindow);	}
void rcurso_ ()
{	XUndefineCursor(theDisplay, theRootWindow);	}
/************** Enter "xfd -center -fn cursor" to get all cursors available */
void pcurso ()
{	pcurso_ ();	}
void pcurso_ ()
{	XDefineCursor(theDisplay, theRootWindow, thePauseCursor);
}
/************** Enter "xfd -center -fn cursor" to get all cursors available */
void pcurxy (x,y)	int *x,*y;
{	pcurxy_ (x,y);	}
void pcurxy_ (x,y)	int *x,*y;
{	Window	theRoot, theChild;
	int	rootX, rootY;
	unsigned  int	buttonState;

	XDefineCursor (theDisplay, theRootWindow, thePauseCursor);
	XQueryPointer (theDisplay, theRootWindow, &theRoot, &theChild,
		&rootX, &rootY, x, y, &buttonState);
	/*printf("Cursor position   x = %d,   y = %d\n",*x,*y);*/
}
/**********************************************************/
/*setlin (style,width) long *style,*width;
{	setlin_ (style,width);	}
setlin_ (style,width) long *style,*width;
{	unsigned	long	theMask;
	XGCValues	theXGCV;
	theXGCV.line_style	=*style;
	theXGCV.line_width	=*width;
	theMask	=GCLineWidth|GCLineStyle;
	XChangeGC(theDisplay,theGCA,theMask,&theXGCV);
}*/
/***********  Line from (x1,y1) to (x2,y2)  *********************/
void drawvm (id,x1,y1,x2,y2)	int *id,*x1,*y1,*x2,*y2;
{	drawvm_ (id,x1,y1,x2,y2);	}
void drawvm_ (id,x1,y1,x2,y2)	int *id,*x1,*y1,*x2,*y2;
{   int		Xx1, Xy1, Xx2, Xy2;
    double	dx1, dy1, dx2, dy2;
    Xx1 =*x1+10;	Xy1 =*y1+10;	Xx2 =*x2+10;	Xy2 =*y2+10;
    if (*id)	XSetForeground(theDisplay,theGCA,~theCurrentColor);
    XDrawLine(theDisplay,Pixmaps[*id],   theGCA,Xx1,Xy1,Xx2,Xy2);
    if (*id)	XSetForeground(theDisplay,theGCA,theCurrentColor);
	if(FlagPSA > 0) { dx1=Xx1*PSsc; dx2=Xx2*PSsc; dy1=Xy1*PSsc;
		dy2=Xy2*PSsc;	PSADrawLine(dx1, dy1, dx2, dy2);
		/*	printf("Drawn\n"); */
			}
}

void drawline_(int *id,int *iXY,int *n)
{
  int j,n2;
  int		Xx1, Xy1, Xx2, Xy2;
  double	dx1, dy1, dx2, dy2;
  if (*id)	XSetForeground(theDisplay,theGCA,~theCurrentColor);

  n2	=2*(*n)-3;
  for(j=0; j < n2; j+=2){
    Xx1 =iXY[j]+10;
    Xy1 =iXY[j+1]+10;
    Xx2 =iXY[j+2]+10;
    Xy2 =iXY[j+3]+10;
    XDrawLine(theDisplay,Pixmaps[*id],theGCA,Xx1,Xy1,Xx2,Xy2);
  }
  if (*id)	XSetForeground(theDisplay,theGCA,theCurrentColor);

  if(FlagPSA > 0){
    n2	=2*(*n);
    for(j=0; j < n2; j+=2){
      dx1 =(iXY[j]+10)*PSsc;
      dy1 =(iXY[j+1]+10)*PSsc;
      fprintf(PSAfile, "%12.5e %12.5e\n",dx1+PS_xA, PS_yA-dy1);
    }
    fprintf(PSAfile, "moveto %d{lineto}repeat\n",(*n)-1);
  }
}
void drawline(int *id,int *iXY,int *n)
{
  drawline_(id,iXY,n);
}


/***********  Draw a line from (x1,y1) to (x2,y2)  ******************/
void d1line (id,x1,y1,x2,y2)	int *id,*x1,*y1,*x2,*y2;
{	d1line_ (id,x1,y1,x2,y2);	}
void d1line_(id,x1,y1,x2,y2)	int *id,*x1,*y1,*x2,*y2;
{   int		Xx1, Xy1, Xx2, Xy2;
    double	dx1, dy1, dx2, dy2;
    Xx1 =*x1+10;	dy1 =*y1/10.+10;	Xy1 =dy1;
    Xx2 =*x2+10;	dy2 =*y2/10.+10;	Xy2 =dy2;
    if (*id)	XSetForeground(theDisplay,theGCA,~theCurrentColor);
    XDrawLine(theDisplay,Pixmaps[*id],   theGCA,Xx1,Xy1,Xx2,Xy2);
    if (*id)	XSetForeground(theDisplay,theGCA,theCurrentColor);
	if(FlagPSA > 0) { dx1=Xx1*PSsc; dx2=Xx2*PSsc;
			  dy1=dy1*PSsc; dy2=dy2*PSsc;
		       /* printf("%5d %d %11.4e %11.4e\n",*y1,*y2,dy1,dy2); */
		       /* PSADrawLine(dx1, dy1, dx2, dy2); */
	fprintf(PSAfile, "%11.4e %11.4e %11.4e %11.4e moveto lineto\n"
		,dx2+PS_xA, PS_yA-dy2, dx1+PS_xA, PS_yA-dy1);
			}
}
/* For g95 compiler use the line:
void d1polyline_(int *id,long *iXY,int *n)
*/
void d1polyline_(int *id,int *iXY,int *n)
{
  int j,n2;
  int		Xx1, Xy1, Xx2, Xy2;
  double	dx1, dy1, dx2, dy2;
  if (*id)	XSetForeground(theDisplay,theGCA,~theCurrentColor);

  n2	=2*(*n)-3;
  for(j=0; j < n2; j+=2){
    Xx1 =iXY[j]+10;
    dy1 =iXY[j+1]/10.+10;
    Xy1 =dy1;
    Xx2 =iXY[j+2]+10;
    dy2 =iXY[j+3]/10.+10;
    Xy2 =dy2;
    //    if(j == 0)printf("%d   %d   %d   %d   %d\n",j,Xx1,Xy1,Xx2,Xy2);
    XDrawLine(theDisplay,Pixmaps[*id],theGCA,Xx1,Xy1,Xx2,Xy2);
  }
  if (*id)	XSetForeground(theDisplay,theGCA,theCurrentColor);

  if(FlagPSA > 0){
    n2	=2*(*n);
    for(j=0; j < n2; j+=2){
      dx1 =(iXY[j]+10)*PSsc;
      dy1 =(iXY[j+1]/10.+10)*PSsc;
      fprintf(PSAfile, "%12.5e %12.5e\n",dx1+PS_xA, PS_yA-dy1);
    }
    fprintf(PSAfile, "moveto %d{lineto}repeat\n",(*n)-1);
#ifdef H
    for(j=0; j < n2; j+=2){
      Xx1 =iXY[j]+10;
      dy1 =iXY[j+1]/10.+10;
      Xy1 =dy1;
      Xx2 =iXY[j+2]+10;
      dy2 =iXY[j+3]/10.+10;
      Xy2 =dy2;

      dx1 =Xx1*PSsc;
      dx2 =Xx2*PSsc;
      dy1 =dy1*PSsc;
      dy2 =dy2*PSsc;
      /* printf("%5d %d %11.4e %11.4e\n",*y1,*y2,dy1,dy2); */
      /* PSADrawLine(dx1, dy1, dx2, dy2); */
      fprintf(PSAfile, "%11.4e %11.4e %11.4e %11.4e moveto lineto\n"
	      ,dx2+PS_xA, PS_yA-dy2, dx1+PS_xA, PS_yA-dy1);
    }
#endif
  }
}
void d1polyline(int *id,int *iXY,int *n)
{
  d1polyline_(id,iXY,n);
}
/***********  Draw a line from (x1,y1) to (x2,y2)  ******************/
void d2line_(id,x1,y1,x2,y2)	int *id,*x1,*y1,*x2,*y2;
{   int		Xx1, Xy1, Xx2, Xy2;
    double	dx1, dy1, dx2, dy2;
    dx1 =*x1/10.+10;	dy1 =*y1/10.+10;    Xx1 =dx1;	Xy1 =dy1;
    dx2 =*x2/10.+10;	dy2 =*y2/10.+10;    Xx2 =dx2;	Xy2 =dy2;
    if (*id)	XSetForeground(theDisplay,theGCA,~theCurrentColor);
    XDrawLine(theDisplay,Pixmaps[*id],   theGCA,Xx1,Xy1,Xx2,Xy2);
    if (*id)	XSetForeground(theDisplay,theGCA,theCurrentColor);
    if(FlagPSA > 0) { dx1=dx1*PSsc; dx2=dx2*PSsc;
		      dy1=dy1*PSsc; dy2=dy2*PSsc;
	fprintf(PSAfile, "%11.4e %11.4e %11.4e %11.4e moveto lineto\n"
		     ,dx2+PS_xA, PS_yA-dy2, dx1+PS_xA, PS_yA-dy1);
	/*	printf("%5d %d %11.4e %11.4e\n",*y1,*y2,dy1,dy2); */
		    }
}
void d2line (id,x1,y1,x2,y2)	int *id,*x1,*y1,*x2,*y2;
{	d2line_ (id,x1,y1,x2,y2);	}

void d2polyline_(int *id,int *iXY,int *n)
{
  int j,n2;
  int		Xx1, Xy1, Xx2, Xy2;
  double	dx1, dy1, dx2, dy2;
  if (*id)	XSetForeground(theDisplay,theGCA,~theCurrentColor);
  n2	=2*(*n)-3;
  for(j=0; j < n2; j+=2){
    dx1 =iXY[j]/10.+10;
    dy1 =iXY[j+1]/10.+10;
    Xx1 =dx1;
    Xy1 =dy1;
    dx2 =iXY[j+2]/10.+10;
    dy2 =iXY[j+3]/10.+10;
    Xx2 =dx2;
    Xy2 =dy2;
    XDrawLine(theDisplay,Pixmaps[*id],theGCA,Xx1,Xy1,Xx2,Xy2);
  }
  if (*id)	XSetForeground(theDisplay,theGCA,theCurrentColor);

  if(FlagPSA > 0){
    fprintf(PSAfile, "newpath\n");
    n2	=2*(*n);
    for(j=0; j < n2; j+=2){
      dx1 =(iXY[j]/10.+10)*PSsc;
      dy1 =(iXY[j+1]/10.+10)*PSsc;
      fprintf(PSAfile, "%12.5e %12.5e\n",dx1+PS_xA, PS_yA-dy1);
    }
    fprintf(PSAfile, "moveto %d{lineto}repeat\n",(*n)-1);

#ifdef H
    for(j=0; j < n2; j+=2){
      dx1 =iXY[j]/10.+10;
      dy1 =iXY[j+1]/10.+10;
      Xx1 =dx1;
      Xy1 =dy1;
      dx2 =iXY[j+2]/10.+10;
      dy2 =iXY[j+3]/10.+10;
      Xx2 =dx2;
      Xy2 =dy2;
      dx1 *=PSsc;
      dx2 *=PSsc;
      dy1 *=PSsc;
      dy2 *=PSsc;
      /* printf("%5d %d %11.4e %11.4e\n",*y1,*y2,dy1,dy2); */
      /* PSADrawLine(dx1, dy1, dx2, dy2); */
      fprintf(PSAfile, "%11.4e %11.4e %11.4e %11.4e moveto lineto\n"
	      ,dx2+PS_xA, PS_yA-dy2, dx1+PS_xA, PS_yA-dy1);
    }
#endif
  }
}
void d2polyline(int *id,int *iXY,int *n)
{
  d2polyline_(id,iXY,n);
}

/********************************************************************/
void rectvm (id,x,y,width,height) int *id,*x,*y,*width,*height;
{	rectvm_ (id,x,y,width,height);	}
void rectvm_ (id,x,y,width,height) int *id,*x,*y,*width,*height;
{	int	Xx,Xy;	unsigned int	Xwidth,Xheight;
	double x1, y1, w1, h1;
	Xx=*x;	Xy=*y;	Xwidth=*width;	Xheight=*height;
    if (*id)	XSetForeground(theDisplay,theGCA,~theCurrentColor);
    XDrawRectangle(theDisplay,Pixmaps[*id],theGCA,Xx,Xy,Xwidth,Xheight);
    if (*id)	XSetForeground(theDisplay,theGCA,theCurrentColor);
    if (FlagPSA > 0)  {	x1=Xx*PSsc; y1=Xy*PSsc; 
			w1=Xwidth*PSsc; h1=(Xheight-80)*PSsc;
			PSADrawRectangle(x1, y1, w1, h1); }
}
/***********************************************************************/
void colovm (clnumb) int *clnumb;
{	colovm_ (clnumb);	}
void colovm_(clnumb) int *clnumb;
{	int	inum;
	if (*clnumb < 32) inum=2*(*clnumb);	else inum=62;
/* if (*clnumb >4 && *clnumb < 9)	printf("%d\n",*clnumb); */
	Change_Color (theGCA, AstraColorNum[inum], AstraColorNum[inum+1]);
	theCurrentColorNo = *clnumb;
	if(FlagPSA) { PSASetForeground(inum); }
}
/***********************************************************************/
void colorb (clnumb) int *clnumb;
{	colorb_ (clnumb);	}
void colorb_(clnumb) int *clnumb;
{	int	inum;
	if (*clnumb < 32) inum=2*(*clnumb);	else inum=62;
/* if (*clnumb >4 && *clnumb < 9)	printf("%d\n",*clnumb); */
	Change_Color (hghGC, AstraColorNum[inum], AstraColorNum[inum+1]);
	theCurrentColorNo = *clnumb;
	if(FlagPSA) { PSASetForeground(inum); }
}
/****** ! Note ! the function cannot be called from FORTRAN *******/
void textbf (x,y,str,str_len) int *x,*y,*str_len; char str[];
{    textbf_(x,y,str,str_len);	}
void textbf_(x,y,str,str_len) int *x,*y,*str_len; char str[];
{	int	Xx,Xy,Xstlen;
	double x1, y1, psth;
	Xx =*x+10;	Xy =*y+10;	Xstlen=*str_len;
	XDrawImageString(theDisplay,theRootWindow,hghGC,Xx,Xy,str,Xstlen);
	if(FlagPSA > 0) { x1=Xx*PSsc; y1=Xy*PSsc;
		psth=F1sh*PSsc;	PSADrawLString(x1,y1,str,Xstlen,psth); }
}
/*******************************************************************/
void textvm (x,y,str,str_len) int *x,*y,*str_len; char str[];
{    textvm_(x,y,str,str_len);	}
void textvm_(x,y,str,str_len) int *x,*y,*str_len; char str[];
{	int	Xx,Xy,Xstlen;
	double x1, y1, psth;
	Xx =*x+10;	Xy =*y+10;	Xstlen=*str_len;
	XDrawImageString(theDisplay,theRootWindow,theGCA,Xx,Xy,str,Xstlen);
	if(FlagPSA > 0) { x1=Xx*PSsc; y1=Xy*PSsc;
		psth=F1sh*PSsc;PSADrawLString(x1,y1,str,Xstlen,psth);
	}
}
/*******************************************************************/
void textnb (x,y,str,str_len) int *x,*y,*str_len; char str[];
{    textnb_(x,y,str,str_len);	}
void textnb_(x,y,str,str_len) int *x,*y,*str_len; char str[];
{	int	Xx,Xy,Xstlen;
	double x1, y1, psth;
	Xx =*x+10;	Xy =*y+10;	Xstlen=*str_len;
	XSetForeground(theDisplay, theGCA,~theCurrentColor);
	XDrawString(theDisplay,Pixmaps[2],theGCA,Xx,Xy,str,Xstlen);
	XSetForeground(theDisplay, theGCA, theCurrentColor);
	if(FlagPSA > 0) { x1=Xx*PSsc; y1=Xy*PSsc;
		psth=F1sh*PSsc;PSADrawLString(x1,y1,str,Xstlen,psth);
	}
}
/*******************************************************************/
void pscom_(str,str_len)
     int *str_len; char str[];
{
  int	j;
  if(FlagPSA > 0) {
    fprintf(PSAfile,"%% ");
    for (j=0; j < *str_len; j++) fprintf(PSAfile,"%c",str[j]);
    fprintf(PSAfile,"\n");
  }
  return;
}
void pscom(str,str_len)
     int *str_len; char str[];
{
  pscom_(str,str_len);
}

/**************************************************************************/
/* Function px2pnm added, 
/**************************************************************************/
static int icount=0, kPnm=0;

int px2pnm_(name)
     char *name;
{
  int i,L;
  char str[0x100];
  FILE *lF;

  sprintf(str,"(xwd -silent -id %ld|xwdtopnm > %s%3.3d.ppm) >& /dev/null",
	  theRootWindow,name,kPnm);
  system(str);
  /* printf("%s\n",str); */
  /* printf("%s%3.3d.pnm  %d %d\n",name,kPnm,XWW,XWH); */
  kPnm++;
  return(0);
}

void PSnumber(name,count)
     char *name; int *count;
{
  int	i;

  icount	=*count;
  for (i=0; i<128; i++){
    if( name[i] == '\0' ) break;
  }
 con1:
  if(*count<10){
    sprintf( name+i, "-%1d.ps\0", *count); goto con2;
  }
  if(*count<100){
    sprintf( name+i, "-%2d.ps\0", *count);
    goto con2;
  }
  if(*count<1000){
    sprintf( name+i, "-%3d.ps\0", *count);
    goto con2;
  }
  sprintf( name+i, ".ps\0");
  return;
 con2:
  PSAfile = fopen(name,"r");
  if(PSAfile){
    fclose(PSAfile);
    *count=*count+1;
    goto con1;
  }
}
/**************************************************************************/
void psopen (PSname, sty, iret) char *PSname; int *sty, *iret;
{
  psopen_(PSname, sty, iret);
}

void psopen_(PSname, sty, iret)
     char *PSname; int *sty, *iret;
{
  if(FlagPSA){
    *iret=-1;
    return;
  }

  PSnumber(PSname,&FigAcount);
  PSAfile = fopen(PSname, "w");

  if(PSAfile != NULL){
    FlagPSA = 1;
    *iret=0;
    fprintf(PSAfile,"%%!PS-Adobe-2.0\n");
    fprintf(PSAfile,"/Lshow {exch dup scale show dup scale}def\n");
    fprintf(PSAfile,"/Rshow {dup stringwidth pop 3 -1 roll\n");
    fprintf(PSAfile,"dup 3 -1 roll mul neg 0 rmoveto\n");
    fprintf(PSAfile,"dup scale show dup scale}def\n");
    fprintf(PSAfile,"/Mshow {dup stringwidth pop 3 -1 roll\n");
    fprintf(PSAfile,"dup 3 -1 roll -0.5 mul mul 0 rmoveto\n");
    fprintf(PSAfile,"dup scale show dup scale}def\n");
    fprintf(PSAfile,"%%Page: %d %d\n",icount,icount);
    switch(*sty){
    case 0:
      PSsc=.7;
      fprintf(PSAfile,"%%Figure %2d\n %g setlinewidth\n",FigAcount,0.4);
      fprintf(PSAfile, "/Courier-Bold findfont 1.02 scalefont setfont\n");
      break;
    case 1:
      PSsc=.98;
      fprintf(PSAfile,"%%Figure %2d\n %g setlinewidth\n",FigAcount,0.5);
      fprintf(PSAfile,"756 0 translate\n");
      fprintf(PSAfile,"90 rotate\n");
      fprintf(PSAfile, "/Courier-Bold findfont 1.02 scalefont setfont\n");
      break;
    default:
      printf(" >>> PSOPEN >>> Unknown style option: %d\n",*sty);
    }
    PSASetForeground(0);
    FigAcount++;
    return;
  }
  else{
    *iret=1;
    return;
  }
}
/**************************************************************************/
void psclos ()
{	psclos_();	}
void psclos_ ()
{	if (FlagPSA) {	fprintf(PSAfile, "stroke\nshowpage\n\n");
				FlagPSA = 0;		fclose(PSAfile);	}
}
/**************************************************************************/
PSASetForeground(i) int i;
{   XColor theRGBColor;		float	r, g, b;
/* Zakharov's colors:
    static char    *CNm[24] = {
	"0 0 0",	"0 0 .6",	"0 .5 0",	"0 .6 .6",
	".6 0 0",	".6 0 .6",	".7 .2 .2",	".8 .8 .8",
	".7 .7 .7",	"0 0 1",	"0 1 0",	"0 1 1",
	"1 0 0",	"1 0 1",	"1 1 0",	"1 1 1",
	"0 0 0",	"0 0 0",	"0 0 0",	"0 0 0",
	"0 0 0",	"0 0 0",	"0 0 0",	"0 0 0"	};*/
	fprintf(PSAfile, "stroke\n");
	theRGBColor.pixel = thePixels[AstraColorNum[i]];
	XQueryColor (theDisplay, theColormap, &theRGBColor);
	if ((theRGBColor.red & theRGBColor.green & theRGBColor.blue) != 0xffff)
	FlagPSA = 1;	else	{FlagPSA = -1;	return; }
	r = theRGBColor.red; g = theRGBColor.green; b = theRGBColor.blue; 
	r = r/0xffff;	     g = g/0xffff;	    b = b/0xffff;
	fprintf(PSAfile, "%5.3f %5.3f %5.3f setrgbcolor\n", r,g,b);
/*	printf("ACN[i] = %d,   RGB  %x  %x  %x,   %5.3f %5.3f %5.3f\n",
 AstraColorNum[i],theRGBColor.red,theRGBColor.green,theRGBColor.blue,r,g,b); */
}
/**************************************************************************/
PSADrawRectangle(x, y, w, h) double x, y, w, h;
{	fprintf(PSAfile, "newpath %%Rectangle\n");
	fprintf(PSAfile, "%10.3e %10.3e  moveto\n", x+PS_xA, PS_yA-y);
   fprintf(PSAfile,"%10.3e %10.3e %10.3e %10.3e rlineto rlineto\n",w,0.,0., -h);
   fprintf(PSAfile,"%10.3e %10.3e %10.3e %10.3e rlineto rlineto\n",-w,0.,0., h);
	fprintf(PSAfile, "stroke\n");
}
/**************************************************************************/
PSADrawLine(x1, y1, x2, y2)
double x1, y1, x2, y2;
{
  fprintf(PSAfile, "%11.4e %11.4e %11.4e %11.4e moveto lineto\n"
	  ,x2+PS_xA, PS_yA-y2, x1+PS_xA, PS_yA-y1);
}

/**************************************************************************/
void PSAMove(x, y) double x, y;
{	fprintf(PSAfile, "%10.3e %10.3e moveto\n", x+PS_xA, PS_yA-y);
}
/**************************************************************************/
void PSALine(x, y) double x, y;
{	fprintf(PSAfile, "%10.3e %10.3e lineto\n", x+PS_xA, PS_yA-y);
}
/**************************************************************************/
PSADrawLString(x, y, ln, n, psfsc) double x, y, psfsc; char *ln; int n;
{	char	tlin[256];	double	psfsc1;
	int	i;
	psfsc1=1./psfsc;
	for(i=0; i <= n-1 ; ++i) *(tlin+i) = *(ln+i);	*(tlin+n) = '\0';
	fprintf(PSAfile, "%10.3e %10.3e moveto\n", x+PS_xA, PS_yA-y);
	fprintf(PSAfile, "%15.8e %15.8e (%s) Lshow\n", psfsc1, psfsc, tlin);
}
/**************************************************************************/
PSADrawRString(x, y, ln, n, psfsc) double x, y, psfsc; char *ln; int n;
{	char	tlin[256];	double	psfsc1;
	int	i;
	psfsc1=1./psfsc;
	for(i=0; i <= n-1 ; ++i) *(tlin+i) = *(ln+i);	*(tlin+n) = '\0';
	fprintf(PSAfile, "%10.3e %10.3e moveto\n", x+PS_xA, PS_yA-y);
	fprintf(PSAfile, "%15.8e %15.8e (%s) Rshow\n", psfsc1, psfsc, tlin);
}
/**************************************************************************/
PSADrawMString(x, y, ln, n, psfsc) double x, y, psfsc; char *ln; int n;
{	char	tlin[256];	double	psfsc1;
	int	i;
	psfsc1=1./psfsc;
	for(i=0; i <= n-1 ; ++i) *(tlin+i) = *(ln+i);	*(tlin+n) = '\0';
	fprintf(PSAfile, "%10.3e %10.3e moveto\n", x+PS_xA, PS_yA-y);
	fprintf(PSAfile, "%15.8e %15.8e (%s) Mshow\n", psfsc1, psfsc, tlin);
}
/**************************************************************************/
