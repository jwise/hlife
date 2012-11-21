/*
 *   This program displays very large life universes, allowing you to
 *   both zoom in (more than one screen pixel per cell) and out (more
 *   than one cell per screen pixel).  It works with (and indeed
 *   integrates) hlife.  It does *not* animate, it only allows you
 *   to pan and scan and zoom a static universe.
 *
 *   It works on X11 and Win32.  No menus; just keyboard shortcuts
 *   of arrow keys and + (to zoom in) and - (to zoom out).  Also,
 *   you can click the mouse on a particular area to zoom in on that
 *   area.
 *
 *   Usage:
 *      hdraw pattern
 *
 *   where pattern is in RLE, picture, or macrocell format.  (For
 *   pictures larger than about 2^32 on a side, only macrocell format
 *   suffices).
 *
 *   More options are included, since this is just a display-enabled
 *   hlife; see the options to hlife for information.  The main
 *   difference is if neither -i nor -m are given, it just displays
 *   the base generation rather than computing anything.
 *
 *   Windows compilation (vc++):
 *      vcvars
 *      cl hdraw.c user32.lib gdi32.lib
 *   X11 compilation (gcc/Linux):
 *      cc -o hdraw hdraw.c -L /usr/X11R6/lib -lX11
 *   X11 compilation (MacOS):
 *      cc -O -o hdraw hdraw.c -L/usr/X11R6/lib -lX11 -I/usr/X11R6/include
 *
 *   This program endeavors to be reasonably quick, and it does not
 *   use explicit multiprecision arithmetic (although the currently
 *   displayed location needs to be expressed in multiprecision
 *   arithmetic).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
 *   Some macros and etc. that are used by the system-specific routines.
 */
#define TITLE "hlife 0.95 Radical Eye Software"
#define UI(a) ui(a)
int vieww, viewh ;          /* view height and width */
struct node *drawroot ;     /* the universe to draw */
unsigned int ibigbuf[128] ; /* a shared buffer for up to 32x32 pixels */
unsigned char *bigbuf = (unsigned char *)ibigbuf ;
int pmag = 1 ;              /* the magnification when zoomed in */
int redraw = 1 ;            /* do we still need to redraw? */
void handlecmd(int cmd, char *arg) ; /* handle an arg */
void draw(struct node *) ;           /* generalized draw */
void ui() ;                          /* overall UI called from within hlife */
void drawpixel(int x, int y) {
   bigbuf[((63-y) << 3) + (x >> 3)] &= (~128) >> (x & 7) ;
}
/*
 *   The Windows code is bare-bones Win32 API calls.  No fancy MFC
 *   here (I can't afford it).
 */
#ifdef _WIN32
#include <windows.h>
HDC globalhdc, memhdc ;
HBRUSH hbr ;
HWND hwnd;
HBITMAP hbm ;
void openscreen() {}
void renderbm(int x, int y) {
   SetBitmapBits(hbm, 512, bigbuf) ;
   if (pmag == 1) {
      BitBlt(globalhdc, x, viewh-y-64, 64, 64, memhdc, 0, 0, SRCCOPY) ;
   } else {
      StretchBlt(globalhdc, x*pmag, (viewh-y-64)*pmag, 64*pmag, 64*pmag,
                                               memhdc, 0, 0, 64, 64, SRCCOPY) ;
   }
   memset(bigbuf, -1, 512) ;
}
void clearrect(int x1, int y1, int w, int h) {
   PatBlt(globalhdc, x1 * pmag, (viewh-y1-h) * pmag, w * pmag, h * pmag,
                                                                   WHITENESS) ;
}
void win32_redraw(struct node *u) {
   RECT rect ;
   GetClientRect(hwnd, &rect);
   vieww = rect.right = (rect.right + pmag - 1) / pmag ;
   viewh = rect.bottom = (rect.bottom + pmag - 1) / pmag ;
   draw(u) ;
}
void domsgs(struct node *u) {
   drawroot = u ;
   while (1) {
      MSG msg;
      int msgcnt = 0 ;
      while (1) {
         while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
            msgcnt++ ;
            if (msg.message == WM_QUIT)
               exit(0) ;
            TranslateMessage(&msg);
            DispatchMessage(&msg);
         }
         if (msgcnt > 0)
            break ;
         WaitMessage() ;
      }
      if (redraw)
         win32_redraw(u) ;
   }
}
long APIENTRY WndProc(HWND hwnd, UINT msg, UINT wParam, LONG lParam) {
   PAINTSTRUCT ps;
   int thistime ;
   char mbuf[80] ;
   switch (msg) {
      case WM_CHAR:
         switch(wParam) {
case '+':
case '=':
            handlecmd('+', 0) ;
            break ;
case '-':
            handlecmd('-', 0) ;
            break ;
         }
         break ;
      case WM_KEYDOWN:
         switch(wParam) {
case VK_LEFT:
            handlecmd('h', 0) ;
            break ;
case VK_RIGHT:
            handlecmd('l', 0) ;
            break ;
case VK_UP:
            handlecmd('j', 0) ;
            break ;
case VK_DOWN:
            handlecmd('k', 0) ;
            break ;
         }
         break ;
      case WM_LBUTTONDOWN:
         sprintf(mbuf, "%d %d", LOWORD(lParam)/pmag, HIWORD(lParam)/pmag) ;
         handlecmd('z', mbuf) ;
         break ;
      case WM_PAINT:
         BeginPaint(hwnd, &ps) ;
         if (drawroot)
            win32_redraw(drawroot) ;
         EndPaint(hwnd, &ps) ;
         return 0;
      case WM_DESTROY:
         PostQuitMessage(0);
         break;
      default:
         return(DefWindowProc(hwnd, msg, wParam, lParam));
   }
   return TRUE;
}
char cmdbuf[1000] ;
char *argv[100] ;
int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
                     LPSTR lpszCmdParam, int nCmdShow) {
   static char szAppName[] = TITLE ;
   WNDCLASS    wc;
   if(! hPrevInstance) {
      wc.style         = CS_HREDRAW | CS_VREDRAW;
      wc.lpfnWndProc   = WndProc;
      wc.cbClsExtra    = 0;
      wc.cbWndExtra    = 0;
      wc.hInstance     = hInstance;
      wc.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
      wc.hCursor       = LoadCursor(NULL, IDC_ARROW);
      wc.hbrBackground = GetStockObject(WHITE_BRUSH);
      wc.lpszMenuName  = NULL;
      wc.lpszClassName = szAppName;
      RegisterClass(&wc);
   }
   hwnd = CreateWindow(szAppName, TITLE,
      WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT,
      400, 400, NULL, NULL, hInstance, NULL);
   vieww = viewh = 400 ;
   globalhdc = GetDC(hwnd) ;
   hbr = GetStockObject(WHITE_BRUSH) ;
   memhdc = CreateCompatibleDC(globalhdc) ;
   hbm = CreateBitmap(64, 64, 1, 1, bigbuf) ;
   SelectObject(memhdc, hbm) ;
   strcpy(cmdbuf, lpszCmdParam) ;
   ShowWindow(hwnd, SW_SHOWNORMAL) ;
   {
      int i, argc = 0, nargc = 1 ;
      char *p = cmdbuf ;
      argv[argc++] = "hlife" ;
      while (*p) {
         while (*p && *p <= ' ')
            p++ ;
         if (*p == 0)
            break ;
         argv[argc++] = p ;
         while (*p > ' ')
            p++ ;
         if (*p == 0)
            break ;
         *p++ = 0 ;
      }
      argv[argc] = 0 ;
      for (i=1; i<argc; i++)
         argv[nargc++] = argv[i] ;
      argv[nargc] = 0 ;
      life_main(nargc, argv) ;
   }
}
#define main life_main
#else
/*
 *   The X11 code is bare-bones X11.
 */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
int whitePixel, blackPixel ;
char *progname = TITLE ;
char *geometry = "400x400+0+0" ;
Display *display ;
GC gc ;
void error(char *s) {
   fprintf(stderr, "%s\n", s) ;
   if (*s == '!')
      exit(10) ;
}
void init_X() {
   if (display == 0) {
      if ((display = XOpenDisplay(0))==0)
         error("! can't open display") ;
      whitePixel = WhitePixel(display, DefaultScreen(display)) ;
      blackPixel = BlackPixel(display, DefaultScreen(display)) ;
      gc = DefaultGC(display, DefaultScreen(display)) ;
      XSetGraphicsExposures(display, gc, 1) ;
   }
}
Window window ;
XSetWindowAttributes xsetwindowattributes ;
void openscreen() {
   XSizeHints xsizehints ;
   int x, y, w, h ;

   init_X() ;
   XParseGeometry(geometry, &x, &y, &w, &h) ;
   vieww = w ;
   viewh = h ;
   xsetwindowattributes.background_pixel = whitePixel ;
   xsetwindowattributes.border_pixel = blackPixel ;
   xsetwindowattributes.backing_store = Always ;
   window = XCreateWindow(display, DefaultRootWindow(display), x, y,
      w, h, 2, DefaultDepth(display, DefaultScreen(display)),
      InputOutput, DefaultVisual(display, DefaultScreen(display)),
      CWBackPixel|CWBorderPixel, &xsetwindowattributes) ;
   if (window == 0)
      error("! couldn't open window") ;
   xsizehints.x = x ;
   xsizehints.y = y ;
   xsizehints.width = w ;
   xsizehints.height = h ;
   xsizehints.flags = USSize | USPosition ;
   XSetStandardProperties(display, window, progname, progname, None,
      0, 0, &xsizehints) ;
   XSelectInput(display, window,
              ButtonPressMask|ExposureMask|KeyPressMask|StructureNotifyMask) ;
   XMapWindow(display, window) ;
}
XImage ximage =
    { 0, 0, 0, XYBitmap, 0, MSBFirst, 16, MSBFirst, 16, 1, 0, 1 } ;
/*
 *   When rendering zoomed in, X11 doesn't have any stretch-blits call.
 *   So we draw a bunch of rectangles.  We try to combine adjacent
 *   rectangles of the same color.  This isn't amazingly fast at the
 *   zoomed-in size of 2x2, which is pretty common.  Oh well.
 */
XRectangle xrs[64] ; /* used when we need to scale things. */
int rect0, rect1 ;
void clearrects() {
   rect0 = 0 ;
   rect1 = 63 ;
}
void flushrects() {
   if (rect0 > 0) {
      XSetForeground(display, gc, blackPixel) ;
      XFillRectangles(display, window, gc, &xrs[0], rect0) ;
      rect0 = 0 ;
   }
   if (rect1 < 63) {
      XSetForeground(display, gc, whitePixel) ;
      XFillRectangles(display, window, gc, &xrs[rect1+1], 63-rect1) ;
      rect1 = 63 ;
   }
}
void outrect(int val, int x1, int x2, int y) {
   if (rect1 > rect0)
      flushrects() ;
   if (val)
      val = rect1-- ;
   else
      val = rect0++ ;
   xrs[val].x = x1 * pmag ;
   xrs[val].y = (viewh - y - 1) * pmag ;
   xrs[val].width = (x2 - x1) * pmag ;
   xrs[val].height = pmag ;
}
void renderbm(int x, int y) {
   XSetForeground(display, gc, whitePixel) ;
   XSetBackground(display, gc, blackPixel) ;
   ximage.width = 64 ;
   ximage.height = 64 ;
   ximage.bytes_per_line = 8 ;
   ximage.data = (char *)bigbuf ;
   if (pmag == 1) {
      XPutImage(display, window, gc, &ximage, 0, 0, x, viewh-y-64, 64, 64) ;
   } else { /* yuck!  build rectangles by hand. */
      int xx, yy ;
      clearrects() ;
      for (yy=0; yy<64; yy++) {
         char *p = (char *)bigbuf + (yy << 3) ;
         for (xx=0; xx<64; ) {
            int run = xx ;
            if (p[xx >> 3] & (128 >> (xx & 7))) {
               xx++ ;
               while (xx < 64 && p[xx >> 3] & (128 >> (xx & 7)))
                  xx++ ;
               outrect(1, run + x, xx + x, 63 - yy + y) ;
            } else {
               xx++ ;
               while (xx < 64 && 0 == (p[xx >> 3] & (128 >> (xx & 7))))
                  xx++ ;
               outrect(0, run + x, xx + x, 63 - yy + y) ;
            }
         }
      }
      flushrects() ;
   }
   memset(bigbuf, -1, 512) ;
}
void domsgs(struct node *u) {
   XEvent event ;
   XKeyEvent *kev ;
   XExposeEvent *xev ;
   XConfigureEvent *cev ;
   XButtonEvent *bev ;
   char mbuf[30] ;
   int len, key = 0 ;
   long arg = 0 ;
   while (1) {
      int first = 1 ;
      redraw = 0 ;
      while (first || XPending(display)) {
         first = 0 ;
         XNextEvent(display, &event) ;
         switch(event.type) {
case KeyPress:
            kev = (XKeyEvent *)&event ;
            len = XLookupString(kev, mbuf, 30, 0, 0) ;
            if (len != 1) {
               key = XKeycodeToKeysym(display, kev->keycode, 0) ;
               len = 1 ;
               switch(key) {
case 65361:       key = 'h' ; break ;
case 65362:       key = 'j' ; break ;
case 65363:       key = 'l' ; break ;
case 65364:       key = 'k' ; break ;
default:          break ;
               }
            } else
               key = mbuf[0] ;
            if ('A' <= key && key <= 'Z')
               key += 32 ;
            else if (1 <= key && key <= 27 && key != 8)
               key += 96 ;
            if ('0' <= key && key <= '9')
               arg = 10 * arg + key - '0' ;
            else {
               switch(key) {
case 'l':         handlecmd('l', 0) ;
                  break ;
case 'h':         handlecmd('h', 0) ;
                  break ;
case 'j':         handlecmd('j', 0) ;
                  break ;
case 'k':         handlecmd('k', 0) ;
                  break ;
case 'q':         exit(0) ;
                  break ;
case '+':
case '=':         handlecmd('+', 0) ;
                  break ;
case '-':         handlecmd('-', 0) ;
                  break ;
default:          break ;
               }
            }
case Expose: case GraphicsExpose:
            xev = (XExposeEvent *)&event ;
            redraw++ ;
            break ;
case ConfigureNotify:
            cev = (XConfigureEvent *)&event ;
            vieww = (cev->width + pmag - 1) / pmag ;
            viewh = (cev->height + pmag - 1) / pmag ;
            redraw++ ;
            break ;
case ButtonPress:
            bev = (XButtonEvent *)&event ;
            sprintf(mbuf, "%d %d", bev->x/pmag, bev->y/pmag) ;
            handlecmd('z', mbuf) ;
            break ;
         }
      }
      if (redraw)
         draw(u) ;
   }
}
void clearrect(int x1, int y1, int w, int h) {
   XSetForeground(display, gc, whitePixel) ;
   XFillRectangle(display, window, gc, x1 * pmag, (viewh-y1-h) * pmag,
                                       w * pmag, h * pmag) ;
}
#endif
/*
 *   Now the shared code.  We incorporate the whole hlife program.
 */
#include "hlife.c"
/*
 *   The llxb/llyb hold the offset of the screen as a high-precision
 *   number, one bit per char.  We use this to select the portion of
 *   the macrocell to display at the very high levels (when we are
 *   well above the screen resolution).  We also store all the lower
 *   bits here too, and generate a conventional 32-bit value for the
 *   lower bits as appropriate.  This way the recursive routines
 *   that draw only need use standard 32-bit math when doing clipping
 *   and the like, but we still retain the high-order bits to for
 *   positioning within the universe.
 */
int mag ; /* magnification (0 means 1:1, 4 means 16:1, and so on) */
int llbits ; /* how many bits in ll. */
char *llxb, *llyb ; /* llx and lly in bit-per-char format, least sig first */
/*
 *   Draw a 4x4 area yielding 1x1, 2x2, or 4x4 pixels.
 */
void draw4x4(unsigned short bits, int llx, int lly, int sw) {
   int i, j ;
   if (sw == 1) {
     if (bits)
        drawpixel(-llx, -lly) ;
   } else if (sw == 2) {
      for (i=0; i<sw; i++)
         for (j=0; j<sw; j++)
	    if (bits & (0x33 << (2 + 8 * j - 2 * i)))
               drawpixel(i-llx, j-lly) ;
   } else if (sw == 4) {
      for (i=0; i<sw; i++)
         for (j=0; j<sw; j++)
            if (bits & (1 << (3 + 4 * j - i)))
               drawpixel(i-llx, j-lly) ;
   }
}
/*
 *   Here, llx and lly are coordinates in screen pixels describing
 *   where the lower left pixel of the screen is.  Draw one node.
 *   This is our main recursive routine.
 */
void drawnode(struct node *n, int llx, int lly, int depth) {
   int sw = 1 << (depth - mag + 1) ;
   if (sw >= 64 &&
       (llx + vieww <= 0 || lly + viewh <= 0 || llx >= sw || lly >= sw))
      return ;
   if (n == zeronode(depth)) {
      if (sw >= 64)
         clearrect(-llx, -lly, sw, sw) ;
   } else if (depth > 2 && sw > 1) {
      if (sw == 64) {
         sw >>= 1 ;
         depth-- ;
         drawnode(n->sw, 0, 0, depth) ;
         drawnode(n->se, -32, 0, depth) ;
         drawnode(n->nw, 0, -32, depth) ;
         drawnode(n->ne, -32, -32, depth) ;
         renderbm(-llx, -lly) ;
      } else {
         sw >>= 1 ;
         depth-- ;
         drawnode(n->sw, llx, lly, depth) ;
         drawnode(n->se, llx-sw, lly, depth) ;
         drawnode(n->nw, llx, lly-sw, depth) ;
         drawnode(n->ne, llx-sw, lly-sw, depth) ;
      }
   } else if (sw == 1) {
      drawpixel(-llx, -lly) ;
   } else {
      struct leaf *l = (struct leaf *)n ;
      sw >>= 1 ;
      draw4x4(l->sw, llx, lly, sw) ;
      draw4x4(l->se, llx-sw, lly, sw) ;
      draw4x4(l->nw, llx, lly-sw, sw) ;
      draw4x4(l->ne, llx-sw, lly-sw, sw) ;
   }
}
/*
 *   This is the top-level draw routine that takes the root node.
 *   It maintains four nodes onto which the screen fits and uses the
 *   high bits of llx/lly to project those four nodes as far down
 *   the tree as possible, so we know we can get away with just
 *   32-bit arithmetic in the above recursive routine.  This way
 *   we don't need any high-precision addition or subtraction to
 *   display an image.
 */
void draw(struct node *root) {
   int d = node_depth(root) ;
   int maxd = vieww ;
   int i, llx=0, lly=0 ;
   struct node *z = zeronode(d) ;
   struct node *sw = root, *nw = z, *ne = z, *se = z ;
   redraw = 0 ;
   if (viewh > maxd)
      maxd = viewh ;
   if (llxb[llbits-1])
      llx = -1 ;
   if (llyb[llbits-1])
      lly = -1 ;
/*   Skip down to top of tree. */
   for (i=llbits-2; i>d && i>=mag; i--) { /* go down to d, but not further than mag */
      llx = (llx << 1) + llxb[i] ;
      lly = (lly << 1) + llyb[i] ;
      if (llx > 2*maxd || lly > 2*maxd || llx < -2*maxd || lly < -2*maxd) {
         clearrect(0, 0, vieww, viewh) ;
         return ;
      }
   }
   /*  Find the lowest four we need to examine */
   while (d > 2 && (d - mag > 28 || (1 << (d - mag)) > 2 * maxd)) {
      if (llxb[d] != llxb[llbits-1]) {
         if (llyb[d] != llyb[llbits-1]) {
            ne = ne->sw ;
            nw = nw->se ;
            se = se->nw ;
            sw = sw->ne ;
         } else {
            ne = se->nw ;
            nw = sw->ne ;
            se = se->sw ;
            sw = sw->se ;
         }
      } else {
         if (llyb[d] != llyb[llbits-1]) {
            ne = nw->se ;
            nw = nw->sw ;
            se = sw->ne ;
            sw = sw->nw ;
         } else {
            ne = sw->ne ;
            nw = sw->nw ;
            se = sw->se ;
            sw = sw->sw ;
         }
      }
      llx = (llx << 1) ;
      lly = (lly << 1) ;
      if (llx < 0)
         llx++ ;
      if (lly < 0)
         lly++ ;
      if (llx > 2*maxd || lly > 2*maxd || llx < -2*maxd || lly < -2*maxd) {
         clearrect(0, 0, vieww, viewh) ;
         return ;
      }
      d-- ;
   }
   /*  At this point we know we can use 32-bit arithmetic. */
   for (i=d; i>=mag; i--) {
      llx = (llx << 1) + llxb[i] ;
      lly = (lly << 1) + llyb[i] ;
   }
   /* clear the border *around* the universe if necessary */
   if (d + 1 <= mag) {
      struct node *z = zeronode(d) ;
      if (llx > 0 || lly > 0 || llx+vieww <= 0 || lly + viewh <= 0 ||
          (sw == z && se == z && nw == z && ne == z)) {
         clearrect(0, 0, vieww, viewh) ;
      } else {
         clearrect(0, 1-lly, vieww, viewh-1+lly) ;
         clearrect(0, 0, vieww, -lly) ;
         clearrect(0, -lly, -llx, 1) ;
         clearrect(1-llx, -lly, vieww-1+llx, 1) ;
         drawpixel(-llx, -lly) ;
      }
   } else {
      maxd = 1 << (d - mag + 2) ;
      clearrect(0, maxd-lly, vieww, viewh-maxd+lly) ;
      clearrect(0, 0, vieww, -lly) ;
      clearrect(0, -lly, -llx, maxd) ;
      clearrect(maxd-llx, -lly, vieww-maxd+llx, maxd) ;
      if (maxd < 128) {
         maxd >>= 1 ;
         drawnode(sw, 0, 0, d) ;
         drawnode(se, -maxd, 0, d) ;
         drawnode(nw, 0, -maxd, d) ;
         drawnode(ne, -maxd, -maxd, d) ;
         renderbm(-llx, -lly) ;
      } else {
         maxd >>= 1 ;
         drawnode(sw, llx, lly, d) ;
         drawnode(se, llx-maxd, lly, d) ;
         drawnode(nw, llx, lly-maxd, d) ;
         drawnode(ne, llx-maxd, lly-maxd, d) ;
      }
   }
}
/*
 *   Add to:  move the screen.  Multiply delta by mag.
 */
void lladdto(char *what, int delta) {
   int i = mag ;
   int carry = 0 ;
   while (i < llbits) {
      int nsum = (delta & 1) + what[i] + carry ;
      carry = nsum >> 1 ;
      what[i] = nsum & 1 ;
      delta >>= 1 ;
      i++ ;
   }
}
int picdepth ;              /*   How deep is the pic? */
/*
 *   Our common command routine.
 */
void handlecmd(int cmd, char *arg) {
   double mult ;
   int x, y ;
   redraw++ ;
   if (arg == 0)
      arg = "" ;
   switch(cmd) {
case '-':
      if (pmag == 1) {
         if (mag >= picdepth)
            break ; /* don't do it! */
         mag++ ;
      } else {
         pmag >>= 1 ;
         vieww *= 2 ;
         viewh *= 2 ;
      }
      lladdto(llxb, -vieww/4) ;
      lladdto(llyb, -viewh/4) ;
      break ;
case '+':
      if (mag > 0) {
         lladdto(llxb, vieww/4) ;
         lladdto(llyb, viewh/4) ;
         mag-- ;
      } else if (pmag < 128) {
         lladdto(llxb, vieww/4) ;
         lladdto(llyb, viewh/4) ;
         pmag *= 2 ;
         vieww = (vieww + 1) / 2 ;
         viewh = (viewh + 1) / 2 ;
      }
      break ;
case 'h':
      if (sscanf(arg, "%lg", &mult) != 1)
         mult = 0.4 ;
      lladdto(llxb, (int)(-mult*vieww)) ;
      break ;
case 'j':
      if (sscanf(arg, "%lg", &mult) != 1)
         mult = 0.4 ;
      lladdto(llyb, (int)(mult*viewh)) ;
      break ;
case 'k':
      if (sscanf(arg, "%lg", &mult) != 1)
         mult = 0.4 ;
      lladdto(llyb, (int)(-mult*viewh)) ;
      break ;
case 'l':
      if (sscanf(arg, "%lg", &mult) != 1)
         mult = 0.4 ;
      lladdto(llxb, (int)(mult*vieww)) ;
      break ;
case 'z':
      if (sscanf(arg, "%d %d", &x, &y) != 2) {
         x = vieww / 2 ;
         y = viewh / 2 ;
      } else
	 y = viewh - y - 1 ;
      lladdto(llxb, x/2) ;
      lladdto(llyb, y/2) ;
      if (mag > 0) {
         mag-- ;
      } else {
         pmag *= 2 ;
         vieww = (vieww + 1) / 2 ;
         viewh = (viewh + 1) / 2 ;
      }
      break ;
   }
}
/*
 *   Finally the main ui call called by hlife after it has computed the
 *   generations.
 */
void ui(struct node *root) {
   struct node *n = root ;
   int d = node_depth(n) ;
   openscreen() ;
   picdepth = d + 1 ;
   llbits = 2 * d + 3 ;
   if (llbits < 32)
      llbits = 32 ;
   llxb = (char *)calloc(llbits, sizeof(char)) ;
   llyb = (char *)calloc(llbits, sizeof(char)) ;
   llxb[d] = llyb[d] = 1 ; /* start in center; ll is excess 2^d */
   lladdto(llxb, -vieww/2) ;
   lladdto(llyb, -viewh/2) ;
   memset(bigbuf, -1, 512) ;
   domsgs(n) ;
}
