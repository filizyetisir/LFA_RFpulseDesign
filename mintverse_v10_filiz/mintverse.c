
/* ----------------------------------------------------------------------
	mintverse.c -- Minimum-time VERSE functions.

	Copyright (c)  Brian A. Hargreaves and Charles H. Cunningham

	The main function of interest is mintverse().  Other support
	functions are included here.

	Please use this software freely.  
	Appropriate credit or acknowledgement 
	is certainly appreciated.
   ---------------------------------------------------------------------- */


/* ========================================================================
	This file is maintained in CVS.  Log entries follow.
   ---------------------------------------------------------------------- 

	$Log: mintverse.c,v $
	Revision 1.14  2005/02/22 23:53:01  brian
	Added and tested functionality for complex B (rf) pulses.
	
	Revision 1.13  2005/02/22 23:19:42  brian
	Added complex B support -> yet to test
	
	Revision 1.12  2004/10/05 19:46:29  brian
	Changes to the mintverse() function to allow a maximum-energy
	constraint to be included, then iterate old scheme, changing amplitude
	limit to converge on the energy limit.
	
	Changed mex and .m functions accordingly.
	Added an energy-constrained example to exampleverse.m
	
	Revision 1.11  2004/10/01 06:41:45  brian
	Added verse() function and a bunch
	of error checking.
	There was a conflict at this stage that accidentally
	got committed - big oops!.
	
	Revision 1.10  2004/09/24 16:18:00  brian
	-Fixed compressmax() to ignore zero-gradient samples
	-Fixed verse() to properly scale waveforms.
	
	Revision 1.9  2004/09/22 20:01:08  brian
	minor edits
	
	Revision 1.8  2004/09/22 19:58:12  brian
	First revision on web.  Added error checks with warnings.
	

   ---------------------------------------------------------------------- 
	End of CVS log entries.
   ======================================================================== */

/*
#define DEBUG
*/


#include <math.h>
#include <stdio.h>
#include <string.h>


double dabs(double num)
{
double absnum;
if (num > 0)	
	absnum=num;
else
	absnum=-num;

return(absnum);
}


double sumarray(double *f, long n)
/* ------------------------------------------------------------
	Function returns the sum of the first n points in
	the array f.
   ------------------------------------------------------------ */
{
long count;
double asum=0;

for (count = 0; count < n; count++)
	asum += *f++;

return asum;
}



int isnonnegative(double *f, long n)
/* ------------------------------------------------------------
	Function checks the first n points of the array f.
	If they are all non-negative, it returns 1, otherwise
	returns 0.
   ------------------------------------------------------------ */
{
int isnonneg=1;
long count;

for (count = 0; count < n; count++)
	if (*f++ < 0)
		isnonneg=0;

return isnonneg;
}


int ispositive(double *f, long n)
/* ------------------------------------------------------------
	Function checks the first n points of the array f.
	If they are all positive, it returns 1, otherwise
	returns 0.

   ------------------------------------------------------------ */
{
int ispos=1;
long count;

for (count = 0; count < n; count++)
	if (*f++ <= 0)
		ispos=0;		

return ispos;
}




long uniformresamplesize(double *dt, long n, double dtr)
/* ------------------------------------------------------------
	Function returns the number of points that will
	are required to uniformly resample an array with
	block-widths in *dt, with a uniform time block width dtr.

	INPUT:
		dt -	n-point array of block widths for input.
		n -	number of points in input.
		dtr -	block width for each output block.

	RETURN VALUE:
		number of output points/blocks that 
		uniformresample() will generate.

   ------------------------------------------------------------ */
{
long nfr;
long fcount;
double tend;


  /* Check some input parameters. */

if (dtr <= 0)
	printf("Warning:  uniformresamplesize() dtr not positive-valued.\n");
if (n <= 0)
	printf("Warning:  uniformresamplesize() n not positive-valued. \n");
if (!(isnonnegative(dt,n)))
	printf("Warning:  uniformresamplesize() dt must be non-negative. \n");


  /* Find number of output points. */

tend = 0;
for (fcount = 0; fcount < n; fcount++) 
	{
	tend+= dt[fcount];
	}
nfr = (long)( (tend+0.99999*dtr) / dtr);	/* Round up end time. */

return nfr;
}



int resample(double *f, double *dt, long n, double *fr, double *dtr,
				long nr)
/* ------------------------------------------------------------
	Function resamples f onto a different time step.
	The "input function" is linearly interpolated between 
	the centers of blocks, and the resampled function is
	generated by accumulating the input area across the
	output blocks.

	The "output" function is essentially blocks with width
	dtr, and have area equal to the area of the input 
	function over that block.

	INPUT:
		f - 	n-point array of heights of block-centers
				in input function to resample.
		dt -	n-point array of block widths for input.
		n -	number of points in input.
		dtr -	nr-point array of block widths for each output block.
		nr = 	number of blocks in output.

	OUTPUT:
		fr - 	array of heights of output blocks.  
			 	**This must be adequate size!.

	NOTES:
		Again, note that the array fr must be of
		adequate size.  Use uniformresamplesize() to determine
		this size.		
   ------------------------------------------------------------ */


{
double ts;	/* Start time of input interval. */
double tf;	/* End time of input interval. */
double fs;	/* Start level of input interval. */
double ff;	/* End level of input interval. */

double tr;	/* Time of next block end in output. */
double ftr;	/* Interpolated level at tr, when needed. */
double farea;	/* Cumulative area over output block. */

long frcount=0; 	/* Output point number. */
long fcount=0;		/* Input point number. */


  /* Check some input parameters. */

if (n <= 0)
	printf("Warning:  resample() n must be greater than zero.\n");
if (nr <= 0)
	printf("Warning:  resample() nr must be greater than zero.\n");
if (!(isnonnegative(dt,n)))
	printf("Warning:  resample() dt values must be non-negative. \n");
if (!(isnonnegative(dtr,nr)))
	printf("Warning:  resample() dtr values must be non-negative. \n");
if (dabs(sumarray(dt,n) - sumarray(dtr,nr)) > .001 * sumarray(dt,n))
	printf("Warning:  resample(): dt and dtr differ by more than 0.1%%\n");


tr = *dtr;		/* First output interval end. */


/* Setup what is known for first input interval. */

ts = 0;			/* Start at 0. */
fs = *f++;		/* First input block height. */
tf = dt[0]/2.0;		/* Finish at the center of the time-block */
ff = fs;		/* Flat first interval. */
farea = 0;		/* Accumulated area of output. */



while (frcount < nr)	/* Continue until all output points assigned. */
	{
	if (tf < tr)	/* Input interval ends before output interval. */	
			/* Simply add area, advance input pointers */
		{
		farea += (fs+ff)*(tf-ts)/2.0;	/* Trapezoidal area. */
		ts = tf;			/* Advance start time.*/
		fs = ff;			/* Update start height. */

		while (tf <= ts)	/* Ignore "zero-length" intervals */
		    {
		    if (fcount < n)	/* Never increment past dt[n-1] */
			{
			fcount++;		   /* Adv. input block count.*/
			tf += (*dt++ + *dt)/2.0;   /* Block end midway 
							between centers */
			ff = *f++;		   /* Block height at center.*/
			}
		    else		/* Make interval flat, ends well past
					what is needed. */
			{
			tf += (tf+tr);	/* Nice and big! */
			}
		    if (tf <= ts)	
			printf("Warning:  Zero-length block in resample input.\n");
			
		    }
		}

	else 	/* Output interval ends before input interval. */
		/* Interpolate to output-end, add area, and advance 
			 output interval. */

		{
		ftr = fs + (ff-fs)/(tf-ts)*(tr-ts);  /* Interp'd level at tr*/
		farea += (fs+ftr) * (tr-ts)/2.0;	/* Add area. */	

		ts = tr;		/* Input interval now starts at tr. */
		fs = ftr;		/* New fs. */
					/* Leave ff and tf the same... */

	
		while ((tr <= ts) && (frcount < nr))	
				/* Go to next non-zero-length interval*/
			{	
			*fr++ = farea/(*dtr++);	   /* Assign output value. */
			tr += *dtr;	/* Increment output interval-end-time*/
			farea = 0;	/* Reset accumulated area. */
			frcount++;	/* Increment output counter. */

		    	if (tr <= ts)
			    printf("Warning:  Zero-length block in resample output.\n");
			}
		}

	}

if ((fcount > n) || (frcount > nr))
	{
	printf("Warning:  resample() over-ran input or output array.\n");
	printf("  Blocks - input: %d (%d).  output: %d (%d).\n",
			fcount,n,frcount,nr);
	}
}



int uniformresample(double *f, double *dt, long n, double *fr, double dtr)
/* ------------------------------------------------------------
	Function resamples f onto a uniformly-sampled time-step.
	The "input function" is linearly interpolated between 
	the centers of blocks.

	The "output" function is essentially blocks that end
	at multiples of dtr, and have area equal to the area
	of the input function over that block.

	INPUT:
		f - 	n-point array of heights of block-centers
				in input function to resample.
		dt -	n-point array of block widths for input.
		n -	number of points in input.
		dtr -	block width for each output block.

	OUTPUT:
		fr - 	array of heights of output blocks.  
			 	**This must be adequate size!.

	NOTES:
		Again, note that the array fr must be of
		adequate size.  Use uniformresamplesize() to determine
		this size.		
   ------------------------------------------------------------ */

{
double ts;	/* Start time of input interval. */
double tf;	/* End time of input interval. */
double fs;	/* Start level of input interval. */
double ff;	/* End level of input interval. */

double tr;	/* Time in output. */
double ftr;	/* Interpolated level at tr, when needed. */
double farea;	/* Cumulative area over output block. */

long nfr;		/* Number of resampled points. */
long frcount=0; 	/* Output point number. */
long fcount=0;		/* Input point number. */


  /* Check some input parameters. */

if (*dt <= 0)	/* just check 1st sample */
	printf("Warning:  uniformresample() dt values must be positive.\n");
if (!(isnonnegative(dt,n)))
	printf("Warning:  uniformresample() dt values must be non-negative.\n");
if (dtr <= 0)
	printf("Warning:  uniformresample() dtr must be greater than zero.\n");
if (n <= 0)
	printf("Warning:  uniformresample() n must be greater than zero.\n");



nfr = uniformresamplesize(dt, n, dtr); 		/* Number of output points. */
tr = dtr;					/* First output interval end. */


/* Setup what is known for first input interval. */

ts = 0;			/* Start at 0. */
fs = *f++;		/* First input block height. */
tf = dt[0]/2.0;		/* Finish at the center of the time-block */
ff = fs;		/* Flat first interval. */
farea = 0;		/* Accumulated area of output. */


while (frcount < nfr)	/* Continue until all output points assigned. */
	{
	if (tf < tr)	/* Input interval ends before output interval. */	
			/* Simply add area, advance input pointers */
		{
		farea += (fs+ff)*(tf-ts)/2.0;	/* Trapezoidal area. */
		ts = tf;			/* Advance start time.*/
		fs = ff;			/* Update start height. */

		while (tf <= ts)
		    {
		    if (fcount < n)	 /* Never increment past dt[n-1] */
			{
			fcount++;		   /* Adv. input block count.*/
			tf += (*dt++ + *dt)/2.0;   /* Block end midway 
							between centers */
			ff = *f++;		   /* Block height at center.*/
			}
		    else		/* Make interval flat, ends well past
					what is needed. */
			{
			tf += *dt/2.0 + 5*dtr;	   /* 5 is "big enough." */
			}
		    if (tf <= ts)
			{
			printf("Warning:  Zero-length sample at %d of %d.\n",
					fcount,n);
			}
		    }
		}

	else 	/* Output interval ends before input interval. */
		/* Interpolate to output-end, add area, and advance 
			 output interval. */

		{
		ftr = fs + (ff-fs)/(tf-ts)*(tr-ts);  /* Interp'd level at tr*/
		farea += (fs+ftr) * (tr-ts)/2.0;	/* Add area. */	

		ts = tr;		/* Input interval now starts at tr. */
		fs = ftr;		/* New fs. */
					/* Leave ff and tf the same... */

		tr += dtr;		/* Increment output interval-end-time*/
		*fr++ = farea/dtr;	/* Assign output value. */
		farea = 0;		/* Reset accumulated area. */
		frcount++;		/* Increment output counter. */
		}

	}
}



double stretchslew(double gh, double gl, double dth, double dtl, double smax)
/* ------------------------------------------------------------
	Function finds the appropriate stretch in time
	of a gradient block of height gh and width dth that
	is adjacent to a gradient block of height gl and width dtl,
	so that the slew rate smax exactly met, while the area
	of the block is preserved.


	INPUT:
		gh = higher gradient (assumed non-negative).
		gl = lower gradient (assumed non-negative).
		dth= width of higher gradient segment.	
		dtl= width of lower gradient segment.	
		smax = maximum slew rate.

	Units:  Arbitrary, but smax must be in gradient units per
		time unit.

	OUTPUT:
		(return value) is stretch factor to stretch
		the higher gradient block width in TIME!

	NOTES:
		Let k be the stretch factor.  The desired k is the
		solution of the quadratic equation:

			gh/k-gl = smax*(k*dth+dtl)/2

		or      (smax*dth) k^2  + (smax*dtl+2gl) k  - 2gh = 0    

		The function tries to solve this efficiently.	
   ---------------------------------------------------------------------- */

{
double stretch;
double r1, r2;
double a,b,c,d,f,g;


	/* Do some error checking - hopefully this will save 
		some problems as this will inevitably be called
		for negative gradients!				*/

if (gh < gl)
	printf("Warning:  stretchslew() - gh and gl are reversed. \n");
if ((gh < 0) || (gl < 0))
	printf("Warning:  stretchslew() - gradient is negative. \n");
if ((dth <= 0) || (dtl <= 0))
	printf("Warning:  stretchslew() - block width is not positive. \n");
if (smax <= 0)
	printf("Warning:  stretchslew() - smax is not positive. \n");




  /* Define a,b,c for quadratic equation a*r^2+b*r+c = 0 */

a = smax * dth;
b = smax * dtl + 2*gl;
c = -2*gh;

  /* Now express roots as r1 = b/2a +- sqrt(b*b-4*a*c)/2a. 	*/
  /*							 	*/
  /* This is done as r1 = d+g, r2=d-g,  d and g intermediate 	*/
  /*	variables to check for errors and minimize operations. 	*/	

d = b/2/a;
f = d*d - c/a;		/* Inside of square root sign for QDF */


if (f < 0)		/* Check - should be real! */
	{
	printf("Error - quadratic roots are complex. \n");
	stretch = -1;	
	}
else
	{
	g = sqrt( f );		
	stretch = -d+g;		/* Greater root. */

	/* Roots will be positive and negative, representing a
	   stretch that keeps dth positive, and the non-realistic
	   stretch that makes dth negative.  Thus we only care
	   about the greater root.				*/

	}

return(stretch);
}
	
		


int adjustslew(double *br,double *bi, double *g,double *dt,long n, 
			double smax, long loc, long recursecount)
/* ------------------------------------------------------------
	Fix slew rate violations, if they exist, by stretching
	the block width so that the slew rate is exactly met.

	INPUT:
		br = real part of RF waveform (arbitrary units.)
		bi = imag part of RF waveform (arbitrary units.)
		g = gradient waveform (arbitrary units.)
		dt= time-delta waveform (arbitrary units.)
		n = number of points in b, g and dt.
		smax = max slew rate in gradient units per time unit.
		loc = current point to adjust (1 <= loc < n).

	OUTPUT:
		g,br, bi,dt are all adjusted.
	
	RETURN: 1 if successful, 0 if not.

	NOTE:	recursive!

   ------------------------------------------------------------ */

{
double stretch;
double slew;
int success;


  /* Do some error checking. */

if (loc <= 0)
	printf("Warning:  adjustslew():  loc must be positive. \n");
if (n <= loc)
	printf("Warning:  adjustslew():  n must be greater than loc. \n");
if (smax <= 0)
	printf("Warning:  adjustslew():  smax must be greater than zero. \n");
if (!(ispositive(dt,n)))
	printf("Warning:  adjustslew() dt values must be positive.\n");



if (g[loc-1] > g[loc])	 	/* Downward slope */
	{
	stretch = stretchslew(g[loc-1],g[loc],dt[loc-1],dt[loc],smax);
	g[loc-1]  = g[loc-1]/stretch;
	br[loc-1]  = br[loc-1]/stretch;
	bi[loc-1]  = bi[loc-1]/stretch;
	dt[loc-1] = dt[loc-1]*stretch;

	if ((stretch > 0) && (loc > 1) )
		{	
		slew = (g[loc-1]-g[loc-2]) / (dt[loc-1]+dt[loc-2]) / 0.5;
		if (dabs(slew) > smax)
			{
			success =adjustslew(br,bi,g,dt,n,smax,loc-1,
						recursecount+1);
			stretch = (double)success+0.5;
			}
		}
	}
else		/* Upward slope. Just reduce height of higher point. */
	{
	stretch = stretchslew(g[loc],g[loc-1],dt[loc],dt[loc-1],smax);
	g[loc]  = g[loc]/stretch;
	br[loc]  = br[loc]/stretch;
	bi[loc]  = bi[loc]/stretch;
	dt[loc] = dt[loc]*stretch;
	}


return (stretch > 0);	/* negative stretch indicates error. */
}




int slewcheck(double *br, double *bi, double *g, double *dt, long n, 
			double gmax, double smax)
/* ------------------------------------------------------------
	Adjust pulse so that slew rates are not violated

	INPUT:
		br = real part of RF waveform (arbitrary units.)
		bi = imag part RF waveform (arbitrary units.)
		g = gradient waveform (arbitrary units.)
		dt= time-delta waveform (arbitrary units.)
		n = number of points in b, g and dt.
		gmax = max gradient value (units of g).
		smax = max slew rate in gradient units per time unit.

	OUTPUT:
		br,bi,g,dt are updted.

	RETURN VALUE:
		1 if successful, 0 otherwise.


   ------------------------------------------------------------ */
{
int count;
double slew;
int success = 1;



  /* Do some error checking. */

if (n <= 0)
	printf("Warning:  slewcheck():  n must be greater than zero. \n");
if (!(ispositive(dt,n)))
	printf("Warning:  slewcheck() dt values must be positive.\n");
if (!(isnonnegative(g,n)))
	printf("Warning:  slewcheck() g values must be non-negative.\n");
if (gmax <= 0)
	printf("Warning:  slewcheck():  gmax must be greater than zero. \n");
if (smax <= 0)
	printf("Warning:  slewcheck():  smax must be greater than zero. \n");



for (count = 1; count < n; count++)
	{
	slew = (g[count]-g[count-1]) / (dt[count]+dt[count-1]) / 0.5;
	if (dabs(slew) > smax)
		success = adjustslew(br,bi,g,dt,n,smax,count,0);
	if (success != 1)
		count = n;	/* break out of loop */
	}

}


double calcenergy(double *br, double *bi, double *dt, long n)
/* ------------------------------------------------------------
	Calculate the energy in a waveform b, with time steps
	in dt and n points.  This is just sum(b_i^2*dt_i).

	INPUT:
		br = real part of waveform, nominally RF (arb. units)
		bi = imag part of waveform, nominally RF (arb. units)
		dt= Time-widths of sample blocks in b. (arb units.)
		n = Number of points in b and dt.

	OUTPUT:
		(return value):  Energy in waveform.

  ------------------------------------------------------------ */
{
long count;
double benergy = 0.0;

if (n < 0)
	printf("Warning:  calcenergy() - length must be non-negative.\n");	

for (count = 0; count < n; count++)
	{
	benergy += ( (*br)*(*br)+(*bi)*(*bi) ) * (*dt);
	br++;
	bi++;
	dt++;
	}


return (benergy);
}
	

int compressmax(double *br, double *bi, double *g, double *dt,
			long n, double bmax, double gmax)
/* ------------------------------------------------------------
	Stretch or shrink the n "infinitesimal" block widths (dt)
	such that at each time step either |b| or g reaches its
	maximum.  Does not change the length of |b|, g or dt.  
	br, bi and bmax must be the same units, while g and gmax must
	also be the same units.

	INPUT:
		br = real part of RF waveform (arbitrary units.)
		bi = imag part of RF waveform (arbitrary units.)
		g = gradient waveform (arbitrary units.)
		dt= time-delta waveform (arbitrary units.)
		n = number of points in b, g and dt.
		bmax = max RF value (units of b).
		gmax = max gradient value (units of g).
	
	OUTPUT:
		br,bi,g,dt are all updated.

	NOTES:

   ------------------------------------------------------------ */

{
long count;
double bfrac, gfrac;



  /* Do some error checking. */

if (n <= 0)
	printf("Warning:  compressmax(): n must be greater than zero. \n");
if (!(ispositive(dt,n)))
	printf("Warning:  compressmax(): dt values must be positive.\n");
if (!(isnonnegative(g,n)))
	printf("Warning:  compressmax(): g values must be non-negative.\n");
if (gmax <= 0)
	printf("Warning:  compressmax(): gmax must be greater than zero. \n");
if (bmax <= 0)
	printf("Warning:  compressmax(): bmax must be greater than zero. \n");




for (count = 0; count < n; count++)
	{

	bfrac = dabs( sqrt((*br)*(*br)+(*bi)*(*bi)) / bmax);
	gfrac = dabs(*g / gmax);

	if ((bfrac > gfrac) && (gfrac > 0))	/* RF will be max'd */
		{
		*br = *br / bfrac;
		*bi = *bi / bfrac;
		*g = *g / bfrac;
		*dt = *dt * bfrac;
		}
	else if (gfrac > 0)			/* Gradient will be max'd */
		{
		*br = *br / gfrac;
		*bi = *bi / gfrac;
		*g = *g / gfrac;
		*dt = *dt * gfrac;
		}
	else
		printf("compressmax():  Zero-gradient (%d), ignored.\n",count);

	br++;
	bi++;
	g++;
	dt++;
	}
}



int verse(double *br, double *bi, double *g, long n, double *gv, 
		double *bvr, double *bvi, long nv)
/* ------------------------------------------------------------
	Convert an RF/gradient pair to a VERSE
	equivalent pair, specified by the given gv/nv.

	INPUT:
		br = real part of n-point RF waveform (arbitrary units.)
		bi = imaginary part of n-point RF waveform (arbitrary units.)
		g = n-point gradient waveform (arbitrary units.)
		n = number of points in b, g and dt.
		gv= nv-point VERSE gradient waveform (same units as g).
		nv= number of points in VERSE waveforms.

	OUTPUT:
		bvr = real part of nv-point VERSE RF waveform (same units as b).
		bvi = imag part of nv-point VERSE RF waveform (same units as b).

	NOTES:
		For the case where waveforms are discrete
		blocks, the gradient waveform can be viewed
		as the time step (block-width).  Thus VERSE is
		nothing more than resampling the initial RF/gradient
		using the VERSE gradient as the time-step.

		No error checking is done.  g and gv should be
		positive, and the sum of g and of gv should be equal.
 
   ------------------------------------------------------------ */

{
double *brwork;
double *biwork;
int returnval;
int count;

if (n <= 0)	
	printf("Warning:  verse(): n must be greater than zero.\n");
if (!(isnonnegative(g,n)))
	printf("Warning:  verse(): g values must be non-negative.\n");
if (nv <= 0)	
	printf("Warning:  verse(): n must be greater than zero.\n");
if (!(isnonnegative(gv,nv)))
	printf("Warning:  verse(): gv values must be non-negative.\n");
if (dabs(sumarray(g,n) - sumarray(gv,nv)) > .001)
	printf("Warning:  verse(): g and gv sums differ by more than .001.\n");


printf("verse()   %d point input -> %d point output. \n",n,nv);

brwork = (double *) malloc(n*sizeof(double));
biwork = (double *) malloc(n*sizeof(double));
/*
bcopy (br,brwork,n*sizeof(double));
bcopy (bi,biwork,n*sizeof(double));
 */

/* Bastien */
memcpy (brwork,br,n*sizeof(double));
memcpy (biwork,bi,n*sizeof(double));


for (count = 0; count < n; count++)
	{
	if (g[count]>0)
		{
		brwork[count] /= g[count];
		biwork[count] /= g[count];
		}
	else if ((g[count]==0) && ((brwork[count]!=0) || (biwork[count]!=0)) )
		{
		printf("Warning:  verse(): Zero-gradient point (%d)",count);
		printf("with non-zero B1.\n");
		}
	}
		

returnval = resample(brwork, g, n, bvr, gv, nv);
returnval = resample(biwork, g, n, bvi, gv, nv);

for (count = 0; count < nv; count++)
	{
	bvr[count] *= gv[count];
	bvi[count] *= gv[count];
	}

free(brwork);
free(biwork);
return (returnval);

}


int mintverse(double *breal, double *bimag, double *g, double *dt, long n, 
			double bmax, double gmax, double smax, double tout, 
			double emax, long *nout, double **bvr, 
			double **bvi, double **gv)

/* ------------------------------------------------------------
	Convert an RF/gradient pair to the minimum-time VERSE
	equivalent pair.  This is the "Do-all" function that:
	
		1) Compresses B1/gradient so one is always maximized.
		2) Removes gradient slew rate violations.
		3) Resamples waveforms uniformly.
		
	INPUT:
		breal = n-point real part of RF waveform (arbitrary units.)
		bimag = n-point imaginary part of RF waveform (a.u.)
		g = n-point gradient waveform (arbitrary units.)
		dt= n-point time-delta waveform (arbitrary units.)
		n = number of points in b, g and dt.
		bmax = max RF value (units of b).
		gmax = max gradient value (units of g).
		smax = max gradient slew rate (units of g per unit dt)
		tout = output sampling period (units of dt)
		emax = max RF energy (units of b*b*dt)  (-1 to not constrain)
	
	OUTPUT:
		nout = number of output samples.
		bvr = address of output real B1 pulse array (allocated here)
		bvi = address of output imaginary B1 pulse array (alloc. here)
		gv = address of output gradient array (allocated here)

	NOTES:
		-Outputs are allocated with malloc.
		-No error checking is done:  breal,bimag,g,dt should 
			be length n and g,dt should be non-negative

		If the gradient and RF are both zero at the ends,
		then the output will also start and end at zero.  This
		is often desirable.  
		
   ------------------------------------------------------------ */


{
double *brwork;		/* Working array for real part of B */
double *biwork;		/* Working array for imaginary part of B */
double *gwork;		/* Working array for gradient */
double *dtwork;		/* Working array for dt */
long count;
double *tv;

double emaxratio=0.98;	/* 	minumum energy/emax for stopping iteration */
double bmaxh;		/* 	bmax that gives energy higher than max. */
double bmaxl;		/*	bmax that gives energy lower than max. */
double bmaxc;		/*	current bmax. */
double benergy;		/*	current energy */
short doneiter=0;	/*	0 to start.  1 when done iterating. */
long niter = 0;		/* 	Number of iterations 	*/
long maxiter = 1000;	/* 	Max number of iterations */

if (n <= 0)
	printf("Warning:  mintverse(): n must be greater than zero. \n");
if (!(ispositive(dt,n)))
	printf("Warning:  mintverse(): dt values must be positive.\n");
if (!(isnonnegative(g,n)))
	printf("Warning:  mintverse(): g values must be non-negative.\n");
if (bmax <= 0)
	printf("Warning:  mintverse(): bmax must be greater than zero. \n");
if (gmax <= 0)
	printf("Warning:  mintverse(): gmax must be greater than zero. \n");
if (smax <= 0)
	printf("Warning:  mintverse(): smax must be greater than zero. \n");





/*	Allocate working arrays and copy inputs to them. */

brwork = (double *) malloc(n*sizeof(double));
biwork = (double *) malloc(n*sizeof(double));
gwork = (double *) malloc(n*sizeof(double));
dtwork = (double *) malloc(n*sizeof(double));

bmaxc = bmax;	 
bmaxh = bmax;	 /* Maximum bmax for energy OR b-amplitude */
bmaxl = 0.0;	


while (doneiter == 0)
	{
	niter++;

    /*
	bcopy(breal,brwork,n*sizeof(double));
	bcopy(bimag,biwork,n*sizeof(double));
	bcopy(g,gwork,n*sizeof(double));
	bcopy(dt,dtwork,n*sizeof(double));
     */
    
    /* Bastien */
	memcpy(brwork,breal,n*sizeof(double));
	memcpy(biwork,bimag,n*sizeof(double));
	memcpy(gwork,g,n*sizeof(double));
	memcpy(dtwork,dt,n*sizeof(double));
	
	/*      Call C-functions Here */

	#ifdef DEBUG
		printf("Calling compressmax(). \n");
	#endif
	compressmax(brwork,biwork,gwork,dtwork,n,bmaxc,gmax);
	
	#ifdef DEBUG
		printf("Calling slewcheck(). \n");
	#endif
	slewcheck(brwork,biwork,gwork,dtwork,n,gmax,smax);

	if (emax > 0)	
		{
		benergy = calcenergy(brwork, biwork, dtwork, n);
		#ifdef DEBUG
			if (!(isnonnegative(dtwork,n)))
				printf("oops! dtwork is not nonnegative. \n");
			printf("Iteration %d:  energy = %g, max is %g \n",
				niter,benergy,emax);
		#endif

		if (benergy > emax)
			{
			bmaxh = bmaxc;
			bmaxc = bmaxl + (bmaxh - bmaxl)/2.0;
			}
		else
			{
			bmaxl = bmaxc;
			if (benergy/emax > emaxratio) 
				doneiter = 1;
			if (niter==1)		/* Energy constraint inactive*/
				doneiter = 1;	
			bmaxc = bmaxl + (bmaxh - bmaxl)/2.0;
			}

		if (bmaxl >= bmaxh)
			{
			printf("Warning:  mintverse() low-E limit > high-E limit.\n");
			printf("  Exiting after %d iterations.\n",niter);
			doneiter = 1;
			}
		if (niter >= maxiter)
			{
			printf("Warning:  mintverse() iteration limit reacheed.\n");
			printf("  Exiting after %d iterations.\n",niter);
			doneiter = 1;
			}
			
		#ifdef DEBUG
			printf("Iteration bmax:  high/low/next = %g,%g,%g\n",
				bmaxh,bmaxl,bmaxc);
		#endif
		}
	else
		doneiter = 1;	/* no energy constraint, so done! */

	}	
	
/*      Calculate length of, and allocate output vectors. 	*/
	
*nout = uniformresamplesize(dtwork,n,tout);
*bvr = (double *) malloc(*nout * sizeof(double));
*bvi = (double *) malloc(*nout * sizeof(double));
*gv = (double *) malloc(*nout * sizeof(double));


tv = (double *) malloc(*nout * sizeof(double));
for (count = 0; count < *nout; count++)
	tv[count] = tout;


/*	Resample with uniform sample step. 	*/

#ifdef DEBUG
	printf("Calling resample(). \n");
#endif

/* 	Either resample() or uniformresample() can be used here... */

/*
resample(brwork,dtwork,n, *bvr, tv, *nout);
resample(biwork,dtwork,n, *bvi, tv, *nout);
resample(gwork,dtwork,n, *gv, tv, *nout);
*/

uniformresample(brwork,dtwork,n, *bvr, tout);
uniformresample(biwork,dtwork,n, *bvi, tout);
uniformresample(gwork,dtwork,n, *gv, tout);


#ifdef DEBUG
	printf("mintverse() freeing memory. \n");
#endif

free(tv);

free(brwork);
free(biwork);
free(gwork);
free(dtwork);
}




