#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"


// For prototype-making (makeprots) if using 256x256 raw images as inputs for the prototypes:
//#define IMAGEXSIZE 256
//#define IMAGEYSIZE 256

// For actual images in the experiment, which are all 320x320:
#define IMAGEXSIZE 320
#define IMAGEYSIZE 320

#define HORIZ 0.0
#define VERTIC (M_PI / 2.0)

//#define SIGMA .5  
#define SIGMA .1  
#define SIGMAS1 0.0  

#define S1 1
#define C1 2
#define S2 3
#define C2 4
#define S3 5
#define C3 6
#define S2b_5 7
#define S2b_7 8
#define S2b_9 9
#define S2b_11 10
#define S2b_13 11

#define MAXGROUP 12

#define MAXSCALES 16


#define DENSE 121
#define RADIUS 122
#define HALFRADIUS 123
#define EVERYOTHER 124
#define POGGIO 125
#define EVERYOTHER4 126


#define MERGE 505
#define BLUR 506
#define KEEP 507

#define HANDLEC1SCALES KEEP // Merge = Poggio-like - merge contiguous scales from S1 to C1


#define EPSS1 HALFRADIUS
#define EPSS2 RADIUS

#define EPSC1 EVERYOTHER4



#define NBSCALES_S1  12


int NBSCALES_S2;
int NBSCALES_C1;
int NBSCALES_C2;






#define MAXPROTS 600
#define NBS2PROTS 600 // Number of S2/S2b prototypes / patches / templates
//#define MAXPROTS 1200
//#define NBS2PROTS 1200 // Number of S2/S2b prototypes / patches / templates
#define NBS3PROTS 2 // Number of S3 prototypes / patches / templates

//#define MAXAFF (MAXPROTS * 3 * 3) // IF we have S3/S4 layers, then MAXAFF is  MAXPROTS possible source layers at each position, 3*3 RF (number of afferents for S3 layer, hopefully the maximum)
#define MAXAFF (4 * 13 * 13) // BUT if we only have S2 layers, then MAXAFF is 4 orientations times 13x13 (maximum RF size)

#define MAXLAYERS 90000
#define MAXNEUR 150000000


#define NBCONNS2PROT  4 * 3 * 3 //  4 possible orientations at each position, 3x3 RF
#define NBCONNS3PROT  NBS2PROTS * 3 * 3 //  NBS2PROTS possible source layers at each position, 3x3 RF



int eps[MAXLAYERS]; // Downsampling factor in terms of cells from the previous layer
int epsreal[MAXLAYERS];  // Same, in terms of actual pixels
int rad[MAXLAYERS];
int offsetreal[MAXLAYERS]; // Offset between the of the image and the layer's first cell, in pixels

int xsize[MAXLAYERS];  // The x- and y-size of the layers, in cells (not in pixels!)
int ysize[MAXLAYERS];
int size[MAXLAYERS];

int start[MAXLAYERS];  // The index of the layer's first neuron, among all neurons

int nbaff[MAXLAYERS];  // only meaningful for S layers
int sourcelayer[MAXLAYERS][10];   // only used for C layers (S layers have many, unpredictable source layers)
int nbsourcelayers[MAXLAYERS];   // only used for C layers (S layers have many, unpredictable source layers)
int nbneurl[MAXLAYERS]; 
int scale[MAXLAYERS]; // The scale of this layer
int sourcescale[MAXLAYERS]; // The scale of the layers (or of the smallest layer) from which this layer takes connections
int group[MAXLAYERS]; // The group of this layer
int sourcegroup[MAXLAYERS]; // The group of the layers from which this layer takes connections
int numprot[MAXLAYERS]; // The prototype associated with this layer 
int diamRF_S1[MAXLAYERS]; // The diameter of the RF for each S1 layer
int radius[MAXLAYERS]; // The diameter of the RF for each S1 layer


//int realx[MAXNEUR];  // The x and y position of each cell, in pixels
//int realy[MAXNEUR];

int *realx ;
int *realy ;


//IMPLEMENT THE SCALE ARRAY


//int layerneur[MAXNEUR];
int *layerneur;


	int numlayer[MAXGROUP][MAXSCALES][MAXPROTS];
	
	int nblayers = 0;
	int nblayersgr[MAXGROUP];
	int first[MAXGROUP];
	int nbprots[MAXGROUP]; 




int nbneur = 0;

//float v[MAXNEUR];
//float *v = (float*) calloc(MAXNEUR, sizeof(float));
float *v;



float pic[IMAGEXSIZE*IMAGEYSIZE];



// These ones are where the prototypes themselves are actually stored:
float StoredProts [MAXGROUP][MAXPROTS][4 * MAXAFF]; // 4 floats per connection
int StoredProtsInt [MAXGROUP][MAXPROTS][4 * MAXAFF]; // 4 floats per connection
float *weightsS1[MAXSCALES][4]; 


float *ProtOfLayer[MAXLAYERS]; // Pointer to the weight set associated to each particular layer; only used for V1, should be merged with the rest... 



void mydie(const char *s)
{
	printf("%s", s); printf("\nBye!");
	fflush(stdout);
	exit(-1);
}


// Build the structure of layer nl by modelling it on the "source" layer nl: 
// Each nl neuron is located on top of a sl neuron (i.e. has the same global coordinates).
// The number of neurons in nl, and their real world positions, are calculated
// based on the structure of sl, the diameter of the nl RF (which determines
// the size of the edges of sl which are not covered by nl) and the eps
// sampling factor (there is one nl neuron evey eps sl neuron).

void buildS1Layer(int nl, int sl, int diam, int myeps)
{

	int myrad = (int) floor (diam / 2);
	radius[nl] = myrad; eps[nl]= myeps;
	nbaff[nl] = diam*diam;
//	if (nbaff[nl] >= MAXAFF) mydie("Too many afferents in S1 layer!");
	diamRF_S1[nl] = diam;

	xsize[nl] = ceil ((double)(xsize[sl] - 2 * myrad) / (double)myeps);
	ysize[nl] = ceil ((double)(ysize[sl] - 2 * myrad) / (double)myeps);
	size[nl] = xsize[nl] * ysize[nl];
	nbneurl[nl] = 0;

	offsetreal[nl] = offsetreal[sl] + myrad * epsreal[sl];
	epsreal[nl] = myeps * epsreal[sl];

	start[nl] = nbneur;
	for (int j=0; j < ysize[nl] ; j++)
		for (int i=0; i < xsize[nl] ; i++)
		{
			int targetx = i * myeps + myrad;
			int targety = j * myeps + myrad;
			int targetn = start[sl] + targetx + targety * xsize[sl];
			realx[nbneur] = realx[targetn];
			realy[nbneur] = realy[targetn];
			layerneur[nbneur] = nl;
			nbneur ++; nbneurl[nl] ++;
			if (nbneur > MAXNEUR) { printf("L %d (group %d), neur %d (total %d)\n", nl, group[nl], nbneurl[nl], nbneur ); mydie("Just too many neurons in total!\n");}
		}
}


int buildSLayer(int gr, int sc, int sourcegr, int sourcesc, int myprot, int myrad, int mynbaff)
{


	int nl = nblayers;

	int myeps = -1;
	if (EPSS2 == DENSE) myeps = 1;
	else if (EPSS2 == RADIUS) myeps = myrad;
	else if (EPSS2 == HALFRADIUS) myeps = myrad/2;
	else mydie("Wrong EPSS2 value!");

	if (myeps == 0) myeps = 1;

	numlayer[gr][sc][myprot] = nblayers;
	numprot[nl] = myprot;
	group[nl] = gr;
	scale[nl] = sc;
	sourcegroup[nl] = sourcegr;
	sourcescale[nl] = sourcesc;
	radius[nl] = myrad; 
	eps[nl]= myeps;
	nbaff[nl] = mynbaff;
	if (nbaff[nl] >= MAXAFF) mydie("Too many afferents in S layer!");
	numprot[nl] = myprot;

//	printf("Layer %d has %d affs\n", nl, mynbaff);

	int sl = numlayer[sourcegr][sourcesc][0]; // Just for the sizing of the layer

	xsize[nl] = ceil ((double)(xsize[sl] - 2 * myrad) / (double)myeps);
	ysize[nl] = ceil ((double)(ysize[sl] - 2 * myrad) / (double)myeps);
	size[nl] = xsize[nl] * ysize[nl];

	if ((xsize[nl] < 1) || (ysize[nl] < 1)) {printf("Layer %d (group %d, scale %d):\n", nl, gr, sc); mydie("Zero-sized layer!\n");}
	//printf("Size layer %d (group %d): %d neurs (%d x %d)\n", nl, group[nl], size[nl], xsize[nl], ysize[nl]);
	nbneurl[nl] = 0;

	offsetreal[nl] = offsetreal[sl] + myrad * epsreal[sl];
	epsreal[nl] = myeps * epsreal[sl];

	start[nl] = nbneur;
	for (int j=0; j < ysize[nl] ; j++)
		for (int i=0; i < xsize[nl] ; i++)
		{
			int targetx = i * myeps + myrad;
			int targety = j * myeps + myrad;
			int targetn = start[sl] + targetx + targety * xsize[sl];
			realx[nbneur] = realx[targetn];
			realy[nbneur] = realy[targetn];
			layerneur[nbneur] = nl;
			nbneur ++; nbneurl[nl] ++;
			if (nbneur > MAXNEUR) { printf("L %d (group %d), neur %d (total %d)\n", nl, group[nl], nbneurl[nl], nbneur ); mydie("Just too many neurons in total!\n");}
		}
	

	nblayersgr[gr] ++;
	if (nblayersgr[gr] == 1) first[gr] = nblayers;

	nblayers ++;
	if (nblayers >= MAXLAYERS) { printf("%d layers - \n", nblayers); mydie("Too many layers!");}

	return nl;
}

int buildCLayer(int gr, int sc, int sourcegr, int sourcesc, int ori, int myrad, int myeps)
{

	int nl = nblayers;


	numlayer[gr][sc][ori] = nblayers;
	numprot[nl] = ori;
	group[nl] = gr;
	scale[nl] = sc;
	sourcegroup[nl] = sourcegr;
	sourcescale[nl] = sourcesc;
	radius[nl] = myrad; 
	eps[nl]= myeps;

	nbaff[nl] = 0; // Will remain 0 for C layers
	nbsourcelayers[nl] = 0; 

	int sl = numlayer[sourcegr][sourcesc][0]; // Just for the sizing of the layer

	xsize[nl] = ceil ((double)(xsize[sl] - 2 * myrad) / (double)myeps);
	ysize[nl] = ceil ((double)(ysize[sl] - 2 * myrad) / (double)myeps);
	size[nl] = xsize[nl] * ysize[nl];
	nbneurl[nl] = 0;

	offsetreal[nl] = offsetreal[sl] + myrad * epsreal[sl];
	epsreal[nl] = myeps * epsreal[sl];

	start[nl] = nbneur;
	for (int j=0; j < ysize[nl] ; j++)
		for (int i=0; i < xsize[nl] ; i++)
		{
			int targetx = i * myeps + myrad;
			int targety = j * myeps + myrad;
			int targetn = start[sl] + targetx + targety * xsize[sl];
			realx[nbneur] = realx[targetn];
			realy[nbneur] = realy[targetn];
			layerneur[nbneur] = nl;
			nbneur ++; nbneurl[nl] ++;
			if (nbneur > MAXNEUR) { printf("L %d, neur %d (total %d)\n", nl, nbneurl[nl], nbneur ); mydie("Just too many neurons in total!\n");}
		}
			
	nblayersgr[gr] ++;
	nblayers ++;
	if (nblayers >= MAXLAYERS) { printf("%d layers - \n", nblayers); mydie("Too many layers!");}

	return nl;
	
}

void runLayerC(int nl) 
{
	int RFctrx, RFctry, xsizesl, ysizesl,  myi, myj, tgtx, tgty, tgtn, sl;

	printf("Running C layer %d, size %d, %d source layers (%d & %d)...\n", nl, size[nl], nbsourcelayers[nl], sourcelayer[nl][0], sourcelayer[nl][1]);

	for (int ns=0; ns < nbsourcelayers[nl]; ns++)
	{
		sl = sourcelayer[nl][ns];
			xsizesl = xsize[sl];
			ysizesl = xsize[sl];
		for (int n = start[nl]; n < start[nl] + size[nl]; n++)
		{
			v[n] = 0;  // I could slap myself sometimes...

			RFctrx = floor( (double) (realx[n] - offsetreal[sl]) / (double) epsreal[sl]); 
			RFctry = floor( (double) (realy[n] - offsetreal[sl]) / (double) epsreal[sl]);
			for (myi = -radius[nl]; myi <= radius[nl]; myi++)
				for (myj = -radius[nl]; myj <= radius[nl]; myj++)
				{
						tgtx = RFctrx + myi; tgty = RFctry + myj;
						tgtn = start[sl] + tgty * xsizesl + tgtx;
						
						if ((tgtx < 0) || (tgtx > xsizesl) || (tgty < 0) || (tgty > ysizesl))
						{
							printf("Layer %d, gr %d, sc %d, from source layer %d, sc %d, size %d (ns = %d)\n", nl, group[nl], scale[nl], sl, scale[sl], size[sl], ns);
							printf("tgtx %d, tgty %d, tgtn %d\n", tgtx, tgty, tgtn);
							mydie("Attempting to pick outside the source layer!");
						}

						if (v[n] < v[tgtn]) v[n] = v[tgtn];
				}
		}
	}
}


void runLayerS_rbf(int nl) // RBF activation function
{
	float myw, sumsqdiff;
	int RFctrx, RFctry, xsizesl, ysizesl, k, myi, myj, tgtx, tgty, tgtn;

	int p = numprot[nl];  int gr = group[nl];
	int sourcegr = sourcegroup[nl]; int sourcesc = sourcescale[nl];

	//printf("%d %d %d \n", nl, size[nl], start[nl]);
	for (int n = start[nl]; n < start[nl] + size[nl]; n++)
	{

		// Just for the size / offset / eps values... Of course we are assuming that
		// all the source layers feeding afferents to a particular S layer have the
		// same scale, and thus the same size(s)...
		// The following can probably be put outside the loop!
		int sl = numlayer[sourcegr][sourcesc][0]; 

		RFctrx = floor( (double) (realx[n] - offsetreal[sl]) / (double) epsreal[sl]); 
		RFctry = floor( (double) (realy[n] - offsetreal[sl]) / (double) epsreal[sl]);
		xsizesl = xsize[sl];
		ysizesl = xsize[sl];


		sumsqdiff=0;
		for (int numconn=0; numconn < nbaff[nl]; numconn++)
		{
			k = numconn * 4; 
			myi = (int)StoredProts[gr][p][k+1];
			//printf("%d|", myi);
			myj = (int)StoredProts[gr][p][k+2];
			myw = StoredProts[gr][p][k+3];
			tgtx = RFctrx + myi; tgty = RFctry + myj;
						if ((tgtx < 0) || (tgtx > xsizesl) || (tgty < 0) || (tgty > ysizesl))
							mydie("Attempting to pick outside the source layer!");
			tgtn = start[numlayer[sourcegr][sourcesc][(int)StoredProts[gr][p][k]]] + tgty * xsizesl + tgtx;
			 sumsqdiff += (v[tgtn] - myw) * (v[tgtn] - myw) ;
		}
		v[n] = expf(-sqrt(sumsqdiff) / (radius[nl]));
//		v[n] = expf(-sqrt(sumsqdiff) / (2* (sqrt(radius[nl])) ));
//		v[n] = expf(-sumsqdiff / (2* (radius[nl]) * (radius[nl]) ));
//		v[n] = expf(-(sumsqdiff) / (2 * 5 * 5));

	}
}

void runLayerS_fndp(int nl) // Fully normalized dot product - no sigma! (Additive constant in the denominator)
{
	float res, norm, normw, myw; 
	int RFctrx, RFctry, xsizesl, ysizesl, k, myi, myj, tgtx, tgty, tgtn;

	int p = numprot[nl];  int gr = group[nl];
	int sourcegr = sourcegroup[nl]; int sourcesc = sourcescale[nl];

	//printf("%d %d %d \n", nl, size[nl], start[nl]);
	for (int n = start[nl]; n < start[nl] + size[nl]; n++)
	{

		// Just for the size / offset / eps values... Of course we are assuming that
		// all the source layers feeding afferents to a particular S layer have the
		// same scale, and thus the same size(s)...
		// The following can probably be put outside the loop!
		int sl = numlayer[sourcegr][sourcesc][0]; 

		RFctrx = floor( (double) (realx[n] - offsetreal[sl]) / (double) epsreal[sl]); 
		RFctry = floor( (double) (realy[n] - offsetreal[sl]) / (double) epsreal[sl]);
		xsizesl = xsize[sl];
		ysizesl = ysize[sl];


		res=0; norm=0; normw=0;
		for (int numconn=0; numconn < nbaff[nl]; numconn++)
		{
			k = numconn * 4; 
			myi = (int)StoredProts[gr][p][k+1];
			//printf("%d|", myi);
			myj = (int)StoredProts[gr][p][k+2];
			myw = StoredProts[gr][p][k+3];
			tgtx = RFctrx + myi; tgty = RFctry + myj;
						if ((tgtx < 0) || (tgtx > xsizesl) || (tgty < 0) || (tgty > ysizesl))
							mydie("Attempting to pick outside the source layer!");
			tgtn = start[numlayer[sourcegr][sourcesc][(int)StoredProts[gr][p][k]]] + tgty * xsizesl + tgtx;
			res += myw * v[tgtn];
			norm += v[tgtn] * v[tgtn];
			normw += myw * myw;   
		}
		/*
		// Faster by ~20%, but less flexible:
		k = 0;
		for (int myp=0; myp < nbprots[sourcegr]; myp++)
		for (int myj = -radius[nl]; myj <= radius[nl]; myj++)
		for (int myi = -radius[nl]; myi <= radius[nl]; myi++)
		{
			myw = StoredProts[gr][p][k*4+3];
				tgtn = start[numlayer[sourcegr][sourcesc][(int)StoredProts[gr][p][k*4]]] + RFctrx + myi + (RFctry + myj) * xsizesl;
			res += myw * v[tgtn];
			norm += v[tgtn] * v[tgtn];
			normw += myw * myw;   
			k ++;
		}
		if (k != nbaff[nl]) { printf("k %d, nbaff[%d] %d\n", k, nl, nbaff[nl]);  mydie("Hmmm, not seen all afferents!\n");}
		*/

//		v[n] = res / (sqrtf(norm) * sqrtf(normw) + SIGMA); 
		v[n] = res / (sqrtf(norm) * sqrtf(normw)); 

		//v[n] = res / norm; 
	}
}

void runLayerS_ndp_old(int nl) // Normalized dot product - now with sigma! (Additive constant in the denominator)
{
	float res, norm, normw, myw; 
	int RFctrx, RFctry, xsizesl, ysizesl, k, myi, myj, tgtx, tgty, tgtn;

	int p = numprot[nl];  int gr = group[nl];
	int sourcegr = sourcegroup[nl]; int sourcesc = sourcescale[nl];

	//printf("%d %d %d \n", nl, size[nl], start[nl]);
	for (int n = start[nl]; n < start[nl] + size[nl]; n++)
	{

		// Just for the size / offset / eps values... Of course we are assuming that
		// all the source layers feeding afferents to a particular S layer have the
		// same scale, and thus the same size(s)...
		// The following can probably be put outside the loop!
		int sl = numlayer[sourcegr][sourcesc][0]; 

		RFctrx = floor( (double) (realx[n] - offsetreal[sl]) / (double) epsreal[sl]); 
		RFctry = floor( (double) (realy[n] - offsetreal[sl]) / (double) epsreal[sl]);
		xsizesl = xsize[sl];
		ysizesl = ysize[sl];


		res=0; norm=0; normw=0;
		for (int numconn=0; numconn < nbaff[nl]; numconn++)
		{
			k = numconn * 4; 
			myi = (int)StoredProts[gr][p][k+1];
			//printf("%d|", myi);
			myj = (int)StoredProts[gr][p][k+2];
			myw = StoredProts[gr][p][k+3];
			tgtx = RFctrx + myi; tgty = RFctry + myj;
						if ((tgtx < 0) || (tgtx > xsizesl) || (tgty < 0) || (tgty > ysizesl))
							mydie("Attempting to pick outside the source layer!");
			tgtn = start[numlayer[sourcegr][sourcesc][(int)StoredProts[gr][p][k]]] + tgty * xsizesl + tgtx;
			
			res += myw * v[tgtn];
			norm += v[tgtn] * v[tgtn];
			normw += myw * myw;   
		}

		v[n] = res / (sqrtf(norm) * sqrtf(normw) + SIGMA); 

		//v[n] = res / norm; 
	}
}


void runLayerS_ndp(int nl) // Normalized dot product - now with sigma! (Additive constant in the denominator)
{
	float res, norm, normw, myw, myv;
	int RFctrx, RFctry, xsizesl, ysizesl, k, myi, myj, tgtx, tgty, tgtn;

	int p = numprot[nl];  int gr = group[nl];
	int sourcegr = sourcegroup[nl]; int sourcesc = sourcescale[nl];

	float sp[4*MAXAFF];
	int spint[4*MAXAFF];

	normw = 0;
	for (int i=0; i < 4*nbaff[nl]; i++)
	{
		sp[i] = StoredProts[gr][p][i];
		if (i % 4 == 3)
			normw += sp[i] * sp[i];
		else
			spint[i] = (int)sp[i];
	}
	int startneurons[MAXAFF];
	for (k=0; k < nbaff[nl]; k++)
	{
		startneurons[k] = start[numlayer[sourcegr][sourcesc][spint[k*4]]];
	}

		// Just for the size / offset / eps values... Of course we are assuming that
		// all the source layers feeding afferents to a particular S layer have the
		// same scale, and thus the same size(s)...
		int sl = numlayer[sourcegr][sourcesc][0]; 
		xsizesl = xsize[sl];
		ysizesl = ysize[sl];

	//printf("%d %d %d \n", nl, size[nl], start[nl]);
	for (int n = start[nl]; n < start[nl] + size[nl]; n++)
	{


		RFctrx = floor( (double) (realx[n] - offsetreal[sl]) / (double) epsreal[sl]); 
		RFctry = floor( (double) (realy[n] - offsetreal[sl]) / (double) epsreal[sl]);


		res=0; norm=0; 
		//normw = 0;

		for (int numconn=0; numconn < nbaff[nl]; numconn++)
		{
			k = numconn * 4; 
			
			myi = spint[k+1];
			myj = spint[k+2];
			myw = sp[k+3];
			
			/*myi = (int)StoredProts[gr][p][k+1];
			myj = (int)StoredProts[gr][p][k+2];
			myw = StoredProts[gr][p][k+3];*/
			
			tgtx = RFctrx + myi; 
			tgty = RFctry + myj;
			tgtn = startneurons[numconn] + tgty * xsizesl + tgtx;
			if ((tgtx < 0) || (tgtx > xsizesl) || (tgty < 0) || (tgty > ysizesl))
						{
							printf("S Layer %d, gr %d, sc %d, from source layer %d, sc %d, size %d \n", nl, group[nl], scale[nl], sl, scale[sl], size[sl]);
							printf("tgtx %d, tgty %d (myi %d myj %d), tgtn %d\n", tgtx, tgty, myi, myj, tgtn);
							mydie("Attempting to pick outside the source layer!");
						}

			//tgtn = start[numlayer[sourcegr][sourcesc][(int)StoredProts[gr][p][k]]] + tgty * xsizesl + tgtx;


			myv = v[tgtn];
			res += myw * myv;
			norm+= myv * myv;
			//normw += myw * myw;   
		}

		//printf("%d\n", resint);
		v[n] =  res / (sqrtf(norm) * sqrtf(normw) + SIGMA); 

		//v[n] = res / norm; 
	}
}

void runS1Layer(int nl)  // CONTRAST INVARIANT (divided by norm of the inputs, as in the hmin code)!
{
	int k, myi, myj, tgtn, RFctrx, RFctry, myrad, xs0; float myw, res, sumsq;

	xs0 = xsize[0];
	for (int n=start[nl]; n < start[nl] + size[nl]; n++)
	{
		// We find the centre of the RF in the source layer, which is layer 0. This makes things much simpler!
		// The RF centre is simply the real (global) coordinates of the cell!
		RFctrx =  realx[n] ;
		RFctry =  realy[n];
		res=0; sumsq=0;
/*		int totalnbconns = diamRF_S1[nl]*diamRF_S1[nl];
//		for (int numconn=0; numconn < totalnbconns; numconn++)
		for (int k=0; k < 4*totalnbconns; k+=4)
		{
			tgtn = (RFctrx + ProtOfLayer[nl][k+1] + (RFctry + ProtOfLayer[nl][k+2]) * xs0);
			myw = ProtOfLayer[nl][k+3];

			//			printf("Conn %d (position %d in the descriptor of connections): Trying to connect from cell %d in layer %d (which is %d, %d away from us)\n", numconn, k, tgtn, sl, myi, myj);
			res += myw * v[tgtn];
			sumsq += v[tgtn] * v[tgtn];
		}*/ // ~ 30% slower!
		k=0;
		myrad=(int)floor(diamRF_S1[nl]/2.0);
		for (myj=-myrad; myj <= myrad; myj++)
			for (myi=-myrad; myi <= myrad; myi++)
			{
				tgtn = RFctrx+myi + (RFctry + myj) * xsize[0];
//			tgtn = (RFctrx + ProtOfLayer[nl][k*4+1] + (RFctry + ProtOfLayer[nl][k*4+2]) * xs0); // ~30% slower! And it's not the *4!
				myw = ProtOfLayer[nl][k*4+3];
				k++;
				res += myw * v[tgtn];
				sumsq += v[tgtn] * v[tgtn];
			}
		if (res < 0) res = -res;
		v[n] = res  / (sqrtf(sumsq) + SIGMAS1);  
	}
}


// Get the activity values of all the cells in the layers specified in layers[] that fall within
// a diameter diam from cell n, in global coordinates
// Store them in resdescr, which has almost the same format as conndescr, except that the orientation
// is stored instead of the actual source layer number
// This is the format for the S2 prototypes.

// this is one of the ill-designed aspects of the code: normally, just giving
// the S2 cell should automatically give you the scale of the C1 layers that
// compose its RF.

void getActivRF_S2(int n, int sc, int diam, float *resdescr)
{
	int posdescr=0;
	int myrad = (int)floor(diam / 2); 
	int tgtx, tgty, tgtn;
	for (int ori=0; ori < 4; ori ++)
	{
		int sl = numlayer[C1][sc][ori];
		for (int myj = -myrad; myj <= myrad; myj++)
			for (int myi = -myrad; myi <= myrad; myi++)
			{
				tgtx = floor( (double) (realx[n] - offsetreal[sl]) / (double) epsreal[sl]);
				tgty = floor( (double) (realy[n] - offsetreal[sl]) / (double) epsreal[sl]);
				tgtx = tgtx + myi; tgty = tgty + myj;
				tgtn = start[sl] + tgty * xsize[sl] + tgtx;
				// We store the orientation rather than the layer number; this will then be
				// translated into various appropriate layer numbers at various scales
				resdescr[posdescr] = ori; 
				resdescr[posdescr+1] = myi;
				resdescr[posdescr+2] = myj;
				resdescr[posdescr+3] = v[tgtn];
				//				printf("tgtn %d from layer %d has activity %.3f.\n", tgtn, sl, resdescr[posdescr+3]);
				posdescr += 4;
			}
	}
}

void getActivRF_S3(int n, int sc, int diam, float *resdescr)
{
	int posdescr=0;
	int myrad = (int)floor(diam / 2); 
	int tgtx, tgty, tgtn;
	for (int p=0; p < NBS2PROTS; p++)
	{
		int sl = numlayer[C2][sc][p];
		for (int myj = -myrad; myj <= myrad; myj++)
			for (int myi = -myrad; myi <= myrad; myi++)
			{
				tgtx = floor( (double) (realx[n] - offsetreal[sl]) / (double) epsreal[sl]);
				tgty = floor( (double) (realy[n] - offsetreal[sl]) / (double) epsreal[sl]);
				tgtx = tgtx + myi; tgty = tgty + myj;
				tgtn = start[sl] + tgty * xsize[sl] + tgtx;
				// We store the prototype rather than the layer number; this will then be
				// translated into various appropriate layer numbers at various scales
				resdescr[posdescr] = p; 
				resdescr[posdescr+1] = myi;
				resdescr[posdescr+2] = myj;
				resdescr[posdescr+3] = v[tgtn];
				//				printf("tgtn %d from layer %d has activity %.3f.\n", tgtn, sl, resdescr[posdescr+3]);
				posdescr += 4;
			}
	}
}


// Get the activities of all the cells that fall within the possible RF of a
// given (S-layer!) cell. Put them in resdescr, in the usual format (a series of floats,
// giving the source layer prototype/orientation, relative X and Y from the
// centre of the RF, and activity value).
// Cell n must be from a S-layer, lest undef occur. Also, sufficient memory
// must be allocated in resdescr beforehand (not a problem, since resdescr is
// usually a pointer to the apparopriate point in StoredProts).

void getRFValsStest(int n, float *resdescr)
{
	int nl = layerneur[n];
	int posdescr=0;
	int myrad = radius[nl];
	int sourcegr = sourcegroup[nl];
	int sourcesc = sourcescale[nl];
	int tgtx, tgty, tgtn;
	//printf("Extracting RF of neuron %d (layer %d, sourcegroup %d, sourcescale %d, radius %d)\n",
	//		n, nl, sourcegr, sourcesc, myrad);
	for (int p=0; p < nbprots[sourcegr]; p++)
	{
		int sl = numlayer[sourcegr][sourcesc][p];
		printf("%d %d %d: %d\n", sourcegr, sourcesc, p, sl); 
		for (int myj = -myrad; myj <= myrad; myj++)
			for (int myi = -myrad; myi <= myrad; myi++)
			{
				tgtx = floor( (double) (realx[n] - offsetreal[sl]) / (double) epsreal[sl]);
				tgty = floor( (double) (realy[n] - offsetreal[sl]) / (double) epsreal[sl]);
				tgtx = tgtx + myi; tgty = tgty + myj;
				tgtn = start[sl] + tgty * xsize[sl] + tgtx;
				resdescr[posdescr] = p; 
				resdescr[posdescr+1] = myi;
				resdescr[posdescr+2] = myj;
				resdescr[posdescr+3] = v[tgtn];
				if (posdescr == 0)
								//printf("tgtn %d from layer %d, has activity %.3f.\n", tgtn, sl, resdescr[posdescr+3]);
				posdescr += 4;
			}

	}
}
void getRFValsS(int n, float *resdescr)
{
	int nl = layerneur[n];
	int posdescr=0;
	int myrad = radius[nl];
	int sourcegr = sourcegroup[nl];
	int sourcesc = sourcescale[nl];
	int tgtx, tgty, tgtn;
	for (int p=0; p < nbprots[sourcegr]; p++)
	{
		int sl = numlayer[sourcegr][sourcesc][p];
		printf("%d %d %d: %d\n", sourcegr, sourcesc, p, sl); 
		for (int myj = -myrad; myj <= myrad; myj++)
			for (int myi = -myrad; myi <= myrad; myi++)
			{
				// NOTE: A "Floating point exception" around here usually
				// implies the existence of zero-sized layers (e.g. if you have
				// too large scales for the size of your input images)
				tgtx = floor( (double) (realx[n] - offsetreal[sl]) / (double) epsreal[sl]);
				tgty = floor( (double) (realy[n] - offsetreal[sl]) / (double) epsreal[sl]);
				tgtx = tgtx + myi; tgty = tgty + myj;
				tgtn = start[sl] + tgty * xsize[sl] + tgtx;
				resdescr[posdescr] = p; 
				resdescr[posdescr+1] = myi;
				resdescr[posdescr+2] = myj;
				resdescr[posdescr+3] = v[tgtn];
								//printf("tgtn %d from layer %d has activity %.3f.\n", tgtn, sl, resdescr[posdescr+3]);
				posdescr += 4;
			}

	}
}


void sortStoredProts(int gr)
{
	// This sorts the different afferents/connections of a given
	// prototype from highest to lowest weight. 
	// Gnome sort! Simple implementation and we don't care much about efficiency.
	for (int p=0; p < nbprots[gr]; p++)
	{
		printf("Sorting afferents for prot %d, group %d\n", p, gr);
		int diam = 2*radius[first[gr]] + 1;
		int maxaff= diam*diam*nbprots[sourcegroup[first[gr]]]; // Ideally, sourcegroup should be a property of groups rather than layers?...
		int numaff=0; 
		while (numaff < maxaff)
		{
			//printf("%d-%f-%f\n", numaff, StoredProts[gr][p][numaff*4 + 3] , StoredProts[gr][p][(numaff-1)*4 + 3]);
			if (numaff == 0 || StoredProts[gr][p][numaff*4 + 3] <= StoredProts[gr][p][(numaff-1)*4 + 3])
				numaff ++;
			else{
				for (int k=0; k < 4; k++)
				{
					float tmp = StoredProts[gr][p][numaff*4+k];
					StoredProts[gr][p][numaff*4 + k] = StoredProts[gr][p][(numaff-1)*4 + k];
					StoredProts[gr][p][(numaff-1)*4 + k] = tmp;
				}
				numaff--;
			}
		}

	}
}

void randStoredProts(int gr)
{
	// This randomizes the order in which the different afferents/connections of a given
	// prototype are specified. The prototype extractor stores all the possible
	// afferents, which is more than needed. So we randomize them before we
	// store them, and then we only use the first nbaffs[gr] of them.
	for (int p=0; p < nbprots[gr]; p++)
	{
		int diam = 2*radius[first[gr]] + 1;
		for (int numaff=0; numaff < diam*diam*nbprots[sourcegroup[first[gr]]]; numaff++)
		{
			int otheraff = random() % (diam*diam*nbprots[sourcegroup[first[gr]]]);
			for (int k=0; k < 4; k++)
			{
				float tmp = StoredProts[gr][p][numaff*4+k];
				StoredProts[gr][p][numaff*4 + k] = StoredProts[gr][p][otheraff*4 + k];
				StoredProts[gr][p][otheraff*4 + k] = tmp;
				//if (k == 0) printf("o%d o%d (numaff %d otheraff %d) | ", (int)StoredProts[gr][p][numaff*4 + k], (int)StoredProts[gr][p][otheraff*4 + k], numaff, otheraff);
			}
		}

	}
}

void writeProt(int gr, int p, const char *fname)
{
	FILE *f = fopen(fname, "w");
//	int diam = 2*radius[first[gr]] + 1;
	for (int i=0; i< MAXAFF; i++)
	{
		fprintf(f, "i%d j%d: or%d w%f |\n", (int)StoredProts[gr][p][i * 4+1], (int)StoredProts[gr][p][i * 4+2], (int)StoredProts[gr][p][i * 4], StoredProts[gr][p][i * 4+3]);
	}
	fclose(f);

}

void addStoC(int sl, int nl)
{
	printf("Adding source layer %d to C layer %d...\n", sl, nl);
	sourcelayer[nl][nbsourcelayers[nl]] = sl;
	nbsourcelayers[nl] ++;
}


void loadpic(const char *filename)
{

	// INPUT FILES TO LOADPIC MUST BE RAW FILES!!

	unsigned char realpic[IMAGEXSIZE*IMAGEYSIZE], tmpc;
	float max = -9999;
	for (int n=0; n < MAXNEUR; n++)
		v[n] = 0;
	if (!strstr(filename, "raw")) {printf("Filename: %s\n", filename); mydie("Image files must be in the raw format (see prepimages.m, utils.h)!\n");}
	FILE *f = fopen(filename, "r");
	if (!f) { printf("Couldn't open the image file [%s]!", filename); mydie("Bye!"); }
	//printf("Reading image %s\n", filename);
	for (int i=0; i < IMAGEXSIZE * IMAGEYSIZE; i++)
		if( !fread(&realpic[i], 1, 1, f)) {printf("%d %d %d - \n", i, IMAGEXSIZE, IMAGEYSIZE); fflush(stdout); mydie("Image file does not contain enough pixels!\n!");}
	if ( fread(&tmpc, 1, 1, f)) mydie("Image file has too many pixels!\n");
	fclose(f);
	for (int i=0; i < IMAGEXSIZE*IMAGEYSIZE; i++)
	{
		pic[i] = .00001 + (float)realpic[i] / 257.0;
		v[i] = .00001 + (float)realpic[i] / 257.0;
		if (max < v[i]) max = v[i];
	}

}


void makeV1Weights(float *weights, int diam, float dir)
{
	float min = 9999; 
	float sum = 0;
	for (int n=0; n < diam * diam; n++)
	{
		int i = n*4;
		//printf("%d ", i);
		weights[i] = 0; // V1 simple cells only takes connections from source layer 0
		int myrad = (int)floor(diam / 2);
		int myx = n % diam - myrad;
		int myy = n / diam - myrad;
		float myw;
		weights[i+1] = (float)myx;
		weights[i+2] = (float)myy;
		weights[i+3] = 0;


		float theta = (float)dir;
		float sigma = 0.0036 * (float)diam * (float)diam + 0.35 * (float)diam + 0.18; 

		//		sigma = sigma * 2.0; 

		float lambda = sigma / .8;
		float i2 = myx * cos(theta) + myy * sin(theta);
		float j2 = -myx * sin(theta) + myy * cos(theta);

		myw = expf ( -(i2 * i2 + .3*.3 * j2 * j2) / (2 * sigma * sigma)) * cos(2 * M_PI * i2 / lambda);

//		if (myw<0) myw = 0;

		// circular cropping:
		if (sqrtf(myx*myx + myy*myy) > myrad)
			myw = 0;

		if (myw < min) min = myw;
		
		sum += myw;

		weights[i+3] = myw;

	}
	float sumsq=0;
		for (int n=0; n < diam * diam; n++) {
			weights[4 * n + 3] -= sum / (float)(diam*diam);
			sumsq += weights[4 * n + 3] * weights[4 * n + 3];
		}
		for (int n=0; n < diam * diam; n++) 
			weights[4 * n + 3] /= sqrtf(sumsq);
}



void readS2Prots(const char *fname)
{
	FILE *f = fopen(fname, "r");
	if (!f) mydie("Couldn't open S2 prototypes file!\n");
	int thrw;
	for (int n=0; n < NBS2PROTS; n++)
		for (int i=0; i < 4* NBCONNS2PROT; i++)
			thrw = fread(&StoredProts[S2][n][i], sizeof(float), 1, f);
	fclose(f);
}
void readS3Prots(const char *fname)
{
	FILE *f = fopen(fname, "r");
	if (!f) mydie("Couldn't open S2 prototypes file!\n");
	int thrw;
	for (int n=0; n < NBS3PROTS; n++)
		for (int i=0; i < 4* NBCONNS3PROT; i++)
		{
			thrw = fread(&StoredProts[S3][n][i], sizeof(float), 1, f);
			//printf("%f ", StoredProts[S3][n][i]);
		}
	fclose(f);
}


void readProts(const char *fname)
{
	FILE *f = fopen(fname, "r");
	if (!f) mydie("Couldn't read prototypes file!");
			int thrw = fread(&StoredProts, sizeof(StoredProts), 1, f);
			thrw = 0;
	fclose(f);

}
void saveProts(const char *fname)
{
	FILE *f = fopen(fname, "w");
		int thrw = fwrite(&StoredProts, sizeof(StoredProts), 1, f);
		thrw = 0;
	fclose(f);

}

void writeS2Prots(const char *fname)
{
//	FILE *f = fopen("s2prots.dat", "w");
	FILE *f = fopen(fname, "w");
	for (int n=0; n < nbprots[S2]; n++)
	{
		// Writing down the prototype in its own file, in ASCII, so Matlab can read it
		char str[80]; 
		sprintf(str, "s2prot%d.txt", n);
		/*FILE *f2 = fopen(str, "w");
		for (int c=0; c < NBCONNS2PROT; c++)
			fprintf(f2, "%f %f %f %f\n", StoredProts[S2][n][4*c], StoredProts[S2][n][4*c+1], StoredProts[S2][n][4*c+2], StoredProts[S2][n][4*c+3]);
		fclose(f2);*/

		for (int i=0; i < 4* NBCONNS2PROT; i++)
			fwrite(&StoredProts[S2][n][i], sizeof(float), 1, f);
	}
	fclose(f);

}
void writeS3Prots(const char *fname)
{
//	FILE *f = fopen("s3prots.dat", "w");
	FILE *f = fopen(fname, "w");
	for (int n=0; n < NBS3PROTS; n++)
		for (int i=0; i < 4* NBCONNS3PROT; i++)
			fwrite(&StoredProts[S3][n][i], sizeof(float), 1, f);
	fclose(f);
}
	

/*void convertstoredprots()
{
//float StoredProts [MAXGROUP][MAXPROTS][4 * MAXAFF]; 
	for (int gr=0; gr < MAXGROUP; gr++)
		for (int p=0; p < MAXPROTS; p++)
			for (int k=0; k < 4*MAXAFF; k++)
				StoredProtsInt[gr][p][k] = (int)(100.0*StoredProts[gr][p][k]);
}
void convertv()
{ for (int i=0; i < MAXNEUR; i++)	vint[i] = (int)(100.0*v[i]);
}*/


void buildNetwork()
{


	nbneur = 0;
	realx = (int*)calloc(MAXNEUR, sizeof(int));
	realy = (int*)calloc(MAXNEUR, sizeof(int));
	v = (float*) calloc(MAXNEUR, sizeof(float));

	if((!realx) || (!realy) || (!v)) mydie("Memory allocation failed!");

	layerneur = (int*)calloc(MAXNEUR, sizeof(int));

	if((!layerneur)) mydie("Memory allocation failed!");

	for (int gr=0; gr < MAXGROUP; gr++)
		nblayersgr[gr] = 0;


	if (HANDLEC1SCALES == MERGE)
		NBSCALES_C1 = NBSCALES_S1 / 2;
	else if (HANDLEC1SCALES == BLUR)
		NBSCALES_C1 = NBSCALES_S1 -1;
	else if (HANDLEC1SCALES == KEEP)
		NBSCALES_C1 = NBSCALES_S1;
	else mydie("How should I handle C1 scales?...");

	NBSCALES_C2 = NBSCALES_C1 / 2;
	NBSCALES_S2 = NBSCALES_C1;

	nbprots[S1] = 4; // 4 orientations
	nbprots[C1] = 4; // 4 orientations
	nbprots[S2] = NBS2PROTS;
	nbprots[C2] = NBS2PROTS;
	nbprots[S2b_5] = NBS2PROTS;
	nbprots[S2b_7] = NBS2PROTS;
	nbprots[S2b_9] = NBS2PROTS;
	nbprots[S2b_11] = NBS2PROTS;
	nbprots[S2b_13] = NBS2PROTS;
	nbprots[S3] = NBS3PROTS;
	nbprots[C3] = NBS3PROTS;

	// Making Layer 0, which will contain the image data. 
	eps[0] = 1;
	rad[0] = 0;
	offsetreal[0] = 0;
	epsreal[0] = 1;
	xsize[0] = IMAGEXSIZE;
	ysize[0] = IMAGEYSIZE;
	size[0] = IMAGEXSIZE*IMAGEYSIZE;
	start[0] = 0;
	for (int j=0; j < ysize[0]; j++) // columns first, then lines
		for (int i=0; i < xsize[0]; i++)
		{
			realx[nbneur] = i; 
			realy[nbneur] = j;
			layerneur[nbneur] = 0;
			nbneur ++;
		}

	nblayers = 1;


	printf("\nS0 layer filled!\n"); fflush(stdout);

	// Now to the actual processing layers
	// First, the S1 layers, at NBSCALES_S1 scales (SHOULD BE AN EVEN NUMBER!)  and 4 orientations

	first[S1] = nblayers;
	for (int sc=0; sc < NBSCALES_S1; sc++)
	{
		int diam = 7 + 2 * sc;
		for (int ori=0; ori < 4; ori++)
		{
			printf("Scale %d, ori %d...\n", sc, ori);
			weightsS1[sc][ori] = (float *) malloc(sizeof(float) * 4 * diam * diam);
			if (weightsS1[sc][ori] == NULL)  mydie("Error in memory allocation!\n"); 
			makeV1Weights(weightsS1[sc][ori], diam, ori * M_PI / 4);

			int dd=-1;
			if (EPSS1 == DENSE) dd = 1;
			else if (EPSS1 == RADIUS) dd = ( (int) floor(diam / 2) );
			else if (EPSS1 == HALFRADIUS) dd = ( (int) floor(diam / 4) );
			else if (EPSS1 == EVERYOTHER) dd = 2;
			else mydie("Wrong EPSS1 value!");
			if (dd < 1) dd = 1;

			buildS1Layer(nblayers, 0, diam, dd);

			ProtOfLayer[nblayers] = weightsS1[sc][ori];
			numlayer[S1][sc][ori] = nblayers;
			scale[nblayers] = sc;
			numprot[nblayers] = ori;
			nblayers++;
			nblayersgr[S1] ++ ;
			//printf("num S1 layer (sc %d ori %d): %d\n", sc, ori, nblayersgr[S1]);
		}
	}


	FILE *f=fopen("v1_0.txt", "w");
	for (int i=0; i< 13 * 13; i++)
	{
		fprintf(f, "%f ", weightsS1[3][0][i*4 + 3]);
		if (i % 13 == 12) fprintf(f, "\n");
	}
	fclose(f);
	f=fopen("v1_1.txt", "w");
	for (int i=0; i< 13 * 13; i++)
	{
		fprintf(f, "%f ", weightsS1[3][1][i*4 + 3]);
		if (i % 13 == 12) fprintf(f, "\n");
	}
	fclose(f);


	printf("\nS1 layers built!\n"); fflush(stdout);

	printf("Size last layer: %d\n", size[nblayers-1]);

	// Now for the C1 layers, NBSCALES_C1 = NBSCALES_S1 /2, at all 4 orientations
	// A C layer should be built based on the smallest source S layer (i.e. the one with the largest scale)

	first[C1] = nblayers;
	for (int ori=0; ori < 4; ori++) 
		for (int sc=0; sc < NBSCALES_C1; sc++)
		{
			int sourcesc;
			int nl=0;

			// The "source scale" is used for sizing this layer - it must be
			// the smaller of the source scales (i.e. the one with the highest
			// scale)
			if (HANDLEC1SCALES == MERGE)
				sourcesc = 2 * sc + 1;
			else if (HANDLEC1SCALES == BLUR)
				sourcesc = sc + 1; 
			else if (HANDLEC1SCALES == KEEP)
				sourcesc = sc ;
			else mydie("How should I handle scales?");


			if (EPSC1 == DENSE)
				nl = buildCLayer(C1, sc, S1, sourcesc, ori, 1, 1);
			else if (EPSC1 == POGGIO)
				nl = buildCLayer(C1, sc, S1, sourcesc, ori, 3 +  sc, 3 + 2*sc);
			else if (EPSC1 == HALFRADIUS)
				nl = buildCLayer(C1, sc, S1, sourcesc, ori, 2 + sc/2, 2 + sc);
			else if (EPSC1 == RADIUS)
				nl = buildCLayer(C1, sc, S1, sourcesc, ori, 3 + sc, 3 + sc);
			else if (EPSC1 == EVERYOTHER)
				nl = buildCLayer(C1, sc, S1, sourcesc, ori, 2, 2);  
			else if (EPSC1 == EVERYOTHER4)
				nl = buildCLayer(C1, sc, S1, sourcesc, ori, 4, 2);  // Radius can be 2 or 4
			else mydie("Wrong EPSC1 value!");

			if (HANDLEC1SCALES == MERGE)
			{
				addStoC(numlayer[S1][2*sc][ori], nl);
				addStoC(numlayer[S1][2*sc+1][ori], nl);
			}
			else if (HANDLEC1SCALES == BLUR)
			{
				addStoC(numlayer[S1][sc][ori], nl);
				addStoC(numlayer[S1][sc+1][ori], nl);
			}
			else if (HANDLEC1SCALES == KEEP)
			{
				addStoC(numlayer[S1][sc][ori], nl);
			}
		}

	printf("\nS1 & C1 layers built!\n"); fflush(stdout);
	printf("Size last layer: %d\n", size[nblayers-1]);

	// S2 layers, one per C1 scale, taking input from C1 layers of all orientations at that particular band

	for (int sc=0; sc < NBSCALES_C1; sc++)
		for (int p=0; p < nbprots[S2]; p++)
			buildSLayer(S2, sc, C1, sc, p, 1, 10); // Radius is always 1 (3x3 RF), and there are 10 afferents

	printf("\nS2 layers built!\n"); fflush(stdout);

	for (int sc=0; sc < NBSCALES_C1; sc++)
		for (int p=0; p < nbprots[S2b_5]; p++)
			buildSLayer(S2b_5, sc, C1, sc, p, 2, 100); // Radius is 2 (5x5 RF), and there are 100 afferents

	printf("\nS2b_5 layers built!\n"); fflush(stdout);

	for (int sc=0; sc < NBSCALES_C1; sc++)
		for (int p=0; p < nbprots[S2b_7]; p++)
			buildSLayer(S2b_7, sc, C1, sc, p, 3, 100); // Radius is 3 (7x7 RF), and there are 100 afferents
	
	printf("\nS2b_7 layers built!\n"); fflush(stdout);

	for (int sc=0; sc < NBSCALES_C1; sc++)
		for (int p=0; p < nbprots[S2b_9]; p++)
			buildSLayer(S2b_9, sc, C1, sc, p, 4, 100); // Radius is 4 (9x9 RF), and there are 100 afferents
	
	printf("\nS2b_9 layers built!\n"); fflush(stdout);

	for (int sc=0; sc < NBSCALES_C1; sc++)
		for (int p=0; p < nbprots[S2b_11]; p++)
			buildSLayer(S2b_11, sc, C1, sc, p, 5, 100); // Radius is 5 (11x11 RF), and there are 100 afferents
	
	printf("\nS2b_11 layers built!\n"); fflush(stdout);

	for (int sc=0; sc < NBSCALES_C1; sc++)
		for (int p=0; p < nbprots[S2b_13]; p++)
			buildSLayer(S2b_13, sc, C1, sc, p, 6, 100); // Radius is 6 (13x13 RF), and there are 100 afferents

	printf("\nS2b_13 layers built!\n"); fflush(stdout);


}


void saveC1(char *filename)
{
	FILE *f = fopen(filename, "w");
	if (!f) {
		printf("Couldn't open file '%s' for writing out C1 !\n", filename);
		mydie("Bye!");
	} 
	printf("Saving C1 outputs in file %s\n", filename);
	fflush(stdout);
	for (int nl = first[C1]; nl < first[C1] + nblayersgr[C1]; nl++)
		for (int n=start[nl]; n < start[nl] + size[nl]; n++)
		{
			//printf("Saving cell %d of layer %d\n", n - start[nl], nl);
			fprintf(f, "%.15f \n", v[n]);
		}
	fclose(f);
}

void saveC2(char *filename)
{
	float C2out[nbprots[C2]];
	FILE *f = fopen(filename, "w");
	if (!f) {
		printf("Couldn't open file '%s' for writing out C2 !\n", filename);
		mydie("Bye!");
	}
	printf("Saving C2 outputs in file %s\n", filename);fflush(stdout);
	for (int p=0; p < nbprots[C2]; p++)
	{
		float max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2!
		{
			int nl = numlayer[S2][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		C2out[p] = max;
		fprintf(f, "%.15f \n", C2out[p]);
	}
	fclose(f);

}

void saveS2bShort(char *filename)
{
//	float C2out[NBS2PROTS]; // S2b has the same number of prots as S2
	FILE *f = fopen(filename, "w");
	if (!f) {
		printf("Couldn't open file '%s' for writing out S2b !\n", filename);
		mydie("Bye!");
	}
	printf("Saving S2b outputs (binary format) in file %s\n", filename);fflush(stdout);
	for (int p=0; p < NBS2PROTS; p++)
	{
		float max = -9999;
		/*for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_5][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_7][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);*/
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_9][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
			{
			unsigned short int myp = p;
			unsigned short int mysc = sc;
			unsigned short int myrx = realx[n];
			unsigned short int myry = realy[n];
			unsigned short int myv = (unsigned short int) (v[n] * 65535.0);
			//	fprintf(f, "%d %d %d %d %.6f \n", p, sc, realx[n], realy[n], v[n]);
				fwrite(&myp, sizeof(unsigned short int), 1, f);
				fwrite(&mysc, sizeof(unsigned short int), 1, f);
				fwrite(&myrx, sizeof(unsigned short int), 1, f);
				fwrite(&myry, sizeof(unsigned short int), 1, f);
				fwrite(&myv, sizeof(unsigned short int), 1, f);
				//if (p == 248)
				//	printf("%d\n", myv);
			}
		}
		/*max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_11][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_13][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);*/
	}

	fclose(f);
}

void saveS2b(char *filename)
{
//	float C2out[NBS2PROTS]; // S2b has the same number of prots as S2
	FILE *f = fopen(filename, "w");
	if (!f) {
		printf("Couldn't open file '%s' for writing out C2b !\n", filename);
		mydie("Bye!");
	}
	printf("Saving S2b outputs in file %s\n", filename);fflush(stdout);
	for (int p=0; p < NBS2PROTS; p++)
	{
		float max = -9999;
		/*for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_5][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_7][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);*/
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_9][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				fprintf(f, "%d %d %d %d %.6f \n", p, sc, realx[n], realy[n], v[n]);
		}
		/*max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_11][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_13][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);*/
	}

	fclose(f);
}

void saveC2b(char *filename)
{
//	float C2out[NBS2PROTS]; // S2b has the same number of prots as S2
	FILE *f = fopen(filename, "w");
	if (!f) {
		printf("Couldn't open file '%s' for writing out C2b !\n", filename);
		mydie("Bye!");
	}
	printf("Saving C2b outputs in file %s\n", filename);fflush(stdout);
	for (int p=0; p < NBS2PROTS; p++)
	{
		float max = -9999;
		/*for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_5][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_7][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);*/
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_9][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);
		/*max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_11][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_13][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);*/
	}

	fclose(f);

/*

	f = fopen(strcat(filename, ".halfscales.out"), "w");
	if (!f) {
		printf("Couldn't open file '%s' for writing out C2b !\n", filename);
		mydie("Bye!");
	}
	printf("Saving C2b outputs in file %s\n", filename);fflush(stdout);
	for (int p=0; p < NBS2PROTS; p++)
	{
		float max = -9999;
		for (int sc=0; sc < NBSCALES_C1/2; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_5][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1/2; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_7][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1/2; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_9][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1/2; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_11][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);
		max = -9999;
		for (int sc=0; sc < NBSCALES_C1/2; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_13][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		fprintf(f, "%.15f \n", max);
	}
	fclose(f);*/
}

void saveC2bMean(char *filename)
{
//	float C2out[NBS2PROTS]; // S2b has the same number of prots as S2
	FILE *f = fopen(filename, "w");
	if (!f) {
		printf("Couldn't open file '%s' for writing out C2b !\n", filename);
		mydie("Bye!");
	}
	printf("Saving C2b outputs in file %s\n", filename);fflush(stdout);
	for (int p=0; p < NBS2PROTS; p++)
	{
		float max = -9999;
		max = -9999;
			int nbvals=0; float meanvals=0;
		for (int sc=0; sc < NBSCALES_C1; sc++) // which is also the number of scales of S2 and S2b!
		{
			int nl = numlayer[S2b_13][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > 0)
				{
					nbvals++;
				   meanvals += v[n];
				}
		}
		float myresult = (float)meanvals / (float)nbvals;
		fprintf(f, "%.15f \n", myresult);
	}

	fclose(f);

}

void saveC3(const char *filename)
{
	float C3out[NBS3PROTS];
	FILE *f = fopen(filename, "w");
	if (!f) {
		printf("Couldn't open file '%s' for writing out C3 !\n", filename);
		mydie("Bye!");
	}
	printf("Saving C3 outputs in file %s\n", filename);fflush(stdout);
	for (int p=0; p < NBS3PROTS; p++)
	{
		float max = -9999;
		for (int sc=0; sc < NBSCALES_C2; sc++)
		{
			int nl = numlayer[S3][sc][p];
			for (int n=start[nl]; n < start[nl] + size[nl]; n++)
				if (v[n] > max)
					max = v[n];
		}
		C3out[p] = max;
		fprintf(f, "%.15f \n", C3out[p]);
	}
	fclose(f);
}


void printInfo()
{

	for (int gr=0; gr < MAXGROUP; gr++)
	{
		int mysize = size[first[gr] + nblayersgr[gr] - 1];
		printf("%d layers in group %d (size last layer %d)\n ", nblayersgr[gr], gr, mysize);
	}
	printf("%d neurons  and %d layers in total.\n", nbneur, nblayers);
	printf("Index of (S2) layer for S2 prototype 0 at scale 0: %d\n", numlayer[S2][0][0]);
	printf("Index of (S3) layer for S3 prototype 0 at scale 0: %d\n", numlayer[S3][0][0]);

}


void saveLayer(int nl, const char* fn)
{
	FILE *f = fopen(fn, "w");
	for (int j=0; j < ysize[nl]; j++)
	{
		for (int i=0; i < xsize[nl]; i++)
			fprintf(f, "%f ", v[start[nl] + i + j * xsize[nl]]);
		fprintf(f,"\n");
	}

	fclose(f);
}

// Memory allocation functions

float **  alloc2Dfloatarray(int nbli, int nbco)
{

float **v;
	v = (float**) malloc(nbli * sizeof(float*));
	for (int i=0; i < nbli; i++)
		v[i] = (float*) malloc(nbco*sizeof(float));

	return v;
}
int **  alloc2Dintarray(int nbli, int nbco)
{

int **v;
	v = (int**) malloc(nbli * sizeof(int*));
	for (int i=0; i < nbli; i++)
		v[i] = (int*) malloc(nbco*sizeof(int));

	return v;
}
