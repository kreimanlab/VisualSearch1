#include "utils.h"
#include "time.h"
#include "string.h"




int main(int argc, char *argv[])
{
	FILE *f;
	char fn[200];
	char listinputimagesfilename[200];

	strcpy(listinputimagesfilename, "listinputfiles.txt");
	int startfile=0; int endfile= 9999999;
	for (int a=1; a < argc; a++)
	{
		if (strstr(argv[a], "-f"))
		{
				strcpy(listinputimagesfilename, &argv[a][2]);
		}
		if (strstr(argv[a], "-s"))
		{
				startfile = atoi(&argv[a][2]);
				if (startfile <= 0) mydie("Invalid starting line number (must be integer >0)!\n");
		}
		if (strstr(argv[a], "-e"))
		{
				endfile = atoi(&argv[a][2]);
				if (endfile <= 0) mydie("Invalid ending line number (must be integer >0)!\n");
		}

	}

	printf("Taking paths of input images from lines %d to %d of file %s...\n", startfile, endfile, listinputimagesfilename);

	// First, making a list of input pictures....
# define MAXNBFILENAMES 11000
	int NBFILENAMES = 0;
	char filenames[MAXNBFILENAMES][200];
	f = fopen(listinputimagesfilename, "r");
	char* thrw; 
	for (int i=0; i < MAXNBFILENAMES; i++)
	{
		thrw = fgets(fn, 199, f);
		if (!thrw) break;
		if (i >= endfile) break;
		if (i+1 < startfile) continue;
		strncpy(filenames[NBFILENAMES], fn, 199);
		filenames[NBFILENAMES][strlen(filenames[NBFILENAMES])-1] = 0; // remove trailing newline
		puts(filenames[NBFILENAMES]);
		NBFILENAMES ++;  
	}
/*	for (int i=0; i< NBFILENAMES; i++)
	{
		thrw = fgets(filenames[i], 99, f);
		filenames[i][strlen(filenames[i])-1] = 0; // remove trailing newline
	}*/
	fclose(f);
	printf("%d files read (first is %s, last is %s)\n", NBFILENAMES, filenames[0], filenames[NBFILENAMES-1]);


	// OK, let's build the layers of the network before we build the connections:


	//	loadpic("/home/thomas/images/natural_images/img_2.pgm.raw");  
	//loadpic("/home/thomas/images/hungetal2005/cph00003.mat.raw");  

	//readProts("protsSort.dat");
	readProts("../prots.dat");

	
	buildNetwork();

	// Now the network is fully ready and operational.
	// Let's apply it to some pics!


	for (int numpic=0; numpic < NBFILENAMES; numpic++)
	{
	clock_t TimeStart = clock();
	printf("Starting chrono...\n");
	

	for (int n=0; n < MAXNEUR; n++) v[n] = 0; // you never know.
	
	loadpic(filenames[numpic]);

		// First, let's run the S1 and C1 layers!

		for (int l = first[S1]; l < first[S1] + nblayersgr[S1]; l++)
		{
		//	printf("Running S1 layer %d..\n", l);
			runS1Layer(l);
		}

		printf("S1 layers run!\n"); fflush(stdout);
		for (int l = first[C1]; l < first[C1] + nblayersgr[C1]; l++)
			runLayerC(l);
		printf("C1 layers run!\n"); fflush(stdout);

		

	//	for (int l = first[S2]; l < first[S2] + nblayersgr[S2]; l++)
	//		runLayerS_ndp(l);
	//	printf("S2 layers run!\n"); fflush(stdout);
		
//		for (int l = first[S2b_5]; l < first[S2b_5] + nblayersgr[S2b_5]; l++)
//			runLayerS_ndp(l);
//		for (int l = first[S2b_7]; l < first[S2b_7] + nblayersgr[S2b_7]; l++)
//			runLayerS_ndp(l);
		for (int l = first[S2b_9]; l < first[S2b_9] + nblayersgr[S2b_9]; l++)
			runLayerS_ndp(l);
//		for (int l = first[S2b_11]; l < first[S2b_11] + nblayersgr[S2b_11]; l++)
//			runLayerS_ndp(l);
//		for (int l = first[S2b_13]; l < first[S2b_13] + nblayersgr[S2b_13]; l++)
//			runLayerS_ndp(l);
		printf("S2b layers run!\n"); fflush(stdout);

		// The saveC2 function actually calculates the C2 outputs
		//sprintf(fn, "./out/C2out_%s_%s.out", "ndp", rindex(filenames[numpic],'/')+1); saveC2(fn);
		sprintf(fn, "./out/C2bout_%s_%s.out", "ndp", rindex(filenames[numpic],'/')+1); saveC2b(fn);
//		sprintf(fn, "./out/C2bout_%s_%s.out", "ndp", rindex(filenames[numpic],'/')+1); saveC2bMean(fn);
//		sprintf(fn, "./out/S2bout_%s_%s.out", "ndp", rindex(filenames[numpic],'/')+1); saveS2b(fn);
		sprintf(fn, "./out/S2bout_%s_%s.out.short", "ndp", rindex(filenames[numpic],'/')+1); saveS2bShort(fn);
	clock_t TimeEnd= clock();
	printf("Time elapsed since chrono started: %f secs.\n", (double)(TimeEnd-TimeStart) / (double)CLOCKS_PER_SEC);
	}	




	//for (int nl=first[S2]; nl < nblayersgr[S2] ; nl++)	
	//for (int nl=0; nl < nblayers ; nl++)	
	//int nl = first[S2b_7];
	int nl = numlayer[S2][0][nbprots[S2] - 1];
//	int nl = numlayer[S2][0][1];
	{
		char s[50]; sprintf(s, "outputlayer%d.txt", nl); 
				printf("Writing output of layer %d, start %d, size %d (xsize %d, ysize %d)..\n", nl, start[nl], size[nl], xsize[nl], ysize[nl]);
		f = fopen(s, "w");
		for (int j=0; j < ysize[nl]; j++)
		{
			for (int i=0; i < xsize[nl]; i++)
			{
				int tgtn = start[nl] + j * xsize[nl] + i;
				fprintf(f, "%.15f ", v[tgtn]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}


	/*f = fopen("protlast.txt", "w");
	for (int i=0; i < MAXAFF; i++)
	{
		int k=nbprots[S2] - 1;
		fprintf(f,"i%d j%d sp%d w%.4f |\n ", (int)StoredProts[S2][k][i*4+1],(int)StoredProts[S2][k][i*4+2],(int)StoredProts[S2][k][i*4+0],   StoredProts[S2][k][i*4+3]); 
	}
	fclose(f);
	f = fopen("prot7.txt", "w");
	for (int i=0; i < MAXAFF; i++)
	{
		int k=7;
		fprintf(f,"i%d j%d sp%d w%.4f |\n ", (int)StoredProts[S2][k][i*4+1],(int)StoredProts[S2][k][i*4+2],(int)StoredProts[S2][k][i*4+0],   StoredProts[S2][k][i*4+3]); 
	}
	fclose(f);*/

	printInfo();
	printf("Nb affs(1st S2 layer) = %d; Nb affs(1st S3 layer) = %d.\n", nbaff[first[S2]], nbaff[first[S3]]);

	return 0;
}
