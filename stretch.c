
/*
 * 	File:		stretch.c
 * 	Programmer:	BForsyth
 * 	Date:		20120928
 * 	Description:	Implement contrast stretch
 * 			Usage ./stretch A.dat B.dat fire_mask.dat 
 * 			Takes either one or 2 input files
 * 			(uint16_t) and a mask file (uint16_t, zero if not fire)
 * 			Stretches the file acording to the specified
 * 			sigma and offset.  If two files given
 * 			it creates a false color with A tied to red
 * 			and B tied to green, and high values of B
 * 			tied to blue
 */
#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#define 	MAX_DIGITS	50
#define		K 		6.5	
#define		BUFF_SIZE	16
#define		DSP_RANGE	16000
#define 	SIGMA_A		.020
#define		OFFSET_A	-15//bigger makes it darker
#define		MASKED_SIGMA_A	.025
#define		MASKED_OFFSET_A	 -127
#define		MASKED_SIGMA_B	4	
#define		MASKED_OFFSET_B	40 
#define		SIGMA_B		.60
#define		OFFSET_B	20
//These constants are determined by plotting the min A for each B on GNUPLOT
//then curve fitting


int 	fillBuff (uint16_t buf[], FILE * inFile);
int 	makeHisto(uint16_t buf[], int histo[], int *max, int *min);
int 	clearHisto (int histo[]);
int	negateMask(uint16_t buf[]);
int 	makeMaskedHisto (uint16_t buf[], uint16_t maskBuf[], int histo[],
	       			int *max, int *min);
float	averageHisto (int histo[]);
float 	stdDev(int histo[], float average);
int 	buildLUT (int histo[], float average, float stdDev, 
		int offset, float sigma);
int	contrastStretch(uint16_t buf[], int histo[], int *clippedBlack,
			int *clippedWhite, FILE *outFile);

int main (int argc, char *argv[])
{
	FILE 	*inAFile, *inBFile, *inMaskFile;
	FILE	*outAFile, *outBFile;
	uint16_t   	bufA[BUFF_SIZE];
	uint16_t 	bufB[BUFF_SIZE];
	uint16_t	bufMask[BUFF_SIZE];
	int 	histoA[DSP_RANGE];
	int 	histoB[DSP_RANGE];
	int 	maskedHistoA[DSP_RANGE];
	int 	maskedHistoB[DSP_RANGE]; //negative of A mask
	int 	count = 0;
	int 	maxA = 0;
	int 	maxB = 0;
	int 	maskedMaxA = 0;
	int 	maskedMaxB = 0;
	int 	minA = 60000;
	int 	minB = 60000;
	int 	maskedMinA = 60000;
	int 	maskedMinB = 60000;
	int 	blackPointA = 0;
	int 	blackPointB = 0;
	int 	whitePointA = 0;
	int 	whitePointB = 0;
	int 	clippedBlackA = 0;
	int 	clippedWhiteA = 0;
	int 	clippedBlackB = 0;
	int 	clippedWhiteB = 0;
	float 	averageA = 0;
	float 	maskedAverageA = 0;
	float 	maskedAverageB = 0;
	float 	stdDevA = 0;
	float 	maskedStdDevA = 0;
	float 	maskedStdDevB = 0;
	float 	averageB = 0;
	float 	stdDevB = 0;
	int	offsetA = OFFSET_A;
	int	offsetB = OFFSET_B;
	int 	maskedOffsetA = MASKED_OFFSET_A;
	int 	maskedOffsetB = MASKED_OFFSET_B;
	float	sigmaA = SIGMA_A;
	float	sigmaB = SIGMA_B;
	float	maskedSigmaA = MASKED_SIGMA_A;
	float	maskedSigmaB = MASKED_SIGMA_B;


	if (argc != 4)	{
		printf("\nError: need infile and outfile");
		exit(1);
	}
	
	outAFile = fopen("outA.dat", "w");
	if (outAFile == 0) {
		printf("\nError: could not open shit");
		exit (1);
	}
	outBFile = fopen("outB.dat", "w");
	if (outBFile == 0) {
		printf("\nError: could not open shit");
		exit (1);
	}
	//inAFile = fopen("rawB.dat", "r");

	inAFile = fopen(argv[1], "r");
	if (inAFile == 0)	{
		printf("\nError: could not open inAfile");
		exit(1);
	}
printf("\n opening the file");	
	
	inBFile = fopen(argv[2], "r");
	if (inBFile == 0)	{
		printf("\nError: could not open inBfile");
		exit(1);
	}
	
//printf("\n opening the last file");	
	inMaskFile = fopen(argv[3], "r");
//printf("\n opened it ");	

	if (inMaskFile == 0)	{
		printf("\nError: could not open outfile");
		exit(1);
	}
	printf("\ngot the right arguments all files opened");
	

//	int k = 0;
	//while ((count = fillBuff(bufA, inAFile)) == fillBuff(bufB, inBFile) &&
	//		(count == fillBuff(bufMask, inMaskFile)) &&
	//	       count != 0)
///////////////////////////////////////////////////dddddddddddddddddddd

printf("\nstarting hard shit");
	clearHisto(histoA);
	clearHisto(histoB);
	clearHisto(maskedHistoA);
	clearHisto(maskedHistoB);
	while (((count = fillBuff(bufA, inAFile))  != 0 ) &&
			(fillBuff(bufMask, inMaskFile) != 0))
	//while ((count = fillBuff(bufA, inAFile)) == fillBuff(bufB, inBFile) &&
	//		(count == fillBuff(bufMask, inMaskFile)) &&
	//	       count != 0)
	{	
		makeHisto (bufA, histoA, &maxA, &minA);
		makeMaskedHisto(bufA, bufMask, 
				maskedHistoA, &maskedMaxA, &maskedMinA);

			
	}
	rewind(inMaskFile);
	while (((count = fillBuff(bufB, inBFile))  != 0 ) &&
			(fillBuff(bufMask, inMaskFile) != 0))
	//while ((count = fillBuff(bufA, inAFile)) == fillBuff(bufB, inBFile) &&
	//		(count == fillBuff(bufMask, inMaskFile)) &&
	//	       count != 0)
	{	
		makeHisto (bufB, histoB, &maxB, &minB);
		negateMask(bufMask); //mask fire
		makeMaskedHisto(bufB, bufMask, 
				maskedHistoB, &maskedMaxB, &maskedMinB);
		//makeHisto (bufB, histoB, &maxB, &minB);
		//printf("\nread = %d", bufA[2]);
		//printf("\n maxA = %d, minA = %d", maxA, minA);
		//makeHisto (bufB, histoB, &maxB, &minB);
		//makeMaskedHisto(bufA, bufMask, maskedHistoA, &maxA, &minA);
		//makeMaskedHisto(bufB, bufMask, maskedHistoB, &maxB, &minB);
	
	}
	count = 0;
	int i;
	for (i=0;i<DSP_RANGE;i++)
		count += histoA[i];
	printf("\nread %d pixels", count);
	averageA = averageHisto(histoA);
	maskedAverageA = averageHisto(maskedHistoA);
	maskedAverageB = averageHisto(maskedHistoB);
	averageB = averageHisto(histoB);
	stdDevA = stdDev(histoA, averageA);
	maskedStdDevA = stdDev(maskedHistoA, maskedAverageA);
	maskedStdDevB = stdDev(maskedHistoB, maskedAverageB);
	stdDevB = stdDev(histoB, averageB);
	printf("\nmaxA = %d minA = %d", maxA, minA);
	printf("\naverage A = %f", averageA);
	printf("\nstddev A = %f, sigmaA = %f, offsetA = %d",
		       	stdDevA, sigmaA, offsetA);
	printf("\nmaskedmaxA = %d maskedminA = %d", maskedMaxA, maskedMinA);
	printf("\nmasked average A = %f", maskedAverageA);
	printf("\nsmasked tddev A = %f, masked sigma = %f, masked offset = %d",
		       	maskedStdDevA,	maskedSigmaA, maskedOffsetA);
	printf("\nmaxB = %d minB = %d", maxB, minB);
	printf("\naverage B = %f", averageB);
	printf("\nstddev B = %f, sigmab = %f, offsetb = %d",
		       	stdDevB, sigmaB, offsetB);
	printf("\nmaskedmaxB = %d maskedminB = %d", maskedMaxB, maskedMinB);
	printf("\nmasked average B = %f", maskedAverageB);
	printf("\nsmasked tddev B = %f, masked sigma = %f, masked offset = %d",
		       	maskedStdDevB,	maskedSigmaB, maskedOffsetB);
//	buildLUT(histoA, averageA, stdDevA, offsetA, sigmaA);
	buildLUT(maskedHistoA, maskedAverageA, maskedStdDevA,
		       	maskedOffsetA, maskedSigmaA);
	//buildLUT(histoB, averageB, stdDevB, offsetB, sigmaB);
	buildLUT(maskedHistoB, maskedAverageB, maskedStdDevB,
		       	maskedOffsetB, maskedSigmaB);
	averageA = averageHisto(histoA);
	maskedAverageA = averageHisto(maskedHistoA);
	averageB = averageHisto(histoB);
	maskedAverageB = averageHisto(maskedHistoB);
	printf("\naverage A pixel value = %f", averageA);
	printf("\nmasked average A pixel value = %f", maskedAverageA);
	printf("\naverage B pixel value = %f", averageB);
//	printf("\nblackpoint A = %d, whitepoint A = %d", blackPointA,
		       					//whitePointA);	
	//stretch the data
	rewind(inAFile);
	rewind(inBFile);
//	rewind(inMaskFile);
	while (((count = fillBuff(bufA, inAFile)) != 0 ) &&
		(fillBuff(bufB, inBFile) !=0))
	{	
		//contrastStretch(bufA, histoA,&clippedBlackA,
		//	       &clippedWhiteA,	outAFile);
		contrastStretch(bufA, maskedHistoA,&clippedBlackA,
			       &clippedWhiteA,	outAFile);
		contrastStretch(bufB, maskedHistoB,&clippedBlackB,
			       &clippedWhiteB,	outBFile);
	}
		printf("\nContrast stretch complete");
	printf("\nclippedblackA %d, clippedwhiteA %d",
			clippedBlackA, clippedWhiteA);
	printf("\nclippedblackB %d, clippedwhiteB %d",
			clippedBlackB, clippedWhiteB);
/*
	fclose(inAFile);
	fclose(inBFile);
	fclose(outAFile);
	*/
	return 0;
}

//fillBuff() reads BUFF_SIZE numbers from inFile into buf[] returns
//# of bytes read
int	fillBuff (uint16_t buf[], FILE * inFile)
{
	return fread(buf, sizeof(uint16_t), BUFF_SIZE, inFile);
}	

int 	makeHisto(uint16_t buf[], int histo[], int *max, int *min)
{
	int i;
	for (i = 0; i< BUFF_SIZE; i++)
	{
		if ((buf[i] < 0) || (buf[i] > 65535))
		{
		       printf("\nout of range input %d", buf[i]);	
		       exit(1);
		}
	       histo[(buf[i])]++;
	       if (*min > buf[i])
		       *min = buf[i];
	       if (*max < buf[i])
		       *max = buf[i];
	}
}
int 	makeMaskedHisto (uint16_t buf[], uint16_t maskBuf[], int histo[],
	       			int *max, int *min)
{ int count = 0;
	int i;
	for (i = 0; i< BUFF_SIZE; i++)
	{
		if ((buf[i] < 0) || (buf[i] > 65535))
		{
		       printf("\nout of range input %d", buf[i]);	
		       exit(1);
		}
		if (maskBuf[i] == 0)
		{
			continue;
		}
	       histo[(buf[i])]++;
	       if (*min > buf[i])
		       *min = buf[i];
	       if (*max < buf[i])
		       *max = buf[i];
	}
	//printf("\nfound Fire = %d", count);
	return 0;

}


int 	clearHisto (int histo[])
{
	int i;
	for (i = 0; i <DSP_RANGE; i++)
		histo[i] = 0;
	return 0;
}
		

float	averageHisto (int histo[])
{
	int i;
	int average = 0;
	float sum = 0;
	float count = 0;
	for (i=0;i<DSP_RANGE; i++)
	{
		sum += (float)histo[i] * (float)i;//this will overflow an int
		count += (float)histo[i];	
	}
	return sum/count;
}

float 	stdDev(int histo[], float average)
{
	int i;
	float count = 0;
	float numerator = 0;
	for (i=0;i<DSP_RANGE;i++)
	{
		numerator += pow(i - average, 2) * (float)histo[i];
		count += (float)histo[i];
	}
	return sqrt((float)numerator / (float)count);
}


int 	buildLUT (int histo[], float average, float stdDev,
	       	int offset, float sigma)
{
	int i;
	float b;
	//float sigma = SIGMA_A;
	float slope;
	b = (float )(65535)/(float)2 * ((float)1 
			- (average + offset)/(sigma * stdDev));
	slope = (float)65535/(2 * sigma * stdDev);

	for (i=0;i<DSP_RANGE;i++)
	{
		if (i <= average - stdDev * sigma + offset)
		{
			histo[i] = 0;
		}
		else if (i >= average + stdDev * sigma + offset)
		{
			histo[i] = 65535;
		}
		else 
		{
		//histo[i] = (int)((float)65535/(float)(2*sigma*stdDev) * i
		//			+ (float)(65535/2)
		//	 	* (1 - (float)average/(float)(sigma*stdDev))
		//		- OFFSET_A);
			histo[i] = (int) (slope * i + b);
			if ((histo[i]<0) || (histo[i]>65535))
			{
				printf("\noutof bounds lut value value %d", histo[i]);
				printf("\nindex i = %d", i);
				exit(1);
			}
		}
	}
}

int	contrastStretch(uint16_t buf[], int histo[],
	       	int *clippedBlack, int *clippedWhite,	FILE *outFile)
{
	int i = 0; 
	int dummy = 0;
	for (i=0;i<BUFF_SIZE;i++)
	{
		dummy = histo[buf[i]];
		buf[i] = (uint16_t)dummy;

		if ((dummy < 0) || (dummy > 65535))
		       printf("\ndummy is out of bounds %d", dummy);	
		if (dummy == 0)
			*clippedBlack = *clippedBlack + 1;
		if (dummy >= 65535)
		{		*clippedWhite = *clippedWhite + 1;
		}
		fwrite(&buf[i], sizeof(uint16_t), 1, outFile);
	}

	return 0;

}

int	negateMask(uint16_t buf[])
{ 
	int i;
	for (i =0;i<BUFF_SIZE;i++)
	{
		if (buf[i] == 0)
			buf[i] = 65535;
		else
			buf[i] = 0;
	}
	return 0;
}
				
//int 	makeMaskedHisto (int buf[], int maskBuf[], int histo[]);
//int 	stdDev(int histo[], int average);
//int 	buildLUT (int histo[], int average, int stdDev);
//int	contrastStretch(int buf[], int histo[], FILE outFile);
/*int 	fireDetect (int bufA[], int bufB[], FILE * outFile, int count)
{
	float compute;
	uint16_t fireval;
	uint16_t fireval16;
	int k = 0;
	for (k = 0; k < count; k++)
	{
		//fireval = K *  bufA[k] -  bufB[k];
//		fireval = fireval + 1000;
		//This implements curve fit
		compute = bufA[k] - A2*pow(bufB[k], 2) - A1*bufB[k] - A0;
//		if (bufB[k] > -4130)
		if (bufB[k] < BMIN)
				fireval = 0;
		else if (bufB[k] > BMAX)
			fireval = bufA[k] + 10000;
		else if (compute > 0)
			fireval = bufA[k]+10000;
		else
			fireval = 0;
		if (max < bufA[k] + 10000)
			max = bufA[k]+10000;
		if (min > bufA[k] + 10000)
			min = bufA[k] + 10000;

		fwrite(&fireval, sizeof(fireval), 1, outFile);
	}
	return k;
}
*/
