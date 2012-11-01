#include "XSBench_header.h"

void load_mats( int ** _mats, double ** _concs, int * _num_nucs)
{
	int *mats[12];
	int mats0[] =  { 58, 59, 60, 61, 40, 42, 43, 44, 45, 46, 1, 2, 3, 7,
	                 8, 9, 10, 29, 57, 47, 48, 0, 62, 15, 33, 34, 52, 53,
									 54, 55, 56, 18, 23, 41, -1 }; //fuel
	int mats1[] =  { 63, 64, 65, 66, 67, -1 }; // cladding
	int mats2[] =  { 24, 41, 4, 5, -1 }; // cold borated water
	int mats3[] =  { 24, 41, 4, 5, -1 }; // hot borated water
	int mats4[] =  { 19, 20, 21, 22, 35, 36, 37, 38, 39, 25, 27, 28, 29,
	                 30, 31, 32, 26, 49, 50, 51, 11, 12, 13, 14, 6, 16,
									 17, -1 }; // RPV
	int mats5[] =  { 24, 41, 4, 5, 19, 20, 21, 22, 35, 36, 37, 38, 39, 25,
	                 49, 50, 51, 11, 12, 13, 14, -1 }; // lower radial reflector
	int mats6[] =  { 24, 41, 4, 5, 19, 20, 21, 22, 35, 36, 37, 38, 39, 25,
	                 49, 50, 51, 11, 12, 13, 14, -1 }; // top reflector / plate
	int mats7[] =  { 24, 41, 4, 5, 19, 20, 21, 22, 35, 36, 37, 38, 39, 25,
	                 49, 50, 51, 11, 12, 13, 14, -1 }; // bottom plate
	int mats8[] =  { 24, 41, 4, 5, 19, 20, 21, 22, 35, 36, 37, 38, 39, 25,
	                 49, 50, 51, 11, 12, 13, 14, -1 }; // bottom nozzle
	int mats9[] =  { 24, 41, 4, 5, 19, 20, 21, 22, 35, 36, 37, 38, 39, 25,
	                 49, 50, 51, 11, 12, 13, 14, -1 }; // top nozzle
	int mats10[] = { 24, 41, 4, 5, 63, 64, 65, 66, 67, -1 }; // top of fuel assemblies
	int mats11[] = { 24, 41, 4, 5, 63, 64, 65, 66, 67, -1 }; // bottom of fuel assemblies

	mats[0] = mats0;
	mats[1] = mats1;
	mats[2] = mats2;
	mats[3] = mats3;
	mats[4] = mats4;
	mats[5] = mats5;
	mats[6] = mats6;
	mats[7] = mats7;
	mats[8] = mats8;
	mats[9] = mats9;
	mats[10] = mats10;
	mats[11] = mats11;

	int num_nucs[12];
	for( int i = 0; i < 12; i++ )
		for( int j = 0; ; j++ )
			if( mats[i][j] == -1 )
			{	
				num_nucs[i] = j;
				break;
			}
	
	double ** concs = (double **)malloc( 12 * sizeof( double *) );
	for( int i = 0; i < 12; i++ )
		concs[i] = (double *)malloc( num_nucs[i] * sizeof(double) );
	for( int i = 0; i < 12; i++ )
		for( int j = 0; j < num_nucs[i]; j++ )
			concs[i][j] = (double) rand() / (double) RAND_MAX;

	// So - now we have 3 arrays
	// 1) 2-D array of material compositions, by index in nuclide ID grid
	// 2) 2-D array of nuclide concentrations, corresponding to above array
	// 3) 1-D array of number of nuclides per material.
	_mats = mats;
	_concs = concs;
	_num_nucs = num_nucs;
}

int pick_mat(void)
{
	// I have a nice spreadsheet supporting these numbers. They are
	// the fractions (by volume) of material in the core. Not a 
	// *perfect* approximation of where XS lookups are going to occur,
	// but this will do a good job of biasing the system nonetheless.

	double dist[12];
	dist[0]  = 0.140;	// fuel
	dist[1]  = 0.052;	// cladding
	dist[2]  = 0.275;	// cold, borated water
	dist[3]  = 0.134;	// hot, borated water
	dist[4]  = 0.154;	// RPV
	dist[5]  = 0.064; // Lower, radial reflector
	dist[6]  = 0.066;	// Upper reflector / top plate
	dist[7]  = 0.055;	// bottom plate
	dist[8]  = 0.008;	// bottom nozzle
	dist[9]  = 0.015;	// top nozzle
	dist[10] = 0.025;	// top of fuel assemblies
	dist[11] = 0.013;	// bottom of fuel assemblies
	
	double roll = (double) rand() / (double) RAND_MAX;

	// makes a pick based on the distro
	for( int i = 0; i < 12; i++ )
	{
		double running = 0;
		for( int j = i; j > 0; j-- )
			running += dist[j];
		if( roll < running )
			return i;
	}

	return 0;
}
