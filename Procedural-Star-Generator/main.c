/*
	DISCLAIMER:
	THIS MAIN DRIVER IS NOT PERMANENT, SUBJECT TO CHANGES IN THE FUTURE

	Currently this is tuned for running mass simulations and data collection. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <windows.h>
#include "star_standard.h"

int main(int argc, char* argv[])
{
	unsigned long int seed = (unsigned int)(time(NULL) ^ GetTickCount64());

	printf("SEED: %lu\n\n", seed);
	srand(seed);

	for (int i = 0; i < 3; i++)
	{
		STAR hStar;

		printf("Star %d:\n", i + 1);
		hStar = star_init_default();
		star_generate_random(hStar);
		star_print_details(hStar);
		star_destroy(&hStar);
	}

	return 0;
}