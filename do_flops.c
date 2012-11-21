#include "XSbench_header.h"

// Note that this file needs to be compiled separately with the -O0 flag
// to prevent the compiler from optimizing away the flops. The included
// makefile takes this into account.
void do_flops( void )
{
	double a = 1.33;
	double b = 2.34;
	int i;
	
	for( i = 0; i < EXTRA_FLOPS; i++ )
	{
		a = a * b;
	}
}
