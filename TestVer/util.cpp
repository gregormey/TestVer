#include "util.h"
#include "output.h"

NTL_CLIENT

long floor_log2(unsigned long x)
{
    int l = -1; // floor_log2(0) will return -1
    while (x != 0u)
    {
        x = x >> 1u;
        ++l;
    }
    return l;
}

RR log2(const RR& x)
{
	return log(x) / log(2);
}
#include <math.h>  
// Calculates log2 of number.  
double log2( double n )  
{  
    // log(n)/log(2) is log2.  
    return log( n ) / log( 2 );  
}
