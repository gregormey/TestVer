#ifndef _UTIL_H
#define _UTIL_H
#include <NTL/RR.h>
//---Constants
#define PI 3.14159265
//---Functions
long floor_log2(unsigned long x);
NTL::RR log2(const NTL::RR& x);
double log2( double x );
#endif
