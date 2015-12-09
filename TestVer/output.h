#ifndef OUTPUT_H
#define OUTPUT_H

#include<fstream>
using namespace std;

//#define PRINT_DEBUG 1
#define PRINT_TIME 1
//#define PRINT_VERBOSE 1
//#define PRINT_MEMORY 1
//#define PRINT_TIME_SPARSE 1
//#define PRINT_TO_SCREEN 1
//#define PRINT_STATS 1
//#define PRINT_PERFORMANCE

#if defined(PRINT_DEBUG) || defined(PRINT_TIME) || defined(PRINT_TIME_SPARSE) || defined(PRINT_STATS) || defined(PRINT_MEMORY) || defined(PRINT_PERFORMANCE)
    #define PRINT
	extern ofstream fout;
    #include <iomanip>
#endif
#if defined(PRINT_VERBOSE)
	extern ofstream fout_verbose;
    #include <iomanip>
#endif
#if defined(PRINT_TIME_SPARSE) || defined(PRINT_TIME)  || defined(PRINT_STATS) || defined(PRINT_PERFORMANCE) || defined(PRINT_VERBOSE)
	#include <time.h>
	struct timespec diff(struct timespec start, struct timespec end);
	double to_seconds(struct timespec t);
#endif

#if defined(PRINT_MEMORY) || defined(PRINT_PERFORMANCE)
    #include <windows.h>
    #include <stdio.h>
    #include "psapi.h"
#endif

#endif
