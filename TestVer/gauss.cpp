#include <fstream>
#include <cstring>
#include "gauss.h"

NTL_CLIENT

//-----------------------------------------------------------
//-------------- Rho_precomp_table (see class def in gauss.h)
//-----------------------------------------------------------
Rho_precomp_table::Rho_precomp_table()
{
	rho_vals = NULL;
}
Rho_precomp_table::Rho_precomp_table(const double v_sigma, const long v_precision, const long range)
{		
	// precompute 2^precision
	const RR max_precs = power2_RR(v_precision);
	// request memory to store table of values 
	rho_vals = new ZZ[range];
	for (long x = 0; x < range; x++)
	{
		rho_vals[x] = TruncToZZ(max_precs * exp(-power(to_RR(x) /to_RR(v_sigma),2)/to_RR(2)));
	}
	rho_vals[0] = rho_vals[0] / 2;
}
void Rho_precomp_table::set_table(const double v_sigma, const long v_precision, const long range)
{
	if(rho_vals)
		delete [] rho_vals;
	// precompute 2^precision
	const RR max_precs = power2_RR(v_precision);
	// request memory to store table of values 
	rho_vals = new ZZ[range];
	for (long x = 0; x < range; x++)
	{
		rho_vals[x] = TruncToZZ(max_precs * exp(-power(to_RR(x) /to_RR(v_sigma),2)/to_RR(2)));
	}
	rho_vals[0] = rho_vals[0] / 2;
}
/*Rho_precomp_table::~Rho_precomp_table()
{
	delete [] rho_vals;
}*/
ZZ Rho_precomp_table::operator()(long *x) const
{
	return  rho_vals[*x];
}
