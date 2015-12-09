#ifndef _GAUSS_GNUTH_YAO_H
#define _GAUSS_GNUTH_YAO_H

#include <NTL/ZZX.h>
#include <NTL/RR.h>

class Gauss_Sampler_KY 
{
public:
    Gauss_Sampler_KY();
    Gauss_Sampler_KY(const double v_sigma, const long v_range, const long v_precision);
    //~Gauss_Sampler_KY();
	//
    void sample(NTL::ZZX *poly, const NTL::ZZ *degr);
    void sample(NTL::ZZ *res);
private:
    long precision; // precision
    long *n_inter; // number of internal nodes at each level
    long **leaves; // leaves at each level
};

#endif
