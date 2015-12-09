#ifndef _GAUSS_H
#define _GAUSS_H

#include <NTL/ZZX.h>
#include <NTL/RR.h>

/***
 * compute a polynomial of degree degr-1 whose coefficients are ZZs
 * distributed according to a discrete Gaussian distribution.
 * There are 2 algorithms here 1) rejection sampling 2) inversion of CDF
 ***/
/*** discrete Gaussian distribution by rejection sampling.
 * The function can be used with a precomputed table of the Gaussian function or computing
 * values on demand. This is determined by how rho is instantiated
 * output:	res = a polynomial of degree degr-1
 *
 * input:	sigma = standard deviation
 * 			range = range where elements are sampled [-range, range]. It is assumed to be a power of two
 * 			degr = max. degree of resulting polynomials + 1
 *			precision = number of bits of precision used in rejection step
 *			rdm_bits1 = a function that takes an integer k as input and returns a Range_type element
 *						chosen uniformly at random in [0,2^k-1]
 *			rdm_bits2 = a function that takes an integer k as input and returns a ZZ element
 *						chosen uniformly at random in [0,2^k-1]
 *			rho = a function that takes a Range_type element x as input and returns the Gaussian function evaluated at x
 *
 * Usage:	Example using precomputed table. This use is prefered when sigma is small
 *			ZZ n = 16;
 *			double sigma = 8;
 *			ZZ range = to_ZZ(2*sigma);
 *			Rho_precomp_table rho_table;
 *			rho_table.set_table(sigma, 32, to_long(range));
 *			ZZX poly;
 *			gauss_sample_rej(&poly, &n, &range, 32, RandomBits_long, RandomBits_ZZ, &rho_table);
 *
 * 			Example on the spot computation of rho. This use is slower but is prefered when sigma is large.
 *			double sigma = to_double(power2_ZZ(30));
 *			ZZ range = to_ZZ(512*sigma);
 *			Rho_basic rho(sigma, 32);
 *			gauss_sample_rej(&poly, &n, &range, 32, RandomBits_ZZ, RandomBits_ZZ, &rho);
 ***/
template <class Range_type, class Rho_fnc>
void gauss_sample_rej(NTL::ZZX *res,
                      const NTL::ZZ *degr,
                      const NTL::ZZ *range,
                      const long precision,
                      Range_type (*rdm_bits1)(const long),
                      NTL::ZZ (*rdm_bits2)(const long),
                      Rho_fnc *rho)
{
    long log_range = NumBits(*range)-1;
    /*	x <- [-range,..,range)
     y <- [0,..,M-1] where M = 2^{precision}
     z = M*\rho_s(x) */
    Range_type x;
    NTL::ZZ y,z;
    // Initialize result to zero
    *res = 0;
    // for each coefficient
    for (int i = 0; i < *degr; i++)
    {
        do
        {
            x = rdm_bits1(log_range);
            //			x = RandomBnd(*range);
            y = rdm_bits2(precision);
            z = (*rho)(&x);
            //with probability \rho_s(x)\in(0,1], output x, otherwise repeat
        }while(y >= z);
        //set coefficient with random sign
        SetCoeff(*res, i, (2*rdm_bits1(1)-1)*x);
    }
}

/** A function object that computes rho
 * Basic implementation. Computes the function on function call
 */
class Rho_basic
{
private:
    const double sigma;
    /// precompute 2^precision
    const NTL::RR max_precs;
public:
    Rho_basic(const double v_sigma, const long v_precision) :
    sigma(v_sigma),
    max_precs(NTL::power2_RR(v_precision)){}
    NTL::ZZ operator()(NTL::ZZ *x) const
    {
        return  ((*x) == 0) ? NTL::TruncToZZ(max_precs * exp(-power(NTL::to_RR(*x) /NTL::to_RR(sigma),2)/NTL::to_RR(2)) / 2) :
        NTL::TruncToZZ(max_precs * NTL::exp(-power(NTL::to_RR(*x) /NTL::to_RR(sigma),2)/NTL::to_RR(2)));
    }
};
/** A function object that computes rho
 *
 * Precomputes a table with all posible values of rho, and simply access the table on function call.
 */
class Rho_precomp_table
{
private:
    NTL::ZZ *rho_vals;
public:
    Rho_precomp_table();
    Rho_precomp_table(const double v_sigma, const long v_precision, const long range);
    void set_table(const double v_sigma, const long v_precision, const long range);
    //~Rho_precomp_table();
    NTL::ZZ operator()(long *x) const;
};
/** Gaussian sampling by inverting CDF as in Peikert (2011)
 * Usage:
 ZZ n;
 long precision;
 long sigma_rlwe, rlwe_range;
 ZZVec cdf_table;
 construct_CDF_table(&cdf_table, sigma_rlwe, rlwe_range, precision);
 ZZX poly;
 gauss_sample_inv_CDF(&poly, &n, precision, RandomBits_long, RandomBits_ZZ, &cdf_table);
 */
void gauss_sample_inv_CDF(NTL::ZZX *res,
                          const NTL::ZZ *degr,
                          const long precision,
                          long (*rdm_bits1)(const long),
                          NTL::ZZ (*rdm_bits2)(const long),
                          NTL::ZZVec* prob_table);
void construct_CDF_table(NTL::ZZVec* cdf_table,
                         const double sigma,
                         const long range,
                         const long precision);

#endif
