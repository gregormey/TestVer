#include <fstream>
#include <cstring>
#include <NTL/vec_RR.h>
#include "gauss.h"

NTL_CLIENT

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
void gauss_sample_inv_CDF(ZZX *res,
							const ZZ *degr,
							const long precision,
							long (*rdm_bits1)(const long),
							ZZ (*rdm_bits2)(const long),
							ZZVec* cdf_table)
{
	ZZ x;
	// Initialize result to zero
	*res = 0;
	// for each coefficient
	for (int i = 0; i < *degr; i++)
	{
		x = rdm_bits2(precision);
		//binary search x in cdf_table
		long a = 0, //< lower bound
			b = cdf_table->length() - 1; //< upper bound
		long c;
		if(x < (*cdf_table)[a])
		{
			c = a;
		}
		else
		{
			do
			{
				//middle of the interval [a,b]
				c = (a+b)/2;
				if(a+1 >= b)
				{
					c = b;
					break;
				}
				//update a or b
				if(x < (*cdf_table)[c])
				{
					//if x is between c-1 and c we're done
					if((*cdf_table)[c-1] <= x)
						break;
					//otherwise update lower bound a
					b = c;
				}
				//otherwise update upper bound b
				else
					a = c;
			}while(true);
		}
		//random sign
		long sign = power_long((long)(-1), rdm_bits1(1));
		//set coefficient
		SetCoeff(*res, i, sign*c);
	}
}
void construct_CDF_table(ZZVec* cdf_table,
							const double sigma,
							const long range,
							const long precision)
{
	cdf_table->SetSize(range, precision);	
	vec_RR vec;
	vec.SetLength(range);
	vec[0] = 0.5;
	for(long x = 1; x < vec.length(); x++)
	{
		vec[x] = vec[x-1] + exp(-power(to_RR(x) /to_RR(sigma),2)/to_RR(2));
	}
	//normalize
	const RR norm_factor = power2_RR(precision) / vec[vec.length()-1];
	for(long x = 0; x < vec.length(); x++)
	{
		(*cdf_table)[x] = to_ZZ(vec[x] * norm_factor);
	}
}
