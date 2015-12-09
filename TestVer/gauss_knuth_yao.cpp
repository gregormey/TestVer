#include <bitset>
#include "util.h"
#include "gauss_knuth_yao.h"

NTL_CLIENT
using namespace std;

Gauss_Sampler_KY::Gauss_Sampler_KY()
{
}
Gauss_Sampler_KY::Gauss_Sampler_KY(const double v_sigma, const long v_range, const long v_precision)
{
    this->precision = v_precision; // it is assumed that this is the precision of the NTL::RR arithmetic
    //compute Gaussian values
    RR *val = new RR[v_range]; // value of Gaussian function as a floating point
    RR sum_gauss = to_RR(0.); // sum of Gaussian values
	for (long x = 0; x < v_range; x++)
	{
        val[x] = x==0 ? exp(-power(to_RR(x) /to_RR(v_sigma),2)/to_RR(2))/2
                        : exp(-power(to_RR(x) /to_RR(v_sigma),2)/to_RR(2));
        sum_gauss += val[x];
    }
    // compute probabilities
    long ex; // exponent
    ZZ *mt = new ZZ[v_range]; // mantissa
    ZZ sum_prob = to_ZZ(0); // sum of probabilities
    // --- some debug variables to be removed
    // ZZ count_bits = to_ZZ(0);
    // long n_nz[v_range];
    ///---
	for (long x = 0; x < v_range; x++)
	{
        //compute probabilities
        RR prob = val[x] / sum_gauss;
        //extract bit expansion
        mt[x] = prob.mantissa();//mantissa
        ex = prob.exponent();
        //shift mantissa to make all exponents equal to precision
        if(-ex < this->precision)
            mt[x]<<=(this->precision+ex);
        else if(-ex > this->precision)
            mt[x]>>=(-ex-this->precision);            
        sum_prob += mt[x];
        // --- some debug computations to be removed
        // n_nz[x] = 0;
        // bool flag1 = true;
        // cout<<endl;
        // for (long clm = this->precision-1; clm >= 0 ; clm--)
        // {
            // cout<<bit(mt[x], clm);
            // count_bits += bit(mt[x], clm);
            // if(flag1 && bit(mt[x], clm))
            // {
                // n_nz[x] = clm;
                // flag1 = false;
            // }
        // }
        // ---
	}
    // --- some debug printing to be removed
    // cout<<endl<<endl;
    // for (long clm = this->precision-1; clm >= 0 ; clm--)
        // cout<<bit(sum_prob, clm);
    // cout<<endl<<endl<<"count_bits="<<count_bits;
    // cout<<endl<<"n_nz=[";
    // for (long x = 0; x < v_range; x++)
        // cout<<n_nz[x]<<",";
    // ---
    // fix the sum of the probabilities to be equal to 1
    ZZ diff = power(to_ZZ(2),this->precision)-sum_prob;
    for (long x = 0; x < diff; x++)
    {
        mt[v_range-1-x]++;
    }
    //request memory for graph description
    this->n_inter = new long[this->precision];
    this->leaves = new long*[this->precision];
    for(long clm = 0; clm < this->precision; clm++)
        leaves[clm] = new long[v_range];
    //construct graph table    
    this->n_inter[0] = 0; //number of internal n_nodess
    long n_nodes = 0; // numnber of nodes
    for (long clm = 0; clm < this->precision; clm++) //start filling last column (lsb)(clm=0) first 
    {
        //compute number of internal leaves
        this->n_inter[clm] = n_nodes/2;
        //add leaves
        long n_leaves = 0; //number of leaves
        for (long x = 0; x < v_range; x++)
        {            
            if (bit(mt[x], clm))
            {
                this->leaves[clm][n_leaves++] = x;
            }
        }
        //compute number of nodes for this column
        n_nodes = n_leaves + this->n_inter[clm];
    }
}
/*Gauss_Sampler_KY::~Gauss_Sampler_KY()
{
    for(long clm = 0; clm < this->precision; clm++)
        delete leaves[clm];
    delete leaves;
    delete n_inter;
}*/
void Gauss_Sampler_KY::sample(NTL::ZZX *poly, const NTL::ZZ *degr)
{    
	for (int i = 0; i < *degr; i++)
	{
        long node = RandomBits_long(1),
            clm = this->precision-1; //start at first column (clm=this->precision-1)
        while(true)
        {
            if(node < this->n_inter[clm])
                node = 2*node + RandomBits_long(1);
            else
                break;
            clm--;
        }
        SetCoeff(*poly, i, (2*RandomBits_long(1)-1)*this->leaves[clm][node - this->n_inter[clm]]);
    }
}

void Gauss_Sampler_KY::sample(NTL::ZZ *res)
{    
	long node = RandomBits_long(1),
            clm = this->precision-1; //start at first column (clm=this->precision-1)
        while(true)
        {
            if(node < this->n_inter[clm])
                node = 2*node + RandomBits_long(1);
            else
                break;
            clm--;
        }
        *res = (2*RandomBits_long(1)-1)*this->leaves[clm][node - this->n_inter[clm]];
}
