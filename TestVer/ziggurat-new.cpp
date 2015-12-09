#include "ziggurat-new.h"


//< helper functions
void rho(RR* res, RR* x, RR* sigma)
{
	*res = exp(-power(*x / *sigma, 2)/to_RR(2));
}

RR rho(RR x, RR sigma)
{
	RR res;
	rho(&res, &x, &sigma);
	return res;
}

RR rho(ZZ x, RR sigma)
{
	return rho(to_RR(x), sigma);
}

void set_rhovals(ZZ* rhovals, int* omega, long* tailcut, RR* sigma)
{
	RR yfactor = power2_RR(*omega);
	if (rhovals)
		delete[] rhovals;
	rhovals = new ZZ[*tailcut+1];
	for (int i = 0; i <= *tailcut; i++)
	{
		rhovals[i] = TruncToZZ(yfactor * rho(to_RR(i), *sigma));
	}
}


//< sample a ZZ in {0,...,m-1}
void rejSample(ZZ* res, ZZ* m)
{
	*res = RandomBnd(*m);
}

void rejSample(RR* res, ZZ* m)
{
	*res = to_RR(RandomBnd(*m));
}

ZZ rejSampleZZ(ZZ m)
{
	ZZ res;
	rejSample(&res, &m);
	return res;
}

RR rejSampleRR(ZZ m)
{
	ZZ res;
	rejSample(&res, &m);
	return to_RR(res);
}

void rejSample(int* res, int* m)
{
	*res = to_int(RandomBnd(to_ZZ(*m)));
}

int rejSample(int m)
{
	int res;
	rejSample(&res, &m);
	return res;
}


//< needed for new variant of Ziggurat
void straightLine(RR *res, int* i, ZZ* x, ZZ* rectxs, RR* rectys)
{
	if (*i == 1) *res = to_RR(rectys[*i]-1)/to_RR(rectxs[*i]-rectxs[*i-1])*to_RR(*x-rectxs[*i]);
	else *res = to_RR(rectys[*i]-rectys[*i-1])/to_RR(rectxs[*i]-rectxs[*i-1])*to_RR(*x-rectxs[*i]);
}

RR straightLine(int i, ZZ x, ZZ* rectxs, RR* rectys)
{
	if (i == 1) return (rectys[i]-1)/to_RR(rectxs[i]-rectxs[i-1])*to_RR(x-rectxs[i]);
	return (rectys[i]-rectys[i-1])/to_RR(rectxs[i]-rectxs[i-1])*to_RR(x-rectxs[i]);
}


/************************************\
|             ZIGGURAT               |
\************************************/
void ziggurat(ZZ* res, int* m, ZZ* rectxs, RR* rectys, RR* sigma, int* omega)
{
	int i; //< rectangle i <- {1,...,m}
	int s = 1 - 2 * rejSample(2); //< sign <- +/- 1;
	
	ZZ xurb, yub;
	RR y, yfactor = power2_RR(*omega);
	
	while (true)
	{
		 rejSample(&i, m);i++; //< sample rectangle uniformly
		 
		 xurb = rectxs[i] + 1;
		 rejSample(res, &xurb); //< res <- [0, \floor{x_{i}}]
		
		if (*res != 0 && *res <= rectxs[i-1]) break;
		else
		{
			// the case x=0 is special due to 0=-0, therefore we have to
			// halve the prob. for 0 which results in 1/2; so with p=1/2
			// p gets accepted and with p=1/2 rejected -> "coin toss"
			if (*res == 0) { if (rejSample(2) == 0) break; }
			
			else {
				yub = TruncToZZ(yfactor * (rectys[i-1] - rectys[i])) + 1;
				rejSample(&y, &yub);
					// y <- [0,\floor{2^{\omega} (\rho_{\sigma}(x_{i-1}) - \rho_{\sigma}(x_{i}))}]
					// for implementation +1 at end of yub because of definition of rejSample
					
				if (y <= yfactor * (rho(*res, *sigma) - rectys[i])) break;
			}
		}
	}
	*res *= s;
}


ZZ ziggurat(int m, ZZ* rectxs, RR* rectys, RR sigma, int omega)
{
	ZZ res;
	ziggurat(&res, &m, rectxs, rectys, &sigma, &omega);
	return res;
}

/************************************\
|             ZIGGURAT               |
|       with precomputation          |
\************************************/
void ziggurat_precomp(ZZ* res, int* m, ZZ* rectxs, RR* rectys, ZZ* rhovals, RR* sigma, int* omega)
{
	int i; //< rectangle i <- {1,...,m}
	int s = 1 - 2 * rejSample(2); //< sign <- +/- 1;
	
	ZZ xurb, yub;
	RR y, yfactor = power2_RR(*omega);
	
	while (true)
	{
		 rejSample(&i, m);i++; //< sample rectangle uniformly
		 
		 xurb = rectxs[i] + 1;
		 rejSample(res, &xurb); //< res <- [0, \floor{x_{i}}]
		
		if (*res != 0 && *res <= rectxs[i-1]) break;
		else
		{
			// the case x=0 is special due to 0=-0, therefore we have to
			// halve the prob. for 0 which results in 1/2; so with p=1/2
			// p gets accepted and with p=1/2 rejected -> "coin toss"
			if (*res == 0) { if (rejSample(2) == 0) break; }
			
			else {
				yub = TruncToZZ(yfactor * (rectys[i-1] - rectys[i])) + 1;
				rejSample(&y, &yub);
					// y <- [0,\floor{2^{\omega} (\rho_{\sigma}(x_{i-1}) - \rho_{\sigma}(x_{i}))}]
					// for implementation +1 at end of yub because of definition of rejSample
					
				if (y <= to_RR(rhovals[to_int(*res)]) - yfactor * rectys[i]) break;
			}
		}
	}
	*res *= s;
}


/************************************\
|             ZIGGURAT               |
|       with precomputation          |
\************************************/
void ziggurat_new(ZZ* res, int* m, ZZ* rectxs, RR* rectys, RR* sigma, int* omega)
{
	int i; //< rectangle i <- {1,...,m}
	int s = 1 - 2 * rejSample(2); //< sign <- +/- 1;
	//cout << "i:" << i << " sign:" << sign << endl;
	
	ZZ xurb, yub;
	RR y, sl, yfactor = power2_RR(*omega);
	
	//int finish = 0;
	
	while (true)
	{
		 rejSample(&i, m);i++; //< sample rectangle uniformly
		 
		 xurb = rectxs[i] + 1;
		 rejSample(res, &xurb); //< res <- [0, \floor{x_{i}}]
		
		if (*res != 0 && *res <= rectxs[i-1]) break;
		else
		{
			// the case x=0 is special due to 0=-0, therefore we have to
			// halve the prob. for 0 which results in 1/2; so with p=1/2
			// p gets accepted and with p=1/2 rejected -> "coin toss"
			if (*res == 0) { if (rejSample(2) == 0) break; }
			else {
				yub = TruncToZZ(yfactor * (rectys[i-1] - rectys[i])) + 1;
				rejSample(&y, &yub);
					// y <- [0,\floor{2^{\omega} (\rho_{\sigma}(x_{i-1}) - \rho_{\sigma}(x_{i}))}]
					// for implementation +1 at end of yub because of definition of rejSample
				
				if (*sigma <= to_RR(rectxs[i-1]))
				{
					straightLine(&sl, &i, res, rectxs, rectys); sl *= yfactor;
					if (y >= sl || y > yfactor * (rho(*res, *sigma) - rectys[i])) continue;
					else break;
				}
				else if (*sigma >= to_RR(rectxs[i]+1))
				{
					straightLine(&sl, &i, res, rectxs, rectys); sl *= yfactor;
					if (y <= sl || y <= yfactor * (rho(*res, *sigma) - rectys[i])) break;
				}
				else { if (y <= yfactor * (rho(*res, *sigma) - rectys[i])) break; }
			}
		}
	}
	*res *= s;
}


/************************************\
|             ZIGGURAT               |
|            NEW VARIANT             |
|        with precomputation         |
\************************************/
void ziggurat_new_precomp(ZZ* res, int* m, ZZ* rectxs, RR* rectys, ZZ* rhovals, RR* sigma, int* omega)
{
	int i; //< rectangle i <- {1,...,m}
	int x;
	int s = 1 - 2 * rejSample(2); //< sign <- +/- 1;
	
	ZZ xurb, yub;
	RR y, sl, yfactor = power2_RR(*omega);
	
	while (true)
	{
		 rejSample(&i, m);i++; //< sample rectangle uniformly
		 
		 xurb = rectxs[i] + 1;
		 rejSample(res, &xurb); //< res <- [0, \floor{x_{i}}]
		 x = to_int(*res);
		
		if (*res != 0 && *res <= rectxs[i-1]) break;
		else
		{
			// the case x=0 is special due to 0=-0, therefore we have to
			// halve the prob. for 0 which results in 1/2; so with p=1/2
			// 0 gets accepted and with p=1/2 rejected -> "coin toss"
			if (*res == 0) { if (rejSample(2) == 0) break; }
			else {
				yub = TruncToZZ(yfactor * (rectys[i-1] - rectys[i])) + 1;
				rejSample(&y, &yub);
					// y <- [0,\floor{2^{\omega} (\rho_{\sigma}(x_{i-1}) - \rho_{\sigma}(x_{i}))}]
					// for implementation +1 at end of yub because of definition of rejSample
				
				if (*sigma <= to_RR(rectxs[i-1]))
				{
					straightLine(&sl, &i, res, rectxs, rectys); sl *= yfactor;
					if (y >= sl || y > to_RR(rhovals[x]) - yfactor * rectys[i]) continue;
					else break;
				}
				else if (*sigma >= to_RR(rectxs[i]+1))
				{
					straightLine(&sl, &i, res, rectxs, rectys); sl *= yfactor;
					if (y <= sl || y <= to_RR(rhovals[x]) - yfactor * rectys[i]) break;
				}
				else { if (y <= to_RR(rhovals[x]) - yfactor * rectys[i]) break; }
			}
		}
	}
	*res *= s;
}
