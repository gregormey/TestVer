#ifndef _ZIGGURAT_NEW_H
#define _ZIGGURAT_NEW_H
#include <fstream>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
NTL_CLIENT

//< helper functions
void rho(RR* res, RR* x, RR* sigma);
RR rho(RR x, RR sigma);
RR rho(ZZ x, RR sigma);
void set_rhovals(ZZ* rhovals, int* omega, long* tailcut, RR* sigma);

//< sample a ZZ in {0,...,m-1}
void rejSample(ZZ* res, ZZ* m);
void rejSample(RR* res, ZZ* m);
ZZ rejSampleZZ(ZZ m);
RR rejSampleRR(ZZ m);
void rejSample(int* res, int* m);
int rejSample(int m);

//< needed for new variant of Ziggurat
void straightLine(RR *res, int* i, ZZ* x, ZZ* rectxs, RR* rectys);
RR straightLine(int i, ZZ x, ZZ* rectxs, RR* rectys);

//< (old) Ziggurat
void ziggurat(ZZ* res, int* m, ZZ* rectxs, RR* rectys, RR* sigma, int* omega);
ZZ ziggurat(int m, ZZ* rectxs, RR* rectys, RR sigma, int omega);
void ziggurat_precomp(ZZ* res, int* m, ZZ* rectxs, RR* rectys, ZZ* rhovals, RR* sigma, int* omega);

//< new variant of Ziggurat
void ziggurat_new(ZZ* res, int* m, ZZ* rectxs, RR* rectys, RR* sigma, int* omega);
void ziggurat_new_precomp(ZZ* res, int* m, ZZ* rectxs, RR* rectys, ZZ* rhovals, RR* sigma, int* omega);

#endif
