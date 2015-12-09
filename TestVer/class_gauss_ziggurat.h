//
//  class_gauss_ziggurat.h
//  Encryption schemes
//
//  Created by Vanessa Erbenich on 11.09.15.
//  Copyright (c) 2015 Vanessa Erbenich. All rights reserved.
//

#ifndef __Encryption_schemes__class_gauss_ziggurat__
#define __Encryption_schemes__class_gauss_ziggurat__

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>
#include "ziggurat-new.h"

class class_gauss_ziggurat{
    
private:
    int m, precision, omega;
    long tailcut;
    RR sigma;
    
    ZZ* rectxs;
    ZZ* rhovals;
    RR* rectys;
    
    NTL::RR rho_new(NTL::RR sigma, NTL::RR x);
    NTL::RR compute_recurrence(NTL::RR* xis, NTL::RR* rhos, NTL::RR c, NTL::RR m, NTL::RR sigma);
    bool create_partitions();
    
public:
    class_gauss_ziggurat();
    class_gauss_ziggurat(int m, int precision, int sigma);
    NTL::vec_ZZ sample_vec(int n);
    NTL::ZZX sample_poly(int n);
    NTL::mat_ZZ sample_mat(int n, int m);
    
};

#endif /* defined(__Encryption_schemes__class_gauss_ziggurat__) */