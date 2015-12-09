//
//  class_gauss_knuth_yao.h
//  Encryption schemes
//
//  Created by Vanessa Erbenich on 14.09.15.
//  Copyright (c) 2015 Vanessa Erbenich. All rights reserved.
//

#ifndef __Encryption_schemes__class_gauss_knuth_yao__
#define __Encryption_schemes__class_gauss_knuth_yao__

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZVec.h>
#include "gauss_knuth_yao.h"

#include <iostream>

class class_gauss_knuth_yao{
    
private:
    double sigma;
    long range, precision;
    Gauss_Sampler_KY gauss_sampler_ky;
    
public:
    class_gauss_knuth_yao();
    class_gauss_knuth_yao(double sigma, long range, long precision);
    NTL::vec_ZZ sample_vec(int n);
    NTL::ZZX sample_poly(int n);
    NTL::mat_ZZ sample_mat(int n, int m);
    
};

#endif /* defined(__Encryption_schemes__class_gauss_knuth_yao__) */