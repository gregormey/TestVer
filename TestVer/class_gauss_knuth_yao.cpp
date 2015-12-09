//
//  class_gauss_knuth_yao.cpp
//  Encryption schemes
//
//  Created by Vanessa Erbenich on 14.09.15.
//  Copyright (c) 2015 Vanessa Erbenich. All rights reserved.
//

#include "class_gauss_knuth_yao.h"
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include "gauss_knuth_yao.h"


#include <iostream>

using namespace NTL;
using namespace std;

class_gauss_knuth_yao::class_gauss_knuth_yao(){}

class_gauss_knuth_yao::class_gauss_knuth_yao(double sigma, long range, long precision){
    this->sigma = sigma;
    this->range = range;
    this->precision = precision;
    this->gauss_sampler_ky = Gauss_Sampler_KY(this->sigma,this->range,this->precision);
}

vec_ZZ class_gauss_knuth_yao::sample_vec(int n) {
    
    ZZX res;
    ZZ n_ZZ = ZZ(n);
    this->gauss_sampler_ky.sample(&res,&n_ZZ);
    
    vec_ZZ res_vec = VectorCopy(res, to_long(n));
    
    return res_vec;
}

ZZX class_gauss_knuth_yao::sample_poly(int n) {
    
    ZZX res;
    ZZ n_ZZ = ZZ(n+1); // Polynom hat Grad n
    this->gauss_sampler_ky.sample(&res,&n_ZZ);
    
    return res;
}

mat_ZZ class_gauss_knuth_yao::sample_mat(int n, int m) {
    
    ZZX res;
    ZZ nm = ZZ(n*m);
    this->gauss_sampler_ky.sample(&res,&nm);
    
    vec_ZZ res_vec = VectorCopy(res, to_long(n*m));
    mat_ZZ res_mat;
    res_mat.SetDims(n,m);
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            res_mat(i+1,j+1) = res_vec[i*m+j];
    
    return res_mat;
}
