//
//  class inverse_cdf.cpp
//  Encryption schemes
//
//  Created by Vanessa Erbenich on 10.09.15.
//  Copyright (c) 2015 Vanessa Erbenich. All rights reserved.
//

#include "class_gauss_inverse_cdf.h"
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include "gauss.h"

#include <iostream>

using namespace NTL;
using namespace std;

class_gauss_inverse_cdf::class_gauss_inverse_cdf(){}

class_gauss_inverse_cdf::class_gauss_inverse_cdf(long sigma, long range, long precision){
    this->sigma = sigma;
    this->range = range;
    this->precision = precision;
    construct_CDF_table(&cdf_table, sigma, range, precision);
}

vec_ZZ class_gauss_inverse_cdf::sample_vec(int n) {
    
    ZZX res;
    ZZ n_ZZ = ZZ(n);
    gauss_sample_inv_CDF(&res, &n_ZZ, precision, RandomBits_long, RandomBits_ZZ, &cdf_table);
    
    vec_ZZ res_vec = VectorCopy(res, to_long(n));
    return res_vec;
}

ZZX class_gauss_inverse_cdf::sample_poly(int n) {
    
    ZZX res;
    ZZ n_ZZ = ZZ(n+1); // Polynom hat Grad n
    gauss_sample_inv_CDF(&res, &n_ZZ, precision, RandomBits_long, RandomBits_ZZ, &cdf_table);
    
    return res;
}

mat_ZZ class_gauss_inverse_cdf::sample_mat(int n, int m) {
    
    ZZX res;
    ZZ nm = ZZ(n*m);
    gauss_sample_inv_CDF(&res, &nm, precision, RandomBits_long, RandomBits_ZZ, &cdf_table);
    
    vec_ZZ res_vec = VectorCopy(res, to_long(n*m));
    mat_ZZ res_mat;
    res_mat.SetDims(n,m);
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            res_mat(i+1,j+1) = res_vec[i*m+j];
    return res_mat;
}
