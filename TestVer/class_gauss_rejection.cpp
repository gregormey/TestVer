//
//  class_gauss_rejection.cpp
//  Encryption schemes
//
//  Created by Vanessa Erbenich on 14.09.15.
//  Copyright (c) 2015 Vanessa Erbenich. All rights reserved.
//

#include "class_gauss_rejection.h"
#include "gauss.h"
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>


#include <iostream>

using namespace NTL;
using namespace std;

class_gauss_rejection::class_gauss_rejection(){}

class_gauss_rejection::class_gauss_rejection(double sigma, ZZ range, long precision){
    this->sigma = sigma;
    this->range = range;
    this->precision = precision;
    rho_table.set_table(sigma, precision, to_long(range));
    
    /*for(long x=0;x<100;x++)
        cout << "rho_table(" << x <<")=" << rho_table(&x) << endl;*/
}

vec_ZZ class_gauss_rejection::sample_vec(int n) {
    ZZX res;
    ZZ n_ZZ = ZZ(n);
    
    gauss_sample_rej(&res, &n_ZZ, &range, precision, RandomBits_long, RandomBits_ZZ, &rho_table);
    
    vec_ZZ res_vec = VectorCopy(res, to_long(n));
    
    return res_vec;
}

ZZX class_gauss_rejection::sample_poly(int n) {
    ZZX res;
    ZZ n_ZZ = ZZ(n+1); // Polynom hat Grad n
    
    gauss_sample_rej(&res, &n_ZZ, &range, precision, RandomBits_long, RandomBits_ZZ, &rho_table);
    
    return res;
}

mat_ZZ class_gauss_rejection::sample_mat(int n, int m) {
    
    ZZX res;
    ZZ nm = ZZ(n*m);
    gauss_sample_rej(&res, &nm, &range, precision, RandomBits_long, RandomBits_ZZ, &rho_table);
    
    vec_ZZ res_vec = VectorCopy(res, to_long(n*m));
    mat_ZZ res_mat;
    res_mat.SetDims(n,m);
    for (int i=0; i<n; i++)
        for (int j=0; j<m; j++)
            res_mat(i+1,j+1) = res_vec[i*m+j];
    return res_mat;
}

