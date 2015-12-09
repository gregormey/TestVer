//
//  class inverse_cdf.h
//  Encryption schemes
//
//  Created by Vanessa Erbenich on 10.09.15.
//  Copyright (c) 2015 Vanessa Erbenich. All rights reserved.
//

#ifndef __Encryption_schemes__class_inverse_cdf__
#define __Encryption_schemes__class_inverse_cdf__

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZVec.h>

#include <iostream>

class class_gauss_inverse_cdf{
    
private:
    long sigma, range, precision;
    NTL::ZZVec cdf_table;
    
public:
    class_gauss_inverse_cdf();
    class_gauss_inverse_cdf(long sigma, long range, long precision);
    NTL::vec_ZZ sample_vec(int n);
    NTL::ZZX sample_poly(int n);
    NTL::mat_ZZ sample_mat(int n, int m);
    
};

#endif /* defined(__Encryption_schemes__class_inverse_cdf__) */