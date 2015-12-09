//
//  class_gauss_rejection.h
//  Encryption schemes
//
//  Created by Vanessa Erbenich on 14.09.15.
//  Copyright (c) 2015 Vanessa Erbenich. All rights reserved.
//

#ifndef __Encryption_schemes__class_gauss_rejection__
#define __Encryption_schemes__class_gauss_rejection__

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZVec.h>
#include "gauss.h"


#include <iostream>

class class_gauss_rejection{
    
private:
    double sigma;
    long precision;
    NTL::ZZ range;
    Rho_precomp_table rho_table;
    
public:
    class_gauss_rejection();
    class_gauss_rejection(double sigma, NTL::ZZ range, long precision);
    NTL::vec_ZZ sample_vec(int n);
    NTL::ZZX sample_poly(int n);
    NTL::mat_ZZ sample_mat(int n, int m);
    
};

#endif /* defined(__Encryption_schemes__class_gauss_rejection__) */