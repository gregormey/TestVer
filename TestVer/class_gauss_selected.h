//
//  class_gauss_selected.hpp
//  Encryption schemes
//
//  Created by Vanessa Erbenich on 18.10.15.
//  Copyright © 2015 Vanessa Erbenich. All rights reserved.
//

#ifndef class_gauss_selected_h
#define class_gauss_selected_h

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZVec.h>

#include "gauss.h"
#include "gauss_knuth_yao.h"
#include "ziggurat-new.h"
#include "class_gauss_inverse_cdf.h"
#include "class_gauss_ziggurat.h"
#include "class_gauss_rejection.h"
#include "class_gauss_knuth_yao.h"

#include <stdio.h>

class class_gauss_selected{
private:
    const int gauss_scheme = 3;
    //const int m = 8192;
    const int m = 100;
    class_gauss_rejection gauss_rejection;
    class_gauss_inverse_cdf gauss_inverse_cdf;
    class_gauss_knuth_yao gauss_knuth_yao;
    class_gauss_ziggurat gauss_ziggurat;
    
public:
    class_gauss_selected();
    class_gauss_selected(double sigma, long factor, long precision);
    NTL::vec_ZZ sample_vec(int n);
    NTL::ZZX sample_poly(int n);
    NTL::mat_ZZ sample_mat(int n, int m);
};

#endif /* class_gauss_selected_h */
