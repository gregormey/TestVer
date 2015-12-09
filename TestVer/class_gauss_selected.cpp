//
//  class_gauss_selected.cpp
//  Encryption schemes
//
//  Created by Vanessa Erbenich on 18.10.15.
//  Copyright © 2015 Vanessa Erbenich. All rights reserved.
//

#include "class_gauss_selected.h"

class_gauss_selected::class_gauss_selected() {}

class_gauss_selected::class_gauss_selected(double sigma, long factor, long precision){
    if(gauss_scheme==1){
        gauss_rejection = class_gauss_rejection(sigma,to_ZZ(factor*sigma),precision);
    } else if(gauss_scheme==2){
        gauss_inverse_cdf = class_gauss_inverse_cdf(to_long(sigma),factor*sigma,precision);
    } else if(gauss_scheme==3){
        gauss_knuth_yao = class_gauss_knuth_yao(sigma, factor*sigma, precision);
    } else {
        gauss_ziggurat = class_gauss_ziggurat(m, to_int(precision), to_int(sigma));
    }
}

vec_ZZ class_gauss_selected::sample_vec(int n) {
    if(gauss_scheme==1){
        return gauss_rejection.sample_vec(n);
    } else if(gauss_scheme==2){
        return gauss_inverse_cdf.sample_vec(n);
    } else if(gauss_scheme==3){
        return gauss_knuth_yao.sample_vec(n);
    } else {
        return gauss_ziggurat.sample_vec(n);
    }
}

ZZX class_gauss_selected::sample_poly(int n) {
    if(gauss_scheme==1){
        return gauss_rejection.sample_poly(n);
    } else if(gauss_scheme==2){
        return gauss_inverse_cdf.sample_poly(n);
    } else if(gauss_scheme==3){
        return gauss_knuth_yao.sample_poly(n);
    } else {
        return gauss_ziggurat.sample_poly(n);
    }
}

mat_ZZ class_gauss_selected::sample_mat(int n, int m) {
    if(gauss_scheme==1){
        return gauss_rejection.sample_mat(n,m);
    } else if(gauss_scheme==2){
        return gauss_inverse_cdf.sample_mat(n,m);
    } else if(gauss_scheme==3){
        return gauss_knuth_yao.sample_mat(n,m);
    } else {
        return gauss_ziggurat.sample_mat(n,m);
    }
}