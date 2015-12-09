//
//  class_gauss_ziggurat.cpp
//  Encryption schemes
//
//  Created by Vanessa Erbenich on 11.09.15.
//  Copyright (c) 2015 Vanessa Erbenich. All rights reserved.
//

#include "class_gauss_ziggurat.h"
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>
#include "ziggurat-new.h"

#include <iostream>

using namespace NTL;
using namespace std;


RR class_gauss_ziggurat::rho_new(RR sigma, RR x) {
    return exp(-power(to_RR(x)/to_RR(sigma), 2)/to_RR(2));
}


RR class_gauss_ziggurat::compute_recurrence(RR* xis, RR* rhos, RR c, RR m, RR sigma) {
    int em = to_int(m)+1;
    RR om = to_RR(1)/to_RR(m), m2 = to_RR(-2), o2 = to_RR(1)/to_RR(2);
    RR area = sigma*om*sqrt(ComputePi_RR()/to_RR(2))*to_RR(c);
    xis[em-1] = to_RR(13*sigma);
    rhos[em-1] = rho_new(sigma, to_RR(TruncToZZ(xis[em-1])+to_ZZ(1)));
    RR sqrtv = m2*log(area/to_RR(TruncToZZ(xis[em-1])+to_ZZ(1)));
    if (sqrtv < to_RR(0)) return to_RR(-1);
    xis[em-2] = to_RR(to_RR(sigma)*sqrt(sqrtv));
    rhos[em-2] = rho_new(sigma, xis[em-2]);
    
    for (int i = em-3; i > 0; i--)
    {
        sqrtv = m2*log(area/to_RR(TruncToZZ(xis[i+1])+to_ZZ(1))+rho_new(sigma, xis[i+1]));
        if (sqrtv < to_RR(0)) return to_RR(-1);
        xis[i] = to_RR(to_RR(sigma)*sqrt(sqrtv));
        rhos[i] = exp(-o2*power(to_RR(xis[i]/sigma), 2));
    }
    rhos[0] = area/to_RR(TruncToZZ(xis[1])+to_ZZ(1)) + rho_new(sigma, xis[1]);
    return rhos[0];
}


bool class_gauss_ziggurat::create_partitions() {
    
    //< set precisions
    RR prec = power2_RR(-precision); //< max. approximation error
    
    //< initialize different variables
    RR* xis = new RR[m+1];
    RR* bestxis = new RR[m+1];
    for (int i = 0; i <= m; i++)
        bestxis[i] = -1;
    RR* rhos = new RR[m+1];
    RR y0;
    RR m_RR = ConvPrec(m, precision);
    RR c = 1 + 1/m_RR;
    RR tailcut_RR = ConvPrec(tailcut, precision);
    
    if (m_RR == 1) {
        bestxis[0] = 0;
        bestxis[1] = tailcut_RR;
        y0 = to_RR(1); // "to_RR(1) + bestdiff" mit "bestdiff=0"
    }
    
    RR bestdiff = to_RR(3);
    
    while (tailcut_RR < to_RR(14*sigma)) {
        xis[m] = tailcut_RR;
        RR cu = to_RR(0), cl = to_RR(1), cc;
        RR difference = to_RR(-1), lastdiff = to_RR(-2);
        
        while (difference < 0 || (abs(difference) > prec && abs(difference-lastdiff) > prec)) {
            cc = c;
            lastdiff = difference;
            difference = compute_recurrence(xis, rhos, c, m_RR, sigma) - to_RR(1);
            if (difference == -2)
                break;
            if (difference >= 0) {
                for (int i = 0; i <= m; i++)
                    bestxis[i] = xis[i];
                cc = c;
                cu = c;
                bestdiff = difference;
            }
            else
                cl = c;
            if (cu < cl)
                c += to_RR(1)/m_RR;
            else
                c = (cu+cl)/to_RR(2);
            if (c >= 11)
                break;
        }
        
        if (difference < 0 || (abs(difference) > prec && abs(difference-lastdiff) > prec))
            tailcut_RR++;
        else
            break;
    }
    
    if (bestxis[m] != -1) {
        // Partition gefunden
        //cout << "Partition gefunden" << endl;
        y0 = to_RR(1) + bestdiff;
    } else {
        // Keine Partition gefunden!
        cout << "Keine Partition gefunden!" << endl;
        return false;
    }
    
    // -------------
    // set partition
    // -------------
    
    rectys[0] = y0;
    for (int i = 0; i <= m; i++) {
        rectxs[i] = TruncToZZ(bestxis[i]);
        if (i != 0)
            rectys[i] = rho_new(bestxis[i], sigma);
    }
    
    return true;
}


class_gauss_ziggurat::class_gauss_ziggurat(){}

class_gauss_ziggurat::class_gauss_ziggurat(int m, int precision, int sigma) {
    
    RR::SetPrecision(precision); //< max. approximation error
    RR::SetOutputPrecision(precision); //< output-precision of numbers
    
    this->m = m;
    this->precision = precision;
    this->sigma = ConvPrec(sigma, precision);
    omega = precision; // Wurde von Patrick so gewählt
    this->tailcut = sigma*13;
    
    rectxs = new ZZ[m+1];
    rectys = new RR[m+1];
    rhovals = new ZZ[tailcut+1];
    
    set_rhovals(rhovals, &omega, &this->tailcut, &this->sigma);
    
    create_partitions();
    
}

vec_ZZ class_gauss_ziggurat::sample_vec(int n) {

    vec_ZZ res;
    res.SetLength(n);
    for (int i=1; i<=n; i++)
        ziggurat_new(&res(i), &m, rectxs, rectys, &sigma, &omega);
    
    return res;
}

ZZX class_gauss_ziggurat::sample_poly(int n) {
    
    ZZX res;
    res.SetLength(n+1);
    for (int i=0; i<=n; i++)
        ziggurat_new(&res[i], &m, rectxs, rectys, &sigma, &omega);
    
    return res;
}

mat_ZZ class_gauss_ziggurat::sample_mat(int n, int m) {
    
    mat_ZZ res;
    res.SetDims(n,m);
    for (int i=1; i<=n; i++)
        for (int j=1; j<=m; j++)
            ziggurat_new(&res(i,j), &m, rectxs, rectys, &sigma, &omega);
    
    return res;
}


