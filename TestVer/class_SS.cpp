/*
 * class_SS.cpp
 *
 *  Created on: 7 dec. 2015
 *  Author:Susanne Riess, following Vanessa Erbenich
 *
 *  Implementation of the PKE scheme by Stehlé and Steinfeld
 *
 *  The scheme is split into:
 *  1.) key generation,
 *  2.) encryption and
 *  3.) decryption
 *
 *  Before every function you find a short description and a statement wether the time of this function
 *  will be measured or not.
 *
 *  The time is measured in main.cpp
 *
 *  Which Gaussian sampling method is used, is stated in class_gauss_selected.h and saved in the parameter *gauss
 */

#include "class_SS.h"
#include <iostream>
#include <random>						// for uniformly random input
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ.h>
#include "class_gauss_selected.h"		// decides the gaussian sampling method

using namespace std;
using namespace NTL;

class_SS::class_SS(long n, uint64_t p, double alpha, double sigma, long factor, long precision) {
	this->n=n, this->p=p;
	ZZ p_new = ZZ(p);
	ZZ_p::init(p_new);

	this->gauss = new class_gauss_selected(alpha,factor,precision);
	this->gauss_sigma = new class_gauss_selected(sigma,factor,precision);

}

//---------------------- PREPARATION---- GENERATE POLYNOMIAL f AND a---------------------------------------------

// Generates the polynomial f = x^n+1
// which defines R_p=Z_p[X]/f
// Uses n
// No time will be measured

void class_SS::generating_f() {

    this->f = ZZ_pX();
    this->f.SetLength(n);
    SetCoeff(this->f,0,1);
    SetCoeff(this->f,n,1);
}


//-------------------------------------- KeyGeneration --------------------------------------------------------//

void class_SS::generating_part_of_Gaussian_secret_key_f() {
    this->f_key = to_ZZ_pX(gauss_sigma->sample_poly(to_int(n)))% this->f;
}

void class_SS::generating_mult_part_of_Gaussian_secret_key_f() {
	this->f_key= ZZ_pX();
	this->f_key.SetLength(n);
	this->f_key = (2 * this->f_key)%this->f;
}

void class_SS::generating_add_part_of_Gaussian_secret_key_f() {
	this->f_key = (1 + this->f_key)%this->f;
}

void class_SS::generating_part_of_Gaussian_public_key_g() {
    this->g = to_ZZ_pX(gauss_sigma->sample_poly(to_int(n)))%this->f;
}

void class_SS::generating_mult_part_of_Gaussian_public_key_h() {
	this->h = ZZ_pX();
	this->h.SetLength(n);
	this->h = (2 * this->g)%this->f;
}

void class_SS::generating_public_key_h() {
	this->h = (this-> h / this->f_key)%this->f;
}

//-------------------------------------- Encryption -----------------------------------------------------------//

void class_SS::generating_Gaussian_s() {
    this->s = to_ZZ_pX(gauss->sample_poly(to_int(n)))% this->f;
}

void class_SS::generating_Gaussian_e() {
    this->e = to_ZZ_pX(gauss->sample_poly(to_int(n)))% this->f;
}

void class_SS::generating_first_mult_part_of_ciphertext_c() {
	this->c_1 = ZZ_pX();
	this->c_1.SetLength(n);
	this->c_1 = (this->h * this->s)%this->f;
}

void class_SS::generating_second_mult_part_of_ciphertext_c() {
	this->c = ZZ_pX();
	this->c.SetLength(n);
	this->c = (2 * this->e)%this->f;
}

void class_SS::generating_add_part_of_ciphertext_c() {
	this->c = (this-> c_1 + this-> c)%this->f;
}

// Generates uniformly random key (lenght n)
// Uses uniform_int_distribution
// Time will be measured

void class_SS::generating_uniform_key() {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 1);

    this->key = vec_ZZ_p();
    this->key.SetLength(n);

    for (long i=1; i<=n; i++)
        this->key(i) = ZZ_p(dis(gen));
}

void class_SS::generating_ciphertext_c() {
	this-> key_mod=conv<ZZ_pX>(this-> key);
	this->c = (this-> c + this-> key_mod)%this->f;
}

//-------------------------------------- Decryption -----------------------------------------------------------//

void class_SS::generating_decrypted_c() {
	this->d = ZZ_pX();
	this->d.SetLength(n);
	this->d = (this->c * this->f_key)%this->f;
}

void class_SS::generating_decrypted_key() {
	this-> d_mod=conv<ZZX>(this-> d);

	for(int i=1; i<=n; i++){
		this->d_mod[i] = (this->d_mod[i])%2;
	}
}
