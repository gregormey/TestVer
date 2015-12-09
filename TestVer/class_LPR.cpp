/*
 * class_LPR.cpp
 *
 *  Created on: 3 dec. 2015
 *      Author: Susanne Riess, following Vanessa Erbenich
 *
 *  Implementation of the PKE scheme described in LPR12 and LPR13
 *
 *  The scheme is split into:
 *  1.) the preparation (computation of f and a),
 *  2.) key generation,
 *  3.) encryption and
 *  4.) decryption
 *
 *  Before every function you find a short description and a statement wether the time of this function
 *  will be measured or not.
 *
 *  The time is measured in main.cpp
 *
 *  Which Gaussian sampling method is used, is stated in class_gauss_selected.h and saved in the parameter *gauss
 */

#include "class_LPR.h"
#include <iostream>
#include <random>						// for uniformly random input
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ.h>
#include "class_gauss_selected.h"		// decides the gaussian sampling method
#include <algorithm>

using namespace std;
using namespace NTL;

class_LPR::class_LPR(long n, uint64_t p, double sigma, long factor, long precision) {
    this->n=n, this->p=p;
    ZZ p_new = ZZ(p);
    ZZ_p::init(p_new);

    this->gauss = new class_gauss_selected(sigma,factor,precision);
}


//---------------------- PREPARATION---- GENERATE POLYNOMIAL f AND a---------------------------------------------

// Generates the polynomial f = x^n+1
// which defines R_p=Z_p[X]/f
// Uses n
// No time will be measured

void class_LPR::generating_f() {

    this->f = ZZ_pX();
    this->f.SetLength(n);
    SetCoeff(this->f,0,1);
    SetCoeff(this->f,n,1);
}


// Generates uniformly random polynomial a over Z_p of degree < n
// The polynomial a is a public parameter of the scheme
// Uses ZZ_pX.random
// No time will be measured

void class_LPR::generating_uniform_polynomial_a() {
      this->a = ZZ_pX();
      random(a, n);
}

//-------------------------------------- KeyGeneration --------------------------------------------------------//


// Generates polynomials e_A and s_A over Z_p mod f with Gauss
// Time will be measured indivually

void class_LPR::generating_Gaussian_error_e_A() {
	this->e_A_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
	this->e_A = e_A_old % this->f;
	//ZZ x=conv<ZZ>(e_A.maxCoeff());
	//ZZ x=conv<ZZ>(max_element(this->e_A,this->e_A+n));

	//while(x>16){
    	//this->e_A_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    	//this->e_A = e_A_old % this->f;
    	//max_element(e_A.begin(),e_A.end());
    //}
}

void class_LPR::generating_Gaussian_secret_key_s_A() {
    this->s_A_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->s_A = s_A_old % this->f;}

// Computes the multiplicative part of the public key p_A mod f (named p_A_mult)
// Uses s_A, a and f
// Time will be measured

void class_LPR::computing_mult_part_of_public_key_p_A() {

    this->p_A_mult = (this->a * this->s_A) % this->f;
}

// Computes the additive part of the public key p_A mod f (named p_ALICE)
// Uses e_A, p_A_mult and f
// Time will be measured

void class_LPR::computing_add_part_of_public_key_p_A() {

    this->p_A = (this->p_A_mult + this->e_A) % this->f;
}


//-------------------------------------- Encryption -----------------------------------------------------------//

// Generates uniformly random key (lenght n)
// Uses uniform_int_distribution
// Time will be measured

void class_LPR::generating_uniform_key() {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 1);

    this->key = vec_ZZ_p();
    this->key.SetLength(n);

    for (long i=0; i<n; i++)
        this->key[i] = ZZ_p(dis(gen));
}

// Encodes key (convert to polynomial and multiply with floor p/2 (lenght n)
// Uses key and p
// Time will be measured

void class_LPR::encoding_key() {

    this->key_preenc = to_ZZ_pX(this->key);	// convertion to polynomial

    this->key_enc = (key_preenc * round(p/2)) % this->f;
}

// Generates polynomials e_B, e_B_1, s_B over Z_p mof f with Gauss
// Uses Gauss
// Time will be measured individually

void class_LPR::generating_Gaussian_error_s_B() {
    this->s_B = to_ZZ_pX(gauss->sample_poly(to_int(n)));
}

void class_LPR::generating_Gaussian_error_e_B() {
    this->e_B = to_ZZ_pX(gauss->sample_poly(to_int(n)));
}

void class_LPR::generating_Gaussian_error_e_B_1() {
    this->e_B_1 = to_ZZ_pX(gauss->sample_poly(to_int(n)));
}

// Computes multiplicative part mod f of ciphertext C=(p_B,c)
// Uses s_B, a, p_A and f
// Time will be measured individually

void class_LPR::computing_mult_part_of_ciphertext_p_B() {

    this->p_B_mult = ZZ_pX();
    this->p_B_mult = (this->a * this->s_B) % this->f;
}

void class_LPR::computing_mult_part_of_ciphertext_c() {

    this->c_mult = ZZ_pX();
    this->c_mult = (this->p_A * this->s_B) % this->f;
}

// Computes additive part mod f of ciphertext C=(p_B,c)
// Uses p_B_mult, c_mult, e_B, e_B_1, key_enc and f
// Time will be measured

void class_LPR::computing_add_part_of_ciphertext_p_B() {

    this->p_B = ZZ_pX();
    this->p_B = (this->p_B_mult + this->e_B) % this->f;
}

void class_LPR::computing_add_part_of_ciphertext_c() {

    this->c = ZZ_pX();
    this->c = (this->c_mult + this->e_B_1 + this->key_enc) % this->f;
}


//-------------------------------------- Decryption -----------------------------------------------------------//

// Computes multiplicative part of decrypt_mult = p_B*s_A mod f
// Uses p_B, s_A and f
// Time will be measured

void class_LPR::computing_mult_part_of_decrypt() {

    this->decrypt_mult = ZZ_pX();
    this->decrypt_mult = (this->p_B * this->s_A) % this->f;
}

// Computes additive part of decrypt = c-decrypt_mult mod f
// Uses decrypt_mult, c and f
// Time will be measured

void class_LPR::computing_add_part_of_decrypt() {

    this->decrypt = ZZ_pX();
    this->decrypt = (this->c - this->decrypt_mult ) % this->f;
}

// Convert decrypt, save in decrypt_int, multiply by p/2 and round, then take mod 2
// Time will be measured

void class_LPR::converting_decrypt() {

    this->decrypt_int = ZZX();
    this->decrypt_int = to_ZZX(this->decrypt);
}


// Determines c_dec by decoding c_new
// Uses polynomial c_new and decoding function
// time

void class_LPR::decoding_decrypt() {

	this->key_dec = ZZX();
    this->key_dec.SetLength(n);

    this->decrypt_int=2*this->decrypt_int;


    for (int i=0; i<n; i++){
    	this->key_dec[i] =  (this-> decrypt_int[i]/ p) %2;

    }
}



