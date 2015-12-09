/*
 * class_BCNS.cpp
 *
 *  Created on: 30 nov. 2015
 *      Author: Susanne Riess, following Vanessa Erbenich
 *
 *  Implementation of the key exchange protocol described in BCNS14, originally designed by Peikert as a KEM
 *
 *  The scheme is split into:
 *  1.) the preparation (computation of f and a),
 *  2.) key generation of Alice,
 *  3.) key generation of Bob,
 *  4.) protocol computations of Bob and
 *  5.) protocol computations of Alice
 *  In the protocol computations of Bob some constants (p/2, p/4 etc) are computed.
 *
 *  Before every function you find a short description and a statedment wether the time of this function
 *  will be measured or not.
 *
 *  The time is measured in main.cpp
 *
 *  Which Gaussian sampling method is used, is stated in class_gauss_selected.h and saved in the parameter *gauss
 */

#include "class_BCNS.h"
#include <iostream>
#include "class_gauss_selected.h"	// decides the gaussian sampling method
#include <random>					// for uniformly random input
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include "math.h"					// for functions floor and ceil

using namespace std;
using namespace NTL;

class_BCNS::class_BCNS(long n, uint64_t p, double sigma, long factor, long precision) {
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

void class_BCNS::generating_f() {

    this->f = ZZ_pX();
    this->f.SetLength(n);
    SetCoeff(this->f,0,1);
    SetCoeff(this->f,n,1);
}


// Generates uniformly random polynomial a over Z_p of degree < n
// The polynomial a is a public parameter of the scheme
// Uses ZZ_pX.random
// No time will be measured

void class_BCNS::generating_uniform_polynomial_a() {
	this->a = ZZ_pX();
	random(a, n);
	cout<<"a"<<endl;
	cout<< a<< endl;
}



//------------------------ KEY GENERATION for ALICE--------------------------------------------------------------


// Generates polynomials e_ALICE and s_ALICE over Z_p mod f with Gauss
// Time will be measured indivually

void class_BCNS::generating_Gaussian_error_e_ALICE() {
    this->e_ALICE_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->e_ALICE = e_ALICE_old % this->f;
}

void class_BCNS::generating_Gaussian_secret_key_s_ALICE() {
    this->s_ALICE_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->s_ALICE = s_ALICE_old % this->f;
}


// Computes the multiplicative part of the public key p_ALICE of Alice mod f (named p_ALICE_mult)
// Uses s_ALICE, a and f
// Time will be measured

void class_BCNS::computing_mult_part_of_public_key_p_ALICE() {

    this->p_ALICE_mult = ZZ_pX();
    this->p_ALICE_mult.SetLength(n);
    this->p_ALICE_mult = (this->a * this->s_ALICE) % this->f;
}

// Computes the additive part of the public key p_ALICE of Alice mod f (named p_ALICE)
// Uses e_ALICE, p_ALICE_mult and f
// Time will be measured

void class_BCNS::computing_add_part_of_public_key_p_ALICE() {

    this->p_ALICE = ZZ_pX();
    this->p_ALICE.SetLength(n);
    this->p_ALICE = (this->p_ALICE_mult + this->e_ALICE) % this->f;
}



//------------------------ KEY GENERATION for BOB----------------------------------------------------------------


// Generates polynomials e_BOB and s_BOB over Z_p mod f with Gauss
// Time will be measured indivually

void class_BCNS::generating_Gaussian_error_e_BOB() {
    this->e_BOB_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->e_BOB = e_BOB_old % this->f;
}

void class_BCNS::generating_Gaussian_secret_key_s_BOB() {
    this->s_BOB_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->s_BOB = s_BOB_old % this->f;
}

// Computes the multiplicative part of the public key p_BOB of BOB mod f (named p_BOB_mult)
// Uses s_BOB, a and f
// Time will be measured

void class_BCNS::computing_mult_part_of_public_key_p_BOB() {

    this->p_BOB_mult = ZZ_pX();
    this->p_BOB_mult.SetLength(n);
    this->p_BOB_mult = (this->a * this->s_BOB) % this->f;
}

// Computes the additive part of the public key p_BOB of BOB mod f (named p_BOB)
// Uses e_BOB, p_BOB_mult and f
// Time will be measured

void class_BCNS::computing_add_part_of_public_key_p_BOB() {

    this->p_BOB = ZZ_pX();
    this->p_BOB.SetLength(n);
    this->p_BOB = (this->p_BOB_mult + this->e_BOB) % this->f;
}



//--------------------------- THE PROTOCOL for BOB---------------------------------------------------------------


// Generates the polynomial e_BOB_2 over Z_p mod f with Gauss
// Time will be measured

void class_BCNS::generating_Gaussian_error_e_BOB_2() {
    this->e_BOB_2_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->e_BOB_2 = e_BOB_2_old % this->f;
}

// Computes the multiplicative part of the polynomial v mod f (named v_mult)
// Uses s_BOB, p_ALICE and f
// Time will be measured

void class_BCNS::computing_mult_part_of_polynomial_v() {

    this->v_mult = ZZ_pX();
    this->v_mult.SetLength(n);
    this->v_mult = (this->p_ALICE * this->s_BOB) % this->f;
}

// Computes the additive part of the polynomial v mod f (named v)
// Uses e_BOB_2, v_mult and f
// Time will be measured

void class_BCNS::computing_add_part_of_polynomial_v() {

    this->v = ZZ_pX();
    this->v.SetLength(n);
    this->v = (this->v_mult + this->e_BOB_2) % this->f;
}

//Creates a uniform random polynomial ran with coefficients in [-1,0,0,1]
// Time will be measured

void class_BCNS::generating_random_ternary_polynomial_ran() {

	int numbers[] = { -1, 0, 0, 1 };
	int length = sizeof(numbers) / sizeof(int);

    this->ran = ZZX();
    this->ran.SetLength(n);

    for (long i=0; i<n; i++){
        this->ran[i] = ZZ(numbers[rand() % length]);
    }
}

// Computes multiplicative part of dbl(v), named dbl_mult
// which is 2*v
// Time will be measured

void class_BCNS::computing_mult_part_of_double_v() {

	this->dbl_mult = ZZX();
	this->dbl_mult.SetLength(n);
	this->dbl_mult=conv<ZZX>(2*this->v);		// for compatibility
	this->f_mod= conv<ZZX>(this->f);			// for compatibility
	this->dbl_mult = this->dbl_mult % this->f_mod;
}

// Computes additive part of dbl(v), named dbl
// which is 2*v+ran
// Time will be measured

void class_BCNS::computing_add_part_of_double_v() {

	this->dbl = ZZX();
	this->dbl.SetLength(n);
	this->dbl = (this->dbl_mult + this->ran) % this->f_mod;
}

// Computes p/2, p/4 and 7p/4
// Time will be measured individually

void class_BCNS::computing_p_half() {
	this->p_half=p_half;
	this-> p_half = this->p/2;
}

void class_BCNS::computing_p_quarter() {
	this->p_quarter=p_quarter;
	this-> p_quarter = this->p/4;
}

void class_BCNS::computing_7p_quarter() {
	this->seven_p_quarter=seven_p_quarter;
	this-> seven_p_quarter = 7*this->p_quarter;
}

// Computes rounded p/2 and rounded 3p/2
// Time will be measured individually

void class_BCNS::computing_rounded_p_half() {
	this->rounded_p_half=rounded_p_half;
	this-> rounded_p_half = (int)round(this->p_half);
}

void class_BCNS::computing_rounded_3p_half() {
	this->rounded_3p_half=rounded_3p_half;
	this-> rounded_3p_half = (int)round(3*this->p_half);
}

// Computes value of the cross rounding function for input dbl(v)
// Time will be measured

void class_BCNS::crossrounding_dbl_v() {

	this->c = ZZ_pX();
	this->c.SetLength(n);

	for (long i=0; i<n; i++) {
			if ((dbl[i]>= this-> rounded_p_half && dbl[i]<=(p-1)) || (dbl[i]>= this-> rounded_3p_half && dbl[i]<2*p)){
				this->c[i]=1;
			}
			else {
				this->c[i]=0;
			}
		}
}

// Computes value of the rounding function for input dbl(v)
// Time will be measured

void class_BCNS::rounding_dbl_v() {

	this->k = ZZ_pX();
	this->k.SetLength(n);

	for (long i=0; i<n; i++) {
			if (dbl[i] >=  this->rounded_p_half && dbl[i]<=(this->rounded_3p_half-1) ) {
				this->k[i]=1;
			}
			else {
				this->k[i]=0;
			}
		}
}



//----------------------------THE PROTOCOL for ALICE-------------------------------------------------------------


// Computes 2*p_BOB
// Time will be measured

void class_BCNS::computing_two_p_BOB() {

	this->two_p_BOB = ZZ_pX();
	this->two_p_BOB.SetLength(n);
	this->two_p_BOB = 2*this->p_BOB % this->f;
}

// Computes 2*p_BOB*s_ALICE, named rec
// Time will be measured

void class_BCNS::computing_two_p_BOB_s_ALICE() {

	this->rec = ZZ_pX();
	this->rec.SetLength(n);
	this->rec = this-> two_p_BOB*this->s_ALICE % this->f;
}

// Computes the value of the reconcilation function (rec) for the input rec
// Time will be measured

void class_BCNS::computing_key_with_rec() {

	this->key = ZZX();
	this->key.SetLength(n);
	this->rec_mod=conv<ZZX>(this->rec); 	// for compatibility

	for (long i=0; i<n; i++) {
		if (this->c[i]==0) {
			if((this->rec_mod[i]>=(this->rounded_p_half+(int)ceil(this->p_quarter)-1)) && (this->rec_mod[i]<(int)floor(this->seven_p_quarter))){
				this->key[i]=1;
			}
			else{
				this->key[i]=0;
			}
		}
		else{
			if((this->rec_mod[i]>=(int)(ceil(this->p_quarter-1))) && (this->rec_mod[i]<=(int)(floor(this->seven_p_quarter)-floor(this->p_half)))){
				this->key[i]=1;
			}
			else{
				this->key[i]=0;
			}
		}
	}
}

