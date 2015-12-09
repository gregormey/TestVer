/*
 * class_JD_ring.cpp
 *
 *  Created on: 1 dec. 2015
 *      Author: Susanne Riess, following Vanessa Erbenich
 *
 *  Implementation of the non-ring-based variant of the key exchange protocol described in JD12
 *
 *  The scheme is split into:
 *  1.) the preparation (computation of f and a),
 *  2.) key generation of Alice,
 *  3.) key generation of Bob and
 *  4.) protocol computations of Bob
 *  In the protocol computations of Bob some constants (p/2, p/4 etc) are computed.
 *
 *  We leave the protocol computation of Alice, as it will just be the same as for Bob without calculating Sig
 *
 *  Before every function you find a short description and a statement wether the time of this function
 *  will be measured or not.
 *
 *  The time is measured in main.cpp
 *
 *  Which Gaussian sampling method is used, is stated in class_gauss_selected.h and saved in the parameter *gauss
 */

#include "class_JD.h"
#include <iostream>
#include "class_gauss_selected.h"	// decides the gaussian sampling method
#include <random>					// for uniformly random input
#include "math.h"					// for functions floor and ceil
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_p.h>

using namespace std;
using namespace NTL;

class_JD::class_JD(long n, uint64_t p, double sigma, long factor, long precision) {
    this->n=n, this->p=p;
    ZZ p_new = ZZ(p);
    ZZ_p::init(p_new);

    this->gauss = new class_gauss_selected(sigma, factor, precision);
}


//---------------------- PREPARATION---- GENERATE MATRIX A-------------------------------------------------------

// Generates uniformly random matrix A (n x n) over Z_p
// Matrix A is a public parameter of the scheme
// Uses uniform_int_distribution
// No time will be measured

void class_JD::generating_uniform_A() {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, p-1);

    this->A = mat_ZZ_p();
    this->A.SetDims(n, n);

    for (long i=1; i<=n; i++)
        for (long j=1; j<=n; j++)
            this->A(i, j) = ZZ_p(dis(gen));
}


//------------------------ KEY GENERATION for ALICE--------------------------------------------------------------

// Generates vectors e_ALICE and s_ALICE over Z_p with Gauss
// Time will be measured individually

void class_JD::generating_Gaussian_error_e_ALICE() {
    this->e_ALICE = to_vec_ZZ_p(gauss->sample_vec(to_int(n)));
}

void class_JD::generating_Gaussian_secret_key_s_ALICE() {
    this->s_ALICE = to_vec_ZZ_p(gauss->sample_vec(to_int(n)));
}

// Computes two times e_ALICE, saves it in two_e_ALICE
// Time will be measured

void class_JD::computing_two_times_e_ALICE() {
	this->two_e_ALICE=vec_ZZ_p();
    this->two_e_ALICE.SetLength(n);
    this->two_e_ALICE = 2 * this->e_ALICE;
}

// Computes the multiplicative part of the public key p_ALICE of Alice (named p_ALICE_mult)
// Uses s_ALICE and A
// Time will be measured

void class_JD::computing_mult_part_of_public_key_p_ALICE() {

    this->p_ALICE_mult = vec_ZZ_p();
    this->p_ALICE_mult.SetLength(n);

    this->p_ALICE_mult = this->A * this->s_ALICE;
}

// Computes the additive part of the public key p_ALICE of Alice (named p_ALICE)
// Uses two_e_ALICE and p_ALICE_mult
// Time will be measured

void class_JD::computing_add_part_of_public_key_p_ALICE() {

    this->p_ALICE = vec_ZZ_p();
    this->p_ALICE.SetLength(n);

    this->p_ALICE = this->p_ALICE_mult + this->two_e_ALICE;
}


//------------------------ KEY GENERATION for BOB--------------------------------------------------------------

// Generates vectors e_BOB and s_BOB over Z_p with Gauss
// Time will be measured individually

void class_JD::generating_Gaussian_error_e_BOB() {
    this->e_BOB = to_vec_ZZ_p(gauss->sample_vec(to_int(n)));
}

void class_JD::generating_Gaussian_secret_key_s_BOB() {
    this->s_BOB = to_vec_ZZ_p(gauss->sample_vec(to_int(n)));
}

// Computes two times e_BOB, saves it in two_e_BOB
// Time will be measured

void class_JD::computing_two_times_e_BOB() {
	this->two_e_BOB=vec_ZZ_p();
    this->two_e_BOB.SetLength(n);
    this->two_e_BOB = 2 * this->e_BOB;
}

// Computes the multiplicative part of the public key p_BOB of BOB (named p_BOB_mult)
// Uses s_BOB and A
// Time will be measured

void class_JD::computing_mult_part_of_public_key_p_BOB() {

    this->p_BOB_mult = vec_ZZ_p();
    this->p_BOB_mult.SetLength(n);

    this->p_BOB_mult = this->A * this->s_BOB;
}

// Computes the additive part of the public key p_BOB of BOB (named p_BOB)
// Uses two_e_BOB and p_BOB_mult
// Time will be measured

void class_JD::computing_add_part_of_public_key_p_BOB() {

    this->p_BOB = vec_ZZ_p();
    this->p_BOB.SetLength(n);

    this->p_BOB = this->p_BOB_mult + this->two_e_BOB;
}


//--------------------------- THE PROTOCOL for BOB---------------------------------------------------------------

// Generates the number e_BOB_2 over Z_p with Gauss
// Time will be measured

void class_JD::generating_Gaussian_error_e_BOB_2() {
	this->e_BOB_2 = to_vec_ZZ_p(gauss->sample_vec(to_int(1)));
}

// Computes two times e_BOB_2, saves it in two_e_BOB_2
// Time will be measured

void class_JD::computing_two_times_e_BOB_2() {
	this->two_e_BOB_2=ZZ_p();
	this->two_e_BOB_2 = 2 * this->e_BOB_2(1);
}

// Computes the multiplicative part of the number k_BOB (named k_BOB_mult)
// Uses s_BOB and p_ALICE
// Time will be measured

void class_JD::computing_mult_part_of_vector_k_BOB() {

    this->k_BOB_mult = ZZ_p();
    InnerProduct(this->k_BOB_mult, this->p_ALICE, this->s_BOB);
}

// Computes the additive part of the number k_BOB (named k_BOB)
// Uses two_e_BOB_2 and k_BOB_mult
// Time will be measured

void class_JD::computing_add_part_of_vector_k_BOB() {

    this->k_BOB = ZZ_p();
    this->k_BOB = this->k_BOB_mult + this->two_e_BOB_2;
}

// Generates uniformly random bit b
// Uses uniform_int_distribution
// Time will be measured

void class_JD::generating_uniform_bit_b() {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 1);

    this->b = ZZ_p();
    this->b = ZZ_p(dis(gen));
}

// Computes floor(p/4)
// Time will be measured individually

void class_JD::computing_floor_p_quarter() {
	this->floor_p_quarter=floor_p_quarter;
	this-> floor_p_quarter = floor(this->p/4);
}

// Computes the value of function S for input k_BOB and for random bit b
// Return value will be saved in Sig
// Time will be measured

void class_JD::computing_S_k_BOB() {

	this->Sig = ZZ();
	this->k_BOB_mod = ZZ();
	this->k_BOB_mod=conv<ZZ>(this->k_BOB); 	// for compatibility

	if (this->b==0) {
		if((this->k_BOB_mod>this->floor_p_quarter) && (this->k_BOB_mod<(this->p-this->floor_p_quarter))){
			this->Sig=1;
		}
		else{
			this->Sig=0;
		}
	}
	else{
		if((this->k_BOB_mod>(this->floor_p_quarter+1)) && (this->k_BOB_mod<(p+this->floor_p_quarter))){
			this->Sig=1;
		}
		else{
			this->Sig=0;
		}
	}

}

// Computes the value of the function E for input k_BOB and Sig
// Return value will be saved in sk_BOB
// Time will be measured

void class_JD::computing_sk_BOB_with_k_BOB_and_Sig(){

	this->sk_BOB = ZZ();

    this->sk_BOB=this->k_BOB_mod+(this->Sig*(int)(p-1)/2);
    this->sk_BOB=this->sk_BOB % 2;
}
