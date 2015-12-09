/*
 * ZZD_one_pass.cpp
 *
 *  Created on: 3 dec. 2015
 *      Author: Susanne Riess, following Vanessa Erbenich
 *
 *  Implementation of the one-pass variant of the key exchange protocol described in ZZD+14
 *
 *  The scheme is split into:
 *  1.) the preparation (computation of f and a),
 *  2.) key generation of Alice,
 *  3.) key generation of Bob
 *  4.) protocol computations of Alice before sending the only message
 *  5.) protocol computations of Alice after sending the only message and
 *  6.) protocol computations of Bob after receiving the only message
 *  In the protocol computations of Bob some constants (p/2, p/4 etc) are computed.
 *
 *  Before every function you find a short description and a statement wether the time of this function
 *  will be measured or not.
 *
 *  The time is measured in main.cpp
 *
 *  Which Gaussian sampling method is used, is stated in class_gauss_selected.h and saved in the parameter *gauss
 */

#include "class_ZZD_one_pass.h"
#include <iostream>
#include "class_gauss_selected.h"	// decides the gaussian sampling method
#include <random>					// for uniformly random input
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include "math.h"					// for functions floor and ceil

using namespace std;
using namespace NTL;

class_ZZD_one_pass::class_ZZD_one_pass(long n, uint64_t p, double alpha, double beta, long factor, long precision) {
    this->n=n, this->p=p;
    ZZ p_new = ZZ(p);
    ZZ_p::init(p_new);

    this->gauss = new class_gauss_selected(alpha, factor, precision);
    this->gauss_beta = new class_gauss_selected(beta, factor, precision);
}


//---------------------- PREPARATION---- GENERATE POLYNOMIAL f AND a---------------------------------------------

// Generates the polynomial f = x^n+1
// which defines R_p=Z_p[X]/f
// Uses n
// No time will be measured

void class_ZZD_one_pass::generating_f() {

    this->f = ZZ_pX();
    this->f.SetLength(n);
    SetCoeff(this->f,0,1);
    SetCoeff(this->f,n,1);
}

// Generates uniformly random polynomial a over Z_p of degree < n
// The polynomial a is a public parameter of the scheme
// Uses ZZ_pX.random
// No time will be measured

void class_ZZD_one_pass::generating_uniform_polynomial_a() {
      this->a = ZZ_pX();
      random(a, n);
      cout<<"a"<<endl;
      cout<<a<<endl;
}


//------------------------ KEY GENERATION for ALICE--------------------------------------------------------------

// Generates polynomials e_ALICE and s_ALICE over Z_p mod f with Gauss
// Time will be measured indivually

void class_ZZD_one_pass::generating_Gaussian_error_e_ALICE() {
    this->e_ALICE_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->e_ALICE = e_ALICE_old % this->f;
}

void class_ZZD_one_pass::generating_Gaussian_secret_key_s_ALICE() {
    this->s_ALICE_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->s_ALICE = s_ALICE_old % this->f;
}

// Computes two times e_ALICE, saves it in two_e_ALICE
// Time will be measured

void class_ZZD_one_pass::computing_two_times_e_ALICE() {
	this->two_e_ALICE=ZZ_pX();
    this->two_e_ALICE.SetLength(n);
    this->two_e_ALICE = (2 * this->e_ALICE) % this->f;
}

// Computes the multiplicative part of the public key p_ALICE of Alice mod f (named p_ALICE_mult)
// Uses s_ALICE, a and f
// Time will be measured

void class_ZZD_one_pass::computing_mult_part_of_public_key_p_ALICE() {

    this->p_ALICE_mult = ZZ_pX();
    this->p_ALICE_mult.SetLength(n);
    this->p_ALICE_mult = (this->a * this->s_ALICE) % this->f;
}

// Computes the additive part of the public key p_ALICE of Alice mod f (named p_ALICE)
// Uses two_e_ALICE, p_ALICE_mult and f
// Time will be measured

void class_ZZD_one_pass::computing_add_part_of_public_key_p_ALICE() {

    this->p_ALICE = ZZ_pX();
    this->p_ALICE.SetLength(n);
    this->p_ALICE = (this->p_ALICE_mult + this->two_e_ALICE) % this->f;
}


//------------------------ KEY GENERATION for BOB----------------------------------------------------------------

// Generates polynomials e_BOB and s_BOB over Z_p mod f with Gauss
// Time will be measured indivually

void class_ZZD_one_pass::generating_Gaussian_error_e_BOB() {
    this->e_BOB_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->e_BOB = e_BOB_old % this->f;
}

void class_ZZD_one_pass::generating_Gaussian_secret_key_s_BOB() {
    this->s_BOB_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->s_BOB = s_BOB_old % this->f;
}

// Computes two times e_BOB, saves it in two_e_BOB
// Time will be measured

void class_ZZD_one_pass::computing_two_times_e_BOB() {
	this->two_e_BOB=ZZ_pX();
    this->two_e_BOB.SetLength(n);
    this->two_e_BOB = (2 * this->e_BOB) % this->f;
}

// Computes the multiplicative part of the public key p_BOB of BOB mod f (named p_BOB_mult)
// Uses s_BOB, a and f
// Time will be measured

void class_ZZD_one_pass::computing_mult_part_of_public_key_p_BOB() {

    this->p_BOB_mult = ZZ_pX();
    this->p_BOB_mult.SetLength(n);
    this->p_BOB_mult = (this->a * this->s_BOB) % this->f;
}

// Computes the additive part of the public key p_BOB of BOB mod f (named p_BOB)
// Uses two_e_BOB, p_BOB_mult and f
// Time will be measured

void class_ZZD_one_pass::computing_add_part_of_public_key_p_BOB() {

    this->p_BOB = ZZ_pX();
    this->p_BOB.SetLength(n);
    this->p_BOB = (this->p_BOB_mult + this->two_e_BOB) % this->f;
}


//-------------------------------------- Computation ALICE before the only mesaage is send PASS -----------------

// Generates polynomials r_ALICE, f_ALICE over Z_p mod f with Gauss (beta)
// Time will be measured individually

void class_ZZD_one_pass::generating_Gaussian_r_ALICE() {
    this->r_ALICE_old = to_ZZ_pX(gauss_beta->sample_poly(to_int(n)));
    this->r_ALICE = r_ALICE_old % this->f;
}

void class_ZZD_one_pass::generating_Gaussian_f_ALICE() {
    this->f_ALICE_old = to_ZZ_pX(gauss_beta->sample_poly(to_int(n)));
    this->f_ALICE = f_ALICE_old % this->f;
}

// Computes two times f_ALICE, saves it in two_f_ALICE
// Time will be measured

void class_ZZD_one_pass::computing_two_times_f_ALICE() {
	this->two_f_ALICE=ZZ_pX();
    this->two_f_ALICE.SetLength(n);
    this->two_f_ALICE = (2 * this->f_ALICE) % this->f;
}

// Computes the multiplicative part of x_ALICE mod f (named x_ALICE_mult)
// Uses r_ALICE, a and f
// Time will be measured

void class_ZZD_one_pass::computing_mult_part_of_x_ALICE() {

    this->x_ALICE_mult = ZZ_pX();
    this->x_ALICE_mult.SetLength(n);
    this->x_ALICE_mult = (this->a * this->r_ALICE) % this->f;
}

// Computes the additive part of the polynomial x_ALICE mod f (named x_ALICE)
// Uses two_f_ALICE, x_ALICE_mult and f
// Time will be measured

void class_ZZD_one_pass::computing_add_part_of_x_ALICE() {

    this->x_ALICE = ZZ_pX();
    this->x_ALICE.SetLength(n);
    this->x_ALICE = (this->two_f_ALICE + this->x_ALICE_mult) % this->f;
}

// Compute polynomial d_A
// Simulate entries as Gaussian distributed with factor alpha
// No time will be measured

void class_ZZD_one_pass::generating_Gaussian_d_A() {
    this->d_A_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->d_A = d_A_old % this->f;
}

// Computes the multiplicative part of hat_r_ALICE mod f (named hat_r_ALICE_mult)
// Uses d_A, s_ALICE and f
// Time will be measured

void class_ZZD_one_pass::computing_mult_part_of_hat_r_ALICE() {

    this->hat_r_ALICE_mult = ZZ_pX();
    this->hat_r_ALICE_mult.SetLength(n);
    this->hat_r_ALICE_mult = (this->s_ALICE * this->d_A) % this->f;
}

// Computes the multiplicative part of hat_f_ALICE mod f (named hat_f_ALICE_mult)
// Uses d_A, e_ALICE and f
// Time will be measured

void class_ZZD_one_pass::computing_mult_part_of_hat_f_ALICE() {

    this->hat_f_ALICE_mult = ZZ_pX();
    this->hat_f_ALICE_mult.SetLength(n);
    this->hat_f_ALICE_mult = (this->e_ALICE * this->d_A) % this->f;
}

// Computes the additive part of the polynomial hat_r_ALICE mod f (named hat_r_ALICE)
// Uses hat_r_ALICE_mult, r_ALICE_mult and f
// Time will be measured

void class_ZZD_one_pass::computing_add_part_of_hat_r_ALICE() {

    this->hat_r_ALICE = ZZ_pX();
    this->hat_r_ALICE.SetLength(n);
    this->hat_r_ALICE = (this->hat_r_ALICE_mult + this->r_ALICE) % this->f;
}

// Computes the additive part of the polynomial hat_f_ALICE mod f (named hat_f_ALICE)
// Uses hat_f_ALICE_mult, f_ALICE_mult and f
// Time will be measured

void class_ZZD_one_pass::computing_add_part_of_hat_f_ALICE() {

    this->hat_f_ALICE = ZZ_pX();
    this->hat_f_ALICE.SetLength(n);
    this->hat_f_ALICE = (this->hat_f_ALICE_mult + this->f_ALICE) % this->f;
}

// Generates the polynomial g_ALICE over Z_p mod f with Gauss (beta)
// Time will be measured

void class_ZZD_one_pass::generating_Gaussian_g_ALICE() {
    this->g_ALICE_old = to_ZZ_pX(gauss_beta->sample_poly(to_int(n)));
    this->g_ALICE = g_ALICE_old % this->f;
}

// Computes two times g_ALICE, saves it in two_g_ALICE
// Time will be measured

void class_ZZD_one_pass::computing_two_times_g_ALICE() {
	this->two_g_ALICE=ZZ_pX();
    this->two_g_ALICE.SetLength(n);
    this->two_g_ALICE = (2 * this->g_ALICE) % this->f;
}

// Computes the multiplicative part of k_ALICE mod f (named k_ALICE_mult)
// Uses p_BOB, hat_r_ALICE and f
// Time will be measured

void class_ZZD_one_pass::computing_mult_part_of_k_ALICE() {

    this->k_ALICE_mult = ZZ_pX();
    this->k_ALICE_mult.SetLength(n);
    this->k_ALICE_mult = (this->p_BOB * this->hat_r_ALICE) % this->f;
}

// Computes the additive part of the polynomial k_ALICE mod f (named k_ALICE)
// Uses k_ALICE_mult, two_g_ALICE and f
// Time will be measured

void class_ZZD_one_pass::computing_add_part_of_k_ALICE() {

    this->k_ALICE = ZZ_pX();
    this->k_ALICE.SetLength(n);
    this->k_ALICE = (this->k_ALICE_mult + this->two_g_ALICE) % this->f;
}

// Comptation of the value of the function Cha for input k_ALICE, saved in Cha
// Time will be measured

void class_ZZD_one_pass::computing_Cha_of_k_ALICE() {

	this->Cha = ZZX();
	this->Cha.SetLength(n);
	this-> p_quarter=p/4;
	this->k_ALICE_mod=conv<ZZX>(this->k_ALICE);

		for (long i=0; i<n; i++) {
			if((this->k_ALICE_mod[i]>(int)(round(this->p_quarter))) && (this->k_ALICE_mod[i]<(p-(int)floor(this->p_quarter)))){
				this->Cha[i]=1;
			}
			else{
				this->Cha[i]=0;
			}
		}
}

// Comptation of the value of the function Mod_2_ALICE for input Cha and k_ALICE, saved in sig_ALICE
// Time will be measured

void class_ZZD_one_pass::evalute_function_Mod2_ALICE(){

	this->sig_ALICE = ZZX();
	this->sig_ALICE.SetLength(n);

	for (long i=0; i<n; i++){
	    this->sig_ALICE[i]=(this->k_ALICE_mod[i]+(this->Cha[i]*(p-1)/2))%this->p;
	    this->sig_ALICE[i]=this->sig_ALICE[i] % 2;
	}
}


//---------------------------------COMPUTATION for BOB after RECEIVING THE ONLY MESSAGE---------------------------

// Generates the polynomial g_BOB over Z_p mod f with Gauss (alpha)
// Time will be measured

void class_ZZD_one_pass::generating_Gaussian_g_BOB() {
    this->g_BOB_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->g_BOB = g_BOB_old % this->f;
}

// Uses d_A from Alice, in reality Bob has to compute it as well

//-------------COMPUTATION OF k_BOB

// First part ((p_ALICE*d_A+x_ALICE)*s_BOB) will be saved in k_BOB_1

// Computes p_ALICE times d_A mod f, saved in k_BOB_1
// Time will be measured

void class_ZZD_one_pass::computing_p_ALICE_times_d_A() {

    this->k_BOB_1 = ZZ_pX();
    this->k_BOB_1.SetLength(n);
    this->k_BOB_1 = (this->p_ALICE * this->d_A) % this->f;
}

// Computes p_ALICE*d_A+x_ALICE mod f which is k_BOB_1+x_ALICE, saved again in k_BOB_1
// Time will be measured

void class_ZZD_one_pass::computing_k_BOB_1_plus_x_ALICE() {

    this->k_BOB_1 = (this->k_BOB_1 + this->x_ALICE) % this->f;
}

// Computes ((p_ALICE*d_A+x_ALICE)*s_BOB) mod f, which is k_BOB_1*s_ALICE, saved again in k_BOB_1
// Time will be measured

void class_ZZD_one_pass::computing_k_BOB_1_times_s_BOB() {

    this->k_BOB_1 = (this->k_BOB_1 * this->s_BOB) % this->f;
}

// Compute now second part of k_BOB, which is 2*d_A*g_BOB, save it in k_BOB_2

// Computes d_A times g_BOB mod f, saved in k_BOB_2
// Time will be measured

void class_ZZD_one_pass::computing_d_A_times_g_BOB() {

    this->k_BOB_2 = ZZ_pX();
    this->k_BOB_2.SetLength(n);
    this->k_BOB_2 = (this->g_BOB * this->d_A) % this->f;
}

// Computes 2*d_A*g_BOB mod f by multiplying k_BOB_2 with 2 mod f, saved in k_BOB_2 again
// Time will be measured

void class_ZZD_one_pass::computing_2_d_A_g_BOB() {

    this->k_BOB_2 = (this->k_BOB_2 * 2) % this->f;
}

// Final computation of k_BOB by adding k_BOB_1 and k_BOB_2
// Time will be measured

void class_ZZD_one_pass::computing_k_BOB_1_plus_k_BOB_2() {

    this->k_BOB = ZZ_pX();
    this->k_BOB.SetLength(n);
    this->k_BOB = (this->k_BOB_1 + this->k_BOB_2) % this->f;
}

// Computation of the value of the function Mod_2 for input Cha and k_BOB, saved in sig_BOB
// Time will be measured

void class_ZZD_one_pass::evalute_function_Mod2_BOB(){

	this->sig_BOB = ZZX();
	this->sig_BOB.SetLength(n);
	this->k_BOB_mod=conv<ZZX>(this->k_BOB);

	for (long i=0; i<n; i++){
	    this->sig_BOB[i]=(this->k_BOB_mod[i]+(this->Cha[i]*(p-1)/2))%this->p;
	    this->sig_BOB[i]=this->sig_BOB[i] % 2;
	}
	cout<< sig_BOB<<endl;
}
