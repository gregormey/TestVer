/*
 * class_JDring.cpp
 *
 *  Created on: 1 dec. 2015
 *      Author: Susanne Riess, following Vanessa Erbenich
 *
 *  Implementation of the ring-based variant of the key exchange protocol described in JD12
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
 *  Before every function you find a short description and a statedment wether the time of this function
 *  will be measured or not.
 *
 *  The time is measured in main.cpp
 *
 *  Which Gaussian sampling method is used, is stated in class_gauss_selected.h and saved in the parameter *gauss
 */

#include "class_JDring.h"
#include <iostream>
#include "class_gauss_selected.h"	// decides the gaussian sampling method
#include <random>					// for uniformly random input
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include "math.h"					// for functions floor and ceil

using namespace std;
using namespace NTL;

class_JDring::class_JDring(long n, int p, double sigma, long factor, long precision) {
    this->n=n, this->p=p;
    ZZ p_new = ZZ(p);
    ZZ_p::init(p_new);

    this->gauss = new class_gauss_selected(sigma, factor, precision);
}


//---------------------- PREPARATION---- GENERATE POLYNOMIAL f AND a---------------------------------------------

// Generates the polynomial f = x^n+1
// which defines R_p=Z_p[X]/f
// Uses n
// No time will be measured

void class_JDring::generating_f() {

    this->f = ZZ_pX();
    this->f.SetLength(n);
    SetCoeff(this->f,0,1);
    SetCoeff(this->f,n,1);
}

// Generates uniformly random polynomial a over Z_p of degree < n
// The polynomial a is a public parameter of the scheme
// Uses ZZ_pX.random
// No time will be measured

void class_JDring::generating_uniform_polynomial_a() {
	this->a = ZZ_pX();
	random(a, n);
      cout << "a" << endl;
      cout << a << endl;
}


//------------------------ KEY GENERATION for ALICE--------------------------------------------------------------

// Generates polynomials e_ALICE and s_ALICE over Z_p mod f with Gauss
// Time will be measured indivually

void class_JDring::generating_Gaussian_error_e_ALICE() {
    this->e_ALICE_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->e_ALICE = e_ALICE_old % this->f;
    cout << "e_ALICE" << endl;
        cout << e_ALICE << endl;
}

void class_JDring::generating_Gaussian_secret_key_s_ALICE() {
    this->s_ALICE_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->s_ALICE = s_ALICE_old % this->f;
    cout << "s_ALICE" << endl;
        cout << s_ALICE << endl;
}

// Computes two times e_ALICE, saves it in two_e_ALICE
// Time will be measured

void class_JDring::computing_two_times_e_ALICE() {
	this->two_e_ALICE=ZZ_pX();
    this->two_e_ALICE.SetLength(n);
    this->two_e_ALICE = (2 * this->e_ALICE) % this->f;
    cout << "two_e_ALICE" << endl;
    cout << two_e_ALICE << endl;
}

// Computes the multiplicative part of the public key p_ALICE of Alice mod f (named p_ALICE_mult)
// Uses s_ALICE, a and f
// Time will be measured

void class_JDring::computing_mult_part_of_public_key_p_ALICE() {

    this->p_ALICE_mult = ZZ_pX();
    this->p_ALICE_mult.SetLength(n);
    this->p_ALICE_mult = (this->a * this->s_ALICE) % this->f;
    cout << "p_ALICE_mult" << endl;
    cout << p_ALICE_mult << endl;
}

// Computes the additive part of the public key p_ALICE of Alice mod f (named p_ALICE)
// Uses two_e_ALICE, p_ALICE_mult and f
// Time will be measured

void class_JDring::computing_add_part_of_public_key_p_ALICE() {

    this->p_ALICE = ZZ_pX();
    this->p_ALICE.SetLength(n);
    this->p_ALICE = (this->p_ALICE_mult + this->two_e_ALICE) % this->f;
    cout << "p_ALICE" << endl;
    cout << p_ALICE << endl;
}


//------------------------ KEY GENERATION for BOB----------------------------------------------------------------

// Generates polynomials e_BOB and s_BOB over Z_p mod f with Gauss
// Time will be measured indivually

void class_JDring::generating_Gaussian_error_e_BOB() {
    this->e_BOB_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->e_BOB = e_BOB_old % this->f;
}

void class_JDring::generating_Gaussian_secret_key_s_BOB() {
    this->s_BOB_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->s_BOB = s_BOB_old % this->f;
}

// Computes two times e_BOB, saves it in two_e_BOB
// Time will be measured

void class_JDring::computing_two_times_e_BOB() {
	this->two_e_BOB=ZZ_pX();
    this->two_e_BOB.SetLength(n);
    this->two_e_BOB = (2 * this->e_BOB) % this->f;
}

// Computes the multiplicative part of the public key p_BOB of BOB mod f (named p_BOB_mult)
// Uses s_BOB, a and f
// Time will be measured

void class_JDring::computing_mult_part_of_public_key_p_BOB() {

    this->p_BOB_mult = ZZ_pX();
    this->p_BOB_mult.SetLength(n);
    this->p_BOB_mult = (this->a * this->s_BOB) % this->f;
}

// Computes the additive part of the public key p_BOB of BOB mod f (named p_BOB)
// Uses two_e_BOB, p_BOB_mult and f
// Time will be measured

void class_JDring::computing_add_part_of_public_key_p_BOB() {

    this->p_BOB = ZZ_pX();
    this->p_BOB.SetLength(n);
    this->p_BOB = (this->p_BOB_mult + this->two_e_BOB) % this->f;
}


//--------------------------- THE PROTOCOL for BOB---------------------------------------------------------------

// Generates the polynomial e_BOB_2 over Z_p mod f with Gauss
// Time will be measured

void class_JDring::generating_Gaussian_error_e_BOB_2() {
    this->e_BOB_2_old = to_ZZ_pX(gauss->sample_poly(to_int(n)));
    this->e_BOB_2 = e_BOB_2_old % this->f;
    cout << "e_BOB_2" << endl;
    cout << e_BOB_2 << endl;
}

// Computes two times e_BOB_2, saves it in two_e_BOB_2
// Time will be measured

void class_JDring::computing_two_times_e_BOB_2() {
	this->two_e_BOB_2=ZZ_pX();
    this->two_e_BOB_2.SetLength(n);
    this->two_e_BOB_2 = (2 * this->e_BOB_2) % this->f;
}

// Computes the multiplicative part of the polynomial k_BOB mod f (named k_BOB_mult)
// Uses s_BOB, p_ALICE and f
// Time will be measured

void class_JDring::computing_mult_part_of_polynomial_k_BOB() {

    this->k_BOB_mult = ZZ_pX();
    this->k_BOB_mult.SetLength(n);
    this->k_BOB_mult = (this->p_ALICE * this->s_BOB) % this->f;
    cout << "k_BOB_mult" << endl;
    cout << k_BOB_mult << endl;
}

// Computes the additive part of the polynomial k_BOB mod f (named k_BOB)
// Uses two_e_BOB_2, k_BOB_mult and f
// Time will be measured

void class_JDring::computing_add_part_of_polynomial_k_BOB() {

    this->k_BOB = ZZ_pX();
    this->k_BOB.SetLength(n);
    this->k_BOB = (this->k_BOB_mult + this->two_e_BOB_2) % this->f;
    cout << "k_BOB" << endl;
    cout << k_BOB << endl;
}

// Generates polynomial b with uniformly random entries 0 or 1
// Uses uniform_int_distribution
// Time will be measured

void class_JDring::generating_uniform_bit_polynomial_b() {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 1);

    this->b = vec_ZZ_p();
    this->b.SetLength(n);

    for (long i=1; i<=n; i++)
        this->b(i) = ZZ_p(dis(gen));
    	cout << b << endl;

}

// Computes floor(p/4)
// Time will be measured individually

void class_JDring::computing_floor_p_quarter() {
	this->floor_p_quarter=floor_p_quarter;
	this-> floor_p_quarter = floor(this->p/4);
}

// Computes the value of function S for input k_BOB and for random bits b
// Return value will be saved in Sig
// Time will be measured

void class_JDring::computing_S_k_BOB() {

	this->Sig = ZZX();
	this->Sig.SetLength(n);
	this->k_BOB_mod = ZZX();
	this->k_BOB_mod.SetLength(n);
	this->k_BOB_mod=conv<ZZX>(this->k_BOB); 	// for compatibility
	cout << Sig << endl;
	cout << k_BOB_mod << endl;


	for (long i=0; i<n; i++) {
		cout << "b[i]" << endl;
		cout << b[i] << endl;
		cout << i << endl;
		if (this->b[i]==0) {
			cout << "k_BOB_mod[i]" << endl;
			cout << k_BOB_mod[i] << endl;
			cout << floor_p_quarter << endl;
			cout << this->p-this->floor_p_quarter<< endl;
			if((this->k_BOB_mod[i]>this->floor_p_quarter) && ((this->k_BOB_mod[i])<(this->p-this->floor_p_quarter))){
				this->Sig[i]=1;
				cout << "Sig[i]" << endl;
				cout << Sig[i] << endl;
			}
			else{
				this->Sig[i]=0;
				cout << "Sig[i]" << endl;
				cout << Sig[i] << endl;
			}
		}
		else{
			if((this->k_BOB_mod[i]>(this->floor_p_quarter+1)) && (this->k_BOB_mod[i]<(p+this->floor_p_quarter))){
				this->Sig[i]=1;
			}
			else{
				this->Sig[i]=0;
			}
		}
	}
	cout<< Sig<<endl;
}

// Computes the value of the function E for input k_BOB and Sig
// Return value will be saved in sk_BOB
// Time will be measured

void class_JDring::computing_sk_BOB_with_k_BOB_and_Sig(){

	this->sk_BOB = ZZX();
	this->sk_BOB.SetLength(n);

	for (long i=0; i<n; i++){
	    this->sk_BOB[i]=this->k_BOB_mod[i]+(this->Sig[i]*(int)(p-1)/2);
	    this->sk_BOB[i]=this->sk_BOB[i] % 2;
	}
	cout << sk_BOB << endl;
}
