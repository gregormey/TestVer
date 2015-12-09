/*
 * class_BCNS.h
 *
 *  Created on: 30 nov. 2015
 *      Author: Susanne Riess, following Vanessa Erbenich
 *
 *  Implementation of the key exchange protocol described in BCNS14, originally designed by Peikert as a KEM
 */

#ifndef CLASS_BCNS_H_
#define CLASS_BCNS_H_

#include <stdio.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include "class_gauss_selected.h"

class class_BCNS {

private:
    long n;							// The degree of all polynomials in this scheme (especially of f)
    uint64_t p;							// The integer modulus for R_p=Z_p[X]/f
    NTL::ZZ_pX f,a;					// Public parameters of the scheme
    NTL::ZZ_pX e_ALICE_old, e_ALICE, s_ALICE_old, s_ALICE, p_ALICE_mult, p_ALICE;
    								// Variables to create the secret (s_ALICE) and public key (p_ALICE) of Alice
    NTL::ZZ_pX e_BOB_old, e_BOB, s_BOB_old, s_BOB, p_BOB_mult, p_BOB;
        							// Variables to create the secret (s_BOB) and public key (p_BOB) of Bob
    NTL::ZZ_pX e_BOB_2_old, e_BOB_2, v_mult, v;
        							// Variables to create the polynomials e_BOB_2 and v for Bob
    NTL::ZZX ran, dbl_mult, dbl, f_mod;
    								// Degree n-1 polynomial with entries in -1, 0, 0, 1 uniformely (called ran) and dbl=2*v+ran
    double p_half, p_quarter, seven_p_quarter;
    								// Stores p/2, p/4 and 7p/4
    int rounded_p_half, rounded_3p_half;
									// Stores rounded (p/2) and (3p/2)
    NTL::ZZ_pX c, k;				// c=cross round(dbl(v)), k=round(dbl(v))
    NTL::ZZ_pX two_p_BOB, rec;		// rec=2*p_BOB*s_ALICE=two_p_BOB*s_ALICE
    NTL::ZZX key, rec_mod;					// final key=output of rec(rec)
    class_gauss_selected *gauss;	// Stated which Gaussian sampling method is used (the decision is made in class_gauss_selected)

public:
	class_BCNS(long n, uint64_t p, double sigma, long factor, long precision);
// Preparation functions
	void generating_f();
	void generating_uniform_polynomial_a();
// Key generation of Alice
	void generating_Gaussian_error_e_ALICE();
	void generating_Gaussian_secret_key_s_ALICE();
	void computing_mult_part_of_public_key_p_ALICE();
	void computing_add_part_of_public_key_p_ALICE();
// Key generation of Bob
	void generating_Gaussian_error_e_BOB();
	void generating_Gaussian_secret_key_s_BOB();
	void computing_mult_part_of_public_key_p_BOB();
	void computing_add_part_of_public_key_p_BOB();
// Protocol for Bob
	void generating_Gaussian_error_e_BOB_2();
	void computing_mult_part_of_polynomial_v();
	void computing_add_part_of_polynomial_v();
	void generating_random_ternary_polynomial_ran();
	void computing_mult_part_of_double_v();
	void computing_add_part_of_double_v();
	//some constants
	void computing_p_half();
	void computing_p_quarter();
	void computing_7p_quarter();
	void computing_rounded_p_half();
	void computing_rounded_3p_half();
	//rounding functions
	void crossrounding_dbl_v();
	void rounding_dbl_v();
// Protocol for Alice
	void computing_two_p_BOB();
	void computing_two_p_BOB_s_ALICE();
	void computing_key_with_rec();
};

#endif /* CLASS_BCNS_H_ */
