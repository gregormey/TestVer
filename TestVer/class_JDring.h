/*
 * class_JDring.h
 *
 *  Created on: 1 dec. 2015
 *      Author: Susanne Riess, following Vanessa Erbenich
 *
 *  Implementation of the ring-based variant of the key exchange protocol described in JD12
 */

#ifndef CLASS_JDring_H_
#define CLASS_JDring_H_

#include <stdio.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include "class_gauss_selected.h"
#include <NTL/vec_ZZ.h>

class class_JDring {

private:
    long n;							// The degree of all polynomials in this scheme (especially of f)
    int p;							// The integer modulus for R_p=Z_p[X]/f
    NTL::ZZ_pX f,a;					// Public parameters of the scheme
    NTL::ZZ_pX e_ALICE_old, e_ALICE, s_ALICE_old, s_ALICE, two_e_ALICE, p_ALICE_mult, p_ALICE;
    								// Variables to create the secret (s_ALICE) and public key (p_ALICE) of Alice
    NTL::ZZ_pX e_BOB_old, e_BOB, s_BOB_old, s_BOB, two_e_BOB, p_BOB_mult, p_BOB;
            						// Variables to create the secret (s_BOB) and public key (p_BOB) of Bob
    NTL::ZZ_pX e_BOB_2_old, e_BOB_2, two_e_BOB_2, k_BOB_mult, k_BOB;
            						// Variables to create the polynomials e_BOB_2 and k_BOB for Bob
    int floor_p_quarter;			// The constant floor(p/4)
    NTL::vec_ZZ_p b;
    NTL::ZZX  Sig, sk_BOB;		// b is a polynomial with coefficients 1 and 0, Sig=S(k_BOB), sk_BOB=E(k_BOB,Sig)
    NTL::ZZX k_BOB_mod;				// =k_BOB, for compatibility reasons
    class_gauss_selected *gauss;	// Stated which Gaussian sampling method is used (the decision is made in class_gauss_selected)

public:
    class_JDring(long n, int p, double sigma, long factor, long precision);
    // Preparation functions
    void generating_f();
    void generating_uniform_polynomial_a();
    // Key generation of Alice
    void generating_Gaussian_error_e_ALICE();
    void generating_Gaussian_secret_key_s_ALICE();
    void computing_two_times_e_ALICE();
    void computing_mult_part_of_public_key_p_ALICE();
    void computing_add_part_of_public_key_p_ALICE();
    // Key generation of Bob
    void generating_Gaussian_error_e_BOB();
    void generating_Gaussian_secret_key_s_BOB();
    void computing_two_times_e_BOB();
    void computing_mult_part_of_public_key_p_BOB();
    void computing_add_part_of_public_key_p_BOB();
    // Protocol for Bob
    void generating_Gaussian_error_e_BOB_2();
    void computing_two_times_e_BOB_2();
    void computing_mult_part_of_polynomial_k_BOB();
    void computing_add_part_of_polynomial_k_BOB();
    void generating_uniform_bit_polynomial_b();
    void computing_floor_p_quarter();
    void computing_S_k_BOB();
    void computing_sk_BOB_with_k_BOB_and_Sig();
};

#endif /* CLASS_JDring_H_ */
