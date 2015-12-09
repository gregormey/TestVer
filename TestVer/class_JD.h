/*
 * class_JD_ring.h
 *
 *  Created on: 1 dec. 2015
 *      Author: Susanne Riess, following Vanessa Erbenich
 *
 *  Implementation of the non-ring-based variant of the key exchange protocol described in JD12
 */

#ifndef CLASS_JD_RING_H_
#define CLASS_JD_RING_H_

#include <stdio.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ_p.h>
#include "class_gauss_selected.h"

class class_JD {

private:
    long n;							// The degree of all polynomials in this scheme (especially of f)
    uint64_t p;							// The integer modulus for R_p=Z_p[X]/f
    NTL::mat_ZZ_p A;				// Matrix A=public parameter of the scheme
    NTL::vec_ZZ_p e_ALICE, s_ALICE, two_e_ALICE, p_ALICE_mult, p_ALICE;
        							// Variables to create the secret key (s_ALICE) and public key (p_ALICE) of Alice
    NTL::vec_ZZ_p e_BOB, s_BOB, two_e_BOB, p_BOB_mult, p_BOB, e_BOB_2;
            						// Variables to create the secret key (s_BOB) and public key (p_BOB) of Bob
    								// and e_BOB_2 is Gaussian number
    NTL::ZZ_p two_e_BOB_2, k_BOB_mult, k_BOB;	//
                					// Variables to create the number k_BOB for Bob
    NTL::ZZ_p b;					// random bit b
    int floor_p_quarter;			// The constant floor(p/4)
    NTL::ZZ Sig, sk_BOB;			// Sig=S(k_BOB), sk_BOB=E(k_BOB,Sig)
    NTL::ZZ k_BOB_mod;				// =k_BOB, for compatibility reasons
    class_gauss_selected *gauss;	// Stated which Gaussian sampling method is used (the decision is made in class_gauss_selected)

public:
    class_JD(long n, uint64_t p, double sigma, long factor, long precision);
    // Preparation functions
    void generating_uniform_A();
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
    void computing_mult_part_of_vector_k_BOB();
    void computing_add_part_of_vector_k_BOB();
    void generating_uniform_bit_b();
    void computing_floor_p_quarter();
    void computing_S_k_BOB();
    void computing_sk_BOB_with_k_BOB_and_Sig();
};

#endif /* CLASS_JD_H_ */
