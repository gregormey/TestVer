/*
 * ZZD_one_pass.h
 *
 *  Created on: 3 dec. 2015
 *      Author: Susanne Riess, following Vanessa Erbenich
 *
 *  Implementation of the one-pass variant of the key exchange protocol described in ZZD+14
 */

#ifndef CLASS_ZZD_one_PASS_H_
#define CLASS_ZZD_one_PASS_H_

#include <stdio.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include "class_gauss_selected.h"

class class_ZZD_one_pass {

private:
	long n;							// The degree of all polynomials in this scheme (especially of f)
	uint64_t p;							// The integer modulus for R_p=Z_p[X]/f
	NTL::ZZ_pX f,a;					// Public parameters of the scheme
	NTL::ZZ_pX e_ALICE_old, e_ALICE, s_ALICE_old, s_ALICE, two_e_ALICE, p_ALICE_mult, p_ALICE;
    								// Variables to create the secret (s_ALICE) and public key (p_ALICE) of Alice
    NTL::ZZ_pX e_BOB_old, e_BOB, s_BOB_old, s_BOB, two_e_BOB, p_BOB_mult, p_BOB;
            						// Variables to create the secret (s_BOB) and public key (p_BOB) of Bob
    NTL::ZZ_pX r_ALICE_old, r_ALICE, f_ALICE_old, f_ALICE, two_f_ALICE, x_ALICE_mult, x_ALICE;
                					// Variables to create r_ALICE, f_ALICE and the message x_ALICE of Alice
    NTL::ZZ_pX d_A_old, d_A, hat_r_ALICE_mult, hat_f_ALICE_mult, hat_r_ALICE, hat_f_ALICE;
                  					// Variables to create hat_r_ALICE, hat_f_ALICE and d_A for Alice
    NTL::ZZ_pX g_ALICE_old, g_ALICE;// Variable for Gaussian g_ALICE
    NTL::ZZ_pX two_g_ALICE,k_ALICE_mult, k_ALICE;
      								// Variable to compute k_ALICE
    int p_quarter;					// constant p/4
    NTL::ZZX Cha, sig_ALICE;		// Cha=cha(k_ALICE), sig_ALICE=Mod2(k_ALICE, Cha)
    NTL::ZZX k_ALICE_mod;			// =k_ALICE, for compatibility reasons
    NTL::ZZ_pX g_BOB_old, g_BOB;	// Variable for Gaussian g_BOB
    NTL::ZZ_pX k_BOB_1, k_BOB_2, k_BOB;
    								// Variable to compute k_BOB
    NTL::ZZX sig_BOB;				// sig_BOB=Mod2(k_BOB, Cha)
    NTL::ZZX k_BOB_mod;				// =k_BOB, for compatibility reasons
    class_gauss_selected *gauss, *gauss_beta;;
									// States which Gaussian sampling method is used (the decision is made in class_gauss_selected)

public:
	class_ZZD_one_pass(long n, uint64_t p, double alpha, double beta, long factor, long precision);
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
    // Computation ALICE before first PASS
    void generating_Gaussian_r_ALICE();
    void generating_Gaussian_f_ALICE();
    void computing_two_times_f_ALICE();
    void computing_mult_part_of_x_ALICE();
    void computing_add_part_of_x_ALICE();
    void generating_Gaussian_d_A();
    void computing_mult_part_of_hat_r_ALICE();
    void computing_add_part_of_hat_r_ALICE();
    void computing_mult_part_of_hat_f_ALICE();
    void computing_add_part_of_hat_f_ALICE();
    	// Computation k_ALICE
    void generating_Gaussian_g_ALICE();
    void computing_two_times_g_ALICE();
    void computing_mult_part_of_k_ALICE();
    void computing_add_part_of_k_ALICE();
    	// Computation of Cha and sig_ALICE
    void computing_Cha_of_k_ALICE();
    void evalute_function_Mod2_ALICE();
    // Computation BOB after receiving the only message
    void generating_Gaussian_g_BOB();
    	// Computation of k_BOB
    void computing_p_ALICE_times_d_A();
    void computing_k_BOB_1_plus_x_ALICE();
    void computing_k_BOB_1_times_s_BOB();
    void computing_d_A_times_g_BOB();
    void computing_2_d_A_g_BOB();
    void computing_k_BOB_1_plus_k_BOB_2();
    	// Computation of sig_BOB
    void evalute_function_Mod2_BOB();
};

#endif /* CLASS_ZZD_one_PASS_H_ */
