/*
 * class_LPR.h
 *
 *  Created on: 3 dec. 2015
 *      Author: Susanne Riess, following Vanessa Erbenich
 */

#ifndef CLASS_LPR_H_
#define CLASS_LPR_H_

#include <stdio.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZX.h>
#include "class_gauss_selected.h"


class class_LPR {

private:
	long n;							// The degree of all polynomials in this scheme (especially of f)
	uint64_t p;							// The integer modulus for R_p=Z_p[X]/f
    NTL::ZZ_pX f,a;					// Public parameters of the scheme
    NTL::ZZ_pX e_A_old, e_A, s_A_old, s_A, p_A_mult, p_A;
    								// Variables to create the secret (s_A) and public key (p_A)
    NTL::vec_ZZ_p key; 				// Variable to create the key
	NTL::ZZ_pX key_preenc, key_enc;
    								// Variables to store the key and its encryption
    NTL::ZZ_pX s_B, e_B, e_B_1;		// Variables for Gaussian errors
    NTL::ZZ_pX p_B_mult, c_mult, p_B, c;
    								// Variables for ciphertexts p_B and c
    NTL::ZZ_pX decrypt_mult, decrypt;// Variable for decryption process
    NTL::ZZX decrypt_int;			// For compatibility
    NTL::ZZX key_dec;			// the decrypted key
    class_gauss_selected *gauss;	// Stated which Gaussian sampling method is used (the decision is made in class_gauss_selected)

public:
    class_LPR(long n, uint64_t p, double sigma, long factor, long precision);
    void generating_f();
    void generating_uniform_polynomial_a();
    void generating_Gaussian_error_e_A();
    void generating_Gaussian_secret_key_s_A();
    void computing_mult_part_of_public_key_p_A();
    void computing_add_part_of_public_key_p_A();
    void generating_uniform_key();
    void encoding_key();
    void generating_Gaussian_error_s_B();
    void generating_Gaussian_error_e_B();
    void generating_Gaussian_error_e_B_1();
    void computing_mult_part_of_ciphertext_p_B();
    void computing_mult_part_of_ciphertext_c();
    void computing_add_part_of_ciphertext_p_B();
    void computing_add_part_of_ciphertext_c();
    void computing_mult_part_of_decrypt();
    void computing_add_part_of_decrypt();
    void converting_decrypt();
    void decoding_decrypt();
};

#endif /* CLASS_LPR_H_ */
