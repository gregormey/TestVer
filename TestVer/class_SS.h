/*
 * class_SS.h
 *
 *  Created on: 7 dec. 2015
 *      Author: Susanne Riess, following Vanessa Erbenich
 */

#ifndef CLASS_SS_H_
#define CLASS_SS_H_

#include <stdio.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZX.h>
#include "class_gauss_selected.h"

class class_SS {

private:
	long n;							// The degree of all polynomials in this scheme (especially of f)
	uint64_t p;						// The integer modulus for R_p=Z_p[X]/f
	NTL::ZZ_pX f;
	NTL::ZZ_pX f_key, g, h;				// Secret and public key
	NTL::ZZ_pX s,e;
	NTL::ZZ_pX c_1,c;
	NTL::vec_ZZ_p key;
	NTL::ZZ_pX key_mod;
	NTL:: ZZ_pX d;
	NTL::ZZX d_mod;
	class_gauss_selected *gauss, *gauss_sigma;	// Stated which Gaussian sampling method is used (the decision is made in class_gauss_selected)


public:
	class_SS(long n, uint64_t p, double alpha, double sigma, long factor, long precision);
	void generating_f();
	// Key Generation
	void generating_part_of_Gaussian_secret_key_f();
	void generating_mult_part_of_Gaussian_secret_key_f();
	void generating_add_part_of_Gaussian_secret_key_f();
	void generating_part_of_Gaussian_public_key_g();
	void generating_mult_part_of_Gaussian_public_key_h();
	void generating_public_key_h();
	// Encryption
	void generating_Gaussian_s();
	void generating_Gaussian_e();
	void generating_first_mult_part_of_ciphertext_c();
	void generating_second_mult_part_of_ciphertext_c();
	void generating_add_part_of_ciphertext_c();
	void generating_uniform_key();
	void generating_ciphertext_c();
	//Decryption
	void generating_decrypted_c();
	void generating_decrypted_key();
};

#endif /* CLASS_SS_H_ */
