#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_RR.h>

#include <iostream>
#include <ctime>
#include <chrono>

#include "class_BCNS.h"
#include "class_JD.h"
#include "class_JDring.h"
#include "class_SS.h"
#include "class_LPR.h"
#include "class_ZZD_three_pass.h"
#include "class_ZZD_one_pass.h"

#include "gauss.h"
#include "gauss_knuth_yao.h"
#include "ziggurat-new.h"
#include "class_gauss_inverse_cdf.h"
#include "class_gauss_ziggurat.h"
#include "class_gauss_rejection.h"
#include "class_gauss_knuth_yao.h"



//hier ist p, was in der Arbeit q ist und p_hat, was in der Arbeit p ist, kann man auch ZZ_qX etc schreiben?

using namespace std;
using namespace NTL;


int main()
{
    int scheme_yes_no;
    
    cout << "Would you like to run all schemes (1) or all samplers(0)?" << endl;
    cin >> scheme_yes_no;
    
    // 0: sampler
    // 1: schemes
    
    if(scheme_yes_no==1){
        
        //--------------------------------------------- Selection of Parameters ---------------------------------------------//
        
        long nr_of_runs;
        
        cout << "number of runs?" << endl;
        cin >> nr_of_runs;

        long n_JD=411,factor_JD=13,precision_JD=128;
        uint64_t p_JD=841741;
        double sigma_JD=3.192;

        long n_JDring=512,factor_JDring=13,precision_JDring=128;
        int p_JDring=841741;
        double sigma_JDring=3.192;

        long n_BCNS=512,factor_BCNS=13,precision_BCNS=128;
        uint64_t p_BCNS=46565383;
        double sigma_BCNS=3.192;

        long n_LPR=512,factor_LPR=13,precision_LPR=128;
        uint64_t p_LPR=46565383;
        double sigma_LPR=3.192;

        long n_SS=1024,factor_SS=13,precision_SS=128;
        uint64_t p_SS=103279992729601;
        double sigma_SS=267445084617, alpha_SS=3.192;

        long n_ZZD_three_pass=1024,factor_ZZD_three_pass=13,precision_ZZD_three_pass=128;
        uint64_t p_ZZD_three_pass=14186338877441;
        double alpha_ZZD_three_pass=3.192, beta_ZZD_three_pass=62914.56;

        long n_ZZD_one_pass=512,factor_ZZD_one_pass=13,precision_ZZD_one_pass=128;
        uint64_t p_ZZD_one_pass=127556609;
        double alpha_ZZD_one_pass=3.192, beta_ZZD_one_pass=31457.28;


        int choose_para;
        cout << "choose parameters (1) or not (0)?" << endl;
        cin >> choose_para;
        
        if (choose_para==1){
            
        	cout << "n_JD?" << endl;
        	cin >> n_JD;
        	cout << "p_JD?" << endl;
        	cin >> p_JD;
        	cout << "sigma_JD?" << endl;
        	cin >> sigma_JD;
        	cout << "factor_JD?" << endl;
        	cin >> factor_JD;
        	cout << "precision_JD?" << endl;
        	cin >> precision_JD;

        	cout << "n_JDring?" << endl;
        	cin >> n_JDring;
        	cout << "p_JDring?" << endl;
        	cin >> p_JDring;
        	cout << "sigma_JDring?" << endl;
        	cin >> sigma_JDring;
        	cout << "factor_JDring?" << endl;
        	cin >> factor_JDring;
        	cout << "precision_JDring?" << endl;
        	cin >> precision_JDring;

            cout << "n_BCNS?" << endl;
            cin >> n_BCNS;
            cout << "p_BCNS?" << endl;
            cin >> p_BCNS;
            cout << "sigma_BCNS?" << endl;
            cin >> sigma_BCNS;
            cout << "factor_BCNS?" << endl;
            cin >> factor_BCNS;
            cout << "precision_BCNS?" << endl;
            cin >> precision_BCNS;

            cout << "n_LPR?" << endl;
            cin >> n_LPR;
            cout << "p_LPR?" << endl;
            cin >> p_LPR;
            cout << "sigma_LPR?" << endl;
            cin >> sigma_LPR;
            cout << "factor_LPR?" << endl;
            cin >> factor_LPR;
            cout << "precision_LPR?" << endl;
            cin >> precision_LPR;

            cout << "n_SS?" << endl;
            cin >> n_SS;
            cout << "p_SS?" << endl;
            cin >> p_SS;
            cout << "sigma_SS?" << endl;
            cin >> sigma_SS;
            cout << "alpha_SS?" << endl;
            cin >> alpha_SS;
            cout << "factor_SS?" << endl;
            cin >> factor_SS;
            cout << "precision_SS?" << endl;
            cin >> precision_SS;

            cout << "n_ZZD_three_pass?" << endl;
            cin >> n_ZZD_three_pass;
            cout << "p_ZZD_three_pass?" << endl;
            cin >> p_ZZD_three_pass;
            cout << "alpha_ZZD_three_pass?" << endl;
            cin >> alpha_ZZD_three_pass;
            cout << "beta_ZZD_three_pass?" << endl;
            cin >> beta_ZZD_three_pass;
            cout << "factor_ZZD_three_pass?" << endl;
            cin >> factor_ZZD_three_pass;
            cout << "precision_ZZD_three_pass?" << endl;
            cin >> precision_ZZD_three_pass;

            cout << "n_ZZD_one_pass?" << endl;
            cin >> n_ZZD_one_pass;
            cout << "p_ZZD_one_pass?" << endl;
            cin >> p_ZZD_one_pass;
            cout << "alpha_ZZD_one_pass?" << endl;
            cin >> alpha_ZZD_one_pass;
            cout << "beta_ZZD_one_pass?" << endl;
            cin >> beta_ZZD_one_pass;
            cout << "factor_ZZD_one_pass?" << endl;
            cin >> factor_ZZD_one_pass;
            cout << "precision_ZZD_one_pass?" << endl;
            cin >> precision_ZZD_one_pass;

        }
        

//------------------------------------KEY EXCHANGE JD----------------------------------------------------------

        class_JD JD = class_JD(n_JD,p_JD,sigma_JD,factor_JD,precision_JD);
/*
        cout << "" << endl;
        cout << "JD" << endl;
        cout << "" << endl;

        double JD_time_generating_Gaussian_error_e_ALICE=0;
        double JD_time_generating_Gaussian_secret_key_s_ALICE=0;
        double JD_time_computing_two_times_e_ALICE=0;
        double JD_time_computing_mult_part_of_public_key_p_ALICE=0;
        double JD_time_computing_add_part_of_public_key_p_ALICE=0;
    // Key generation of Bob
        double JD_time_generating_Gaussian_error_e_BOB=0;
        double JD_time_generating_Gaussian_secret_key_s_BOB=0;
        double JD_time_computing_two_times_e_BOB=0;
        double JD_time_computing_mult_part_of_public_key_p_BOB=0;
        double JD_time_computing_add_part_of_public_key_p_BOB=0;
    // Protocol for Bob
        double JD_time_generating_Gaussian_error_e_BOB_2=0;
        double JD_time_computing_two_times_e_BOB_2=0;
        double JD_time_computing_mult_part_of_vector_k_BOB=0;
        double JD_time_computing_add_part_of_vector_k_BOB=0;
        double JD_time_generating_uniform_bit_b=0;
        double JD_time_computing_floor_p_quarter=0;
        double JD_time_computing_S_k_BOB=0;
        double JD_time_computing_sk_BOB_with_k_BOB_and_Sig=0;

        for (int k=1; k<=nr_of_runs; k++) {

        	// ----------------------------Preparation functions-------------------------------------------------
        	JD.generating_uniform_A();

        	// ----------------------------Key generation of Alice-----------------------------------------------
        	clock_t start = clock();
        	JD.generating_Gaussian_error_e_ALICE();
        	clock_t end = clock();
        	JD_time_generating_Gaussian_error_e_ALICE += double(end)-double(start);

        	start = clock();
        	JD.generating_Gaussian_secret_key_s_ALICE();
        	end = clock();
        	JD_time_generating_Gaussian_secret_key_s_ALICE += double(end)-double(start);

        	start = clock();
        	JD.computing_two_times_e_ALICE();
        	end = clock();
        	JD_time_computing_two_times_e_ALICE += double(end)-double(start);

        	start = clock();
        	JD.computing_mult_part_of_public_key_p_ALICE();
        	end = clock();
        	JD_time_computing_mult_part_of_public_key_p_ALICE += double(end)-double(start);

        	start = clock();
        	JD.computing_add_part_of_public_key_p_ALICE();
        	end = clock();
        	JD_time_computing_add_part_of_public_key_p_ALICE += double(end)-double(start);


        	// ----------------------------Key generation of Bob-------------------------------------------------
        	start = clock();
        	JD.generating_Gaussian_error_e_BOB();
        	end = clock();
        	JD_time_generating_Gaussian_error_e_BOB += double(end)-double(start);

        	start = clock();
        	JD.generating_Gaussian_secret_key_s_BOB();
        	end = clock();
        	JD_time_generating_Gaussian_secret_key_s_BOB += double(end)-double(start);

        	start = clock();
        	JD.computing_two_times_e_BOB();
        	end = clock();
        	JD_time_computing_two_times_e_BOB += double(end)-double(start);


        	start = clock();
        	JD.computing_mult_part_of_public_key_p_BOB();
        	end = clock();
        	JD_time_computing_mult_part_of_public_key_p_BOB += double(end)-double(start);

        	start = clock();
        	JD.computing_add_part_of_public_key_p_BOB();
        	end = clock();
        	JD_time_computing_add_part_of_public_key_p_BOB += double(end)-double(start);


        	// ---------------------------Protocol for Bob-------------------------------------------------------
        	start = clock();
        	JD.generating_Gaussian_error_e_BOB_2();
        	end = clock();
        	JD_time_generating_Gaussian_error_e_BOB_2 += double(end)-double(start);

        	start = clock();
        	JD.computing_two_times_e_BOB_2();
        	end = clock();
        	JD_time_computing_two_times_e_BOB_2 += double(end)-double(start);

        	start = clock();
        	JD.computing_mult_part_of_vector_k_BOB();
        	end = clock();
        	JD_time_computing_mult_part_of_vector_k_BOB += double(end)-double(start);

        	start = clock();
        	JD.computing_add_part_of_vector_k_BOB();
        	end = clock();
        	JD_time_computing_add_part_of_vector_k_BOB += double(end)-double(start);

        	start = clock();
        	JD.generating_uniform_bit_b();
        	end = clock();
        	JD_time_generating_uniform_bit_b += double(end)-double(start);

        	start = clock();
        	JD.computing_floor_p_quarter();
        	end = clock();
        	JD_time_computing_floor_p_quarter += double(end)-double(start);

        	start = clock();
        	JD.computing_S_k_BOB();
        	end = clock();
        	JD_time_computing_S_k_BOB += double(end)-double(start);

        	start = clock();
        	JD.computing_sk_BOB_with_k_BOB_and_Sig();
        	end = clock();
        	JD_time_computing_sk_BOB_with_k_BOB_and_Sig += double(end)-double(start);

        }

        cout << "" << endl;
        cout << "Key Generation Alice:" << endl;
        cout << "JD_time_generating_Gaussian_error_e_ALICE:" << 1000.0 * JD_time_generating_Gaussian_error_e_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_generating_Gaussian_secret_key_s_ALICE:" << 1000.0 * JD_time_generating_Gaussian_secret_key_s_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_computing_two_times_e_ALICE:" << 1000.0 * JD_time_computing_two_times_e_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_computing_mult_part_of_public_key_p_ALICE:" << 1000.0 * JD_time_computing_mult_part_of_public_key_p_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_computing_add_part_of_public_key_p_ALICE:" << 1000.0 * JD_time_computing_add_part_of_public_key_p_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

        cout << "Key Generation Bob:" << endl;
        cout << "JD_time_generating_Gaussian_error_e_BOB:" << 1000.0 * JD_time_generating_Gaussian_error_e_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_generating_Gaussian_secret_key_s_BOB:" << 1000.0 * JD_time_generating_Gaussian_secret_key_s_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_computing_two_times_e_BOB:" << 1000.0 * JD_time_computing_two_times_e_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_computing_mult_part_of_public_key_p_BOB:" << 1000.0 * JD_time_computing_mult_part_of_public_key_p_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_computing_add_part_of_public_key_p_BOB:" << 1000.0 * JD_time_computing_add_part_of_public_key_p_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

        cout << "The Protocol for Bob:" << endl;
        cout << "JD_time_generating_Gaussian_error_e_BOB_2:" << 1000.0 * JD_time_generating_Gaussian_error_e_BOB_2 / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_computing_two_times_e_BOB_2:" << 1000.0 * JD_time_computing_two_times_e_BOB_2 / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_computing_mult_part_of_vector_k_BOB:" << 1000.0 * JD_time_computing_mult_part_of_vector_k_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_computing_add_part_of_vector_k_BOB:" << 1000.0 * JD_time_computing_add_part_of_vector_k_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_generating_uniform_bit_b:" << 1000.0 * JD_time_generating_uniform_bit_b / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_computing_floor_p_quarter:" << 1000.0 * JD_time_computing_floor_p_quarter / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_computing_S_k_BOB:" << 1000.0 * JD_time_computing_S_k_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JD_time_computing_sk_BOB_with_k_BOB_and_Sig:" << 1000.0 * JD_time_computing_sk_BOB_with_k_BOB_and_Sig / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "" << endl;
*/

//------------------------------------KEY EXCHANGE JDring----------------------------------------------------------

        class_JDring JDring= class_JDring(n_JDring, p_JDring, sigma_JDring, factor_JDring, precision_JDring);
/*
        cout << "" << endl;
        cout << "JDring" << endl;
        cout << "" << endl;

        double JDring_time_generating_Gaussian_error_e_ALICE=0;
        double JDring_time_generating_Gaussian_secret_key_s_ALICE=0;
        double JDring_time_computing_two_times_e_ALICE=0;
        double JDring_time_computing_mult_part_of_public_key_p_ALICE=0;
        double JDring_time_computing_add_part_of_public_key_p_ALICE=0;
        // Key generation of Bob
        double JDring_time_generating_Gaussian_error_e_BOB=0;
        double JDring_time_generating_Gaussian_secret_key_s_BOB=0;
        double JDring_time_computing_two_times_e_BOB=0;
        double JDring_time_computing_mult_part_of_public_key_p_BOB=0;
        double JDring_time_computing_add_part_of_public_key_p_BOB=0;
        // Protocol for Bob
        double JDring_time_generating_Gaussian_error_e_BOB_2=0;
        double JDring_time_computing_two_times_e_BOB_2=0;
        double JDring_time_computing_mult_part_of_polynomial_k_BOB=0;
        double JDring_time_computing_add_part_of_polynomial_k_BOB=0;
        double JDring_time_generating_uniform_bit_polynomial_b=0;
        double JDring_time_computing_floor_p_quarter=0;
        double JDring_time_computing_S_k_BOB=0;
        double JDring_time_computing_sk_BOB_with_k_BOB_and_Sig=0;

        for (int k=1; k<=nr_of_runs; k++) {

        	// ----------------------------Preparation functions-------------------------------------------------
            JDring.generating_f();
            JDring.generating_uniform_polynomial_a();

            // ----------------------------Key generation of Alice-----------------------------------------------
            clock_t start = clock();
            JDring.generating_Gaussian_error_e_ALICE();
            clock_t end = clock();
            JDring_time_generating_Gaussian_error_e_ALICE += double(end)-double(start);

            start = clock();
            JDring.generating_Gaussian_secret_key_s_ALICE();
            end = clock();
            JDring_time_generating_Gaussian_secret_key_s_ALICE += double(end)-double(start);

            start = clock();
            JDring.computing_two_times_e_ALICE();
            end = clock();
            JDring_time_computing_two_times_e_ALICE += double(end)-double(start);

            start = clock();
            JDring.computing_mult_part_of_public_key_p_ALICE();
            end = clock();
            JDring_time_computing_mult_part_of_public_key_p_ALICE += double(end)-double(start);

            start = clock();
            JDring.computing_add_part_of_public_key_p_ALICE();
            end = clock();
            JDring_time_computing_add_part_of_public_key_p_ALICE += double(end)-double(start);


            // ----------------------------Key generation of Bob-------------------------------------------------
            start = clock();
            JDring.generating_Gaussian_error_e_BOB();
            end = clock();
            JDring_time_generating_Gaussian_error_e_BOB += double(end)-double(start);

            start = clock();
            JDring.generating_Gaussian_secret_key_s_BOB();
            end = clock();
            JDring_time_generating_Gaussian_secret_key_s_BOB += double(end)-double(start);

            start = clock();
            JDring.computing_two_times_e_BOB();
            end = clock();
            JDring_time_computing_two_times_e_BOB += double(end)-double(start);

            start = clock();
            JDring.computing_mult_part_of_public_key_p_BOB();
            end = clock();
            JDring_time_computing_mult_part_of_public_key_p_BOB += double(end)-double(start);

            start = clock();
            JDring.computing_add_part_of_public_key_p_BOB();
            end = clock();
            JDring_time_computing_add_part_of_public_key_p_BOB += double(end)-double(start);

            cout << "JDring" << endl;
            // ---------------------------Protocol for Bob-------------------------------------------------------
            start = clock();
            JDring.generating_Gaussian_error_e_BOB_2();
            end = clock();
            JDring_time_generating_Gaussian_error_e_BOB_2 += double(end)-double(start);

            start = clock();
            JDring.computing_two_times_e_BOB_2();
            end = clock();
            JDring_time_computing_two_times_e_BOB_2 += double(end)-double(start);

            start = clock();
            JDring.computing_mult_part_of_polynomial_k_BOB();
            end = clock();
            JDring_time_computing_mult_part_of_polynomial_k_BOB += double(end)-double(start);

            start = clock();
            JDring.computing_add_part_of_polynomial_k_BOB();
            end = clock();
            JDring_time_computing_add_part_of_polynomial_k_BOB += double(end)-double(start);

            cout << "JDring" << endl;

            start = clock();
            JDring.generating_uniform_bit_polynomial_b();
            end = clock();
            JDring_time_generating_uniform_bit_polynomial_b += double(end)-double(start);

            cout << "JDring" << endl;

            start = clock();
            JDring.computing_floor_p_quarter();
            end = clock();
            JDring_time_computing_floor_p_quarter += double(end)-double(start);

            cout << "JDring" << endl;

            start = clock();
            JDring.computing_S_k_BOB();
            end = clock();
            JDring_time_computing_S_k_BOB += double(end)-double(start);

            cout << "JDring" << endl;

            start = clock();
            JDring.computing_sk_BOB_with_k_BOB_and_Sig();
            end = clock();
            JDring_time_computing_sk_BOB_with_k_BOB_and_Sig += double(end)-double(start);

            cout << "JDring" << endl;

        }

        cout << "" << endl;
        cout << "Key Generation Alice:" << endl;
        cout << "JDring_time_generating_Gaussian_error_e_ALICE:" << 1000.0 * JDring_time_generating_Gaussian_error_e_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_generating_Gaussian_secret_key_s_ALICE:" << 1000.0 * JDring_time_generating_Gaussian_secret_key_s_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_computing_two_times_e_ALICE:" << 1000.0 * JDring_time_computing_two_times_e_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_computing_mult_part_of_public_key_p_ALICE:" << 1000.0 * JDring_time_computing_mult_part_of_public_key_p_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_computing_add_part_of_public_key_p_ALICE:" << 1000.0 * JDring_time_computing_add_part_of_public_key_p_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

        cout << "Key Generation Bob:" << endl;
        cout << "JDring_time_generating_Gaussian_error_e_BOB:" << 1000.0 * JDring_time_generating_Gaussian_error_e_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_generating_Gaussian_secret_key_s_BOB:" << 1000.0 * JDring_time_generating_Gaussian_secret_key_s_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_computing_two_times_e_BOB:" << 1000.0 * JDring_time_computing_two_times_e_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_computing_mult_part_of_public_key_p_BOB:" << 1000.0 * JDring_time_computing_mult_part_of_public_key_p_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_computing_add_part_of_public_key_p_BOB:" << 1000.0 * JDring_time_computing_add_part_of_public_key_p_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

        cout << "The Protocol for Bob:" << endl;
        cout << "JDring_time_generating_Gaussian_error_e_BOB_2:" << 1000.0 * JDring_time_generating_Gaussian_error_e_BOB_2 / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_computing_two_times_e_BOB_2:" << 1000.0 * JDring_time_computing_two_times_e_BOB_2 / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_computing_mult_part_of_polynomial_k_BOB:" << 1000.0 * JDring_time_computing_mult_part_of_polynomial_k_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_computing_add_part_of_polynomial_k_BOB:" << 1000.0 * JDring_time_computing_add_part_of_polynomial_k_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_generating_uniform_bit_polynomial_b:" << 1000.0 * JDring_time_generating_uniform_bit_polynomial_b / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_computing_floor_p_quarter:" << 1000.0 * JDring_time_computing_floor_p_quarter / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_computing_S_k_BOB:" << 1000.0 * JDring_time_computing_S_k_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "JDring_time_computing_sk_BOB_with_k_BOB_and_Sig:" << 1000.0 * JDring_time_computing_sk_BOB_with_k_BOB_and_Sig / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        cout << "" << endl;

*/

//------------------------------------KEY EXCHANGE BCNS----------------------------------------------------------

            class_BCNS BCNS = class_BCNS(n_BCNS,p_BCNS,sigma_BCNS,factor_BCNS,precision_BCNS);

            cout << "" << endl;
            cout << "BCNS" << endl;
            cout << "" << endl;

            double BCNS_time_generating_Gaussian_error_e_ALICE=0;
            double BCNS_time_generating_Gaussian_secret_key_s_ALICE=0;
            double BCNS_time_computing_mult_part_of_public_key_p_ALICE=0;
            double BCNS_time_computing_add_part_of_public_key_p_ALICE=0;
            // Key generation of Bob
            double BCNS_time_generating_Gaussian_error_e_BOB=0;
            double BCNS_time_generating_Gaussian_secret_key_s_BOB=0;
            double BCNS_time_computing_mult_part_of_public_key_p_BOB=0;
            double BCNS_time_computing_add_part_of_public_key_p_BOB=0;
            // Protocol for Bob
            double BCNS_time_generating_Gaussian_error_e_BOB_2=0;
            double BCNS_time_computing_mult_part_of_polynomial_v=0;
            double BCNS_time_computing_add_part_of_polynomial_v=0;
            double BCNS_time_generating_random_ternary_polynomial_ran=0;
            double BCNS_time_computing_mult_part_of_double_v=0;
            double BCNS_time_computing_add_part_of_double_v=0;
            	//some constants
            double BCNS_time_computing_p_half=0;
            double BCNS_time_computing_p_quarter=0;
            double BCNS_time_computing_7p_quarter=0;
            double BCNS_time_computing_rounded_p_half=0;
            double BCNS_time_computing_rounded_3p_half=0;
            	//rounding functions
            double BCNS_time_crossrounding_dbl_v=0;
            double BCNS_time_rounding_dbl_v=0;
            // Protocol for Alice
            double BCNS_time_computing_two_p_BOB=0;
            double BCNS_time_computing_two_p_BOB_s_ALICE=0;
            double BCNS_time_computing_key_with_rec=0;

            for (int k=1; k<=nr_of_runs; k++) {

               	// ----------------------------Preparation functions-------------------------------------------------
              	BCNS.generating_f();
               	BCNS.generating_uniform_polynomial_a();

               	// ----------------------------Key generation of Alice-----------------------------------------------
               	clock_t start = clock();
               	BCNS.generating_Gaussian_error_e_ALICE();
               	clock_t end = clock();
               	BCNS_time_generating_Gaussian_error_e_ALICE += double(end)-double(start);

               	start = clock();
                BCNS.generating_Gaussian_secret_key_s_ALICE();
                end = clock();
                BCNS_time_generating_Gaussian_secret_key_s_ALICE += double(end)-double(start);

                start = clock();
                BCNS.computing_mult_part_of_public_key_p_ALICE();
                end = clock();
                BCNS_time_computing_mult_part_of_public_key_p_ALICE += double(end)-double(start);

                start = clock();
                BCNS.computing_add_part_of_public_key_p_ALICE();
                end = clock();
                BCNS_time_computing_add_part_of_public_key_p_ALICE += double(end)-double(start);

                	// ----------------------------Key generation of Bob-------------------------------------------------
               	start = clock();
               	BCNS.generating_Gaussian_error_e_BOB();
               	end = clock();
               	BCNS_time_generating_Gaussian_error_e_BOB += double(end)-double(start);

               	start = clock();
               	BCNS.generating_Gaussian_secret_key_s_BOB();
               	end = clock();
               	BCNS_time_generating_Gaussian_secret_key_s_BOB += double(end)-double(start);

               	start = clock();
               	BCNS.computing_mult_part_of_public_key_p_BOB();
               	end = clock();
               	BCNS_time_computing_mult_part_of_public_key_p_BOB += double(end)-double(start);

               	start = clock();
               	BCNS.computing_add_part_of_public_key_p_BOB();
               	end = clock();
               	BCNS_time_computing_add_part_of_public_key_p_BOB += double(end)-double(start);

                	// ---------------------------Protocol for Bob-------------------------------------------------------
               	start = clock();
               	BCNS.generating_Gaussian_error_e_BOB_2();
               	end = clock();
               	BCNS_time_generating_Gaussian_error_e_BOB_2 += double(end)-double(start);

               	start = clock();
               	BCNS.computing_mult_part_of_polynomial_v();
               	end = clock();
               	BCNS_time_computing_mult_part_of_polynomial_v += double(end)-double(start);

               	start = clock();
               	BCNS.computing_add_part_of_polynomial_v();
               	end = clock();
               	BCNS_time_computing_add_part_of_polynomial_v += double(end)-double(start);


               	start = clock();
                BCNS.generating_random_ternary_polynomial_ran();
               	end = clock();
               	BCNS_time_generating_random_ternary_polynomial_ran += double(end)-double(start);

               	start = clock();
               	BCNS.computing_mult_part_of_double_v();
               	end = clock();
               	BCNS_time_computing_mult_part_of_double_v += double(end)-double(start);

               	start = clock();
               	BCNS.computing_add_part_of_double_v();
               	end = clock();
               	BCNS_time_computing_add_part_of_double_v += double(end)-double(start);

               		//-------------------------------some constants
               	start = clock();
               	BCNS.computing_p_half();
               	end = clock();
               	BCNS_time_computing_p_half += double(end)-double(start);

               	start = clock();
               	BCNS.computing_p_quarter();
               	end = clock();
               	BCNS_time_computing_p_quarter += double(end)-double(start);

               	start = clock();
               	BCNS.computing_7p_quarter();
               	end = clock();
               	BCNS_time_computing_7p_quarter += double(end)-double(start);

               	start = clock();
               	BCNS.computing_rounded_p_half();
               	end = clock();
               	BCNS_time_computing_rounded_p_half += double(end)-double(start);

               	start = clock();
               	BCNS.computing_rounded_3p_half();
               	end = clock();
               	BCNS_time_computing_rounded_3p_half += double(end)-double(start);

                		//--------------------------------rounding functions

               	start = clock();
               	BCNS.crossrounding_dbl_v();
               	end = clock();
               	BCNS_time_crossrounding_dbl_v += double(end)-double(start);

               	start = clock();
               	BCNS.rounding_dbl_v();
               	end = clock();
               	BCNS_time_rounding_dbl_v += double(end)-double(start);

               	// -------------------------------------Protocol for Alice-------------------------------------------
               	start = clock();
               	BCNS.computing_two_p_BOB();
               	end = clock();
               	BCNS_time_computing_two_p_BOB += double(end)-double(start);

               	start = clock();
               	BCNS.computing_two_p_BOB_s_ALICE();
               	end = clock();
               	BCNS_time_computing_two_p_BOB_s_ALICE += double(end)-double(start);

               	start = clock();
               	BCNS.computing_key_with_rec();
               	end = clock();
               	BCNS_time_computing_key_with_rec += double(end)-double(start);

                }

            cout << "" << endl;
            cout << "Key Generation Alice:" << endl;
            cout << "BCNS_time_generating_Gaussian_error_e_ALICE:" << 1000.0 * BCNS_time_generating_Gaussian_error_e_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_generating_Gaussian_secret_key_s_ALICE:" << 1000.0 * BCNS_time_generating_Gaussian_secret_key_s_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_mult_part_of_public_key_p_ALICE:" << 1000.0 * BCNS_time_computing_mult_part_of_public_key_p_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_add_part_of_public_key_p_ALICE:" << 1000.0 * BCNS_time_computing_add_part_of_public_key_p_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "Key Generation Bob:" << endl;
            cout << "BCNS_time_generating_Gaussian_error_e_BOB:" << 1000.0 * BCNS_time_generating_Gaussian_error_e_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_generating_Gaussian_secret_key_s_BOB:" << 1000.0 * BCNS_time_generating_Gaussian_secret_key_s_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_mult_part_of_public_key_p_BOB:" << 1000.0 * BCNS_time_computing_mult_part_of_public_key_p_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_add_part_of_public_key_p_BOB:" << 1000.0 * BCNS_time_computing_add_part_of_public_key_p_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "The Protocol for Bob:" << endl;
            cout << "BCNS_time_generating_Gaussian_error_e_BOB_2:" << 1000.0 * BCNS_time_generating_Gaussian_error_e_BOB_2 / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_mult_part_of_polynomial_v:" << 1000.0 * BCNS_time_computing_mult_part_of_polynomial_v / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_add_part_of_polynomial_v:" << 1000.0 * BCNS_time_computing_add_part_of_polynomial_v / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_generating_random_ternary_polynomial_ran:" << 1000.0 * BCNS_time_generating_random_ternary_polynomial_ran / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_mult_part_of_double_v:" << 1000.0 * BCNS_time_computing_mult_part_of_double_v / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_add_part_of_double_v:" << 1000.0 * BCNS_time_computing_add_part_of_double_v / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
             //some constants
            cout << "BCNS_time_computing_p_half:" << 1000.0 * BCNS_time_computing_p_half / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_p_quarter:" << 1000.0 * BCNS_time_computing_p_quarter / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_7p_quarter:" << 1000.0 * BCNS_time_computing_7p_quarter / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_rounded_p_half:" << 1000.0 * BCNS_time_computing_rounded_p_half / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_rounded_3p_half:" << 1000.0 * BCNS_time_computing_rounded_3p_half / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
             //rounding functions
            cout << "BCNS_time_crossrounding_dbl_v:" << 1000.0 * BCNS_time_crossrounding_dbl_v / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_rounding_dbl_v:" << 1000.0 * BCNS_time_rounding_dbl_v / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "The Protocol for Alice:" << endl;
            cout << "BCNS_time_computing_two_p_BOB:" << 1000.0 * BCNS_time_computing_two_p_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_two_p_BOB_s_ALICE:" << 1000.0 * BCNS_time_computing_two_p_BOB_s_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "BCNS_time_computing_key_with_rec:" << 1000.0 * BCNS_time_computing_key_with_rec / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "" << endl;


//------------------------------------KEY EXCHANGE LPR----------------------------------------------------------

            class_LPR LPR = class_LPR(n_LPR,p_LPR,sigma_LPR,factor_LPR,precision_LPR);
/*
            cout << "" << endl;
            cout << "LPR" << endl;
            cout << "" << endl;

            double LPR_time_generating_Gaussian_error_e_A=0;
            double LPR_time_generating_Gaussian_secret_key_s_A=0;
            double LPR_time_computing_mult_part_of_public_key_p_A=0;
            double LPR_time_computing_add_part_of_public_key_p_A=0;
            double LPR_time_generating_uniform_key=0;
            double LPR_time_encoding_key=0;
            double LPR_time_generating_Gaussian_error_s_B=0;
            double LPR_time_generating_Gaussian_error_e_B=0;
            double LPR_time_generating_Gaussian_error_e_B_1=0;
            double LPR_time_computing_mult_part_of_ciphertext_p_B=0;
            double LPR_time_computing_mult_part_of_ciphertext_c=0;
            double LPR_time_computing_add_part_of_ciphertext_p_B=0;
            double LPR_time_computing_add_part_of_ciphertext_c=0;
            double LPR_time_computing_mult_part_of_decrypt=0;
            double LPR_time_computing_add_part_of_decrypt=0;
            double LPR_time_converting_decrypt=0;
            double LPR_time_decoding_decrypt=0;

            for (int k=1; k<=nr_of_runs; k++) {

            	// ----------------------------Preparation functions-------------------------------------------------
            	LPR.generating_f();
            	LPR.generating_uniform_polynomial_a();

            	// ----------------------------Key generation-----------------------------------------------
            	clock_t start = clock();
            	LPR.generating_Gaussian_error_e_A();
            	clock_t end = clock();
            	LPR_time_generating_Gaussian_error_e_A += double(end)-double(start);

            	start = clock();
            	LPR.generating_Gaussian_secret_key_s_A();
            	end = clock();
            	LPR_time_generating_Gaussian_secret_key_s_A += double(end)-double(start);

            	start = clock();
            	LPR.computing_mult_part_of_public_key_p_A();
            	end = clock();
            	LPR_time_computing_mult_part_of_public_key_p_A += double(end)-double(start);

            	start = clock();
            	LPR.computing_add_part_of_public_key_p_A();
            	end = clock();
            	LPR_time_computing_add_part_of_public_key_p_A += double(end)-double(start);

            // ----------------------------Encryption-------------------------------------------------
            	start = clock();
            	LPR.generating_uniform_key();
            	end = clock();
            	LPR_time_generating_uniform_key += double(end)-double(start);

            	start = clock();
            	LPR.encoding_key();
            	end = clock();
            	LPR_time_encoding_key += double(end)-double(start);

            	start = clock();
            	LPR.generating_Gaussian_error_s_B();
            	end = clock();
            	LPR_time_generating_Gaussian_error_s_B += double(end)-double(start);

            	start = clock();
            	LPR.generating_Gaussian_error_e_B();
            	end = clock();
            	LPR_time_generating_Gaussian_error_e_B += double(end)-double(start);

            	start = clock();
            	LPR.generating_Gaussian_error_e_B_1();
            	end = clock();
            	LPR_time_generating_Gaussian_error_e_B_1 += double(end)-double(start);

            	start = clock();
            	LPR.computing_mult_part_of_ciphertext_p_B();
            	end = clock();
            	LPR_time_computing_mult_part_of_ciphertext_p_B += double(end)-double(start);

            	start = clock();
            	LPR.computing_mult_part_of_ciphertext_c();
            	end = clock();
            	LPR_time_computing_mult_part_of_ciphertext_c += double(end)-double(start);

            	start = clock();
            	LPR.computing_add_part_of_ciphertext_p_B();
            	end = clock();
            	LPR_time_computing_add_part_of_ciphertext_p_B += double(end)-double(start);

            	start = clock();
            	LPR.computing_add_part_of_ciphertext_c();
            	end = clock();
            	LPR_time_computing_add_part_of_ciphertext_c += double(end)-double(start);

            	start = clock();
            	LPR.computing_mult_part_of_decrypt();
            	end = clock();
            	LPR_time_computing_mult_part_of_decrypt += double(end)-double(start);

            	start = clock();
            	LPR.computing_add_part_of_decrypt();
            	end = clock();
            	LPR_time_computing_add_part_of_decrypt += double(end)-double(start);

            	start = clock();
            	LPR.converting_decrypt();
            	end = clock();
            	LPR_time_converting_decrypt += double(end)-double(start);

            	start = clock();
            	LPR.decoding_decrypt();
            	end = clock();
            	LPR_time_decoding_decrypt += double(end)-double(start);

            }

            cout << "" << endl;
            cout << "Key Generation:" << endl;
            cout << "LPR_time_generating_Gaussian_error_e_A:" << 1000.0 * LPR_time_generating_Gaussian_error_e_A / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "LPR_time_generating_Gaussian_secret_key_s_A:" << 1000.0 * LPR_time_generating_Gaussian_secret_key_s_A / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "LPR_time_computing_mult_part_of_public_key_p_A:" << 1000.0 * LPR_time_computing_mult_part_of_public_key_p_A / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "LPR_time_computing_add_part_of_public_key_p_A:" << 1000.0 * LPR_time_computing_add_part_of_public_key_p_A / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

            cout << "Encryption:" << endl;
            cout << "LPR_time_generating_uniform_key:" << 1000.0 * LPR_time_generating_uniform_key / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "LPR_time_encoding_key:" << 1000.0 * LPR_time_encoding_key / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "LPR_time_generating_Gaussian_error_s_B:" << 1000.0 * LPR_time_generating_Gaussian_error_s_B / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "LPR_time_generating_Gaussian_error_e_B:" << 1000.0 * LPR_time_generating_Gaussian_error_e_B / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "LPR_time_generating_Gaussian_error_e_B_1:" << 1000.0 * LPR_time_generating_Gaussian_error_e_B_1 / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "LPR_time_computing_mult_part_of_ciphertext_p_B:" << 1000.0 * LPR_time_computing_mult_part_of_ciphertext_p_B / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "LPR_time_computing_mult_part_of_ciphertext_c:" << 1000.0 * LPR_time_computing_mult_part_of_ciphertext_c / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "LPR_time_computing_add_part_of_ciphertext_p_B:" << 1000.0 * LPR_time_computing_add_part_of_ciphertext_p_B / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "LPR_time_computing_add_part_of_ciphertext_c:" << 1000.0 * LPR_time_computing_add_part_of_ciphertext_c / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

        	cout << "Decryption:" << endl;
        	cout << "LPR_time_computing_mult_part_of_decrypt:" << 1000.0 * LPR_time_computing_mult_part_of_decrypt / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "LPR_time_computing_add_part_of_decrypt:" << 1000.0 * LPR_time_computing_add_part_of_decrypt / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "LPR_time_converting_decrypt:" << 1000.0 * LPR_time_converting_decrypt / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "LPR_time_decoding_decrypt:" << 1000.0 * LPR_time_decoding_decrypt / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "" << endl;

*/

//------------------------------------KEY EXCHANGE SS----------------------------------------------------------

        	class_SS SS = class_SS(n_SS, p_SS, alpha_SS, sigma_SS,factor_SS,precision_SS);
/*
        	cout << "" << endl;
        	cout << "SS" << endl;
        	cout << "" << endl;

        	double SS_time_generating_part_of_Gaussian_secret_key_f=0;
        	double SS_time_generating_mult_part_of_Gaussian_secret_key_f=0;
        	double SS_time_generating_add_part_of_Gaussian_secret_key_f=0;
        	double SS_time_generating_part_of_Gaussian_public_key_g=0;
        	double SS_time_generating_mult_part_of_Gaussian_public_key_h=0;
        	double SS_time_generating_public_key_h=0;
        	// Encryption
        	double SS_time_generating_Gaussian_s=0;
        	double SS_time_generating_Gaussian_e=0;
        	double SS_time_generating_first_mult_part_of_ciphertext_c=0;
        	double SS_time_generating_second_mult_part_of_ciphertext_c=0;
        	double SS_time_generating_add_part_of_ciphertext_c=0;
        	double SS_time_generating_uniform_key=0;
        	double SS_time_generating_ciphertext_c=0;
        	//Decryption
        	double SS_time_generating_decrypted_c=0;
        	double SS_time_generating_decrypted_key=0;

        	for (int k=1; k<=nr_of_runs; k++) {

        		// ----------------------------Preparation functions-------------------------------------------------
        	    SS.generating_f();

        	    // ----------------------------Key generation-----------------------------------------------

        	    clock_t start = clock();
        	    SS.generating_part_of_Gaussian_secret_key_f();
        	    clock_t end = clock();
        	    SS_time_generating_part_of_Gaussian_secret_key_f += double(end)-double(start);

        	    start = clock();
        	    SS.generating_mult_part_of_Gaussian_secret_key_f();
        	    end = clock();
        	    SS_time_generating_mult_part_of_Gaussian_secret_key_f += double(end)-double(start);

				start = clock();
        	    SS.generating_add_part_of_Gaussian_secret_key_f();
        	    end = clock();
        	    SS_time_generating_add_part_of_Gaussian_secret_key_f += double(end)-double(start);

				start = clock();
        	    SS.generating_part_of_Gaussian_public_key_g();
        	    end = clock();
        	    SS_time_generating_part_of_Gaussian_public_key_g += double(end)-double(start);

				start = clock();
        	    SS.generating_mult_part_of_Gaussian_public_key_h();
        	    end = clock();
        	    SS_time_generating_mult_part_of_Gaussian_public_key_h += double(end)-double(start);

				start = clock();
        	    SS.generating_public_key_h();
        	    end = clock();
        	    SS_time_generating_public_key_h += double(end)-double(start);


				// ----------------------------Encryption-----------------------------------------------

				start = clock();
        	    SS.generating_Gaussian_s();
        	    end = clock();
        	    SS_time_generating_Gaussian_s += double(end)-double(start);

				start = clock();
        	    SS.generating_Gaussian_e();
        	    end = clock();
        	    SS_time_generating_Gaussian_e += double(end)-double(start);

				start = clock();
        	    SS.generating_first_mult_part_of_ciphertext_c();
        	    end = clock();
        	    SS_time_generating_first_mult_part_of_ciphertext_c += double(end)-double(start);

				start = clock();
        	    SS.generating_second_mult_part_of_ciphertext_c();
        	    end = clock();
        	    SS_time_generating_second_mult_part_of_ciphertext_c += double(end)-double(start);

				start = clock();
        	    SS.generating_add_part_of_ciphertext_c();
        	    end = clock();
        	    SS_time_generating_add_part_of_ciphertext_c += double(end)-double(start);

				start = clock();
        	    SS.generating_uniform_key();
        	    end = clock();
        	    SS_time_generating_uniform_key += double(end)-double(start);

				start = clock();
        	    SS.generating_ciphertext_c();
        	    end = clock();
        	    SS_time_generating_ciphertext_c += double(end)-double(start);

				// ----------------------------Decryption-----------------------------------------------

				start = clock();
        	    SS.generating_decrypted_c();
        	    end = clock();
        	    SS_time_generating_decrypted_c += double(end)-double(start);

				start = clock();
        	    SS.generating_decrypted_key();
        	    end = clock();
        	    SS_time_generating_decrypted_key += double(end)-double(start);

        	}

        	cout << "" << endl;
        	cout << "Key Generation:" << endl;
        	cout << "SS_time_generating_part_of_Gaussian_secret_key_f:" << 1000.0 * SS_time_generating_part_of_Gaussian_secret_key_f / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "SS_time_generating_mult_part_of_Gaussian_secret_key_f:" << 1000.0 * SS_time_generating_mult_part_of_Gaussian_secret_key_f / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "SS_time_generating_add_part_of_Gaussian_secret_key_f:" << 1000.0 * SS_time_generating_add_part_of_Gaussian_secret_key_f / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "SS_time_generating_part_of_Gaussian_public_key_g:" << 1000.0 * SS_time_generating_part_of_Gaussian_public_key_g / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "SS_time_generating_mult_part_of_Gaussian_public_key_h:" << 1000.0 * SS_time_generating_mult_part_of_Gaussian_public_key_h / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "SS_time_generating_public_key_h:" << 1000.0 * SS_time_generating_public_key_h / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

        	cout << "Encryption:" << endl;
        	cout << "SS_time_generating_Gaussian_s:" << 1000.0 * SS_time_generating_Gaussian_s / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "SS_time_generating_Gaussian_e:" << 1000.0 * SS_time_generating_Gaussian_e / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "SS_time_generating_first_mult_part_of_ciphertext_c:" << 1000.0 * SS_time_generating_first_mult_part_of_ciphertext_c / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "SS_time_generating_second_mult_part_of_ciphertext_c:" << 1000.0 * SS_time_generating_second_mult_part_of_ciphertext_c / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "SS_time_generating_add_part_of_ciphertext_c:" << 1000.0 * SS_time_generating_add_part_of_ciphertext_c / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "SS_time_generating_uniform_key:" << 1000.0 * SS_time_generating_uniform_key / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "SS_time_generating_ciphertext_c:" << 1000.0 * SS_time_generating_ciphertext_c / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

        	cout << "Decryption:" << endl;
        	cout << "SS_time_generating_decrypted_c:" << 1000.0 * SS_time_generating_decrypted_c / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "SS_time_generating_decrypted_key:" << 1000.0 * SS_time_generating_decrypted_key / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
        	cout << "" << endl;
*/

//------------------------------------KEY EXCHANGE ZZD_three_pass----------------------------------------------------------

            class_ZZD_three_pass ZZD_three_pass = class_ZZD_three_pass(n_ZZD_three_pass,p_ZZD_three_pass,alpha_ZZD_three_pass,beta_ZZD_three_pass,factor_ZZD_three_pass,precision_ZZD_three_pass);
/*
            cout << "" << endl;
            cout << "ZZD_three_pass" << endl;
            cout << "" << endl;

            double ZZD_three_time_generating_Gaussian_error_e_ALICE=0;
            double ZZD_three_time_generating_Gaussian_secret_key_s_ALICE=0;
            double ZZD_three_time_computing_two_times_e_ALICE=0;
            double ZZD_three_time_computing_mult_part_of_public_key_p_ALICE=0;
            double ZZD_three_time_computing_add_part_of_public_key_p_ALICE=0;
            // Key generation of Bob
            double ZZD_three_time_generating_Gaussian_error_e_BOB=0;
            double ZZD_three_time_generating_Gaussian_secret_key_s_BOB=0;
            double ZZD_three_time_computing_two_times_e_BOB=0;
            double ZZD_three_time_computing_mult_part_of_public_key_p_BOB=0;
            double ZZD_three_time_computing_add_part_of_public_key_p_BOB=0;
            // Computation ALICE before first PASS
            double ZZD_three_time_generating_Gaussian_r_ALICE=0;
            double ZZD_three_time_generating_Gaussian_f_ALICE=0;
            double ZZD_three_time_computing_two_times_f_ALICE=0;
            double ZZD_three_time_computing_mult_part_of_x_ALICE=0;
            double ZZD_three_time_computing_add_part_of_x_ALICE=0;
            double ZZD_three_time_generating_Gaussian_d_A=0;
            double ZZD_three_time_computing_mult_part_of_hat_r_ALICE=0;
            double ZZD_three_time_computing_add_part_of_hat_r_ALICE=0;
            double ZZD_three_time_computing_mult_part_of_hat_f_ALICE=0;
            double ZZD_three_time_computing_add_part_of_hat_f_ALICE=0;
            // Computation BOB before second PASS
            double ZZD_three_time_generating_Gaussian_r_BOB=0;
            double ZZD_three_time_generating_Gaussian_f_BOB=0;
            double ZZD_three_time_computing_two_times_f_BOB=0;
            double ZZD_three_time_computing_mult_part_of_x_BOB=0;
            double ZZD_three_time_computing_add_part_of_x_BOB=0;
            double ZZD_three_time_generating_Gaussian_d_B=0;
            double ZZD_three_time_computing_mult_part_of_hat_r_BOB=0;
            double ZZD_three_time_computing_add_part_of_hat_r_BOB=0;
            double ZZD_three_time_computing_mult_part_of_hat_f_BOB=0;
            double ZZD_three_time_computing_add_part_of_hat_f_BOB=0;
            // Computation BOB after second pass
            double ZZD_three_time_generating_Gaussian_g_BOB=0;
            	// Computation of k_BOB:
            double ZZD_three_time_computing_p_ALICE_times_d_A=0;
            double ZZD_three_time_computing_k_BOB_1_plus_x_ALICE=0;
            double ZZD_three_time_computing_k_BOB_1_times_hat_r_BOB=0;
            double ZZD_three_time_computing_d_A_times_g_BOB=0;
            double ZZD_three_time_computing_2_d_A_g_BOB=0;
            double ZZD_three_time_computing_k_BOB_1_plus_k_BOB_2=0;
            	// rest of computation for BOB
            double ZZD_three_time_computing_Cha_of_k_B=0;
        	double ZZD_three_time_evalute_function_Mod2=0;

            for (int k=1; k<=nr_of_runs; k++) {

                // ----------------------------Preparation functions-------------------------------------------------
                ZZD_three_pass.generating_f();
                ZZD_three_pass.generating_uniform_polynomial_a();

                // ----------------------------Key generation of Alice-----------------------------------------------
                clock_t start = clock();
                ZZD_three_pass.generating_Gaussian_error_e_ALICE();
                clock_t end = clock();
                ZZD_three_time_generating_Gaussian_error_e_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.generating_Gaussian_secret_key_s_ALICE();
                end = clock();
                ZZD_three_time_generating_Gaussian_secret_key_s_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_two_times_e_ALICE();
                end = clock();
                ZZD_three_time_computing_two_times_e_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_mult_part_of_public_key_p_ALICE();
                end = clock();
                ZZD_three_time_computing_mult_part_of_public_key_p_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_add_part_of_public_key_p_ALICE();
                end = clock();
                ZZD_three_time_computing_add_part_of_public_key_p_ALICE += double(end)-double(start);


                // ----------------------------Key generation of Bob-------------------------------------------------
                start = clock();
                ZZD_three_pass.generating_Gaussian_error_e_BOB();
                end = clock();
                ZZD_three_time_generating_Gaussian_error_e_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.generating_Gaussian_secret_key_s_BOB();
                end = clock();
                ZZD_three_time_generating_Gaussian_secret_key_s_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_two_times_e_BOB();
                end = clock();
                ZZD_three_time_computing_two_times_e_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_mult_part_of_public_key_p_BOB();
                end = clock();
                ZZD_three_time_computing_mult_part_of_public_key_p_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_add_part_of_public_key_p_BOB();
                end = clock();
                ZZD_three_time_computing_add_part_of_public_key_p_BOB += double(end)-double(start);


                // ---------------------Computation ALICE before first PASS-------------------------------------------
                start = clock();
                ZZD_three_pass.generating_Gaussian_r_ALICE();
                end = clock();
                ZZD_three_time_generating_Gaussian_r_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.generating_Gaussian_f_ALICE();
                end = clock();
                ZZD_three_time_generating_Gaussian_f_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_two_times_f_ALICE();
                end = clock();
                ZZD_three_time_computing_two_times_f_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_mult_part_of_x_ALICE();
                end = clock();
                ZZD_three_time_computing_mult_part_of_x_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_add_part_of_x_ALICE();
                end = clock();
                ZZD_three_time_computing_add_part_of_x_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.generating_Gaussian_d_A();
                end = clock();
                ZZD_three_time_generating_Gaussian_d_A += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_mult_part_of_hat_r_ALICE();
                end = clock();
                ZZD_three_time_computing_mult_part_of_hat_r_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_add_part_of_hat_r_ALICE();
                end = clock();
                ZZD_three_time_computing_add_part_of_hat_r_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_mult_part_of_hat_f_ALICE();
                end = clock();
                ZZD_three_time_computing_mult_part_of_hat_f_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_add_part_of_hat_f_ALICE();
                end = clock();
                ZZD_three_time_computing_add_part_of_hat_f_ALICE += double(end)-double(start);

                // ---------------------Computation BOB before second PASS-------------------------------------------
                start = clock();
                ZZD_three_pass.generating_Gaussian_r_BOB();
                end = clock();
                ZZD_three_time_generating_Gaussian_r_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.generating_Gaussian_f_BOB();
                end = clock();
                ZZD_three_time_generating_Gaussian_f_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_two_times_f_BOB();
                end = clock();
                ZZD_three_time_computing_two_times_f_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_mult_part_of_x_BOB();
                end = clock();
                ZZD_three_time_computing_mult_part_of_x_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_add_part_of_x_BOB();
                end = clock();
                ZZD_three_time_computing_add_part_of_x_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.generating_Gaussian_d_B();
                end = clock();
                ZZD_three_time_generating_Gaussian_d_B += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_mult_part_of_hat_r_BOB();
                end = clock();
                ZZD_three_time_computing_mult_part_of_hat_r_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_add_part_of_hat_r_BOB();
                end = clock();
                ZZD_three_time_computing_add_part_of_hat_r_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_mult_part_of_hat_f_BOB();
                end = clock();
                ZZD_three_time_computing_mult_part_of_hat_f_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_add_part_of_hat_f_BOB();
                end = clock();
                ZZD_three_time_computing_add_part_of_hat_f_BOB += double(end)-double(start);

                // ---------------------Computation BOB after second PASS-------------------------------------------
                start = clock();
                ZZD_three_pass.generating_Gaussian_g_BOB();
                end = clock();
                ZZD_three_time_generating_Gaussian_g_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_p_ALICE_times_d_A();
                end = clock();
                ZZD_three_time_computing_p_ALICE_times_d_A += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_k_BOB_1_plus_x_ALICE();
                end = clock();
                ZZD_three_time_computing_k_BOB_1_plus_x_ALICE += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_k_BOB_1_times_hat_r_BOB();
                end = clock();
                ZZD_three_time_computing_k_BOB_1_times_hat_r_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_d_A_times_g_BOB();
                end = clock();
                ZZD_three_time_computing_d_A_times_g_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_2_d_A_g_BOB();
                end = clock();
                ZZD_three_time_computing_2_d_A_g_BOB += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_k_BOB_1_plus_k_BOB_2();
                end = clock();
                ZZD_three_time_computing_k_BOB_1_plus_k_BOB_2 += double(end)-double(start);

                start = clock();
                ZZD_three_pass.computing_Cha_of_k_B();
                end = clock();
                ZZD_three_time_computing_Cha_of_k_B += double(end)-double(start);

                start = clock();
                ZZD_three_pass.evalute_function_Mod2();
                end = clock();
                ZZD_three_time_evalute_function_Mod2 += double(end)-double(start);

            }

            cout << "" << endl;
            cout << "Key Generation Alice:" << endl;
            cout << "ZZD_three_time_generating_Gaussian_error_e_ALICE:" << 1000.0 * ZZD_three_time_generating_Gaussian_error_e_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_generating_Gaussian_secret_key_s_ALICE:" << 1000.0 * ZZD_three_time_generating_Gaussian_secret_key_s_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_two_times_e_ALICE:" << 1000.0 * ZZD_three_time_computing_two_times_e_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_mult_part_of_public_key_p_ALICE:" << 1000.0 * ZZD_three_time_computing_mult_part_of_public_key_p_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_add_part_of_public_key_p_ALICE:" << 1000.0 * ZZD_three_time_computing_add_part_of_public_key_p_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

            cout << "Key Generation Bob:" << endl;
            cout << "ZZD_three_time_generating_Gaussian_error_e_BOB:" << 1000.0 * ZZD_three_time_generating_Gaussian_error_e_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_generating_Gaussian_secret_key_s_BOB:" << 1000.0 * ZZD_three_time_generating_Gaussian_secret_key_s_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_two_times_e_BOB:" << 1000.0 * ZZD_three_time_computing_two_times_e_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_mult_part_of_public_key_p_BOB:" << 1000.0 * ZZD_three_time_computing_mult_part_of_public_key_p_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_add_part_of_public_key_p_BOB:" << 1000.0 * ZZD_three_time_computing_add_part_of_public_key_p_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

            cout << "Computation ALICE before first PASS:" << endl;
            cout << "ZZD_three_time_generating_Gaussian_r_ALICE:" << 1000.0 * ZZD_three_time_generating_Gaussian_r_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_generating_Gaussian_f_ALICE:" << 1000.0 * ZZD_three_time_generating_Gaussian_f_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_two_times_f_ALICE:" << 1000.0 * ZZD_three_time_computing_two_times_f_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_mult_part_of_x_ALICE:" << 1000.0 * ZZD_three_time_computing_mult_part_of_x_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_add_part_of_x_ALICE:" << 1000.0 * ZZD_three_time_computing_add_part_of_x_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_generating_Gaussian_d_A:" << 1000.0 * ZZD_three_time_generating_Gaussian_d_A / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_mult_part_of_hat_r_ALICE:" << 1000.0 * ZZD_three_time_computing_mult_part_of_hat_r_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_add_part_of_hat_r_ALICE:" << 1000.0 * ZZD_three_time_computing_add_part_of_hat_r_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_mult_part_of_hat_f_ALICE:" << 1000.0 * ZZD_three_time_computing_mult_part_of_hat_f_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_add_part_of_hat_f_ALICE:" << 1000.0 * ZZD_three_time_computing_add_part_of_hat_f_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

            cout << "Computation BOB before second PASS:" << endl;
            cout << "ZZD_three_time_generating_Gaussian_r_BOB:" << 1000.0 * ZZD_three_time_generating_Gaussian_r_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_generating_Gaussian_f_BOB:" << 1000.0 * ZZD_three_time_generating_Gaussian_f_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_two_times_f_BOB:" << 1000.0 * ZZD_three_time_computing_two_times_f_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_mult_part_of_x_BOB:" << 1000.0 * ZZD_three_time_computing_mult_part_of_x_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_add_part_of_x_BOB:" << 1000.0 * ZZD_three_time_computing_add_part_of_x_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_generating_Gaussian_d_B:" << 1000.0 * ZZD_three_time_generating_Gaussian_d_B / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_mult_part_of_hat_r_BOB:" << 1000.0 * ZZD_three_time_computing_mult_part_of_hat_r_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_add_part_of_hat_r_BOB:" << 1000.0 * ZZD_three_time_computing_add_part_of_hat_r_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_mult_part_of_hat_f_BOB:" << 1000.0 * ZZD_three_time_computing_mult_part_of_hat_f_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_add_part_of_hat_f_BOB:" << 1000.0 * ZZD_three_time_computing_add_part_of_hat_f_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

            cout << "Computation BOB after second PASS:" << endl;
            cout << "ZZD_three_time_generating_Gaussian_g_BOB:" << 1000.0 * ZZD_three_time_generating_Gaussian_g_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_p_ALICE_times_d_A:" << 1000.0 * ZZD_three_time_computing_p_ALICE_times_d_A / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_k_BOB_1_plus_x_ALICE:" << 1000.0 * ZZD_three_time_computing_k_BOB_1_plus_x_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_k_BOB_1_times_hat_r_BOB:" << 1000.0 * ZZD_three_time_computing_k_BOB_1_times_hat_r_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_d_A_times_g_BOB:" << 1000.0 * ZZD_three_time_computing_d_A_times_g_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_2_d_A_g_BOB:" << 1000.0 * ZZD_three_time_computing_2_d_A_g_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_k_BOB_1_plus_k_BOB_2:" << 1000.0 * ZZD_three_time_computing_k_BOB_1_plus_k_BOB_2 / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_computing_Cha_of_k_B:" << 1000.0 * ZZD_three_time_computing_Cha_of_k_B / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "ZZD_three_time_evalute_function_Mod2:" << 1000.0 * ZZD_three_time_evalute_function_Mod2 / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "" << endl;
*/

//------------------------------------KEY EXCHANGE ZZD_one_pass----------------------------------------------------------

            class_ZZD_one_pass ZZD_one_pass = class_ZZD_one_pass(n_ZZD_one_pass,p_ZZD_one_pass,alpha_ZZD_one_pass,beta_ZZD_one_pass,factor_ZZD_one_pass,precision_ZZD_one_pass);
/*
            cout << "" << endl;
            cout << "ZZD_one_pass" << endl;
            cout << "" << endl;

            double time_generating_Gaussian_error_e_ALICE=0;
            double time_generating_Gaussian_secret_key_s_ALICE=0;
            double time_computing_two_times_e_ALICE=0;
            double time_computing_mult_part_of_public_key_p_ALICE=0;
            double time_computing_add_part_of_public_key_p_ALICE=0;
            // Key generation of Bob
            double time_generating_Gaussian_error_e_BOB=0;
            double time_generating_Gaussian_secret_key_s_BOB=0;
            double time_computing_two_times_e_BOB=0;
            double time_computing_mult_part_of_public_key_p_BOB=0;
            double time_computing_add_part_of_public_key_p_BOB=0;
            // Computation ALICE before first PASS
            double time_generating_Gaussian_r_ALICE=0;
            double time_generating_Gaussian_f_ALICE=0;
            double time_computing_two_times_f_ALICE=0;
            double time_computing_mult_part_of_x_ALICE=0;
            double time_computing_add_part_of_x_ALICE=0;
            double time_generating_Gaussian_d_A=0;
            double time_computing_mult_part_of_hat_r_ALICE=0;
            double time_computing_add_part_of_hat_r_ALICE=0;
            double time_computing_mult_part_of_hat_f_ALICE=0;
            double time_computing_add_part_of_hat_f_ALICE=0;
            // Computation k_ALICE
            double time_generating_Gaussian_g_ALICE=0;
            double time_computing_two_times_g_ALICE=0;
            double time_computing_mult_part_of_k_ALICE=0;
            double time_computing_add_part_of_k_ALICE=0;
            	// Computation of Cha and sig_ALICE
            double time_computing_Cha_of_k_ALICE=0;
            double time_evalute_function_Mod2_ALICE=0;
            // Computation BOB after receiving the only message
            double time_generating_Gaussian_g_BOB=0;
            	// Computation of k_BOB
            double time_computing_p_ALICE_times_d_A=0;
            double time_computing_k_BOB_1_plus_x_ALICE=0;
            double time_computing_k_BOB_1_times_s_BOB=0;
            double time_computing_d_A_times_g_BOB=0;
            double time_computing_2_d_A_g_BOB=0;
            double time_computing_k_BOB_1_plus_k_BOB_2=0;
            	// Computation of sig_BOB
            double time_evalute_function_Mod2_BOB=0;


            for (int k=1; k<=nr_of_runs; k++) {

                // ----------------------------Preparation functions-------------------------------------------------
                ZZD_one_pass.generating_f();
                ZZD_one_pass.generating_uniform_polynomial_a();

                // ----------------------------Key generation of Alice-----------------------------------------------
                clock_t start = clock();
                ZZD_one_pass.generating_Gaussian_error_e_ALICE();
                clock_t end = clock();
                time_generating_Gaussian_error_e_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.generating_Gaussian_secret_key_s_ALICE();
                end = clock();
                time_generating_Gaussian_secret_key_s_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_two_times_e_ALICE();
                end = clock();
                time_computing_two_times_e_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_mult_part_of_public_key_p_ALICE();
                end = clock();
                time_computing_mult_part_of_public_key_p_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_add_part_of_public_key_p_ALICE();
                end = clock();
                time_computing_add_part_of_public_key_p_ALICE += double(end)-double(start);


                // ----------------------------Key generation of Bob-------------------------------------------------
                start = clock();
                ZZD_one_pass.generating_Gaussian_error_e_BOB();
                end = clock();
                time_generating_Gaussian_error_e_BOB += double(end)-double(start);

                start = clock();
                ZZD_one_pass.generating_Gaussian_secret_key_s_BOB();
                end = clock();
                time_generating_Gaussian_secret_key_s_BOB += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_two_times_e_BOB();
                end = clock();
                time_computing_two_times_e_BOB += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_mult_part_of_public_key_p_BOB();
                end = clock();
                time_computing_mult_part_of_public_key_p_BOB += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_add_part_of_public_key_p_BOB();
                end = clock();
                time_computing_add_part_of_public_key_p_BOB += double(end)-double(start);


                // ---------------------Computation ALICE before first PASS-------------------------------------------
                start = clock();
                ZZD_one_pass.generating_Gaussian_r_ALICE();
                end = clock();
                time_generating_Gaussian_r_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.generating_Gaussian_f_ALICE();
                end = clock();
                time_generating_Gaussian_f_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_two_times_f_ALICE();
                end = clock();
                time_computing_two_times_f_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_mult_part_of_x_ALICE();
                end = clock();
                time_computing_mult_part_of_x_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_add_part_of_x_ALICE();
                end = clock();
                time_computing_add_part_of_x_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.generating_Gaussian_d_A();
                end = clock();
                time_generating_Gaussian_d_A += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_mult_part_of_hat_r_ALICE();
                end = clock();
                time_computing_mult_part_of_hat_r_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_add_part_of_hat_r_ALICE();
                end = clock();
                time_computing_add_part_of_hat_r_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_mult_part_of_hat_f_ALICE();
                end = clock();
                time_computing_mult_part_of_hat_f_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_add_part_of_hat_f_ALICE();
                end = clock();
                time_computing_add_part_of_hat_f_ALICE += double(end)-double(start);

                // ----------k_ALICE------
                start = clock();
                ZZD_one_pass.generating_Gaussian_g_ALICE();
                end = clock();
                time_generating_Gaussian_g_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_two_times_g_ALICE();
                end = clock();
                time_computing_two_times_g_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_mult_part_of_k_ALICE();
                end = clock();
                time_computing_mult_part_of_k_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_add_part_of_k_ALICE();
                end = clock();
                time_computing_add_part_of_k_ALICE += double(end)-double(start);

                //------  Computation of Cha and sig_ALICE
                start = clock();
                ZZD_one_pass.computing_Cha_of_k_ALICE();
                end = clock();
                time_computing_Cha_of_k_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.evalute_function_Mod2_ALICE();
                end = clock();
                time_evalute_function_Mod2_ALICE += double(end)-double(start);


                // -----------------------Computation BOB after receiving the only message
                start = clock();
                ZZD_one_pass.generating_Gaussian_g_BOB();
                end = clock();
                time_generating_Gaussian_g_BOB += double(end)-double(start);

                // --------- Computation of k_BOB
                start = clock();
                ZZD_one_pass.computing_p_ALICE_times_d_A();
                end = clock();
                time_computing_p_ALICE_times_d_A += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_k_BOB_1_plus_x_ALICE();
                end = clock();
                time_computing_k_BOB_1_plus_x_ALICE += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_k_BOB_1_times_s_BOB();
                end = clock();
                time_computing_k_BOB_1_times_s_BOB += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_d_A_times_g_BOB();
                end = clock();
                time_computing_d_A_times_g_BOB += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_2_d_A_g_BOB();
                end = clock();
                time_computing_2_d_A_g_BOB += double(end)-double(start);

                start = clock();
                ZZD_one_pass.computing_k_BOB_1_plus_k_BOB_2();
                end = clock();
                time_computing_k_BOB_1_plus_k_BOB_2 += double(end)-double(start);

                // ------------- Computation of sig_BOB
                start = clock();
                ZZD_one_pass.evalute_function_Mod2_BOB();
                end = clock();
                time_evalute_function_Mod2_BOB += double(end)-double(start);

            }

            cout << "" << endl;
            cout << "Key Generation Alice:" << endl;
            cout << "time_generating_Gaussian_error_e_ALICE:" << 1000.0 * time_generating_Gaussian_error_e_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_generating_Gaussian_secret_key_s_ALICE:" << 1000.0 * time_generating_Gaussian_secret_key_s_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_two_times_e_ALICE:" << 1000.0 * time_computing_two_times_e_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_mult_part_of_public_key_p_ALICE:" << 1000.0 * time_computing_mult_part_of_public_key_p_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_add_part_of_public_key_p_ALICE:" << 1000.0 * time_computing_add_part_of_public_key_p_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

            cout << "Key Generation Bob:" << endl;
            cout << "time_generating_Gaussian_error_e_BOB:" << 1000.0 * time_generating_Gaussian_error_e_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_generating_Gaussian_secret_key_s_BOB:" << 1000.0 * time_generating_Gaussian_secret_key_s_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_two_times_e_BOB:" << 1000.0 * time_computing_two_times_e_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_mult_part_of_public_key_p_BOB:" << 1000.0 * time_computing_mult_part_of_public_key_p_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_add_part_of_public_key_p_BOB:" << 1000.0 * time_computing_add_part_of_public_key_p_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

            cout << "Computation ALICE before first PASS:" << endl;
            cout << "time_generating_Gaussian_r_ALICE:" << 1000.0 * time_generating_Gaussian_r_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_generating_Gaussian_f_ALICE:" << 1000.0 * time_generating_Gaussian_f_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_two_times_f_ALICE:" << 1000.0 * time_computing_two_times_f_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_mult_part_of_x_ALICE:" << 1000.0 * time_computing_mult_part_of_x_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_add_part_of_x_ALICE:" << 1000.0 * time_computing_add_part_of_x_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_generating_Gaussian_d_A:" << 1000.0 * time_generating_Gaussian_d_A / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_mult_part_of_hat_r_ALICE:" << 1000.0 * time_computing_mult_part_of_hat_r_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_add_part_of_hat_r_ALICE:" << 1000.0 * time_computing_add_part_of_hat_r_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_mult_part_of_hat_f_ALICE:" << 1000.0 * time_computing_mult_part_of_hat_f_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_add_part_of_hat_f_ALICE:" << 1000.0 * time_computing_add_part_of_hat_f_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_generating_Gaussian_g_ALICE:" << 1000.0 * time_generating_Gaussian_g_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_two_times_g_ALICE:" << 1000.0 * time_computing_two_times_g_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_mult_part_of_k_ALICE:" << 1000.0 * time_computing_mult_part_of_k_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_add_part_of_k_ALICE:" << 1000.0 * time_computing_add_part_of_k_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_Cha_of_k_ALICE:" << 1000.0 * time_computing_Cha_of_k_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_evalute_function_Mod2_ALICE:" << 1000.0 * time_evalute_function_Mod2_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;

            cout << "Computation BOB after receiving the only message:" << endl;
            cout << "time_generating_Gaussian_g_BOB:" << 1000.0 * time_generating_Gaussian_g_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_p_ALICE_times_d_A:" << 1000.0 * time_computing_p_ALICE_times_d_A / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_k_BOB_1_plus_x_ALICE:" << 1000.0 * time_computing_k_BOB_1_plus_x_ALICE / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_k_BOB_1_times_s_BOB:" << 1000.0 * time_computing_k_BOB_1_times_s_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "computing_d_A_times_g_BOB:" << 1000.0 * time_computing_d_A_times_g_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_2_d_A_g_BOB:" << 1000.0 * time_computing_2_d_A_g_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_computing_k_BOB_1_plus_k_BOB_2:" << 1000.0 * time_computing_k_BOB_1_plus_k_BOB_2 / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "time_evalute_function_Mod2_BOB:" << 1000.0 * time_evalute_function_Mod2_BOB / CLOCKS_PER_SEC / nr_of_runs << "ms" << endl;
            cout << "" << endl;
*/


    }else if(scheme_yes_no==0){
        
        
        
        //------------------------------------------------------------------------------------------------------------------------------------------------//
        //------------------------------------------------------------------------------------------------------------------------------------------------//
        
        
        cout << "Parameterwahl für Gauss Sampler:" << endl;
        
        double sigma = 36.11;
        long precision = 128;
        ZZ factor = ZZ(13);
        int m = 8192;
        double n = 2048;
        
        int choose_para;
        cout << "choose parameters (1) or not (0)?" << endl;
        cin >> choose_para;
        
        if (choose_para==1){
            
            cout << "sigma?" << endl;
            cin >> sigma;
            cout << "precision?" << endl;
            cin >> precision;
            cout << "factor for range" << endl;
            cin >> factor;
            cout << "m" << endl;
            cin >> m;
            cout << "n" << endl;
            cin >> n;
        }
        
        // Anzahl der Durchläufe
        double nr_of_runs_gauss;
        
        cout << "Wie viele Durchläufe?" << endl;
        cin >> nr_of_runs_gauss;
        
        // ------------------------------------------------------------
        clock_t start_all = clock();
        // ------------------------------------------------------------
        
        // 0: Selected
        
        /*class_gauss_selected gauss_selected = class_gauss_selected(sigma,to_long(factor),precision);
         
         double time_selected = 0;
         
         for (int k=1; k<=nr_of_runs_gauss; k++) {
         clock_t start = clock();
         ZZX res_rej = gauss_selected.sample_poly(0);
         clock_t end = clock();
         time_selected += double(end)-double(start);
         }
         
         cout << "time for selected: " << 1000.0 * time_selected / CLOCKS_PER_SEC / nr_of_runs_gauss << "ms" << endl;
         cout << "time for selected * n: " << 1000.0 * time_selected / CLOCKS_PER_SEC / nr_of_runs_gauss * n << "ms" << endl;
         
         //cout << "time for selected: " << double(time_selected) / double(nr_of_runs_gauss) << "ms" << endl;
         //cout << "time for selected * n: " << double(time_selected) / double(nr_of_runs_gauss) * n << "ms" << endl;
         
         double time_selected_n = 0;
         
         for (int k=1; k<=nr_of_runs_gauss; k++) {
         clock_t start = clock();
         ZZX res_rej = gauss_selected.sample_poly(n-1);
         clock_t end = clock();
         time_selected_n += double(end)-double(start);
         }
         
         cout << "time for selected poly_n: " << 1000.0 * time_selected_n / CLOCKS_PER_SEC / nr_of_runs_gauss << "ms" << endl;*/
        
        // 1: Rejection
        
        class_gauss_rejection gauss_rejection = class_gauss_rejection(sigma,to_ZZ(factor*sigma),precision);
        
        
        double time_rejection_n = 0;
        
        for (int k=1; k<=nr_of_runs_gauss; k++) {
            clock_t start = clock();
            ZZX res_rej = gauss_rejection.sample_poly(n-1);
            clock_t end = clock();
            time_rejection_n += double(end)-double(start);
        }
        
        cout << "time for rejection poly_n-1: " << 1000.0 * time_rejection_n / CLOCKS_PER_SEC / nr_of_runs_gauss << "ms" << endl;
        
        // 2: Inverse CDF
        
        class_gauss_inverse_cdf gauss_inverse_cdf = class_gauss_inverse_cdf(to_long(sigma),to_long(factor)*sigma,precision);
        
        double time_inverse_cdf_n = 0;
        
        for (int k=1; k<=nr_of_runs_gauss; k++) {
            
            clock_t start = clock();
            ZZX res_cdf = gauss_inverse_cdf.sample_poly(n-1);
            clock_t end = clock();
            time_inverse_cdf_n += double(end)-double(start);
        }
        
        cout << "time for inverse cdf poly_n-1: " << 1000.0 * time_inverse_cdf_n / CLOCKS_PER_SEC / nr_of_runs_gauss << "ms" << endl;
        
        
        // 3: Knuth-Yao
        
        class_gauss_knuth_yao gauss_knuth_yao = class_gauss_knuth_yao(sigma, to_long(factor)*sigma, precision);
        
        double time_knuth_yao_n = 0;
        
        for (int k=1; k<=nr_of_runs_gauss; k++) {
            clock_t start = clock();
            ZZX res_cdf = gauss_knuth_yao.sample_poly(n-1);
            clock_t end = clock();
            time_knuth_yao_n += double(end)-double(start);
        }
        
        cout << "time for knuth-yao poly_n-1: " << 1000.0 * time_knuth_yao_n / CLOCKS_PER_SEC / nr_of_runs_gauss << "ms" << endl;
        
        // 4: Ziggurat
        
        class_gauss_ziggurat gauss_ziggurat = class_gauss_ziggurat(m, to_int(precision), to_int(sigma));
        
        double time_ziggurat_n = 0;
        
        for (int k=1; k<=nr_of_runs_gauss; k++) {
            
            clock_t start = clock();
            ZZX res_cdf = gauss_ziggurat.sample_poly(n-1);
            clock_t end = clock();
            time_ziggurat_n += double(end)-double(start);
        }
        
        cout << "time for ziggurat poly_n-1: " << 1000.0 * time_ziggurat_n / CLOCKS_PER_SEC / nr_of_runs_gauss << "ms" << endl;
        
        // ------------------------------------------------------------
        clock_t end_all = clock();
        cout << "CPU time used: " << 1000.0 * double(end_all-start_all) / CLOCKS_PER_SEC << " ms" << endl;
        // ------------------------------------------------------------
        
        getchar();
        
    }
    
    // prevent from terminating
    cout << "END" << endl;
    int pause = 0;
    cin >> pause;
}
