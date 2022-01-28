/*
 * Copyright (c) 2018 XLAB d.o.o.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <gmp.h>
#include "cifer/test.h"

#include "cifer/innerprod/fullysec/damgard.h"
#include "cifer/sample/uniform.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



void readIntFile(cfe_mat* output, char* filename, int row_num, int col_num) {
    FILE* fp = fopen(filename, "r");
    cfe_mat data;
    cfe_mat_init(output, row_num, col_num);
    if (!fp)
        printf("Can't open file\n");
  
    else {
        // Here we have taken size of
        // array 1024 you can modify it
        char buffer[1000000];
  
        int row = 0;
        int column = 0;
        mpz_t val;
        while (fgets(buffer, 1000000, fp)) {
            column = 0;
            char* value = strtok(buffer, " ");

            // printf("row = %d\n", row);
            while (value) {  
                mpz_init(val);
                mpz_set_str(val, value, 10);
                value = strtok(NULL, " ");
                cfe_mat_set(output, val, row, column);
                column++;
                mpz_clear(val);
            }
            row++;  
        }
  
        // Close the file
        fclose(fp);
    }
}

void readModel(cfe_vec* output, char* filename, int feature_num, int precision) {
    FILE* fp = fopen(filename, "r");
    cfe_mat data;
    cfe_vec_init(output, feature_num);
    if (!fp)
        printf("Can't open file\n");
  
    else {
        // Here we have taken size of
        // array 1024 you can modify it
        char buffer[1000000];
  
        int row = 0;
        int column = 0;
        mpz_t val;
        fgets(buffer, 1000000, fp);
        while (fgets(buffer, 1000000, fp)) {
            column = 0;
            char* value = strtok(buffer, ",");

            // printf("row = %d\n", row);
            // printf("row, column = %d, %d\n", row, column);
            while (value) {  
                if (column == 1) {
                    mpz_init(val);
                    float val_f = atof(value);
                    int val_i = (int) (val_f * pow(2, precision));
                    mpz_set_si(val, val_i);
                    // gmp_printf("%s, %.5f, %d, %Zd\n", value, val_f, val_i, val);
                    cfe_vec_set(output, val, row);
                    mpz_clear(val);
                }
                value = strtok(NULL, ",");
                column++;                
            }
            row++;  
        }
  
        // Close the file
        fclose(fp);
    }
}

void saveToCsv(cfe_vec* vec, char* filename, int dim) {
    FILE* fp = fopen(filename, "w+");
    mpz_t val;
    mpz_init(val);
    for (int i = 0; i < dim; i++) {
        cfe_vec_get(val, &vec, i);
        fprintf(fp, "%Zd\n", val);
    }
}

MunitResult test_damgard_end_to_end(const MunitParameter *params, void *data) {
    size_t start, end;
    float time_init, time_enc, time_keygen, time_dec;
    
    size_t l = 5001;
    size_t bdd = 0;
    printf("featrue = %d, weight bound = %d\n", l, bdd);
    mpz_t bound, bound_neg, key1, key2, xy_check, xy;
    mpz_inits(bound, bound_neg, key1, key2, xy_check, xy, NULL);
    mpz_set_ui(bound, 2);
    mpz_pow_ui(bound, bound, 10);
    mpz_neg(bound_neg, bound);
    cfe_damgard s, encryptor, decryptor;
    cfe_error err;

    
    size_t modulus_len;
    const char *precomp = munit_parameters_get(params, "parameters");

    char* filename_data = "../lr_data/dataset/gisette_valid.data";
    char* filename_W = "../lr_data/model/W_df18.csv";
    cfe_mat valid_data;
    cfe_vec W;
    readIntFile(&valid_data, filename_data, 1000, l);
    // cfe_mat_print(&data);
    readModel(&W, filename_W, l, bdd);
    // cfe_vec_print(&W);

    cfe_vec mpk, x, y;
    cfe_vec_inits(l, &x, &y, NULL);

    size_t row = 0;
    cfe_mat_get_row(&x, &valid_data, row);
    cfe_vec_copy(&y, &W);

    // cfe_uniform_sample_range_vec(&x, bound_neg, bound);
    // cfe_uniform_sample_range_vec(&y, bound_neg, bound);
    cfe_vec_dot(xy_check, &x, &y);
    cfe_vec_print(&x);
    cfe_vec_print(&y);
    gmp_printf("xy, %Zd\n", xy_check);

    bool origin = 0;
    if (strcmp(precomp, "precomputed") == 0 || origin == 1) {
        if (strcmp(precomp, "precomputed") == 0) {
            // modulus_len defines the security of the scheme, the higher the better
            modulus_len = 1024;
            err = cfe_damgard_precomp_init(&s, l, modulus_len, bound);
        } else if (strcmp(precomp, "random") == 0) {
            modulus_len = 512;
            err = cfe_damgard_init(&s, l, modulus_len, bound);
        } else {
            err = CFE_ERR_INIT;
        }
        munit_assert(err == 0);
        
        cfe_damgard_sec_key msk;
        cfe_vec ciphertext;

        start = clock();
        cfe_damgard_sec_key_init(&msk, &s);
        cfe_damgard_pub_key_init(&mpk, &s);
        cfe_damgard_generate_master_keys(&msk, &mpk, &s);
        end = clock();
        time_init = (float) (end - start) / CLOCKS_PER_SEC;
        printf("init time = %.7f\n", time_init);

        start = clock();
        cfe_damgard_fe_key key;
        cfe_damgard_fe_key_init(&key);
        err = cfe_damgard_derive_fe_key(&key, &s, &msk, &y);
        munit_assert(err == 0);
        end = clock();
        time_keygen = (float) (end - start) / CLOCKS_PER_SEC;
        printf("keygen time = %.7f\n", time_keygen);

        start = clock();
        cfe_damgard_copy(&encryptor, &s);
        cfe_damgard_ciphertext_init(&ciphertext, &encryptor);
        err = cfe_damgard_encrypt(&ciphertext, &encryptor, &x, &mpk);
        munit_assert(err == 0);
        end = clock();
        time_enc = (float) (end - start) / CLOCKS_PER_SEC;
        printf("enc time = %.7f\n", time_enc);


        start = clock();
        cfe_damgard_copy(&decryptor, &s);
        err = cfe_damgard_decrypt(xy, &decryptor, &ciphertext, &key, &y);
        munit_assert(err == 0);
        end = clock();
        time_dec = (float) (end - start) / CLOCKS_PER_SEC;
        printf("dec time = %.7f\n", time_dec);

        munit_assert(mpz_cmp(xy, xy_check) == 0);


        mpz_clears(bound, bound_neg, key1, key2, xy_check, xy, NULL);
        cfe_vec_frees(&x, &y, &mpk, &ciphertext, NULL);

        cfe_damgard_sec_key_free(&msk);
        cfe_damgard_fe_key_free(&key);
        cfe_damgard_free(&s);
        cfe_damgard_free(&encryptor);
        cfe_damgard_free(&decryptor);

    } else {
        float total_init = 0;
        float total_enc = 0;
        float  total_keygen = 0;
        float total_dec = 0;

        long repeat = 1;
        for (int i = 0; i < repeat; i++)
        {

        printf("run ours...\n");
        modulus_len = 3072;

        // s.g = random
        // s.h = g^r
        // err = cfe_damgard_init_modified(&s, l, modulus_len, bound);
        err = cfe_damgard_precomp_init(&s, l, modulus_len, bound);
        end = clock();
        cfe_damgard_sec_key msk;

        start = clock();
        cfe_damgard_sec_key_init_modified(&msk, &s);
        cfe_damgard_generate_master_keys_modified(&msk, &s);
        end = clock();
        time_init = (float) (end - start) / CLOCKS_PER_SEC;
        printf("init time = %.7f\n", time_init);

        // gmp_printf("p = %Zd\n", s.p);
        // gmp_printf("g = %Zd\n", s.g);
        // gmp_printf("l = %d\n", s.l);

        // printf("s\n");
        // cfe_vec_print(&msk.s);
        // printf("\nt\n");
        // cfe_vec_print(&msk.t);
        // printf("\nc\n");
        // cfe_vec_print(&msk.c);
        // printf("\nc_inv\n");
        // cfe_vec_print(&msk.c_inv);
        // printf("\nk\n");
        // cfe_vec_print(&msk.k);
        // gmp_printf("\na = %Zd\n", msk.a);

        start = clock();
        cfe_vec fe_key;
        cfe_damgard_fe_key_init_modified(&fe_key, &s);
        err = cfe_damgard_derive_fe_key_modified(&fe_key, &s, &msk, &y);
        munit_assert(err == 0);
        end = clock();
        time_keygen = (float) (end - start) / CLOCKS_PER_SEC;
        printf("keygen time = %.7f\n", time_keygen);

        cfe_damgard_ciphertext_fh ciphertext;
        start = clock();
        cfe_damgard_copy(&encryptor, &s);
        cfe_damgard_ciphertext_init_modified(&ciphertext, &encryptor);
        err = cfe_damgard_encrypt_modified(&ciphertext, &encryptor, &x, &msk);
        munit_assert(err == 0);
        end = clock();
        time_enc = (float) (end - start) / CLOCKS_PER_SEC;
        printf("enc time = %.7f\n", time_enc);


        start = clock();
        cfe_damgard_copy(&decryptor, &s);
        err = cfe_damgard_decrypt_modified(xy, &decryptor, &ciphertext, &fe_key, &y);
        printf("decerr = %d\n", err);
        munit_assert(err == 0);
        end = clock();
        time_dec = (float) (end - start) / CLOCKS_PER_SEC;
        printf("dec time = %.7f\n", time_dec);
        printf("i = %d\n", i);

        gmp_printf("real: %Zd, ours: %Zd\n", xy_check, xy);
        // gmp_printf("comp : %Zd\n", mpz_cmp(xy, xy_check));
        // munit_assert(mpz_cmp(xy, xy_check) == 0);
        

        // mpz_clears(bound, bound_neg, key1, key2, xy_check, xy, NULL);
        // cfe_vec_frees(&x, &y, &ciphertext, NULL);

        // cfe_damgard_sec_key_free(&msk);
        // cfe_vec_free(&fe_key);
        // cfe_damgard_free(&s);
        // cfe_damgard_free(&encryptor);
        // cfe_damgard_free(&decryptor);
        total_init += time_init;
        total_enc += time_enc;
        total_keygen += time_keygen;
        total_dec += time_dec;
        }
        printf("avg init = %.7f\n", total_init / repeat);
        printf("avg keygen = %.7f\n", total_keygen/ repeat);
        printf("avg enc = %.7f\n", total_enc/ repeat);
        printf("avg dec = %.7f\n", total_dec/ repeat);
    }

    return MUNIT_OK;
}

// char *damgard_param[] = {
//         (char *) "precomputed", (char *) "random", NULL
// };
char *damgard_param[] = {
        (char *) "random", NULL
};

MunitParameterEnum damgard_params[] = {
        {(char *) "parameters", damgard_param},
        {NULL,                  NULL},
};

MunitTest simple_ip_damgard_tests[] = {
        {(char *) "/end-to-end", test_damgard_end_to_end, NULL, NULL, MUNIT_TEST_OPTION_NONE, damgard_params},
        {NULL,                   NULL,                    NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL}
};

MunitSuite damgard_suite = {
        (char *) "/innerprod/fullysec/damgard", simple_ip_damgard_tests, NULL, 1, MUNIT_SUITE_OPTION_NONE
};
